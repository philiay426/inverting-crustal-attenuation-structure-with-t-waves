#!/usr/bin/env python3
"""多台站预处理 + 五图（可选 1:1）— plotting-standards 规范（2026-09-15 版）

用法（在事件目录的 process_flow/ 下跑）：
    python process_image.py --only SA16,SA17 --figs 1_raw,packet,3_spectrogram \
           --packet --square --xmax 2000 --skip-resp

复用步骤:
  1. 复制到新事件目录的 process_flow/
  2. 改 EVENT / BASE（★ 两处路径必须同时更新：/mnt/e/ 与 /mnt/f/ 不一致是常见 bug）
  3. 改 UTCDateTime 到发震时刻（仅缺 `o` 字段时才会用到）
  4. 确保 stations.csv 与 raw/（含 raw/SACPZ/）已就绪

图型（--figs 取子集）:
  1_raw              原始波形（ZRT，保留 DC）
  2_remove_response  去响应后速度波形（未带通）
  bandpass_2-8Hz     2–8 Hz 带通 + 红线包络（波形为主）
  packet             --packet 启用：2–8 Hz 包络为主角（红填充+轮廓，淡灰波形背景）
  3_spectrogram      0–10 Hz 频谱（必须用去响应数据，否则 0–2 / 8–10 Hz 空白）

常用开关:
  --only A,B         只处理指定台站
  --skip-resp        不重写 remove_response/*.sac
  --square           图件 1:1 —— raw / decon / bandpass / packet / spectrogram 五类都支持
  --xmax <s>         raw / packet 横轴上界（默认 T-S+500）
  --env-smooth <s>   包络平滑窗长；0 = 裸 |Hilbert|（Seizmo / 论文口径）
  --env-mode mean|max|lp   包络后处理方式（lp = 0.25 Hz 低通，论文口径）
  --out-suffix _x    文件名后缀，避免覆盖既有产物
  --tt-both          同一张图并列两套到时口径
  --archive          处理归档的非目标区域台站（输出到 archive_*/sacprocess/）

边缘情况:
  - SAC 头缺少 stla/stlo → 从 SACPZ 解析台站坐标
  - sac.mag 常为 -12345 → 震级需另查
  - 多 location code (00/10/20) → 自动选第一个匹配分量
  - SAC `o` 有效时 `get_t0_offset()` 返回发震时刻在数组中的位置（`o-b`）；缺失时用 UTCDateTime fallback
  - 图1 保留 DC 不 demean；图2–4 用去响应数据（已 demean+detrend）
  - 文件名时间解析: parts[-4]=year, parts[-3]=julday（非反！），parts[-2] 可能带 _proc 需 strip"""

# ============================== ★ 修改点 1: 事件名和路径（改这两行 + 下方 UTCDateTime）
import sys, os, glob, re, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from obspy import read
from scipy.signal import butter, filtfilt, hilbert, spectrogram as _sp_spectrogram
from scipy.ndimage import uniform_filter1d

EVENT = "1997-04-28_mw66-prince-edward"                                    # ★ 事件目录名
BASE = "/mnt/e/sharedate/case/zhou_2021_test/1997-04-28_mw66-prince-edward"  # ★ 事件根目录
RAW_DIR = os.path.join(BASE, "raw")
SACPZ_DIR = os.path.join(RAW_DIR, "SACPZ")

# ★ 步0 数据目录（2026-09-28 规范）：新事件 = `process_flow/0_data_processing/`，
#   旧事件（2026-09-28 前处理） = `<event>/sacprocess/` —— 自动探测，两种布局都能跑
def _detect_stage0(base):
    new = os.path.join(base, "process_flow", "0_data_processing")
    old = os.path.join(base, "sacprocess")
    if os.path.isdir(old) and not os.path.isdir(new):
        return old
    return new
STAGE0 = _detect_stage0(BASE)
OUT_RESP_DIR = os.path.join(STAGE0, "remove_response")
OUT_IMG_DIR = os.path.join(STAGE0, "waveform_image")
OUT_ENV_DIR = os.path.join(STAGE0, "envelope")
CSV_PATH = os.path.join(BASE, "process_flow", "1_station_distribution", "stations.csv")

# --archive 模式：处理归档的非马达加斯加台站（南非大陆等），
# 读 all74 stations.csv，输出到 archive_non_madagascar/sacprocess/（归档目录沿用旧布局）
def _setup_paths(archive=False):
    global RAW_DIR, SACPZ_DIR, OUT_RESP_DIR, OUT_IMG_DIR, OUT_ENV_DIR, CSV_PATH
    if archive:
        ARCH = os.path.join(BASE, "archive_non_madagascar")
        RAW_DIR = os.path.join(ARCH, "raw")
        SACPZ_DIR = os.path.join(RAW_DIR, "SACPZ")
        OUT_RESP_DIR = os.path.join(ARCH, "sacprocess", "remove_response")
        OUT_IMG_DIR = os.path.join(ARCH, "sacprocess", "waveform_image")
        OUT_ENV_DIR = os.path.join(ARCH, "sacprocess", "envelope")
        CSV_PATH = os.path.join(ARCH, "process_flow", "1_station_distribution_all74", "stations.csv")
    os.makedirs(OUT_RESP_DIR, exist_ok=True)
    os.makedirs(OUT_IMG_DIR, exist_ok=True)
    os.makedirs(OUT_ENV_DIR, exist_ok=True)
    print(f'[步0 数据目录] {os.path.dirname(OUT_RESP_DIR)}')

_setup_paths(archive=False)

# ============================== 预处理参数
TAPER_MAX_PERCENTAGE = 0.05
PRE_FILT = (0.001, 0.005, 45.0, 50.0)
WATER_LEVEL = 60.0
BP_FREQMIN, BP_FREQMAX, BP_CORNERS = 2.0, 8.0, 4
ENVELOPE_SMOOTH_S = 6.0      # 包络平滑窗长（秒）；≤0 → 不做平滑（裸 |Hilbert| 包络 = Seizmo/论文口径）
ENV_SUFFIX = ''              # 包络 npz / 图件文件名后缀（--out-suffix；导出试验版本用，避免覆盖既有产物）
FIG_XMAX = None              # 图件 x 轴上界（秒）；None → T-S+500（--xmax，作用于 raw / packet 图）
FIG_SQUARE = False           # 图件宽高比 1:1（--square，作用于 raw / decon / bandpass / packet / spectrogram）
ENV_SMOOTH_MODE = 'mean'     # 包络后处理：'mean'=矩形平均｜'max'=滑动最大值｜'lp'=0.25 Hz 低通（论文口径）
ENVELOPE_LP_HZ = 0.25        # --env-mode lp 时的低通截止频率（Hz）
T_TAIL = 500

# ============================== rcParams — plotting-standards §一
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9,
    'axes.labelsize': 10, 'axes.titlesize': 10,
    'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'legend.fontsize': 7,
    'axes.linewidth': 0.6,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

FIG_W = 19 / 2.54; FIG_H = FIG_W * 0.55
ARRIVAL_COLORS = {'P': '#4daf4a', 'S': '#ff9848', 'T-P': '#4d9fff', 'T-S': '#ff4d4d'}

# ============================== ★ 修改点 2: 时间偏移
def get_t0_offset(tr):
    """返回发震时刻在数据数组中的偏移（秒，通常为+600）。"""
    sac = getattr(tr.stats, 'sac', None)
    o = getattr(sac, 'o', None) if sac is not None else None
    b = getattr(sac, 'b', 0.0) if sac is not None else 0.0
    if o is not None and 0 < abs(o) < 1e6:
        # b-o 是数据起点相对发震的偏移；本函数返回其相反数，供数组切片使用
        return float(o - b)
    # 缺 o 时由事件时间解析回退；模板需按事件修改该时间
    from obspy import UTCDateTime
    origin = UTCDateTime(1997, 4, 28, 12, 7, 38)  # ★ 改为实际发震时刻
    return float(origin - tr.stats.starttime)

# ============================== SACPZ / SAC 查找 / CSV
def parse_sacpz(fp):
    """解析 SACPZ — 只取第一个 ZEROS/POLES 块（部分文件重复写 3 遍，全读会 21 zeros+33 poles → simulate 振幅爆炸）"""
    zeros, poles = [], []
    const = 1.0; a0 = None; sens = None; iz = ip = False
    seen_block = False
    with open(fp) as f:
        for line in f:
            s = line.strip()
            if s.startswith('*'):
                if 'SENSITIVITY' in s.upper():
                    try: sens = float(re.findall(r'[-+]?\d+\.?\d*', s.split(':')[1])[0])
                    except: pass
                if s.upper().lstrip('* ').startswith('A0'):
                    try: a0 = float(s.split('=')[1].split()[0])
                    except: pass
                continue
            if s.startswith('ZEROS'):
                if seen_block: break  # 只取第一个块
                iz = True; ip = False; seen_block = True; continue
            if s.startswith('POLES'): ip = True; iz = False; continue
            if s.startswith('CONSTANT'):
                const = float(s.split()[1]); iz = ip = False; continue
            if iz and s:
                p = s.split()
                if len(p) >= 2: zeros.append(complex(float(p[0]), float(p[1])))
            elif ip and s:
                p = s.split()
                if len(p) >= 2: poles.append(complex(float(p[0]), float(p[1])))
    if a0 and sens: pass
    elif const and sens: a0 = const / sens
    else: a0 = const; sens = 1.0
    return zeros, poles, a0, sens

def remove_response_obs(tr, paz_file):
    z, p, a0, s = parse_sacpz(paz_file)
    paz = {'zeros': z, 'poles': p, 'gain': a0, 'sensitivity': s}
    tr_out = tr.copy()
    tr_out.data = tr_out.data.astype(np.float64)
    tr_out.detrend('demean'); tr_out.detrend('linear')
    tr_out.taper(max_percentage=TAPER_MAX_PERCENTAGE)
    tr_out.simulate(paz_remove=paz, remove_sensitivity=True,
                    water_level=WATER_LEVEL, pre_filt=PRE_FILT)
    return tr_out

def load_stations():
    """读 stations.csv 台站列表。

    ⚠️ 2026-09-02: 必须用 DictReader 按列名解析，勿用固定索引——
       stations.csv 已从 13 列升级为 15 列（新增 snr_r/snr_t），
       旧代码 p[12] 取到的是 lon 而非 baz，R/T 旋转会错误。
    """
    import csv as _csv
    sts = []
    with open(CSV_PATH) as f:
        for r in _csv.DictReader(line for line in f if not line.startswith('#')):
            if 'sta' not in r or not r['sta']:
                continue
            sts.append({'sta': r['sta'], 'az': float(r['az']), 'dist': float(r['dist']),
                        'land_km': float(r['land_km']), 'ocean_km': float(r['ocean_km']),
                        'tt_p': float(r['tt_p']), 'tt_s': float(r['tt_s']),
                        'tt_tp': float(r['tt_tp']), 'tt_ts': float(r['tt_ts']),
                        'baz': float(r['baz']),
                        'lon': float(r['lon']) if r.get('lon') else None,
                        'lat': float(r['lat']) if r.get('lat') else None})
    return sts

def find_sac_triplet(sta):
    r = {}
    for ch in ['BHZ', 'BHN', 'BHE']:
        m = glob.glob(os.path.join(RAW_DIR, f"*{sta}*{ch}*.SAC"))
        if m: r[ch] = m[0]
    return r if len(r) == 3 else None

def find_sacpz(sta, ch):
    for pat in [f"SACPZ.*.{sta}..{ch}", f"SACPZ..{sta}..{ch}",
                f"SACPZ.*.{sta}.*.{ch}", f"SACPZ..{sta}.*.{ch}"]:
        m = glob.glob(os.path.join(SACPZ_DIR, pat))
        if m: return m[0]
    return None

# ============================== 绘图辅助
def setup_axes(fig, title):
    axs = fig.subplots(3, 1, sharex=True)
    fig.subplots_adjust(hspace=0.08, top=0.92, bottom=0.10, left=0.07, right=0.94)
    fig.suptitle(title, fontsize=10, fontweight='bold')
    return axs

def get_zrt(trs, t_min, t_max, baz):
    sr = trs['BHZ'].stats.sampling_rate
    t0 = get_t0_offset(trs['BHZ'])
    i1 = int((t_min + t0) * sr); i2 = int((t_max + t0) * sr)
    i2 = min(i2, len(trs['BHZ'].data))
    dt = np.arange(i2 - i1) / sr + (i1 / sr - t0)
    z = trs['BHZ'].data[i1:i2].astype(np.float64)
    n = trs['BHN'].data[i1:i2].astype(np.float64)
    e = trs['BHE'].data[i1:i2].astype(np.float64)
    ml = min(len(z), len(n), len(e))
    dt, z, n, e = dt[:ml], z[:ml], n[:ml], e[:ml]
    # 2026-09-02: 统一 obspy 约定（R 背离震源），与单台偏振分析一致
    from obspy.signal.rotate import rotate_ne_rt
    r, t = rotate_ne_rt(n, e, float(baz))
    return dt, z, r, t

# ============================== 四图
def plot_raw(trs_raw, arrivals, sta, baz):
    t_min = 0
    t_max = float(FIG_XMAX) if FIG_XMAX else arrivals['T-S'] + T_TAIL
    t, z, r, tc = get_zrt(trs_raw, t_min, t_max, baz)
    p98 = np.percentile(np.abs(np.concatenate([z, r, tc])), 98) or 1.0
    y_max = 1.05 * max(np.max(np.abs(z / p98)), np.max(np.abs(r / p98)), np.max(np.abs(tc / p98)), 1.0)
    figsize = (FIG_W, FIG_W) if FIG_SQUARE else (FIG_W, FIG_H)   # 1:1 方图
    fig = plt.figure(figsize=figsize, dpi=150)
    axs = setup_axes(fig, f'{sta} — Raw Waveforms  (ZRT)')
    for ax, (l, d) in zip(axs, [('Z', z), ('R', r), ('T', tc)]):
        ax.plot(t, d / p98, 'k', lw=0.3)
        ax.text(0.985, 0.88, l, transform=ax.transAxes, fontsize=10,
                fontweight='bold', color='black', va='top', ha='right')
        ax.set_xlim(t_min, t_max); ax.set_ylim(-y_max, y_max)
        ax.tick_params(labelsize=8, length=3)
        for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
        ax.axhline(0, color='gray', lw=0.3)
        ax.axvline(x=0, color='gray', lw=0.5, ls=':', alpha=0.5)
    axs[-1].set_xlabel('Time from origin (s)', fontsize=10, fontweight='bold')
    axs[1].set_ylabel('Norm. Amp.', fontsize=9, fontweight='bold')
    return fig

def plot_decon(trs_d, arrivals, sta, baz):
    t_min, t_max = 0, arrivals['T-S'] + T_TAIL
    t, z, r, tc = get_zrt(trs_d, t_min, t_max, baz)
    p98 = np.percentile(np.abs(np.concatenate([z, r, tc])), 98) or 1.0
    y_max = 1.05 * max(np.max(np.abs(z / p98)), np.max(np.abs(r / p98)), np.max(np.abs(tc / p98)), 1.0)
    figsize = (FIG_W, FIG_W) if FIG_SQUARE else (FIG_W, FIG_H)   # 1:1 方图
    fig = plt.figure(figsize=figsize, dpi=150)
    axs = setup_axes(fig, f'{sta} — Instrument Response Removed  (ZRT)')
    for i, (ax, (l, d)) in enumerate(zip(axs, [('Z', z), ('R', r), ('T', tc)])):
        ax.plot(t, d / p98, 'k', lw=0.3)
        ax.text(0.985, 0.88, l, transform=ax.transAxes, fontsize=10,
                fontweight='bold', color='black', va='top', ha='right')
        ax.set_xlim(t_min, t_max); ax.set_ylim(-y_max, y_max)
        ax.tick_params(labelsize=8, length=3)
        for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
        ax.axhline(0, color='gray', lw=0.3)
        ax.axvline(x=0, color='gray', lw=0.5, ls=':', alpha=0.5)
        for j, (ph, tt, clr) in enumerate(arrivals['lines']):
            ax.axvline(x=tt, color=clr, lw=1.0, ls='--', alpha=0.85,
                       label=ph if i == 0 else None)
        if i == 0:
            ax.legend(loc='upper right', fontsize=7, frameon=False, ncol=4,
                      handlelength=1.2, handletextpad=0.4, columnspacing=0.8,
                      bbox_to_anchor=(0.985, 1.18))
    axs[-1].set_xlabel('Time from origin (s)', fontsize=10, fontweight='bold')
    axs[1].set_ylabel('Norm. Amp.', fontsize=9, fontweight='bold')
    return fig

def plot_bandpass(trs_d, arrivals, sta, baz):
    sr = trs_d['BHZ'].stats.sampling_rate
    t_min, t_max = 0, arrivals['T-S'] + T_TAIL
    z = trs_d['BHZ'].data.astype(np.float64)
    n = trs_d['BHN'].data.astype(np.float64)
    e = trs_d['BHE'].data.astype(np.float64)
    ml = min(len(z), len(n), len(e)); z, n, e = z[:ml], n[:ml], e[:ml]
    # 2026-09-02: 统一 obspy 约定
    from obspy.signal.rotate import rotate_ne_rt
    r, tc = rotate_ne_rt(n, e, float(baz))
    b, a = butter(4, [BP_FREQMIN, BP_FREQMAX], btype='band', fs=sr)
    comps = [('Z', z), ('R', r), ('T', tc)]
    filt = {l: filtfilt(b, a, d) for l, d in comps}
    # 包络 = |解析信号|（Hilbert）。ENVELOPE_SMOOTH_S > 0：再做矩形滑动平均（可视化平滑）；
    # ENVELOPE_SMOOTH_S ≤ 0：保留裸包络 —— 与 Seizmo envelope.m / 论文口径一致（真包络无平滑）。
    if ENVELOPE_SMOOTH_S > 0:
        ws = max(int(ENVELOPE_SMOOTH_S * sr), 1); k = np.ones(ws) / ws
        env = {l: np.convolve(np.abs(hilbert(filt[l])), k, mode='same') for l, d in comps}
    else:
        env = {l: np.abs(hilbert(filt[l])) for l, d in comps}
    t0 = get_t0_offset(trs_d['BHZ'])
    # 保存与预处理图完全一致的 Z/R/T 包络，剖面图直接读取此数组
    env_file = os.path.join(OUT_ENV_DIR, f'{sta}_envelope_2-8Hz{ENV_SUFFIX}.npz')
    time_full = np.arange(ml, dtype=np.float64) / sr - t0
    np.savez_compressed(
        env_file,
        time=time_full,
        envelope_z=env['Z'][:ml], envelope_r=env['R'][:ml], envelope_t=env['T'][:ml],
        waveform_z=filt['Z'][:ml], waveform_r=filt['R'][:ml], waveform_t=filt['T'][:ml],
        delta=float(1.0 / sr), t_off=float(t0),
        band_low=float(BP_FREQMIN), band_high=float(BP_FREQMAX),
        smooth_seconds=float(ENVELOPE_SMOOTH_S),
    )
    i1 = max(0, int((t_min + t0) * sr)); i2 = min(int((t_max + t0) * sr), ml)
    t = np.arange(i2 - i1) / sr + (i1 / sr - t0)
    all_wave = np.concatenate([np.abs(filt[l][i1:i2]) for l, _ in comps])
    ep98 = np.percentile(all_wave, 98) or 1.0
    # 纵轴定标（2026-09-02 更新）: P 波初动常远强于 T 波（AMPY ~12×ep98），
    # 旧版用 99% 分位会把 P 波截断（P_peak 6-12×ep98 >> y轴 1.5×ep98）。
    # 现在以 max(整体99%分位, P/S窗内峰值) 定标，保证 P 波完整显示，
    # 同时避免孤立毛刺撑轴。
    all_norm = np.concatenate(
        [np.abs(filt[l][i1:i2]) / ep98 for l, _ in comps] +
        [env[l][i1:i2] / ep98 for l, _ in comps])
    y_robust = float(np.percentile(all_norm, 99))
    p_peak = 0.0
    if arrivals.get('P') is not None and arrivals.get('S') is not None:
        ip1 = max(i1, int((float(arrivals['P']) - 10 + t0) * sr))
        ip2 = min(i2, int((float(arrivals['S']) + 60 + t0) * sr))
        if ip2 > ip1:
            p_peak = float(np.max(np.concatenate(
                [np.abs(filt[l][ip1:ip2]) / ep98 for l, _ in comps])))
    y_max = 1.05 * max(y_robust, p_peak, 1.0)
    figsize = (FIG_W, FIG_W) if FIG_SQUARE else (FIG_W, FIG_H)   # 1:1 方图
    fig = plt.figure(figsize=figsize, dpi=150)
    axs = setup_axes(fig, f'{sta} — 2–8 Hz Bandpass  (ZRT)')
    for i, (l, d) in enumerate(comps):
        ax = axs[i]; y = filt[l][i1:i2]
        ax.plot(t, y / ep98, color='#555555', lw=0.3, alpha=0.5)
        ax.plot(t, env[l][i1:i2] / ep98, color='#d62728', lw=0.8, alpha=0.8)
        ax.text(0.985, 0.88, l, transform=ax.transAxes, fontsize=10,
                fontweight='bold', color='black', va='top', ha='right')
        ax.set_xlim(t_min, t_max); ax.set_ylim(-y_max, y_max); ax.tick_params(labelsize=8, length=3)
        for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
        ax.axhline(0, color='gray', lw=0.3)
        ax.axvline(x=0, color='gray', lw=0.5, ls=':', alpha=0.5)
        for j, (ph, tt, clr) in enumerate(arrivals['lines']):
            ax.axvline(x=tt, color=clr, lw=0.6, ls='--', alpha=0.85,
                       label=f'{ph}={tt:.0f}s' if i == 0 else None)
        if i == 0:
            ax.legend(loc='upper left', fontsize=7, frameon=False, ncol=4,
                      handlelength=1.2, handletextpad=0.4, columnspacing=0.8,
                      bbox_to_anchor=(0.01, 1.15))
    axs[-1].set_xlabel('Time from origin (s)', fontsize=10, fontweight='bold')
    axs[1].set_ylabel('Norm. Amp.', fontsize=9, fontweight='bold')
    return fig

def plot_packet(trs_d, arrivals, sta, baz):
    """波包显示（--packet）：以 2-8 Hz 包络为主角，突出 T-P / T-S 波包的纺锤形轮廓。
    淡灰 = 2-8 Hz 带通波形（背景参考）；红色填充+轮廓 = 包络（= 波包形状本身）。
    包络口径与 profile npz 一致：ENVELOPE_SMOOTH_S>0 做矩形平均，<=0 为裸 |Hilbert|。
    """
    sr = trs_d['BHZ'].stats.sampling_rate
    t_min = 0
    t_max = float(FIG_XMAX) if FIG_XMAX else arrivals['T-S'] + T_TAIL
    z = trs_d['BHZ'].data.astype(np.float64)
    n = trs_d['BHN'].data.astype(np.float64)
    e = trs_d['BHE'].data.astype(np.float64)
    ml = min(len(z), len(n), len(e)); z, n, e = z[:ml], n[:ml], e[:ml]
    # 2026-09-02: 统一 obspy 约定
    from obspy.signal.rotate import rotate_ne_rt
    r, tc = rotate_ne_rt(n, e, float(baz))
    b, a = butter(4, [BP_FREQMIN, BP_FREQMAX], btype='band', fs=sr)
    comps = [('Z', z), ('R', r), ('T', tc)]
    filt = {l: filtfilt(b, a, d) for l, d in comps}
    env_raw = {l: np.abs(hilbert(filt[l])) for l, d in comps}     # 裸包络（|解析信号|，不平滑）
    ws = max(int(ENVELOPE_SMOOTH_S * sr), 1) if ENVELOPE_SMOOTH_S > 0 else 0
    if ENV_SMOOTH_MODE == 'lp':
        # 论文口径：包络 → 0.25 Hz 零相位低通（F-K 输入用的同一处理）
        bb, aa = butter(4, ENVELOPE_LP_HZ / (0.5 * sr), btype='low')
        env = {l: filtfilt(bb, aa, env_raw[l]) for l in env_raw}
    elif ENV_SMOOTH_MODE == 'max' and ws > 0:
        # 滑动最大值：严格 ≥ |波形|（上包络），不压低波包峰值
        from scipy.ndimage import maximum_filter1d
        env = {l: maximum_filter1d(env_raw[l], size=ws) for l in env_raw}
    elif ws > 0:
        k = np.ones(ws) / ws
        env = {l: np.convolve(env_raw[l], k, mode='same') for l in env_raw}
    else:
        env = env_raw
    t0 = get_t0_offset(trs_d['BHZ'])
    i1 = max(0, int((t_min + t0) * sr)); i2 = min(int((t_max + t0) * sr), ml)
    t = np.arange(i2 - i1) / sr + (i1 / sr - t0)
    ep98 = np.percentile(np.concatenate([np.abs(filt[l][i1:i2]) for l, _ in comps]), 98) or 1.0
    y_max = 1.05 * max(float(np.max(np.concatenate(
        [np.abs(filt[l][i1:i2]) / ep98 for l, _ in comps] +
        [env[l][i1:i2] / ep98 for l, _ in comps]))), 1.0)
    figsize = (FIG_W, FIG_W) if FIG_SQUARE else (FIG_W, FIG_H)   # 1:1 方图
    fig = plt.figure(figsize=figsize, dpi=150)
    axs = setup_axes(fig, f'{sta} — T-wave Packets  (2–8 Hz envelope, ZRT)')
    for i, (l, _) in enumerate(comps):
        ax = axs[i]
        ax.plot(t, filt[l][i1:i2] / ep98, color='#999999', lw=0.3, alpha=0.35)
        # 2026-09-13: 去掉填充阴影，只保留包络轮廓线
        ax.plot(t, env[l][i1:i2] / ep98, color='#d62728', lw=1.0)
        ax.text(0.985, 0.88, l, transform=ax.transAxes, fontsize=10,
                fontweight='bold', color='black', va='top', ha='right')
        ax.set_xlim(t_min, t_max); ax.set_ylim(-y_max, y_max); ax.tick_params(labelsize=8, length=3)
        for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
        ax.axhline(0, color='gray', lw=0.3)
        for j, (ph, tt, clr) in enumerate(arrivals['lines']):
            ax.axvline(x=tt, color=clr, lw=0.6, ls='--', alpha=0.85,
                       label=f'{ph}={tt:.0f}s' if i == 0 else None)
        for ph, tt, clr in arrivals.get('lines2', []):      # 第二套：论文恒速式（点线）
            ax.axvline(x=tt, color=clr, lw=0.9, ls=':', alpha=0.75)
        if i == 0:
            ax.legend(loc='upper left', fontsize=7, frameon=False, ncol=4,
                      handlelength=1.2, handletextpad=0.4, columnspacing=0.8,
                      bbox_to_anchor=(0.01, 1.15))
    l2 = arrivals.get('lines2', []) or []
    if len(l2) == 2:      # 两套到时对照说明（英文，图内不出现中文）
        axs[0].text(0.985, 0.66,
                    'dashed = ray theory (TauP): T-P %.0f / T-S %.0f s\n'
                    'dotted = paper const. vel.:  T-P %.0f / T-S %.0f s'
                    % (arrivals['T-P'], arrivals['T-S'], l2[0][1], l2[1][1]),
                    transform=axs[0].transAxes, ha='right', va='top', fontsize=7)
    axs[-1].set_xlabel('Time from origin (s)', fontsize=10, fontweight='bold')
    axs[1].set_ylabel('Norm. Amp.', fontsize=9, fontweight='bold')
    return fig

def plot_spectrogram(trs_d, arrivals, sta, baz):
    sr = trs_d['BHZ'].stats.sampling_rate; t0 = get_t0_offset(trs_d['BHZ'])
    t_min = 0
    t_max = float(FIG_XMAX) if FIG_XMAX else arrivals['T-S'] + T_TAIL
    z = trs_d['BHZ'].data.astype(np.float64)
    n = trs_d['BHN'].data.astype(np.float64)
    e = trs_d['BHE'].data.astype(np.float64)
    ml = min(len(z), len(n), len(e)); z, n, e = z[:ml], n[:ml], e[:ml]
    # 2026-09-02: 统一 obspy 约定
    from obspy.signal.rotate import rotate_ne_rt
    r, tc = rotate_ne_rt(n, e, float(baz))
    comps = [('Z', z), ('R', r), ('T', tc)]
    nperseg = min(int(sr * 12.8), max(256, ml // 8))
    noverlap = int(nperseg * 0.9375); nfft = nperseg * 2
    sc = {}; av = []
    for l, d in comps:
        fs, ts, Sxx = _sp_spectrogram(d, fs=sr, window='hann', nperseg=nperseg,
                                       noverlap=noverlap, nfft=nfft,
                                       scaling='spectrum', mode='psd')
        Sdb = 10 * np.log10(np.maximum(Sxx, 1e-30)); fm = fs <= 10
        sc[l] = (fs, ts, Sdb, fm); av.append(Sdb[fm, :].ravel())
    av = np.concatenate(av); vmin, vmax = np.percentile(av, 2), np.percentile(av, 95)
    if vmax <= vmin: vmax = vmin + 10
    fw = 19.0 / 2.54; fh = (fw if FIG_SQUARE else fw * 0.62)   # --square → 1:1
    fig = plt.figure(figsize=(fw, fh), dpi=150)
    axs = fig.subplots(3, 1, sharex=True)
    fig.subplots_adjust(hspace=0.12, top=0.90, bottom=0.12, left=0.09, right=0.90)
    fig.suptitle(f'{sta} — Spectrogram  0–10 Hz', fontsize=10, fontweight='bold', y=0.97)
    im = None
    for i, (l, d) in enumerate(comps):
        ax = axs[i]; fs, ts, Sdb, fm = sc[l]
        im = ax.pcolormesh(ts - t0, fs[fm], Sdb[fm, :], cmap='viridis',
                            vmin=vmin, vmax=vmax, shading='auto')
        ax.text(0.985, 0.90, l, transform=ax.transAxes, fontsize=9,
                fontweight='bold', color='white', va='top', ha='right')
        ax.set_xlim(t_min, t_max); ax.set_ylim(0, 10)
        ax.xaxis.set_major_locator(MultipleLocator(500))
        ax.tick_params(labelsize=8, length=3)
        for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
        ax.axvline(x=0, color='gray', lw=0.5, ls=':', alpha=0.5)
        for j, (ph, tt, clr) in enumerate(arrivals['lines']):
            ax.axvline(x=tt, color=clr, lw=1.0, ls='--', alpha=0.85,
                       label=ph if i == 0 else None)
        for ph, tt, clr in arrivals.get('lines2', []):   # 第二套：论文恒速式（点线）
            ax.axvline(x=tt, color=clr, lw=1.2, ls=':', alpha=0.8)
        if i == 0:
            ax.legend(loc='upper right', fontsize=7, frameon=False, ncol=4,
                      handlelength=1.2, handletextpad=0.4, columnspacing=0.8,
                      bbox_to_anchor=(0.985, 1.18))
    axs[-1].set_xlabel('Time from origin (s)', fontsize=9, fontweight='bold')
    axs[1].set_ylabel('Frequency (Hz)', fontsize=9, fontweight='bold')
    ca = fig.add_axes([0.915, 0.12, 0.014, 0.78])
    cb = fig.colorbar(im, cax=ca, orientation='vertical')
    cb.set_label('Power (dB)', fontsize=9, fontweight='bold')
    cb.ax.tick_params(labelsize=8); cb.outline.set_linewidth(0.6)
    return fig

# ============================== 主流程
def main():
    import argparse
    ArgumentParser(description='多台站预处理 + 五图（可选 1:1）')
    ap.add_argument("--archive", action="store_true",
                    help="处理归档的非马达加斯加台站（南非大陆等），输出到 archive_non_madagascar/sacprocess/")
    ap.add_argument("--only", default=None, help="只处理指定台站（逗号分隔），默认全部")
    ap.add_argument("--env-smooth", type=float, default=None,
                    help="包络平滑窗长（秒）；0 = 不做平滑（裸 |Hilbert| 包络）")
    ap.add_argument("--out-suffix", default="",
                    help="输出图件/包络 npz 的文件名后缀，避免覆盖既有产物")
    ap.add_argument("--figs", default=None,
                    help="只出指定图（逗号分隔）：1_raw,2_remove_response,bandpass_2-8Hz,3_spectrogram,packet")
    ap.add_argument("--skip-resp", action="store_true",
                    help="不重写 sacprocess/remove_response/*.sac（只出图/包络）")
    ap.add_argument("--packet", action="store_true",
                    help="额外输出波包显示图（以 2-8 Hz 包络为主角，突出 T-P/T-S 波包）")
    ap.add_argument("--tt-paper", action="store_true",
                    help="T-P/T-S 改用论文式(1)恒速公式：L_ocean/1.48 + L_land/6.0(P)、3.6(S)")
    ap.add_argument("--tt-both", action="store_true",
                    help="同时画两套到时：实线=射线理论(stations.csv)、点线=论文恒速式")
    ap.add_argument("--xmax", type=float, default=None,
                    help="packet 图 x 轴上界（秒），默认 T-S+500")
    ap.add_argument("--square", action="store_true", help="图件宽高比 1:1（raw/decon/bandpass/packet/spectrogram）")
    ap.add_argument("--env-mode", choices=['mean', 'max', 'lp'], default='mean',
                    help="包络后处理：mean=矩形滑动平均｜max=滑动最大值｜lp=0.25 Hz 低通（论文口径）")
    ap.add_argument("--env-lp-hz", type=float, default=0.25,
                    help="--env-mode lp 时的低通截止频率（Hz，默认 0.25）")
    args = ap.parse_args()
    global ENVELOPE_SMOOTH_S, ENV_SUFFIX, FIG_XMAX, FIG_SQUARE, ENV_SMOOTH_MODE, ENVELOPE_LP_HZ
    if args.env_smooth is not None: ENVELOPE_SMOOTH_S = float(args.env_smooth)
    ENV_SUFFIX = args.out_suffix or ''
    FIG_XMAX = args.xmax
    FIG_SQUARE = bool(args.square)
    ENV_SMOOTH_MODE = args.env_mode
    ENVELOPE_LP_HZ = float(args.env_lp_hz)
    only = set(x.strip() for x in args.only.split(',')) if args.only else None
    figs = set(x.strip() for x in args.figs.split(',')) if args.figs else None
    if args.archive:
        _setup_paths(archive=True)
        print("★ archive 模式: 只处理马达加斯加岛范围外的归档台站")
    stations = load_stations()
    print(f"台站数: {len(stations)}")
    for idx, s in enumerate(stations):
        sta = s['sta']
        if only and sta not in only: continue
        # archive 模式只处理岛外台站（all74 CSV 含岛内 25 站，需跳过）
        lon, lat = s.get('lon'), s.get('lat')
        if args.archive and lon is not None and lat is not None:
            if 43.0 <= float(lon) <= 60.0 and -26.0 <= float(lat) <= -11.0:
                continue
        sacs = find_sac_triplet(sta)
        if not sacs: print(f"  [{idx+1}/{len(stations)}] {sta} — ❌ SAC不全"); continue
        baz = s['baz']
        tt_tp, tt_ts = float(s['tt_tp']), float(s['tt_ts'])
        if args.tt_paper:
            # 2026-09-13: 论文式(1) 恒速口径（Zhou et al. 2021, Eq.1）——
            # 海水段 L_ocean/1.48 + 陆地段 L_land/6.0 (T-P) 或 L_land/3.6 (T-S)。
            # 我们默认的 tt_tp/tt_ts 陆地段用 TauP(iasp91)，比论文恒速式早 ~20-25 s。
            tt_tp = float(s['ocean_km']) / 1.48 + float(s['land_km']) / 6.0
            tt_ts = float(s['ocean_km']) / 1.48 + float(s['land_km']) / 3.6
        lines2 = []
        if args.tt_both:
            # 第二套：论文式(1) 恒速（用于与射线理论同图对比）
            c_tp = float(s['ocean_km']) / 1.48 + float(s['land_km']) / 6.0
            c_ts = float(s['ocean_km']) / 1.48 + float(s['land_km']) / 3.6
            lines2 = [('T-P', c_tp, ARRIVAL_COLORS['T-P']),
                      ('T-S', c_ts, ARRIVAL_COLORS['T-S'])]
        arrivals = {'P': s['tt_p'], 'S': s['tt_s'], 'T-P': tt_tp, 'T-S': tt_ts,
                    'lines': [('P', s['tt_p'], ARRIVAL_COLORS['P']),
                              ('S', s['tt_s'], ARRIVAL_COLORS['S']),
                              ('T-P', tt_tp, ARRIVAL_COLORS['T-P']),
                              ('T-S', tt_ts, ARRIVAL_COLORS['T-S'])],
                    'lines2': lines2}
        trs_raw = {ch: read(sacs[ch])[0] for ch in ['BHZ', 'BHN', 'BHE']}
        trs_d = {}
        for ch in ['BHZ', 'BHN', 'BHE']:
            paz = find_sacpz(sta, ch)
            if not paz: print(f"  [{idx+1}/{len(stations)}] {sta} — ❌ SACPZ {ch}"); break
            tr_r = remove_response_obs(trs_raw[ch], paz)
            on = os.path.basename(sacs[ch]).replace('.SAC', '_resp.sac')
            if not args.skip_resp:
                tr_r.write(os.path.join(OUT_RESP_DIR, on), format='SAC')
            trs_d[ch] = tr_r
        else:
            print(f"  [{idx+1}/{len(stations)}] {sta} — ✅ "
                  f"P={arrivals['P']:.0f}s T-P={arrivals['T-P']:.0f}s "
                  f"land={s['land_km']:.0f}km")
            sd = os.path.join(OUT_IMG_DIR, sta); os.makedirs(sd, exist_ok=True)
            fig_list = [(plot_raw, '1_raw'), (plot_decon, '2_remove_response'),
                        (plot_bandpass, 'bandpass_2-8Hz'), (plot_spectrogram, '3_spectrogram')]
            if args.packet:
                fig_list.append((plot_packet, 'packet'))
            for pf, sfx in fig_list:
                if figs and sfx not in figs: continue
                fig = pf(trs_raw if sfx == '1_raw' else trs_d, arrivals, sta, baz)
                fig.savefig(os.path.join(sd, f"{sta}_{sfx}{ENV_SUFFIX}.png"),
                            dpi=300, bbox_inches='tight')
                plt.close(fig)
    print(f"\n完成: {len(stations)} 台站")

if __name__ == '__main__': main()
