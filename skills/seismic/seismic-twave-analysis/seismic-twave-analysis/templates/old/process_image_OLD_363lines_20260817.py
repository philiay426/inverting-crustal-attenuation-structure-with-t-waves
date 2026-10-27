#!/usr/bin/env python3
"""多台站预处理 + 四图 — plotting-standards 规范（2026-08-01）

用法：python process_image.py

复用步骤:
  1. 复制到新事件目录
  2. 改 EVENT/BASE（★ 两处路径必须同时更新：/mnt/e/ 和 /mnt/f/ 不一致常见 bug）
  3. 改 UTCDateTime 到发震时刻
  4. 确保 stations.csv 和 raw/ 已就绪

边缘情况:
  - SAC 头缺少 stla/stlo → 脚本从 SACPZ 解析台站坐标
  - sac.mag 常为 -12345 → 震级从 USGS 目录查
  - 多 location code (00/10/20) → 脚本自动选第一个匹配的分量
  - SAC `o` 有效时，`get_t0_offset()` 返回发震时刻在数组中的位置（`o-b`）；缺失时用 UTCDateTime fallback
  - 原始波形 DC 偏移: 图1 保留不 demean，图2-4 用去响应数据（已 demean+detrend）
  - 频谱图必须用去响应数据（非 BP 滤波），否则 0-2Hz 和 8-10Hz 空白
  - 文件名时间解析: parts[-4]=year, parts[-3]=julday（非反！），parts[-2] 可能带 _proc 需 strip"""

# ============================== ★ 修改点 1: 事件名和路径（改这两行 + 下方 UTCDateTime）
import sys, os, glob, re, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from obspy import read
from scipy.signal import butter, filtfilt, hilbert, spectrogram as _sp_spectrogram
from scipy.ndimage import uniform_filter1d

EVENT = "CHANGE_ME"
BASE = f"/mnt/e/sharedate/case/GROUP/{EVENT}"
RAW_DIR = os.path.join(BASE, "raw")
SACPZ_DIR = os.path.join(RAW_DIR, "SACPZ")
OUT_RESP_DIR = os.path.join(BASE, "sacprocess", "remove_response")
OUT_IMG_DIR = os.path.join(BASE, "sacprocess", "waveform_image")
OUT_ENV_DIR = os.path.join(BASE, "sacprocess", "envelope")
CSV_PATH = os.path.join(BASE, "process_flow", "1_station_distribution", "stations.csv")

os.makedirs(OUT_RESP_DIR, exist_ok=True)
os.makedirs(OUT_IMG_DIR, exist_ok=True)
os.makedirs(OUT_ENV_DIR, exist_ok=True)

# ============================== 预处理参数
TAPER_MAX_PERCENTAGE = 0.05
PRE_FILT = (0.001, 0.005, 45.0, 50.0)
WATER_LEVEL = 60.0
BP_FREQMIN, BP_FREQMAX, BP_CORNERS = 2.0, 8.0, 4
ENVELOPE_SMOOTH_S = 6.0
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
    origin = UTCDateTime(2026, 1, 1, 0, 0, 0)  # ★ 改为实际发震时刻
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
    sts = []
    with open(CSV_PATH) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('sta,'): continue
            p = line.split(',')
            if len(p) < 13: continue
            sts.append({'sta': p[0], 'az': float(p[1]), 'dist': float(p[2]),
                        'land_km': float(p[3]), 'ocean_km': float(p[4]),
                        'tt_p': float(p[5]), 'tt_s': float(p[6]),
                        'tt_tp': float(p[7]), 'tt_ts': float(p[8]),
                        'baz': float(p[12])})
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
    t0 = get_t0_offset(trs['BHZ']); br = np.radians(baz)
    i1 = int((t_min + t0) * sr); i2 = int((t_max + t0) * sr)
    i2 = min(i2, len(trs['BHZ'].data))
    dt = np.arange(i2 - i1) / sr + (i1 / sr - t0)
    z = trs['BHZ'].data[i1:i2].astype(np.float64)
    n = trs['BHN'].data[i1:i2].astype(np.float64)
    e = trs['BHE'].data[i1:i2].astype(np.float64)
    ml = min(len(z), len(n), len(e))
    dt, z, n, e = dt[:ml], z[:ml], n[:ml], e[:ml]
    r = n * np.cos(br) + e * np.sin(br)
    t = -n * np.sin(br) + e * np.cos(br)
    return dt, z, r, t

# ============================== 四图
def plot_raw(trs_raw, arrivals, sta, baz):
    t_min, t_max = 0, arrivals['T-S'] + T_TAIL
    t, z, r, tc = get_zrt(trs_raw, t_min, t_max, baz)
    p98 = np.percentile(np.abs(np.concatenate([z, r, tc])), 98)
    fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=150)
    axs = setup_axes(fig, f'{sta} — Raw Waveforms  (ZRT)')
    for ax, (l, d) in zip(axs, [('Z', z), ('R', r), ('T', tc)]):
        ax.plot(t, d / p98, 'k', lw=0.3)
        ax.text(0.985, 0.88, l, transform=ax.transAxes, fontsize=10,
                fontweight='bold', color='black', va='top', ha='right')
        ax.set_xlim(t_min, t_max); ax.set_ylim(-2.5, 2.5)
        ax.tick_params(labelsize=8, length=3)
        for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
        ax.axhline(0, color='gray', lw=0.3)
        ax.axvline(x=0, color='gray', lw=0.5, ls=':', alpha=0.5)
    axs[-1].set_xlabel('Time from origin (s)', fontsize=10, fontweight='bold')
    axs[1].set_ylabel('Counts', fontsize=9, fontweight='bold')
    return fig

def plot_decon(trs_d, arrivals, sta, baz):
    t_min, t_max = 0, arrivals['T-S'] + T_TAIL
    t, z, r, tc = get_zrt(trs_d, t_min, t_max, baz)
    p98 = np.percentile(np.abs(np.concatenate([z, r, tc])), 98)
    fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=150)
    axs = setup_axes(fig, f'{sta} — Instrument Response Removed  (ZRT)')
    for i, (ax, (l, d)) in enumerate(zip(axs, [('Z', z), ('R', r), ('T', tc)])):
        ax.plot(t, d / p98, 'k', lw=0.3)
        ax.text(0.985, 0.88, l, transform=ax.transAxes, fontsize=10,
                fontweight='bold', color='black', va='top', ha='right')
        ax.set_xlim(t_min, t_max); ax.set_ylim(-2.5, 2.5)
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
    axs[1].set_ylabel('Velocity (m/s)', fontsize=9, fontweight='bold')
    return fig

def plot_bandpass(trs_d, arrivals, sta, baz):
    sr = trs_d['BHZ'].stats.sampling_rate
    t_min, t_max = 0, arrivals['T-S'] + T_TAIL
    br = np.radians(baz)
    z = trs_d['BHZ'].data.astype(np.float64)
    n = trs_d['BHN'].data.astype(np.float64)
    e = trs_d['BHE'].data.astype(np.float64)
    ml = min(len(z), len(n), len(e)); z, n, e = z[:ml], n[:ml], e[:ml]
    r = n * np.cos(br) + e * np.sin(br); tc = -n * np.sin(br) + e * np.cos(br)
    b, a = butter(4, [BP_FREQMIN, BP_FREQMAX], btype='band', fs=sr)
    comps = [('Z', z), ('R', r), ('T', tc)]
    filt = {l: filtfilt(b, a, d) for l, d in comps}
    ws = max(int(ENVELOPE_SMOOTH_S * sr), 1); k = np.ones(ws) / ws
    env = {l: np.convolve(np.abs(hilbert(filt[l])), k, mode='same') for l, d in comps}
    t0 = get_t0_offset(trs_d['BHZ'])
    # 保存与预处理图完全一致的 Z/R/T 包络，剖面图直接读取此数组
    env_file = os.path.join(OUT_ENV_DIR, f'{sta}_envelope_2-8Hz.npz')
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
    all_env = np.concatenate([env[l][i1:i2] for l, _ in comps])
    ep98 = np.percentile(all_env, 98)
    fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=150)
    axs = setup_axes(fig, f'{sta} — 2–8 Hz Bandpass  (ZRT)')
    for i, (l, d) in enumerate(comps):
        ax = axs[i]; y = filt[l][i1:i2]
        ax.plot(t, y / ep98, color='#555555', lw=0.3, alpha=0.5)
        ax.plot(t, env[l][i1:i2] / ep98, color='#d62728', lw=0.8, alpha=0.8)
        ax.text(0.985, 0.88, l, transform=ax.transAxes, fontsize=10,
                fontweight='bold', color='black', va='top', ha='right')
        ax.set_xlim(t_min, t_max); ax.set_ylim(-2.5, 2.5); ax.tick_params(labelsize=8, length=3)
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
    axs[1].set_ylabel('Velocity (m/s)', fontsize=9, fontweight='bold')
    return fig

def plot_spectrogram(trs_d, arrivals, sta, baz):
    sr = trs_d['BHZ'].stats.sampling_rate; t0 = get_t0_offset(trs_d['BHZ'])
    t_min, t_max = 0, arrivals['T-S'] + T_TAIL; br = np.radians(baz)
    z = trs_d['BHZ'].data.astype(np.float64)
    n = trs_d['BHN'].data.astype(np.float64)
    e = trs_d['BHE'].data.astype(np.float64)
    ml = min(len(z), len(n), len(e)); z, n, e = z[:ml], n[:ml], e[:ml]
    r = n * np.cos(br) + e * np.sin(br); tc = -n * np.sin(br) + e * np.cos(br)
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
    fw = 19.0 / 2.54; fh = fw * 0.62
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
    stations = load_stations()
    print(f"台站数: {len(stations)}")
    for idx, s in enumerate(stations):
        sta = s['sta']; sacs = find_sac_triplet(sta)
        if not sacs: print(f"  [{idx+1}/{len(stations)}] {sta} — ❌ SAC不全"); continue
        baz = s['baz']
        arrivals = {'P': s['tt_p'], 'S': s['tt_s'], 'T-P': s['tt_tp'], 'T-S': s['tt_ts'],
                    'lines': [('P', s['tt_p'], ARRIVAL_COLORS['P']),
                              ('S', s['tt_s'], ARRIVAL_COLORS['S']),
                              ('T-P', s['tt_tp'], ARRIVAL_COLORS['T-P']),
                              ('T-S', s['tt_ts'], ARRIVAL_COLORS['T-S'])]}
        trs_raw = {ch: read(sacs[ch])[0] for ch in ['BHZ', 'BHN', 'BHE']}
        trs_d = {}
        for ch in ['BHZ', 'BHN', 'BHE']:
            paz = find_sacpz(sta, ch)
            if not paz: print(f"  [{idx+1}/{len(stations)}] {sta} — ❌ SACPZ {ch}"); break
            tr_r = remove_response_obs(trs_raw[ch], paz)
            on = os.path.basename(sacs[ch]).replace('.SAC', '_resp.sac')
            tr_r.write(os.path.join(OUT_RESP_DIR, on), format='SAC'); trs_d[ch] = tr_r
        else:
            print(f"  [{idx+1}/{len(stations)}] {sta} — ✅ "
                  f"P={arrivals['P']:.0f}s T-P={arrivals['T-P']:.0f}s "
                  f"land={s['land_km']:.0f}km")
            sd = os.path.join(OUT_IMG_DIR, sta); os.makedirs(sd, exist_ok=True)
            for pf, sfx in [(plot_raw, '1_raw'), (plot_decon, '2_remove_response'),
                            (plot_bandpass, 'bandpass_2-8Hz'), (plot_spectrogram, '3_spectrogram')]:
                fig = pf(trs_raw if sfx == '1_raw' else trs_d, arrivals, sta, baz)
                fig.savefig(os.path.join(sd, f"{sta}_{sfx}.png"), dpi=300, bbox_inches='tight')
                plt.close(fig)
    print(f"\n完成: {len(stations)} 台站")

if __name__ == '__main__': main()
