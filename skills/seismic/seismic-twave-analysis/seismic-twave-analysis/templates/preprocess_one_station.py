#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步 0 预处理（单台示例）：原始 SAC -> 去响应(速度) -> 2-8 Hz 带通 -> 包络 -> 5 张 QC 图

对应 skill：seismic/twave-crustal-attenuation  §2 步 0
数据来源：事件目录下的原始 SAC（只读，不复制波形）

产物（默认写本步目录）：
  {STA}_1_raw.png             原始三分量（count），横轴 = 相对发震时刻
  {STA}_2_remove_response.png 去响应后速度（m/s）
  {STA}_3_spectrogram.png     去响应 BHZ 的 0-10 Hz 谱图（看 50 Hz / 窄带干扰）
  {STA}_bandpass_2-8Hz.png    2-8 Hz 带通三分量（整段滤，不切窗）
  {STA}_packet.png            包络 + 信号窗/噪声窗（QC 用）
  {STA}_envelope_2-8Hz.npz    time / waveform_z,n,e / envelope_*（下游 F-K、剖面用）

★ 三个本步必须落实的口径（对应 skill §10 陷阱 1/2/3）：
  1) SAC o 字段：0<|o|<1e6 时 t_off = b - o。本事件 o=600、b=0 -> 数据从 -600 s 起；
     而 starttime 恰好等于发震时刻 —— 用 starttime==origin 判断起点会整体错 600 s。
  2) 去响应取 SACPZ 里的**第一组**零极点（文件常把同一组重复写 3 遍，全读会振幅爆炸），
     反卷积得位移；再频域微分得速度（论文口径）。位移版：OUTPUT=DISP。
  3) 滤波在整段上做、再切窗：先切后滤会因窗边瞬变虚高噪声。

用法：
  /home/sy/miniforge3/envs/obspy/bin/python preprocess_one_station.py
  STA=<台站> /home/sy/miniforge3/envs/obspy/bin/python preprocess_one_station.py
"""
import os
import re
import glob
import numpy as np
from obspy import read
from obspy.signal.filter import bandpass
from scipy.signal import hilbert, spectrogram
from scipy.ndimage import uniform_filter1d
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe          # 到时描边（withStroke）

# ==================== 期刊制图标准（seismic/plotting-standards §一） ====================
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],   # 退化链
    "font.size": 9,
    "axes.labelsize": 11, "axes.labelweight": "bold",
    "axes.titlesize": 12,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "legend.fontsize": 8.5,
    "axes.linewidth": 0.6,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,        # 全框 + 四边朝内刻度
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
# 到时配色（skill §一.5 固定值，不自选）
ARRIVAL_COLORS = {"P": "#4daf4a", "S": "#ff9848", "T-P": "#4d9fff", "T-S": "#ff4d4d"}

# ==================== ★ 参数区 ====================
SRC = os.environ.get("SRC", "").strip()                                      # 事件根（必填）
if not SRC:
    raise SystemExit("请用环境变量 SRC 指定事件根目录（其下有 raw/ 与 raw/SACPZ/）")
STA = os.environ.get("STA", "").strip()                                      # 台站名（必填）
if not STA:
    raise SystemExit("请用环境变量 STA 指定台站名，例如 STA=<台站>")
OUT = os.environ.get("OUT") or os.path.dirname(os.path.abspath(__file__))   # 本步目录本身
# 输出分两个子目录（照案例 sacprocess/ 布局）：图 -> waveform_image/，数据 -> envelope/
FIG_OUT = os.environ.get("FIG_OUT") or os.path.join(OUT, "waveform_image")
ENV_OUT = os.environ.get("ENV_OUT") or os.path.join(OUT, "envelope")
PRE_BAND = (2.0, 8.0)                    # Hz，T 波主频带
PRE_FILT = (0.001, 0.005, 45.0, 50.0)    # 去响应预滤波（process_image.py L59）
WATER_LEVEL = 60.0                       # 去卷积水位（L60）
TAPER_MAX = 0.05                         # 去卷积前 taper（L58）
ENV_SMOOTH = 6.0                         # s，包络平滑（本项目 6 s）
OUTPUT = os.environ.get("OUTPUT", "VEL") # VEL 速度 / DISP 位移
DPI = 300
FIG_WAVE = (7.2, 7.2)    # 波形图图幅 1:1（正方形画布）
FIG_SPEC = (7.2, 7.2)    # 时频图也 1:1（三块面板 Z/R/T，每块是宽扁条带）
# 谱图口径（照 process_image.py L418-427）
SPEC_WIN   = 12.8        # s，nperseg = SPEC_WIN * fs（原脚本 sr*12.8）
SPEC_OVERLAP = 0.9375    # noverlap 比例
SPEC_NFFT_X = 2          # nfft = nperseg * 该倍数
SPEC_VMIN_PCT, SPEC_VMAX_PCT = 2.0, 95.0   # 色标范围 = 0-10 Hz dB 的 P2 / P95
# 波形图横轴（相对发震）：从事件 0 时刻，到理论 T-S + XLIM_TAIL 秒
XLIM_T0   = 0.0          # 横轴起点 (s)
XLIM_TAIL = 500.0        # 横轴终点 = tt_ts + XLIM_TAIL (s)
SNR_MIN   = 6.0          # dB，低于此值该台慎用（项目筛台阈值之一）
SNR_GAP   = 10.0         # s，噪声窗末端距 P 的间隔（避免吃进 P 波初动；0 = 贴 P，与旧表口径一致）
# ===================================================

RAW = os.path.join(SRC, "raw")
SACPZ = os.path.join(RAW, "SACPZ")
STACSV = os.path.join(SRC, "process_flow", "1_station_distribution", "stations.csv")
# 每台一个独立子文件夹：waveform_image/<STA>/<STA>_*.png
FIG_STA = os.environ.get("FIG_STA") or os.path.join(FIG_OUT, STA)
os.makedirs(FIG_STA, exist_ok=True)
os.makedirs(ENV_OUT, exist_ok=True)


# ---------------------------------------------------------------- 输入
def load_raw(sta):
    """读三分量原始 SAC"""
    got = {}
    for c in ("BHZ", "BHN", "BHE"):
        m = sorted(glob.glob(os.path.join(RAW, f"*.{sta}.*{c}*.SAC")))
        if not m:
            raise SystemExit(f"缺分量 {c}: {sta}")
        got[c] = read(m[0])[0]
    return got


def t_off_of(tr):
    """数据起点相对发震时刻的秒数（o 字段口径）"""
    o = getattr(tr.stats.sac, "o", 0.0)
    b = getattr(tr.stats.sac, "b", 0.0)
    if 0 < abs(o) < 1e6:
        return b - o, f"b-o = {b:g} - {o:g}"
    return 0.0, "缺 o 字段"


def read_marks(sta):
    """从 stations.csv 取该台直达 P/S 与 T-P/T-S 到时、陆段/海段（按列名读，禁固定索引）

    四个量都在表里：tt_p / tt_s 是 TauP 直达体波；tt_tp / tt_ts 是 T 波（= L_ocean/1.48 + L_land/v）。
    两者相差 600-700 s，图里要同时标，否则只看 T-P/T-S 会与 300-500 s 的大振幅"对不上"。
    """
    import csv
    with open(STACSV) as f:
        for row in csv.DictReader(l for l in f if not l.startswith("#")):
            if row["sta"].strip() == sta:
                return dict(tt_p=float(row["tt_p"]), tt_s=float(row["tt_s"]),
                            tt_tp=float(row["tt_tp"]), tt_ts=float(row["tt_ts"]),
                            land=float(row["land_km"]), ocean=float(row["ocean_km"]),
                            baz=float(row["baz"]))
    raise SystemExit(f"stations.csv 里没有 {sta}")


def find_sacpz(sta, ch):
    for pat in (f"SACPZ.*.{sta}..{ch}", f"SACPZ..{sta}..{ch}",
                f"SACPZ.*.{sta}.*.{ch}", f"SACPZ..{sta}.*.{ch}"):
        m = glob.glob(os.path.join(SACPZ, pat))
        if m:
            return m[0]
    raise SystemExit(f"缺 SACPZ: {sta} {ch}")


# ---------------------------------------------------------------- 去响应
def parse_sacpz(fp):
    """解析 SACPZ，★ 只取第一组 ZEROS/POLES（源码 process_image.py L102-125）"""
    zeros, poles = [], []
    const = 1.0
    a0 = sens = None
    iz = ip = seen = False
    for line in open(fp):
        t = line.strip()
        if t.startswith("*"):
            if "SENSITIVITY" in t.upper():
                try:
                    sens = float(re.findall(r"[-+]?\d+\.?\d*", t.split(":")[1])[0])
                except Exception:
                    pass
            if t.upper().lstrip("* ").startswith("A0"):
                try:
                    a0 = float(t.split("=")[1].split()[0])
                except Exception:
                    pass
            continue
        if t.startswith("ZEROS"):
            if seen:
                break                      # ★ 只取第一块
            iz, ip, seen = True, False, True
            continue
        if t.startswith("POLES"):
            ip, iz = True, False
            continue
        if t.startswith("CONSTANT"):
            const = float(t.split()[1])
            iz = ip = False
            continue
        if iz and t:
            q = t.split()
            if len(q) >= 2:
                zeros.append(complex(float(q[0]), float(q[1])))
        elif ip and t:
            q = t.split()
            if len(q) >= 2:
                poles.append(complex(float(q[0]), float(q[1])))
    if not (a0 and sens):
        a0 = const / sens if (const and sens) else const
        sens = sens or 1.0
    return zeros, poles, a0, sens


def deconvolve(tr, comp, output=OUTPUT):
    """PAZ 反卷积 -> 位移(m)；output='VEL' 再频域微分 -> 速度(m/s)"""
    z, p, a0, sens = parse_sacpz(find_sacpz(STA, comp))
    paz = {"zeros": z, "poles": p, "gain": a0, "sensitivity": sens}
    out = tr.copy()
    out.data = out.data.astype(np.float64)
    out.detrend("demean")
    out.detrend("linear")
    out.taper(max_percentage=TAPER_MAX)
    out.simulate(paz_remove=paz, remove_sensitivity=True,
                 water_level=WATER_LEVEL, pre_filt=PRE_FILT)
    if output == "VEL":
        n = len(out.data)
        freqs = np.fft.rfftfreq(n, d=out.stats.delta)
        out.data = np.fft.irfft(np.fft.rfft(out.data) * 2j * np.pi * freqs, n=n)
    return out


def snr_welch_proj(x, dt, t_off, tt_p, tt_tp, tt_ts, band=(2.0, 8.0), gap=SNR_GAP):
    """T 波 SNR (dB)，**项目口径** —— 照 step1_map_and_table.py L104-132

    信号窗 [tt_tp - 20, tt_ts + 100]；噪声窗取等长的 P 波之前段 [tt_p - L - gap, tt_p - gap]，
    gap=10 s 让窗尾离 P 初动留一段间隔（gap=0 即旧口径，窗尾贴 P，会吃进一点 P 波）；
    Welch PSD 在 2-8 Hz 频带内求和后取 10*log10(能量比)（= 同频带振幅比的 20*log10）。
    该函数在**位移**波形上算，结果应与 stations.csv 的 snr 列（同口径）一致。
    """
    from scipy.signal import welch
    sr = 1.0 / dt
    t = np.arange(len(x)) * dt + t_off
    s1 = int(max(0.0, tt_tp - 20 - t_off) / dt)
    s2 = int(min(t[-1], tt_ts + 100 - t_off) / dt)
    L = (tt_ts + 100) - (tt_tp - 20)
    n1 = int(max(0.0, tt_p - L - gap - t_off) / dt)
    n2 = int(max(0.0, tt_p - gap - t_off) / dt)
    if s2 <= s1 or n2 <= n1:
        return np.nan, (np.nan, np.nan), (np.nan, np.nan)
    fs_, ps = welch(x[s1:s2].astype(float), sr, nperseg=min(256, s2 - s1))
    fn_, pn = welch(x[n1:n2].astype(float), sr, nperseg=min(256, n2 - n1))
    es = np.sum(ps[(fs_ >= band[0]) & (fs_ <= band[1])])
    en = np.sum(pn[(fn_ >= band[0]) & (fn_ <= band[1])])
    val = 10 * np.log10(es / en) if (es > 0 and en > 0) else np.nan
    return val, (t[s1], t[s2 - 1]), (t[n1], t[n2 - 1])


def smooth_env(x, fs, win_s):
    """Hilbert 包络 + 滑动平均（本项目 6 s）"""
    return uniform_filter1d(np.abs(hilbert(x)), size=max(1, int(round(win_s * fs))),
                            mode="nearest")


# ---------------------------------------------------------------- 主流程
def main():
    print(f"=== 步0 单台预处理 | {STA} | 输出量纲 {OUTPUT} ===")
    raw = load_raw(STA)
    dt = raw["BHZ"].stats.delta
    t0, how = t_off_of(raw["BHZ"])
    print(f"起点: t_off = {t0:+.1f} s ({how}) | npts={raw['BHZ'].stats.npts} "
          f"dt={dt:g} s | starttime={raw['BHZ'].stats.starttime}")

    disp, vel = {}, {}
    for c, tr in raw.items():
        d = deconvolve(tr, c, output="DISP")
        disp[c] = d
        if OUTPUT == "VEL":
            n = len(d.data)
            freqs = np.fft.rfftfreq(n, d=d.stats.delta)
            v = d.copy()
            v.data = np.fft.irfft(np.fft.rfft(d.data) * 2j * np.pi * freqs, n=n)
            vel[c] = v
        else:
            vel[c] = d
        print(f"  {c}: 原始 max|A|={np.abs(tr.data).max():.3g} count -> "
              f"位移 {np.abs(d.data).max():.3g} m -> {OUTPUT} {np.abs(vel[c].data).max():.3g}")

    n = min(len(v.data) for v in vel.values())
    t = t0 + np.arange(n) * dt
    z = vel["BHZ"].data[:n].astype(float)
    nn = vel["BHN"].data[:n].astype(float)
    ee = vel["BHE"].data[:n].astype(float)

    fs = 1.0 / dt
    zf = bandpass(z, PRE_BAND[0], PRE_BAND[1], fs, corners=4, zerophase=True)
    nf = bandpass(nn, PRE_BAND[0], PRE_BAND[1], fs, corners=4, zerophase=True)
    ef = bandpass(ee, PRE_BAND[0], PRE_BAND[1], fs, corners=4, zerophase=True)
    env_z = smooth_env(zf, fs, ENV_SMOOTH)

    mk = read_marks(STA)
    tt_p, tt_s, tt_tp, tt_ts = mk["tt_p"], mk["tt_s"], mk["tt_tp"], mk["tt_ts"]
    land, ocean, baz = mk["land"], mk["ocean"], mk["baz"]
    ARRIVALS = [("P", tt_p, ARRIVAL_COLORS["P"]), ("S", tt_s, ARRIVAL_COLORS["S"]),
                ("T-P", tt_tp, ARRIVAL_COLORS["T-P"]), ("T-S", tt_ts, ARRIVAL_COLORS["T-S"])]
    print(f"理论到时: P={tt_p:.1f} s, S={tt_s:.1f} s, T-P={tt_tp:.1f} s, T-S={tt_ts:.1f} s "
          f"| 陆段 {land:g} km, 海段 {ocean:g} km")

    # 包络峰要看两处：全局最大（常常是 S 波）与 T 波窗内最大（判据只看后者）
    i_pk_glob = int(np.argmax(env_z))
    m_win = (t >= tt_ts - 30) & (t <= tt_ts + 60)
    if m_win.any():
        i_pk = int(np.flatnonzero(m_win)[np.argmax(env_z[m_win])])
        print(f"包络峰: 全局 {t[i_pk_glob]:+.1f} s (多半是 S 波) | "
              f"T 波窗内 {t[i_pk]:+.1f} s（与 T-S 差 {t[i_pk] - tt_ts:+.1f} s, 判据 < 60 s）")
    else:
        i_pk = i_pk_glob
        print(f"包络峰: 全局 {t[i_pk_glob]:+.1f} s（T 波窗内无数据）")

    def rms(x):
        return float(np.sqrt(np.mean(np.square(x)))) if x.size else np.nan

    # ★ SNR 主口径 = 项目口径（位移上的 Welch 2-8 Hz 能量比），与 stations.csv 的 snr 列同源
    snr, sig_win, noi_win = snr_welch_proj(disp["BHZ"].data, dt, t0, tt_p, tt_tp, tt_ts)
    snr_gap0, _, noi_win0 = snr_welch_proj(disp["BHZ"].data, dt, t0, tt_p, tt_tp, tt_ts, gap=0.0)
    snr_vel, _, _ = snr_welch_proj(z, dt, t0, tt_p, tt_tp, tt_ts)
    # 对照：同一台只换到 T-S 附近局部窗（旧写法），保留以便看两种口径差多少
    m_sig = (t >= tt_ts - 10) & (t <= tt_ts + 30)
    m_noi = (t >= tt_ts - 160) & (t <= tt_ts - 60)
    snr_loc = 20 * np.log10(rms(env_z[m_sig]) / rms(env_z[m_noi]))
    print(f"SNR(Welch 2-8 Hz, 位移) = {snr:.1f} dB   信号窗 {sig_win[0]:.0f}-{sig_win[1]:.0f} s, "
          f"噪声窗 {noi_win[0]:.0f}-{noi_win[1]:.0f} s (末距 P {SNR_GAP:g} s)   （< {SNR_MIN:g} dB 慎用该台）")
    print(f"  对照: 噪声窗贴 P(gap=0) = {snr_gap0:.1f} dB  [窗 {noi_win0[0]:.0f}-{noi_win0[1]:.0f} s] "
          f"| 速度波形 = {snr_vel:.1f} dB | T-S 局部窗包络比 = {snr_loc:.1f} dB")

    unit = "m/s" if OUTPUT == "VEL" else "m"

    # 图1 原始
    fig, ax = plt.subplots(3, 1, figsize=FIG_WAVE, sharex=True)   # 1:1
    for a_, c, lab in zip(ax, ("BHZ", "BHN", "BHE"), ("Z", "N", "E")):
        d = raw[c].data
        a_.plot(t[0] + np.arange(len(d)) * dt, d, "k", lw=0.3)
        a_.set_ylabel(f"{lab} (count)")
    ax[0].axvline(0, color="0.5", lw=0.8, ls="--")
    ax[0].set_title(f"{STA} raw SAC   t=0 origin, data start {t0:+.0f} s", fontsize=10)
    ax[-1].set_xlabel("Time relative to origin (s)")
    ax[-1].set_xlim(XLIM_T0, tt_ts + XLIM_TAIL)          # 横轴: 0 ~ T-S+500 s
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_STA, f"{STA}_1_raw.png"), dpi=DPI)
    plt.close(fig)

    # 图2 去响应
    fig, ax = plt.subplots(3, 1, figsize=FIG_WAVE, sharex=True)   # 1:1
    for a_, y, lab in zip(ax, (z, nn, ee), ("Z", "N", "E")):
        a_.plot(t, y, "k", lw=0.4)
        a_.set_ylabel(f"{lab} ({unit})")
        for j, (ph, tt, clr) in enumerate(ARRIVALS):
            a_.axvline(tt, color=clr, lw=0.7, ls="--", alpha=0.85,
                       label=f"{ph} {tt:.0f} s" if a_ is ax[0] else None)
    ax[0].legend(loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=4, frameon=False,
                 handlelength=1.0, handletextpad=0.35, columnspacing=1.0)
    ax[0].set_title(f"{STA}  after remove_response ({OUTPUT})", pad=20)
    ax[-1].set_xlabel("Time relative to origin (s)")
    ax[-1].set_xlim(XLIM_T0, tt_ts + XLIM_TAIL)          # 横轴: 0 ~ T-S+500 s
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_STA, f"{STA}_2_remove_response.png"), dpi=DPI)
    plt.close(fig)

    # 图3 谱图 ★ 三分量 Z/R/T（口径照 process_image.py L402-455）
    #   - 水平分量先按 baz 旋转到 R/T（obspy rotate_ne_rt，与案例脚本同一约定）
    #   - 整段算、不做前后截断；横轴 = t[0] + ts_（ts_ 相对输入段起点，不加就整体错位）
    #   - nperseg=12.8 s、noverlap 93.75%、nfft=2x、mode=psd、hann
    #   - 色标取三分量合并的 P2~P95，三块面板共用一条色标
    from obspy.signal.rotate import rotate_ne_rt
    r_rot, t_rot = rotate_ne_rt(nn, ee, baz)       # R = 径向, T = 切向
    comps = [("Z", z), ("R", r_rot), ("T", t_rot)]
    nperseg = int(fs * SPEC_WIN)
    spc, pool = {}, []
    for lab, d in comps:
        f_, ts_, S = spectrogram(d, fs=fs, window="hann", nperseg=nperseg,
                                 noverlap=int(nperseg * SPEC_OVERLAP),
                                 nfft=nperseg * SPEC_NFFT_X, scaling="spectrum", mode="psd")
        sdb = 10 * np.log10(np.maximum(S, 1e-30))
        fm = f_ <= fs / 2.0
        spc[lab] = (ts_, sdb, fm)
        pool.append(sdb[fm, :].ravel())
    # ★ 色标范围固定 P2~P95（三分量合并统计），否则自动 min-max 会被低频强噪声占满
    vmin, vmax = np.percentile(np.concatenate(pool), [SPEC_VMIN_PCT, SPEC_VMAX_PCT])
    if vmax <= vmin:
        vmax = vmin + 10.0

    fig, axs = plt.subplots(3, 1, figsize=FIG_SPEC, sharex=True)   # 1:1
    # 与其余四张图取同一图幅/同一左边距；右边只让出一条色标条的位置
    fig.subplots_adjust(hspace=0.08, top=0.895, bottom=0.105, left=0.105, right=0.905)
    for i, (lab, _) in enumerate(comps):
        ax = axs[i]
        ts_, sdb, fm = spc[lab]
        tspec = t[0] + ts_                             # ★ 绝对计时
        pc = ax.pcolormesh(tspec, f_[fm], sdb[fm, :], shading="auto", cmap="viridis",
                           vmin=vmin, vmax=vmax)
        ax.set_xlim(XLIM_T0, tt_ts + XLIM_TAIL)        # 与波形图同一横轴
        ax.set_ylim(0, fs / 2.0)
        ax.text(0.995, 0.88, lab, transform=ax.transAxes, ha="right", va="top",
                fontsize=11, fontweight="bold", color="white")
        for ph, tt, clr in ARRIVALS:               # 与其他图同一套配色（全项目统一）
            ln = ax.axvline(tt, color=clr, lw=0.8, ls="--", alpha=0.95,
                            label=f"{ph} {tt:.0f} s" if i == 0 else None)
            # 黑色描边：黑线打底(1.5) + 彩色线(0.8)压中间 -> 两侧各露约 0.35 pt 黑边
            # 注意: 描边宽度必须 > 彩色线宽，否则黑底被完全盖住、描边看不见
            ln.set_path_effects([pe.withStroke(linewidth=1.5, foreground="black")])
    axs[0].legend(loc="lower center", bbox_to_anchor=(0.5, 1.02), ncol=4, frameon=False,
                  handlelength=1.2, handletextpad=0.4, columnspacing=1.4)
    axs[-1].set_xlabel("Time relative to origin (s)")
    axs[1].set_ylabel("Frequency (Hz)")
    ca = fig.add_axes([0.917, 0.105, 0.012, 0.790])    # 单独色标轴（三面板共用）
    cb = fig.colorbar(pc, cax=ca)
    cb.set_label("Power (dB)")                         # 只写物理量，不加任何多余标注
    cb.ax.tick_params(labelsize=9)
    cb.outline.set_linewidth(0.6)
    fig.suptitle(f"{STA}  Z/R/T spectrogram 0-{fs/2:g} Hz (Nyquist {fs:g} sps)", y=0.968)
    fig.savefig(os.path.join(FIG_STA, f"{STA}_3_spectrogram.png"), dpi=DPI)
    plt.close(fig)

    # 图4 包络图 ★ 三分量 Z/R/T（口径照 process_image.py plot_packet L333-405）
    #   - 淡灰 = 2-8 Hz 带通波形（背景参考，lw 0.3）；红 = 包络（|Hilbert| 后 6 s 平滑）
    #   - 归一化：三分量合并的 P98（全局，不是逐台逐分量），纵轴对称 ±y_max
    rf = bandpass(r_rot, PRE_BAND[0], PRE_BAND[1], fs, corners=4, zerophase=True)
    tf = bandpass(t_rot, PRE_BAND[0], PRE_BAND[1], fs, corners=4, zerophase=True)
    env_r = smooth_env(rf, fs, ENV_SMOOTH)
    env_t = smooth_env(tf, fs, ENV_SMOOTH)
    eco = [("Z", zf, env_z), ("R", rf, env_r), ("T", tf, env_t)]
    ep98 = float(np.percentile(np.concatenate([np.abs(d) for _, d, _ in eco]), 98)) or 1.0
    # 只画正半边（不镜像），纵轴从 0 起到全局包络峰值的 1.02 倍 —— 同样的振幅占比，观感放大一倍
    y_top = 1.02 * max(float(np.max(e_ / ep98)) for _, _, e_ in eco)
    fig, axs = plt.subplots(3, 1, figsize=FIG_WAVE, sharex=True)     # 1:1
    for i, (lab, w_, e_) in enumerate(eco):
        ax = axs[i]
        ax.plot(t, e_ / ep98, "k", lw=0.8)                           # 黑色包络线（只要包络，不画波形）
        ax.set_ylim(0, y_top)                                        # 半边显示，y 从 0 起
        ax.set_ylabel("Envelope" if i == 1 else "")
        ax.text(0.99, 0.90, lab, transform=ax.transAxes, ha="right", va="top",
                fontsize=11, fontweight="bold")
        for ph, tt, clr in ARRIVALS:
            ax.axvline(tt, color=clr, lw=0.7, ls="--", alpha=0.85,
                       label=f"{ph} {tt:.0f} s" if i == 0 else None)
    axs[0].legend(loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=4, frameon=False,
                  handlelength=1.0, handletextpad=0.35, columnspacing=1.0)
    axs[0].set_title(f"{STA}  {PRE_BAND[0]:g}-{PRE_BAND[1]:g} Hz envelope", pad=20)
    axs[-1].set_xlabel("Time relative to origin (s)")
    axs[-1].set_xlim(XLIM_T0, tt_ts + XLIM_TAIL)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_STA, f"{STA}_4_envelope.png"), dpi=DPI)
    plt.close(fig)

    # 图5 带通
    fig, ax = plt.subplots(3, 1, figsize=FIG_WAVE, sharex=True)   # 1:1
    for a_, y, lab in zip(ax, (zf, rf, tf), ("Z", "R", "T")):     # 与其他图一致，用旋转后的 Z/R/T
        a_.plot(t, y, "k", lw=0.4)
        a_.text(0.99, 0.90, lab, transform=a_.transAxes, ha="right", va="top",
                fontsize=11, fontweight="bold")                    # 只标分量，不写频带/单位
        for j, (ph, tt, clr) in enumerate(ARRIVALS):
            a_.axvline(tt, color=clr, lw=0.7, ls="--", alpha=0.85,
                       label=f"{ph} {tt:.0f} s" if a_ is ax[0] else None)
    ax[0].legend(loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=4, frameon=False,
                 handlelength=1.0, handletextpad=0.35, columnspacing=1.0)
    ax[0].set_title(f"{STA}  bandpass {PRE_BAND[0]:g}-{PRE_BAND[1]:g} Hz", pad=20)
    ax[-1].set_xlabel("Time relative to origin (s)")
    ax[-1].set_xlim(XLIM_T0, tt_ts + XLIM_TAIL)          # 横轴: 0 ~ T-S+500 s
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_STA, f"{STA}_bandpass_2-8Hz.png"), dpi=DPI)
    plt.close(fig)

    # 图6 包络 + 窗（SNR QC 图；与图4 不同：这张只画 Z 分量并标出算 SNR 的两个窗）
    fig, ax = plt.subplots(figsize=FIG_WAVE)                     # 1:1
    ax.plot(t, env_z, "k", lw=0.9, label="Z envelope (2-8 Hz, 6 s smooth)")
    ax.plot(t[i_pk], env_z[i_pk], "o", ms=5, mfc="none", mec="#d62728", mew=1.2,
            label=f"envelope peak {t[i_pk]:.0f} s")
    # 两个窗 = SNR 计算实际用的窗（项目口径：信号 T-P-20 ~ T-S+100；噪声取等长的 P 之前段）
    ax.axvspan(sig_win[0], sig_win[1], color="#ff7f0e", alpha=0.25,
               label=f"signal window ({sig_win[1]-sig_win[0]:.0f} s)")
    ax.axvspan(noi_win[0], noi_win[1], color="#d62728", alpha=0.18,
               label=f"noise window ({noi_win[1]-noi_win[0]:.0f} s, ends P-{SNR_GAP:g} s)")
    for ph, tt, clr in ARRIVALS:
        ax.axvline(tt, color=clr, lw=1.0, ls="--", alpha=0.85, label=f"{ph} {tt:.0f} s")
    ax.set_xlim(XLIM_T0, tt_ts + XLIM_TAIL)              # 横轴: 0 ~ T-S+500 s
    ax.set_xlabel("Time relative to origin (s)")
    ax.set_ylabel(f"Envelope ({unit})")
    # 图例放轴外上方两行：轴内任何角落都会压到窗色块或 T-P/T-S 竖线
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.015), ncol=4, frameon=False,
              handlelength=1.2, handletextpad=0.4, columnspacing=1.4)
    ax.set_title(f"{STA}  envelope 6 s  |  SNR = {snr:.1f} dB (Welch 2-8 Hz)  |  peak {t[i_pk]:.0f} s",
                 pad=44)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_STA, f"{STA}_packet.png"), dpi=DPI)
    plt.close(fig)

    # 缓存给下游（F-K / 剖面 / 偏振）
    # ★ 同时存 Z/N/E 与旋转后的 Z/R/T：前者对应原始方位，后者是剖面图/衰减脚本的口径
    np.savez_compressed(os.path.join(ENV_OUT, f"{STA}_envelope_2-8Hz.npz"),
                        time=t, waveform_z=zf, waveform_n=nf, waveform_e=ef,
                        waveform_r=rf, waveform_t=tf,
                        envelope_z=env_z,
                        envelope_n=smooth_env(nf, fs, ENV_SMOOTH),
                        envelope_e=smooth_env(ef, fs, ENV_SMOOTH),
                        envelope_r=env_r, envelope_t=env_t,
                        dt=dt, t_off=t0, tt_tp=tt_tp, tt_ts=tt_ts,
                        baz=baz, land_km=land, ocean_km=ocean, output=OUTPUT,
                        band=(PRE_BAND[0], PRE_BAND[1]), smooth_seconds=ENV_SMOOTH)

    print("\n产物：")
    for f in sorted(glob.glob(os.path.join(FIG_STA, f"{STA}_*")) +
                    glob.glob(os.path.join(ENV_OUT, f"{STA}_*"))):
        print("  ", os.path.basename(f))
    print("步0 完成。")


if __name__ == "__main__":
    main()
