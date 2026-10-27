#!/usr/bin/env python3
"""
步 2 定版图件 B：全时段偏振统计（PCA 滑窗，竖排三栏）

  第一栏：径向占比 R/(R+T)，红，0-1
  第二栏：切向占比 T/(R+T)，蓝，0-1     （与第一栏互补，两栏相加恒为 1）
  第三栏：线性度（蓝实线）+ 平面度（红虚线），共用 0-1 纵轴

  三栏竖排、共享时间轴，只有最下面一栏写 x 轴标题；
  时窗阴影 T-P 红 / T-S 蓝 + 同色中心虚线，三栏严格对齐；
  图例放在最上面一栏之外（横排四项），不压曲线；
  平滑：曲线先做滑动均值（默认 9 点 x 步长），线性度/平面度先 3 点中值再均值。

口径与玫瑰图一致：2-8 Hz 带通 -> BAZ 旋转 -> PCA 滑窗（窗 1 s / 步 0.25 s / analytic）
输出：<STA>_polarization_fullscan.png

环境变量
  SRC   必填  去响应后 SAC 目录
  CSV   必填  stations.csv（sta / tt_tp / tt_ts / baz / land_km / dist）
  STA   必填  台站名
  OUT   可选  输出目录（默认当前目录）
  NET   可选  台网代码（默认 XA）
  EVENT_TIME 可选  发震时刻 ISO（SAC 头 o 无效时兜底时间轴）
  FMIN/FMAX/HALF/PWIN/PSTEP/PAD_PRE/PAD_POST/SMOOTH_N 可选
"""
import csv
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from obspy import read, Stream, UTCDateTime
from obspy.signal.filter import bandpass
from obspy.signal.rotate import rotate_ne_rt

for _p in ("/home/sy/codex", "/home/sy/codex/t_phase_analysis"):
    if os.path.isdir(_p):
        sys.path.insert(0, _p)
from t_phase_analysis.pca import pca          # noqa: E402

SRC = os.environ.get("SRC")
CSV = os.environ.get("CSV")
STA = os.environ.get("STA")
if not (SRC and CSV and STA):
    raise SystemExit("必填环境变量：SRC=<去响应SAC目录> CSV=<stations.csv> STA=<台站名>")
if not os.path.isdir(SRC):
    raise SystemExit(f"SRC 不存在：{SRC}")
if not os.path.isfile(CSV):
    raise SystemExit(f"CSV 不存在：{CSV}")
OUT = os.environ.get("OUT", os.getcwd())
NET = os.environ.get("NET", "XA")
EVENT_TIME = os.environ.get("EVENT_TIME")
EV = UTCDateTime(EVENT_TIME) if EVENT_TIME else None

FMIN = float(os.environ.get("FMIN", 2.0))
FMAX = float(os.environ.get("FMAX", 8.0))
PWIN = float(os.environ.get("PWIN", 1.0))
PSTEP = float(os.environ.get("PSTEP", 0.25))
HALF = float(os.environ.get("HALF", 25.0))
PAD_PRE = float(os.environ.get("PAD_PRE", 220.0))
PAD_POST = float(os.environ.get("PAD_POST", 260.0))
SMOOTH_N = int(os.environ.get("SMOOTH_N", 9))
DPI = 300
C_R, C_T = "#c62828", "#1565c0"
C_LIN, C_PLA = "#2b6cb0", "#c0392b"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9.5, "axes.labelsize": 10.5,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 8.5,
    "axes.linewidth": 0.7, "xtick.direction": "in", "ytick.direction": "in",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def sta_row(sta):
    rows = list(csv.DictReader((l for l in open(CSV) if not l.startswith("#"))))
    for r in rows:
        if r["sta"].strip() == sta:
            return r
    raise SystemExit(f"{CSV} 里没有 {sta}")


def main():
    row = sta_row(STA)
    ttp, tts, baz = float(row["tt_tp"]), float(row["tt_ts"]), float(row["baz"])
    land_km, epi_km = float(row["land_km"]), float(row["dist"])

    trs = {}
    for ch in ("BHZ", "BHN", "BHE"):
        m = glob.glob(f"{SRC}/{NET}.{STA}.*.{ch}*_resp.sac") or \
            glob.glob(f"{SRC}/*{STA}*.{ch}*_resp.sac")
        if not m:
            raise SystemExit(f"缺 {STA} {ch}（SRC={SRC}）")
        trs[ch] = read(m[0])[0]
    z = trs["BHZ"]
    sr = 1.0 / z.stats.delta
    sac = z.stats.sac
    o_val = getattr(sac, "o", None)
    if o_val is not None and 0 < abs(o_val) < 1e6:      # 时间轴 b-o 口径
        t_off = float(sac.b) - float(o_val)
    elif EV is not None:
        t_off = float(z.stats.starttime - EV) + float(sac.b)
    else:
        raise SystemExit("SAC 头里 o 无效且未给 EVENT_TIME，无法确定时间轴")
    t = np.arange(z.stats.npts) / sr + t_off

    bp = {ch: bandpass(trs[f"BH{ch}"].data.astype(float), FMIN, FMAX, df=sr, corners=4)
          for ch in ("Z", "N", "E")}
    r_bp, t_bp = rotate_ne_rt(bp["N"], bp["E"], baz)

    stb = Stream()
    for ch in ("Z", "N", "E"):
        trc = trs[f"BH{ch}"].copy()
        trc.data = bp[ch]
        stb += trc
    res = pca(stb, window_length_s=PWIN, window_step_s=PSTEP,
              back_azimuth=baz, covariance_method="analytic")
    tw = np.asarray(res["time_windows_s"], dtype=float) + t_off
    lin = np.asarray(res["linearity"], dtype=float)
    pla = np.asarray(res["planarity"], dtype=float)

    # 逐窗 R/T 能量占比（同一窗长做滑动平均，R+T=1）
    n = max(3, int(round(2 * PWIN * sr)) | 1)
    ex = np.convolve(r_bp ** 2, np.ones(n) / n, mode="same")
    ey = np.convolve(t_bp ** 2, np.ones(n) / n, mode="same")
    with np.errstate(invalid="ignore", divide="ignore"):
        frac_t_trace = ey / (ex + ey)
    frac_t = np.interp(tw, t, frac_t_trace)
    frac_r = 1.0 - frac_t
    ok = np.isfinite(frac_r) & np.isfinite(lin)

    x0, x1 = ttp - PAD_PRE, tts + PAD_POST
    m = ok & (tw >= x0) & (tw <= x1)

    fig, axs = plt.subplots(3, 1, figsize=(9.2, 7.6), sharex=True,
                            gridspec_kw={"hspace": 0.10, "left": 0.105, "right": 0.985,
                                         "top": 0.915, "bottom": 0.085})

    def smooth(x, k=SMOOTH_N):
        pad = k // 2
        pad_ = np.r_[np.full(pad, x[0]), x, np.full(pad, x[-1])]
        return np.convolve(pad_, np.ones(k) / k, mode="valid")

    def med3(x, k=3):
        pad_ = np.r_[np.full(k // 2, x[0]), x, np.full(k // 2, x[-1])]
        return np.nanmedian(np.vstack([pad_[i:i + k] for i in range(len(x))]), axis=1)

    ax_r, ax_t, ax_a = axs
    ax_r.plot(tw[m], smooth(frac_r[m]), color=C_R, lw=1.1)
    ax_r.set_ylim(0, 1.02); ax_r.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax_r.set_ylabel("Radial  R/(R+T)", fontsize=10.5, fontweight="bold")
    ax_r.text(0.04, 0.93, "R", color=C_R, transform=ax_r.transAxes,
              fontsize=11, fontweight="bold")

    ax_t.plot(tw[m], smooth(frac_t[m]), color=C_T, lw=1.1)
    ax_t.set_ylim(0, 1.02); ax_t.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax_t.set_ylabel("Transverse  T/(R+T)", fontsize=10.5, fontweight="bold")
    ax_t.text(0.04, 0.93, "T", color=C_T, transform=ax_t.transAxes,
              fontsize=11, fontweight="bold")

    ax_a.plot(tw[m], smooth(med3(lin[m])), color=C_LIN, lw=1.1)
    ax_a.plot(tw[m], smooth(med3(pla[m])), color=C_PLA, lw=1.1, ls="--")
    ax_a.set_ylim(0, 1.02); ax_a.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax_a.set_ylabel("Linearity / Planarity", fontsize=10.5, fontweight="bold")

    for a_ in axs:
        for cen, col in ((ttp, C_R), (tts, C_T)):
            a_.axvspan(cen - HALF, cen + HALF, color=col, alpha=0.18, lw=0, zorder=0)
            a_.axvline(cen, color=col, ls="--", lw=0.8, zorder=1)
        a_.set_xlim(x0, x1)
    ax_a.set_xlabel("Time from origin (s)", fontsize=10.5)

    ax_r.legend(handles=[Line2D([], [], color=C_LIN, lw=1.1, label="Linearity"),
                         Line2D([], [], color=C_PLA, lw=1.1, ls="--", label="Planarity"),
                         Patch(facecolor=C_R, alpha=0.18, label="T-P window"),
                         Patch(facecolor=C_T, alpha=0.18, label="T-S window")],
                loc="lower left", bbox_to_anchor=(0.0, 1.005), ncol=4,
                frameon=False, fontsize=8.5)

    fig.suptitle(f"{STA}   land {land_km:.0f} km   epi {epi_km:.0f} km",
                 fontsize=11, fontweight="bold", y=0.975)
    os.makedirs(OUT, exist_ok=True)
    fn = f"{OUT}/{STA}_polarization_fullscan.png"
    fig.savefig(fn, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("=>", fn)
    for cen, lab in ((ttp, "T-P"), (tts, "T-S")):
        mm = ok & (np.abs(tw - cen) <= HALF)
        print(f"{lab} 窗 +/-{HALF:.0f} s: R {np.nanmean(frac_r[mm]):.3f} / T {np.nanmean(frac_t[mm]):.3f} "
              f"| bias {np.nanmean(2 * frac_r[mm] - 1):+.3f} "
              f"| 线性度 {np.nanmedian(lin[mm]):.3f} 平面度 {np.nanmedian(pla[mm]):.3f}")


if __name__ == "__main__":
    main()
