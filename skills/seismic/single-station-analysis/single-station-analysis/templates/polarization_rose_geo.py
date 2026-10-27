#!/usr/bin/env python3
"""
步 2 定版图件 A：单台偏振玫瑰（地理 N/E 方位）+ Z/R/T 包络时窗图

版式
  上排：左 T-P 玫瑰、右 T-S 玫瑰；完整整圆，0=N、顺时针，方位刻度 N/E/S/W
        柱：10 度分箱，柱宽 = 分箱宽（相邻柱条无缝），环形 3 点平滑后再画
        黑线：理论偏振方向，**只画半段**（圆心到外圈），且画在偏振样本较多的那一侧
              T-P 沿径向 baz，T-S 沿切向 baz+90
        副标题：T-P / T-S + radial / transverse 扇区占比（30 度半宽，权重 lambda1）
  下排：Z / R / T 三条黑色包络（2-8 Hz、6 s 平滑、三条共用 P98 归一、垂直偏移堆叠）
        T-P / T-S 时窗阴影 + 同色中心虚线；图例右上两项；x 范围 T-P-220 -> T-S+260

口径
  方向样本 = 2-8 Hz 带通 -> BAZ 旋转 ZRT -> PCA 滑窗（窗 1 s / 步 0.25 s / analytic）主轴
  权重 lambda1；地理方位 = 主轴 R-T 角 + BAZ；时窗 = 理论到时 +/- HALF
  PCA 主轴只有轴向含义（v 与 -v 等价），180 度模糊不可避免

环境变量
  SRC   必填  去响应后 SAC 所在目录（文件名含 <NET>.<STA>.*.BH[ZNE]）
  CSV   必填  stations.csv（需列 sta / tt_tp / tt_ts / baz / land_km / dist）
  STA   必填  台站名
  OUT   可选  输出目录（默认当前目录）
  NET   可选  台网代码（默认 XA）
  EVENT_TIME 可选  发震时刻 ISO 串，仅在 SAC 头里 o 无效时用于兜底时间轴
  FMIN/FMAX/HALF/PWIN/PSTEP 可选  默认 2/8/25/1.0/0.25
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
from scipy.signal import hilbert
from obspy import read, Stream, UTCDateTime
from obspy.signal.filter import bandpass
from obspy.signal.rotate import rotate_ne_rt
from obspy.signal.util import smooth

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
HALF = float(os.environ.get("HALF", 25.0))
PWIN = float(os.environ.get("PWIN", 1.0))
PSTEP = float(os.environ.get("PSTEP", 0.25))
SMOOTH_S = 6.0
DPI = 300

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9.5, "axes.labelsize": 10.5, "axes.titlesize": 11,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9,
    "axes.linewidth": 0.7, "xtick.direction": "in", "ytick.direction": "in",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
C_ROSE = "#5B9BD5"
C_TP, C_TS = "#5B9BD5", "#9370DB"


def sta_row(sta):
    rows = list(csv.DictReader((l for l in open(CSV) if not l.startswith("#"))))
    for r in rows:
        if r["sta"].strip() == sta:
            return r
    raise SystemExit(f"{CSV} 里没有 {sta}")


def load_3c(sta, net, src):
    trs = {}
    for ch in ("BHZ", "BHN", "BHE"):
        m = glob.glob(f"{src}/{net}.{sta}.*.{ch}*_resp.sac")
        if not m:
            m = glob.glob(f"{src}/*{sta}*.{ch}*_resp.sac")
        if not m:
            raise SystemExit(f"缺 {sta} {ch}（SRC={src}）")
        trs[ch] = read(m[0])[0]
    z = trs["BHZ"]
    sr = 1.0 / z.stats.delta
    sac = z.stats.sac
    o_val = getattr(sac, "o", None)
    # 时间轴 b-o 口径：o 有效就用 b-o，否则用 starttime 与发震时刻之差兜底
    if o_val is not None and 0 < abs(o_val) < 1e6:
        t_off = float(sac.b) - float(o_val)
    elif EV is not None:
        t_off = float(z.stats.starttime - EV) + float(sac.b)
    else:
        raise SystemExit("SAC 头里 o 无效且未给 EVENT_TIME，无法确定时间轴")
    t = np.arange(z.stats.npts) / sr + t_off
    return trs, z, sr, t, t_off


def main():
    row = sta_row(STA)
    ttp, tts, baz = float(row["tt_tp"]), float(row["tt_ts"]), float(row["baz"])
    land_km, epi_km = float(row["land_km"]), float(row["dist"])

    trs, z, sr, t, t_off = load_3c(STA, NET, SRC)
    ch_bp = {ch: bandpass(trs[f"BH{ch}"].data.astype(float), FMIN, FMAX, df=sr, corners=4)
             for ch in ("Z", "N", "E")}
    r_bp, t_bp = rotate_ne_rt(ch_bp["N"], ch_bp["E"], baz)
    comps = {"Z": ch_bp["Z"], "R": r_bp, "T": t_bp}

    stb = Stream()
    for ch in ("Z", "N", "E"):
        trc = trs[f"BH{ch}"].copy()
        trc.data = ch_bp[ch]
        stb += trc
    res = pca(stb, window_length_s=PWIN, window_step_s=PSTEP,
              back_azimuth=baz, covariance_method="analytic")
    pv = np.asarray(res["principal_vector"], dtype=float)
    lam1 = np.asarray(res["eigenvalues"], dtype=float)[:, 0]
    tw = np.asarray(res["time_windows_s"], dtype=float) + t_off
    ang_rt = np.degrees(np.arctan2(pv[:, 2], pv[:, 1])) % 360.0     # 0 = 径向
    ang_geo = (ang_rt + baz) % 360.0
    ok = np.isfinite(ang_geo) & np.isfinite(lam1)

    def sector(angs, w, half=30.0):
        d = np.abs(((angs - baz + 90.0) % 180.0) - 90.0)            # 到径向轴的最小角
        return (float(np.nansum(w[d <= half]) / np.nansum(w)),
                float(np.nansum(w[d >= 90.0 - half]) / np.nansum(w)))

    fig = plt.figure(figsize=(10.0, 8.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.98, 1.0], hspace=0.26, wspace=0.32,
                          left=0.085, right=0.945, top=0.905, bottom=0.075)

    def rose(gaxi, center, kind):
        m = ok & (np.abs(tw - center) <= HALF)
        angs, w = ang_geo[m], lam1[m]
        frac_r, frac_t = sector(angs, w)
        frac = frac_r if kind == "T-P" else frac_t
        a2 = np.radians(angs) * 2.0
        th0 = (np.degrees(np.arctan2(np.nansum(w * np.sin(a2)),
                                     np.nansum(w * np.cos(a2)))) / 2.0) % 180.0

        ax = fig.add_subplot(gaxi, polar=True)
        ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
        edges = np.arange(0.0, 360.1, 10.0)
        h, _ = np.histogram(angs % 360.0, bins=edges, weights=w)
        h = (np.roll(h, 1) + h + np.roll(h, -1)) / 3.0              # 环形 3 点平滑，消除空档
        h = h / h.max() * 60.0 if h.max() > 0 else h
        ax.bar(np.radians((edges[:-1] + edges[1:]) / 2.0), h, width=np.radians(10.0),
               color=C_ROSE, alpha=0.80, edgecolor=C_ROSE, lw=0.0, bottom=0)

        # 理论偏振方向：取样本多的那一侧，只画半段
        dev = np.abs(((angs - th0 + 180.0) % 360.0) - 180.0)
        w1 = np.nansum(w[dev <= 90.0]); w2 = np.nansum(w[dev > 90.0])
        prefer = th0 if w1 >= w2 else (th0 + 180.0) % 360.0
        d_th = baz if kind == "T-P" else baz + 90.0
        cand = [d_th % 360.0, (d_th + 180.0) % 360.0]
        d_th = min(cand, key=lambda x: abs(((x - prefer + 180.0) % 360.0) - 180.0))
        ax.plot([np.radians(d_th)] * 2, [0, 66], color="black", lw=1.6, zorder=6)

        ax.set_ylim(0, 66)
        ax.set_yticks([20, 40, 60])
        ax.tick_params(axis="y", labelsize=8.5)
        ax.set_rlabel_position(35)
        ax.set_xticks(np.radians([0, 90, 180, 270]))
        ax.set_xticklabels(["N", "E", "S", "W"], fontsize=10)
        ax.grid(color="0.80", lw=0.6)
        ax.spines["polar"].set_color("black"); ax.spines["polar"].set_linewidth(1.0)
        sub = "radial" if kind == "T-P" else "transverse"
        ax.set_title(f"{kind}\n{sub} {frac * 100:.0f}%", fontsize=10, pad=12)
        print(f"[rose] {kind}: {sub} {frac * 100:.0f}%  理论线 {d_th:.1f} deg  样本 {m.sum()}")

    rose(gs[0, 0], ttp, "T-P")
    rose(gs[0, 1], tts, "T-S")

    ax = fig.add_subplot(gs[1, :])
    x0, x1 = ttp - 220.0, tts + 260.0
    mm = (t >= x0) & (t <= x1)
    gs_ = np.percentile(np.abs(np.concatenate([comps[c][mm] for c in ("Z", "R", "T")])), 98) or 1.0
    for ch, off in (("Z", 0.78), ("R", 0.50), ("T", 0.22)):
        env = smooth(np.abs(hilbert(comps[ch])), max(3, int(SMOOTH_S * sr))) / gs_
        ax.plot(t[mm], env[mm] * 0.26 + off, color="black", lw=0.8)
        ax.text(x0 - 3.0, off + 0.01, ch, ha="right", va="center",
                fontsize=10.5, fontweight="bold")
    ax.set_yticks([])
    ax.set_ylim(0.10, 1.10)
    for cen, col, lab in ((ttp, C_TP, "T-P window"), (tts, C_TS, "T-S window")):
        ax.axvspan(cen - HALF, cen + HALF, color=col, alpha=0.25, lw=0)
        ax.axvline(cen, color=col, ls="--", lw=1.0)
    ax.set_xlim(x0, x1)
    ax.set_xlabel("Time from origin (s)")
    ax.legend(handles=[Patch(facecolor=C_TP, alpha=0.20, label="T-P window"),
                       Patch(facecolor=C_TS, alpha=0.20, label="T-S window")],
              loc="upper right", frameon=True, framealpha=0.9, edgecolor="0.7", fontsize=8.5)

    fig.suptitle(f"{STA}   land {land_km:.0f} km   epi {epi_km:.0f} km",
                 fontsize=11, fontweight="bold", y=0.985)
    os.makedirs(OUT, exist_ok=True)
    fn = f"{OUT}/{STA}_rose_geo_TP_TS.png"
    fig.savefig(fn, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print("=>", fn)


if __name__ == "__main__":
    main()
