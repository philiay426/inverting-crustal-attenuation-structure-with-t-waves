#!/usr/bin/env python3
"""复刻 Zhou et al. (2021) Fig.3 —— 三分量记录剖面（3 列 × 2 行，6 面板）。

上排 (a)(b)(c)：纵轴 = 震中距 Epicentral distance (km)，横轴 = 时间 (s)，0 = 发震时刻
下排 (d)(e)(f)：纵轴 = 陆地段距离 Continental path distance (km)，横轴 = 时间 (s)，
                0 = T 波到达大陆坡（t − ocean_km/1.48）
红/蓝短线：T-P_local / T-S_local 理论到时（论文式(1) 恒速口径：L_ocean/1.48 + L_land/6.0 或 3.6）

波形：sacprocess/envelope/{STA}_envelope_2-8Hz.npz 的 waveform_z/r/t（2-8 Hz 零相位带通），
      逐台 P98 归一化（每台振幅尺度一致，便于按台比较波包形态）。
"""
import os
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- 路径：事件根/BASE 由环境变量传入；输出目录默认 = 本脚本所在目录 ----
import os as _os
BASE = _os.environ.get("T_BASE", "").strip()
if not BASE:
    raise SystemExit("请用环境变量 T_BASE 指定事件根目录")
CSV = _os.environ.get("T_CSV", f"{BASE}/process_flow/1_station_distribution/stations.csv")
ENV = _os.environ.get("T_ENV", f"{BASE}/sacprocess/envelope")
_MYOUT = _os.path.dirname(_os.path.abspath(__file__))          # 图直接写本步目录
OUT = _os.environ.get("T_OUT", _os.path.join(_MYOUT, "fig3_record_sections_3x2.png"))

V_OCEAN, V_P, V_S = 1.48, 6.0, 3.6     # km/s，论文式(1)
TT_MODE = _os.environ.get("TT_MODE", "proj")   # ★ 到时口径：'proj'=投影走时(stations.csv, 陆地段 TauP)
                                              #              'const'=论文式(1) 恒速 6.0/3.6
if TT_MODE == "proj":                   # 投影口径输出到 *_proj.png，恒速版保留原名
    OUT = OUT.replace(".png", "_proj.png")
XLIM_TOP = (0.0, 1600.0)               # 相对发震
XLIM_BOT = (-1000.0, 500.0)            # 相对 T 波上岸
AMP_KM_DIST = 50.0                     # ★ 上排（震中距轴）包络显示幅度 ±km —— 调大小改这里
AMP_KM_LAND = 25.0                      # ★ 下排（陆地区距轴）包络显示幅度 ±km
SEG_FRAC = 1.05                        # ★ 到时短竖线长度 = SEG_FRAC × 间距（调短改这里）
COMP_TITLE = {"Z": "Vertical", "R": "Radial", "T": "Tangential"}
PANEL_TAG = {"Z": "abc", "R": "abc", "T": "abc"}

plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 8, "axes.labelsize": 9, "axes.titlesize": 10,
                     "pdf.fonttype": 42, "ps.fonttype": 42})

# ---------------- 读台站 + 波形 ----------------
rows = [r for r in csv.DictReader(l for l in open(CSV) if not l.startswith("#"))]
rows = [r for r in rows if r.get("sta")]
data = {}
for r in rows:
    sta = r["sta"]
    f = f"{ENV}/{sta}_envelope_2-8Hz.npz"
    if not os.path.exists(f):
        continue
    z = np.load(f)
    data[sta] = dict(time=z["time"].astype(float),
                     waves={"Z": z["waveform_z"].astype(float),
                            "R": z["waveform_r"].astype(float),
                            "T": z["waveform_t"].astype(float)},
                     envs={"Z": z["envelope_z"].astype(float),
                           "R": z["envelope_r"].astype(float),
                           "T": z["envelope_t"].astype(float)},
                     dist=float(r["dist"]), land=float(r["land_km"]), ocean=float(r["ocean_km"]),
                     snr=float(r["snr"]))
    data[sta]["tocean"] = data[sta]["ocean"] / V_OCEAN
    if TT_MODE == "proj":
        # 投影走时：stations.csv 已由投影剖面算法给出（海洋段 1.48 + 陆地段 TauP(iasp91)）
        data[sta]["tt_tp"] = float(r["tt_tp"])
        data[sta]["tt_ts"] = float(r["tt_ts"])
    else:
        data[sta]["tt_tp"] = data[sta]["ocean"] / V_OCEAN + data[sta]["land"] / V_P   # 论文式(1)
        data[sta]["tt_ts"] = data[sta]["ocean"] / V_OCEAN + data[sta]["land"] / V_S
print(f"台站数: {len(data)}  (dist {min(d['dist'] for d in data.values()):.0f}–"
      f"{max(d['dist'] for d in data.values()):.0f} km, "
      f"land {min(d['land'] for d in data.values()):.0f}–{max(d['land'] for d in data.values()):.0f} km)")

by_dist = sorted(data, key=lambda s: data[s]["dist"])
by_land = sorted(data, key=lambda s: data[s]["land"])

def spacing(keys, field):
    v = np.array([data[k][field] for k in keys])
    return float(np.median(np.diff(np.sort(v))))

SP_DIST = spacing(by_dist, "dist")
SP_LAND = spacing(by_land, "land")
print(f"中位间距: dist {SP_DIST:.1f} km | land {SP_LAND:.1f} km")

def draw_row(ax, comp, keys, field, xlim, align_land, seg_len, show_ylabel, amp):
    """在 ax 上画一个分量的一排包络剖面（6 s 平滑，非负）。align_land=True 时横轴减上岸走时。"""
    for sta in keys:
        d = data[sta]
        x = d["time"] - (d["tocean"] if align_land else 0.0)
        w = d["envs"][comp]                       # 6 s 平滑包络（= 波包轮廓）
        m = (x >= xlim[0]) & (x <= xlim[1])
        if not np.any(m):
            continue
        y0 = d[field]
        # 逐台峰值归一化：每台包络等比放大到 amp（km），允许相邻台重叠
        norm = np.max(w[m]) or 1.0
        ax.plot(x[m], y0 + w[m] / norm * amp, color="black", lw=0.35, alpha=0.60,   # 变淡
                rasterized=True, zorder=3)
        # 理论到时：短竖线（红 = T-P_local，蓝 = T-S_local）
        if align_land:
            # 下排横轴已按上岸时刻对齐 → 理论线 = 纯陆地段走时（随口径自动正确）
            x_tp, x_ts = d["tt_tp"] - d["tocean"], d["tt_ts"] - d["tocean"]
        else:
            x_tp, x_ts = d["tt_tp"], d["tt_ts"]
        for xx, c in ((x_tp, "red"), (x_ts, "blue")):
            if xlim[0] <= xx <= xlim[1]:
                ax.plot([xx, xx], [y0 - seg_len / 2, y0 + seg_len / 2],
                        color=c, lw=0.8, alpha=0.45, zorder=6)   # 变短、变淡
    ax.set_xlim(*xlim)
    ylo = min(data[k][field] for k in keys)
    yhi = max(data[k][field] for k in keys)
    pad = max(0.6 * (SP_LAND if align_land else SP_DIST), 1.05 * amp)   # 保证放大后不被裁
    ax.set_ylim(ylo - pad, yhi + pad)
    ax.tick_params(direction="in", labelsize=7.5)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    if show_ylabel:
        ax.set_ylabel("Continental path distance (km)" if align_land else "Epicentral distance (km)",
                      fontsize=9, fontweight="bold")

fig, axs = plt.subplots(2, 3, figsize=(10.2, 9.0))
fig.subplots_adjust(left=0.075, right=0.985, top=0.925, bottom=0.065, wspace=0.22, hspace=0.20)

for col, comp in enumerate(["Z", "R", "T"]):
    # 上排：震中距 / 相对发震
    ax = axs[0, col]
    draw_row(ax, comp, by_dist, "dist", XLIM_TOP, False, SEG_FRAC * SP_DIST, col == 0, AMP_KM_DIST)
    ax.set_title(f"({chr(97+col)})", loc="left", fontsize=9)
    ax.text(0.5, 1.115, COMP_TITLE[comp], transform=ax.transAxes, ha="center", va="bottom",
            fontsize=10, fontweight="bold")
    ax.text(0.40, 1.005, "T-P$_{local}$", transform=ax.transAxes, color="red",
            fontsize=8.5, ha="left", va="bottom")
    ax.text(0.72, 1.005, "T-S$_{local}$", transform=ax.transAxes, color="blue",
            fontsize=8.5, ha="left", va="bottom")
    ax.set_xlabel("Time (s)", fontsize=8.5)
    # 下排：陆地段距离 / 相对上岸
    ax = axs[1, col]
    draw_row(ax, comp, by_land, "land", XLIM_BOT, True, SEG_FRAC * SP_LAND, col == 0, AMP_KM_LAND)
    ax.set_title(f"({chr(100+col)})", loc="left", fontsize=9)
    ax.text(0.40, 1.005, "T-P$_{local}$", transform=ax.transAxes, color="red",
            fontsize=8.5, ha="left", va="bottom")
    ax.text(0.72, 1.005, "T-S$_{local}$", transform=ax.transAxes, color="blue",
            fontsize=8.5, ha="left", va="bottom")
    ax.set_xlabel("Time (s)", fontsize=8.5)

os.makedirs(os.path.dirname(OUT), exist_ok=True)   # 输出目录可能还没建
fig.savefig(OUT, dpi=300, bbox_inches="tight")
plt.close(fig)
print("已保存:", OUT)
