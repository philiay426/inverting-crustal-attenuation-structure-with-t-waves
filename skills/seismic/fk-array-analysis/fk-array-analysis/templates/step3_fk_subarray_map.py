#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step3_fk_subarray_map.py —— F-K 两个子阵（SUB1 / SUB2）的台站分布图（单面板，无放大局部）

产物: <T_OUTDIR>/<T_OUTNAME>.png   600 dpi
     · 区域地图（灰白大陆 + -1000 m 等深线虚线）
     · 其余台站灰三角（背景）；SUB1 = 红三角；SUB2 = 蓝三角
     · 沙滩球 = GCMT 042897C（279/88/1）
     · 红色/蓝色大圆虚线 = 震中 → 各子阵质心（真实大圆，geographiclib）
单面板，无放大局部

运行: T_BASE=<事件目录> SUB1=A,B,C [SUB2=D,E] EVENT_LAT=.. EVENT_LON=.. python step3_fk_subarray_map.py
"""
import os
import csv
import sys
import numpy as np
from geographiclib.geodesic import Geodesic

BASE = os.environ.get("T_BASE", "").strip()
if not BASE:
    raise SystemExit("请用环境变量 T_BASE 指定事件根目录")
FK = f"{BASE}/process_flow/4_fk_analysis"
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, f"{FK}/synth_forward")
try:                                     # 出版样式（缺失则用 pygmt 默认）
    from _pubstyle import apply_pub_style
    apply_pub_style()
except Exception as e:
    print(f"⚠ 未加载 _pubstyle（{e}），使用 pygmt 默认样式")

import pygmt  # noqa: E402
from obspy.geodetics import gps2dist_azimuth  # noqa: E402

# ============================== ★ 参数区 ==============================
OUT_NAME = os.environ.get("T_OUTNAME", "subarray_map")
OUT_DIR = os.environ.get("T_OUTDIR") or f"{FK}/figures"
OUT_DPI = 600
SAVE_PNG, SAVE_PDF = True, False

CSV = f"{BASE}/process_flow/1_station_distribution/stations.csv"
# ★ 两个子阵的台站清单按事件填（SUB1=A,B,C / SUB2=D,E 传入），组名用 NAME1 / NAME2 覆盖
G1 = [x for x in os.environ.get("SUB1", "").split(",") if x.strip()]
G2 = [x for x in os.environ.get("SUB2", "").split(",") if x.strip()]
if not G1 and not G2:
    raise SystemExit("请用 SUB1=A,B,C / SUB2=D,E 指定子阵台站")
NAME1 = os.environ.get("NAME1", "subarray 1")   # 图例里的组名
NAME2 = os.environ.get("NAME2", "subarray 2")

# ★ 事件参数（环境变量；默认值只是占位）
EVLA = float(os.environ.get("EVENT_LAT", "0.0"))
EVLO = float(os.environ.get("EVENT_LON", "0.0"))
EV_DEPTH_KM = float(os.environ.get("EVENT_DEP", "10.0"))
EVENT_LABEL = os.environ.get("EVENT_LABEL", "")   # 图上震中标两行文字的第二行
MECA_NP = (float(os.environ.get("CMT_STRIKE", "0")),
           float(os.environ.get("CMT_DIP", "90")),
           float(os.environ.get("CMT_RAKE", "0")))
MECA_MW = float(os.environ.get("EVENT_MW", "0.0"))
MECA_SCALE, MECA_CONVENTION = "0.6c", "aki"

# 地图范围：REGION=lon_min,lon_max,lat_min,lat_max（按事件裁剪）
_REG = [float(v) for v in os.environ.get("REGION", "0,360,-60,60").split(",")]
PROJ, REG = "M15c", _REG
COAST_RES = "h"
LAND_FILL, WATER_FILL = "gray92", "white"
SHORE_PEN = "0.35p,gray25"
BORDER_PEN = "1/0.25p,gray60"
FRAME = ["xa10f5+lLongitude", "ya5f2+lLatitude", "WSne"]

C_G1, C_G2, C_OTH = "red3", "blue3", "gray55"
T_G1, T_G2, T_OTH = 0.34, 0.34, 0.22
PEN_G1, PEN_G2, PEN_OTH = "0.3p,black", "0.3p,black", "0.2p,gray30"

DRAW_ISOBATH = True
ISO_LEVEL = [-1000]
ISO_PEN = "0.7p,black,--"
ISO_GRID = "@earth_relief_01m"

LINE_NPTS = 300
DRAW_VEC = True                          # 震中→两个质心的红色/蓝色大圆虚线

LEG_SPEC = (
    "G 0.12c\n"
    f"S 0.3c t {T_G1}c {C_G1} {PEN_G1} 0.7c {NAME1} ({len(G1)} stations)\n"
    "G 0.15c\n"
    f"S 0.3c t {T_G2}c {C_G2} {PEN_G2} 0.7c {NAME2} ({len(G2)} stations)\n"
    "G 0.15c\n"
    f"S 0.3c t {T_OTH}c {C_OTH} {PEN_OTH} 0.7c other stations\n"
    "G 0.15c\n"
    f"S 0.3c - 0.6c - 0.8p,black,-- 0.7c -1000 m depth\n"
)
# ⚠ 图例框宽要按最长条目的文字宽度定，否则文字会溢出框外：
#   9p 时最长条目（组名 + 台数）实测约 4.8 cm
#   -> 框宽 >= 0.7(dx2) + 4.82 + 0.3(边距) = 5.8 cm，取 6.0c
# ⚠ 字号不要写进 S 记录里（GMT 会报 "9p,Helvetica,blac not a valid number"），
#   改用 pygmt.config 临时把图例文字字体设成 9p（只影响 legend 这一次调用）
LEG_FONT = "9p,Helvetica,black"
LEG_BOX_W = "6.0c"
LEG_POS = f"jBL+w{LEG_BOX_W}+o0.5c/0.5c"

F_FONT = "9p,Helvetica-Bold,black"
BOX = dict(fill="white@20", pen="0.3p,gray45")
# =====================================================================

rows = {}
with open(CSV) as f:
    for r in csv.DictReader(l for l in f if not l.startswith("#")):
        rows[r["sta"]] = r
all_sta = list(rows)
missing = [s for s in G1 + G2 if s not in rows]
if missing:
    raise SystemExit(f"stations.csv 里缺: {missing}")


def centroid(stas):
    la = np.mean([float(rows[s]["lat"]) for s in stas])
    lo = np.mean([float(rows[s]["lon"]) for s in stas])
    return la, lo


def dmin_dmax(stas):
    lon = np.array([float(rows[s]["lon"]) for s in stas])
    lat = np.array([float(rows[s]["lat"]) for s in stas])
    x = (lon - lon.mean()) * 111.195 * np.cos(np.radians(lat.mean()))
    y = (lat - lat.mean()) * 111.195
    D = np.hypot(x[:, None] - x[None, :], y[:, None] - y[None, :])
    iu = np.triu_indices(len(stas), 1)
    return D[iu].min(), D[iu].max()


info = {}
for name, grp in (("G1", G1), ("G2", G2)):
    la, lo = centroid(grp)
    dmn, dmx = dmin_dmax(grp)
    baz = float(np.mean([gps2dist_azimuth(float(rows[s]["lat"]), float(rows[s]["lon"]), EVLA, EVLO)[1]
                         for s in grp]))
    info[name] = dict(lat=la, lon=lo, dmin=dmn, dmax=dmx, baz=baz)
    print(f"{name}: {len(grp)} 台, 质心 ({la:.4f}, {lo:.4f}), d_min {dmn:.1f} km, 孔径 {dmx:.1f} km, "
          f"真实大圆 BAZ 均值 {baz:.2f}°, C=111.195/(2 d_min)={111.195/(2*dmn):.4f}")

fig = pygmt.Figure()
fig.coast(region=REG, projection=PROJ, resolution=COAST_RES,
          land=LAND_FILL, water=WATER_FILL, shorelines=SHORE_PEN,
          borders=[BORDER_PEN], frame=FRAME)
if DRAW_ISOBATH:
    fig.grdcontour(grid=ISO_GRID, region=REG, projection=PROJ, levels=ISO_LEVEL, pen=ISO_PEN)

oth = [s for s in all_sta if s not in G1 and s not in G2]
fig.plot(x=[float(rows[s]["lon"]) for s in oth], y=[float(rows[s]["lat"]) for s in oth],
         style=f"t{T_OTH}c", fill=C_OTH, pen=PEN_OTH)
for name, grp, col in (("G1", G1, C_G1), ("G2", G2, C_G2)):
    fig.plot(x=[float(rows[s]["lon"]) for s in grp], y=[float(rows[s]["lat"]) for s in grp],
             style=f"t{T_G1 if name == 'G1' else T_G2}c", fill=col,
             pen=PEN_G1 if name == "G1" else PEN_G2)

fig.meca(spec=np.array([[EVLO, EVLA, EV_DEPTH_KM, MECA_NP[0], MECA_NP[1], MECA_NP[2], MECA_MW]]),
         scale=MECA_SCALE, convention=MECA_CONVENTION, pen="0.4p,black")

if DRAW_VEC:
    for name, col in (("G1", C_G1), ("G2", C_G2)):
        inv = Geodesic.WGS84.InverseLine(EVLA, EVLO, info[name]["lat"], info[name]["lon"])
        pts = [inv.Position(d) for d in np.linspace(0, inv.s13, LINE_NPTS)]
        fig.plot(x=[p["lon2"] for p in pts], y=[p["lat2"] for p in pts], pen=f"1.2p,{col},-")
        print(f"  {name}: 震中→质心 {inv.s13/1000:.0f} km，质心端来向 BAZ {(inv.Position(inv.s13)['azi2']+180)%360:.2f}°")

# 组标签（放在各子阵旁边，白底框）
fig.text(x=info["G1"]["lon"] - 2.6, y=info["G1"]["lat"] + 1.0, text="G1", font="11p,Helvetica-Bold,red3",
         justify="CM", no_clip=True, **BOX)
fig.text(x=info["G2"]["lon"] + 2.6, y=info["G2"]["lat"] + 1.0, text="G2", font="11p,Helvetica-Bold,blue3",
         justify="CM", no_clip=True, **BOX)

# 事件文字（两行，锚在沙滩球左上，向右会出框）
fig.text(x=EVLO + 0.4, y=EVLA + 2.2, text=f"Epicenter\n{EVENT_LABEL}", font=F_FONT,
         justify="BR", no_clip=True, **BOX)

# 图例：pygmt 只接受「spec 文件名」-> 写临时文件后传路径（同 step1_map_and_table.py 的做法）
import tempfile
lf = tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False)
lf.write(LEG_SPEC)
lf.close()
with pygmt.config(FONT_ANNOT_PRIMARY=LEG_FONT, FONT_LABEL=LEG_FONT):
    fig.legend(spec=lf.name, position=LEG_POS, box="+gwhite+p0.5p,black")
os.unlink(lf.name)

os.makedirs(OUT_DIR, exist_ok=True)
out = os.path.join(OUT_DIR, OUT_NAME)
if SAVE_PNG:
    fig.savefig(out + ".png", dpi=OUT_DPI)
if SAVE_PDF:
    fig.savefig(out + ".pdf")
print("saved:", out + (".png" if SAVE_PNG else "") + (" / .pdf" if SAVE_PDF else ""))
