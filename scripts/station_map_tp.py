#!/usr/bin/env python3
"""
台站分布图 — T-P 理论到时着色
"""
import numpy as np
import json
from geographiclib.geodesic import Geodesic
import pygmt, tempfile, os

OUT_DIR = "/mnt/f/sharedate/case/2014-01-13-mw63-process/process_flow/1_station_distribution"

EVLA, EVLO = 19.0428, -66.8096
EVDP = 20.0

REGION_FULL = [-95, -65, 17, 50]
REGION_ZOOM = [-93, -65, 30, 50]
PROJ_LEFT = 'M14c'
PROJ_RIGHT = 'M17c'

os.makedirs(OUT_DIR, exist_ok=True)
geod = Geodesic.WGS84

# ── 读取理论走时 ──
with open(f"{OUT_DIR}/all_travel_times.json") as f:
    data = json.load(f)

lons, lats, t_tp, stas = [], [], [], []
for d in data:
    if d['tt_tp'] is None or d['tt_tp'] <= 0:
        continue
    lons.append(d['lon'])
    lats.append(d['lat'])
    t_tp.append(d['tt_tp'] / 60.0)  # 转为分钟
    stas.append(d['sta'])

t_min = np.floor(min(t_tp))
t_max = np.ceil(max(t_tp))
print(f"{len(t_tp)} stations | T-P {t_min:.0f}–{t_max:.0f} min")

# ═══════════════════════════════════════════════════
fig = pygmt.Figure()

# ═══ 左栏：全域 ═══
fig.grdimage(grid='@earth_relief_30s', region=REGION_FULL, projection=PROJ_LEFT,
             shading='+a45+nt0.5', cmap='gray', transparency=50)

fig.coast(shorelines='0.2p,gray40',
          borders=['1/0.25p,gray50', '2/0.15p,gray70'],
          area_thresh=500, resolution='f')

fig.basemap(frame=['xa5f1g5', 'ya5f1g5',
                   f'WSne+t2014-01-13 Mw 6.3 Puerto Rico  |  T-P theoretical arrival'])

# 震中距圆弧
for rr in range(500, 4001, 500):
    fig.plot(x=[EVLO], y=[EVLA], style=f'c{rr/111.12:.1f}c',
             pen='0.5p,gray30,--')

# 台站 — T-P 理论到时着色
pygmt.makecpt(cmap='hot', series=[t_min, t_max, 1], reverse=True)
tf_s = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False)
for lo, la, tt in zip(lons, lats, t_tp):
    tf_s.write(f"{lo} {la} {tt}\n")
tf_s.close()
fig.plot(data=tf_s.name, style='t0.18c', cmap=True, pen='0')
os.unlink(tf_s.name)

# 震源机制球
fig.meca(spec={'strike1': 81, 'dip1': 27, 'rake1': 15,
               'strike2': 337, 'dip2': 83, 'rake2': 116,
               'mantissa': 3.97, 'exponent': 25},
         longitude=EVLO, latitude=EVLA, depth=10,
         scale='0.5c', compression_fill='red', extension_fill='white',
         outline='0.5p,black')
fig.text(x=EVLO+0.5, y=EVLA+0.8, text='Mw 6.3',
         font='9p,Helvetica-Bold,red', fill='white@60', pen='0')

# 色标
fig.colorbar(frame='xa2f1+l"T-P arrival time (min)"',
             position='jBL+w6c/0.25c+h+o1c/0.6c+ml')

# 放大区域框
zx = [REGION_ZOOM[0], REGION_ZOOM[1], REGION_ZOOM[1], REGION_ZOOM[0], REGION_ZOOM[0]]
zy = [REGION_ZOOM[2], REGION_ZOOM[2], REGION_ZOOM[3], REGION_ZOOM[3], REGION_ZOOM[2]]
fig.plot(x=zx, y=zy, pen='1.5p,black,-.')

# ═══ 右栏：台网放大 ═══
fig.shift_origin(xshift='16c')

fig.grdimage(grid='@earth_relief_30s', region=REGION_ZOOM, projection=PROJ_RIGHT,
             shading='+a45+nt0.5', cmap='gray', transparency=50)

fig.coast(shorelines='0.3p,gray50',
          borders=['1/0.3p,gray60', '2/0.2p,gray80'],
          area_thresh=200, resolution='f')

fig.basemap(frame=['xa2f1g2', 'ya2f1g2', 'WSne'])

# 台站 (右栏)
pygmt.makecpt(cmap='hot', series=[t_min, t_max, 1], reverse=True)
tf_r = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False)
for lo, la, tt in zip(lons, lats, t_tp):
    tf_r.write(f"{lo} {la} {tt}\n")
tf_r.close()
fig.plot(data=tf_r.name, style='t0.38c', cmap=True, pen='0')
for lo, la, sta in zip(lons, lats, stas):
    fig.text(x=lo+0.28, y=la-0.28, text=sta,
             font='5p,Helvetica,gray30', fill='white@40', pen='0')
os.unlink(tf_r.name)

# 震源机制球 (右栏)
fig.meca(spec={'strike1': 81, 'dip1': 27, 'rake1': 15,
               'strike2': 337, 'dip2': 83, 'rake2': 116,
               'mantissa': 3.97, 'exponent': 25},
         longitude=EVLO, latitude=EVLA, depth=10,
         scale='0.55c', compression_fill='red', extension_fill='white',
         outline='0.5p,black')
fig.text(x=EVLO+1.0, y=EVLA+0.6, text='Mw 6.3',
         font='8p,Helvetica-Bold,red', fill='white@60', pen='0')

fig.text(x=REGION_ZOOM[1], y=REGION_ZOOM[3]-0.5, text='Station Network (Zoom)',
         font='10p,Helvetica-Bold,black', fill='white@60', pen='0',
         justify='TR')

fig.text(x=(REGION_ZOOM[0]+REGION_ZOOM[1])/2, y=REGION_ZOOM[3]+1.5,
         text=f"{EVLA}\u00b0N, {abs(EVLO)}\u00b0W  |  Depth {EVDP} km  |  {len(t_tp)} stations",
         font='14p,Helvetica-Bold,black', justify='CM', no_clip=True)

fig.basemap(map_scale='jBL+o0.5c/0.5c+w200k')

# ═══ 保存 ═══
OUT_PNG = f"{OUT_DIR}/figure_map_tp.png"
fig.savefig(OUT_PNG, dpi=300)
print(f"\nSaved: {OUT_PNG}")
print("Done.")
