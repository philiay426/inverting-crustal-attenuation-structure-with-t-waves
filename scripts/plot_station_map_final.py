#!/usr/bin/env python3
"""
台站分布图（期刊级）
- 15s地形 + geo色标 + SNR着色
- 双色标: 海拔(底部) + SNR(右侧)
- 图例: 台站/震中/-1000m
- 方位角射线 + 角度标注
- 智能标注: land最近+最远+SNR最高
- 文字白描边
- 排除框: 指定矩形区域的台站不显示
"""
import os, csv, numpy as np, pygmt, tempfile, math, sys, yaml
from geographiclib.geodesic import Geodesic

CODE_ROOT = "/home/sy/codex"
CONFIG_PATH = os.path.join(CODE_ROOT, "config", "events.yaml")

STATION_SIZE = 0.5
STATION_PEN  = "0.3p,black"
LABEL_FONT   = "10p,Helvetica-Bold,black"
LABEL_MIN_DIST = 0.8
WHITE_STROKE = 0.03
MAP_WIDTH    = "18c"
DPI = 200
MARGIN = 6
GRID_RES = "@earth_relief_15s"

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('event')
    args = parser.parse_args()

    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f)[args.event]

    EVLA, EVLO = cfg['lat'], cfg['lon']
    EVDP, EVMW = cfg['depth_km'], cfg['mw']
    evdate = cfg['date']
    case_dir = cfg['case_dir']
    BASE = f"/mnt/f/sharedate/case/{case_dir}"
    OUT  = f"{BASE}/process_flow/1_station_distribution"
    CSV  = f"{OUT}/stations.csv"
    os.makedirs(OUT, exist_ok=True)

    # 排除框: [lon_min, lon_max, lat_min, lat_max]，落在框内的台站被剔除
    EXCLUDE_BOXES = [
        [32.0, 37.5, -11.5, -7.0],
    ]

    def _excluded(lon, lat):
        for lo, hi, la, ha in EXCLUDE_BOXES:
            if lo <= lon <= hi and la <= lat <= ha:
                return True
        return False

    stations = []
    n_dropped = 0
    with open(CSV) as f:
        for row in csv.DictReader((l for l in f if not l.startswith("#"))):
            lon, lat = float(row["lon"]), float(row["lat"])
            if _excluded(lon, lat):
                n_dropped += 1
                continue
            stations.append(dict(sta=row["sta"], lon=lon,
                                 lat=lat, snr=float(row["snr"]),
                                 land_km=float(row["land_km"])))
    print(f"台站: {len(stations)} (剔除 {n_dropped})")

    snrs = [s["snr"] for s in stations]
    lons_all = [s["lon"] for s in stations] + [EVLO]
    lats_all = [s["lat"] for s in stations] + [EVLA]
    R = [min(lons_all)-MARGIN, max(lons_all)+MARGIN,
         min(lats_all)-MARGIN, max(lats_all)+MARGIN]

    pygmt.config(MAP_GRID_PEN_PRIMARY="0.3p,gray40,--", PS_MEDIA="A0")
    fig = pygmt.Figure()

    elev_cpt = os.path.join(tempfile.gettempdir(), "elev_relief.cpt")
    pygmt.makecpt(cmap="geo", series=[-8000, 5000, 100], continuous=True, output=elev_cpt)

    fig.grdimage(grid=GRID_RES, region=R, projection=f"M{MAP_WIDTH}",
                 cmap=elev_cpt, shading="+a315+nt0.8+m0")

    fig.coast(shorelines="0.3p,gray30", borders=["1/0.3p,gray50","2/0.2p,gray70"],
              resolution="h", water="lightblue@60", land="white@80")

    fig.basemap(frame=["xa10f5", "ya10f5", "WSne"])

    fig.grdcontour(grid=GRID_RES, region=R, levels=[-1000],
                   pen="0.8p,black,--")

    snr_valid = [s for s in snrs if np.isfinite(s) and s > 0]
    snr_min, snr_max = max(0, min(snr_valid)), max(snr_valid) if snr_valid else 10
    pygmt.makecpt(cmap="viridis", series=[math.floor(snr_min), math.ceil(snr_max), 1])

    tf = tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False)
    for s in stations: tf.write(f"{s['lon']} {s['lat']} {s['snr']}\n")
    tf.close()
    fig.plot(data=tf.name, style=f"t{STATION_SIZE}c", cmap=True, pen=STATION_PEN)
    os.unlink(tf.name)

    labeled = []
    top = sorted(stations, key=lambda x: x["snr"], reverse=True)[:15]
    sel = {}
    by_land = sorted(top, key=lambda x: x["land_km"])
    sel[by_land[0]["sta"]] = by_land[0]
    sel[by_land[-1]["sta"]] = by_land[-1]
    sel[top[0]["sta"]] = top[0]
    label_stations = list(sel.values())

    def _offsets(s):
        prefs = []
        if s["lon"] < EVLO: prefs.extend([(0.8, 0.5), (0.8, -0.5)])
        else:              prefs.extend([(-0.8, 0.5), (-0.8, -0.5)])
        if s["lat"] < EVLA: prefs.extend([(0.5, 0.8), (-0.5, 0.8)])
        else:               prefs.extend([(0.5, -0.8), (-0.5, -0.8)])
        return prefs

    geod = Geodesic.WGS84
    for s in label_stations:
        fixed_offsets = {"AMPY": (-0.8, -0.5)}
        if s["sta"] in fixed_offsets:
            best_dx, best_dy = fixed_offsets[s["sta"]]
        else:
            best_dx, best_dy = 0.8, 0.5
            for dx, dy in _offsets(s):
                lx_test, ly_test = s["lon"]+dx, s["lat"]+dy
                if not any(abs(lx_test-lx)<LABEL_MIN_DIST and abs(ly_test-ly)<LABEL_MIN_DIST
                           for lx, ly in labeled):
                    best_dx, best_dy = dx, dy
                    break
        labeled.append((s["lon"], s["lat"]))
        lx, ly = s["lon"]+best_dx, s["lat"]+best_dy
        fig.plot(x=[s["lon"], lx], y=[s["lat"], ly], pen="0.4p,gray30,--")
        for ddx, ddy in [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]:
            fig.text(x=lx, y=ly, text=s["sta"], font="10p,Helvetica-Bold,white",
                     justify="BL", offset=f"{ddx*WHITE_STROKE}c/{ddy*WHITE_STROKE}c")
        fig.text(x=lx, y=ly, text=s["sta"], font=LABEL_FONT, justify="BL")

        g = geod.Inverse(EVLA, EVLO, s["lat"], s["lon"])
        az = (g["azi1"] + 360) % 360
        g2 = geod.Direct(EVLA, EVLO, az, g["s12"] + 3.0 * 111000)
        fig.plot(x=[EVLO, g2["lon2"]], y=[EVLA, g2["lat2"]],
                 pen="0.5p,black,--")
        fig.text(x=g2["lon2"]+0.5, y=g2["lat2"]+0.3,
                 text=f"{az:.0f}°", font="9p,Helvetica,black",
                 fill="white@70", pen="0")

    fig.plot(x=EVLO, y=EVLA, style="a0.8c", fill="red", pen="0.7p,black")

    if EVLO - 3 > R[0] + 1:
        lx, ly = EVLO - 3, EVLA - 1; jt = "BR"
    elif EVLO + 3 < R[1] - 1:
        lx, ly = EVLO + 3, EVLA - 1; jt = "BL"
    else:
        lx, ly = EVLO, EVLA - 2; jt = "CT"
    etxt = f"{evdate}  Mw {EVMW}"
    for ddx, ddy in [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]:
        fig.text(x=lx, y=ly, text=etxt, font="14p,Helvetica-Bold,white",
                 justify=jt, offset=f"{ddx*WHITE_STROKE}c/{ddy*WHITE_STROKE}c")
    fig.text(x=lx, y=ly, text=etxt, font="14p,Helvetica-Bold,black", justify=jt)

    fig.colorbar(frame="xa10f1+lSNR (dB)", position="JMR+w7c/0.3c+o1c/0c+v")
    fig.colorbar(cmap=elev_cpt, frame="xa2000f500+lElevation (m)",
                 position="JBC+w10c/0.3c+h+o0c/1.6c")

    legend_spec = (
        "G 0.1c\n"
        f"S 0.3c t {STATION_SIZE}c gray30 {STATION_PEN} 0.8c Station (colored by SNR)\n"
        "G 0.15c\n"
        "S 0.3c a 0.4c red 0.7p,black 0.8c Epicenter\n"
        "G 0.15c\n"
        "S 0.3c - 0.6c - 0.8p,black,-- 0.8c -1000 m depth\n"
    )
    lf = tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False)
    lf.write(legend_spec); lf.close()
    fig.legend(spec=lf.name, position="jBR+w4.3c+o0.3c/0.3c", box="+gwhite+p0.5p,black")
    os.unlink(lf.name)

    fig.text(x=EVLO+1, y=EVLA+1, text="-1000 m", font="8p,Helvetica,black", angle=-30)

    out_file = os.path.join(OUT, "figure_map_final.png")
    fig.savefig(out_file, dpi=DPI)
    print(f"✓ {out_file}")

if __name__ == '__main__':
    main()
