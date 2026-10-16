#!/usr/bin/env python3
"""
逐站转换点计算 — 通用事件版 (event-agnostic)
方法: 大圆路径与 ETOPO -1000m 等深线的交点 = 唯一转换点

用法:
  python find_conversion_point.py puerto_rico
  python find_conversion_point.py mindanao_2026

输出:
  {sta}_conversion_point.png — 地图标注
  station_conversion_points.csv — 逐站转换点表 (批量模式)
"""

import os, argparse, yaml, numpy as np, tempfile, csv
from geographiclib.geodesic import Geodesic
import pygmt

# ---------------------------------------------------------------------------
# 参数解析 + 配置读取
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='Per-station conversion point calculation')
parser.add_argument('event', help='Event name key in config/events.yaml (e.g. puerto_rico)')
parser.add_argument('--config', default=os.path.join(os.path.dirname(__file__), '..', 'config', 'events.yaml'),
                    help='Path to events.yaml config file')
args = parser.parse_args()

config_path = os.path.abspath(args.config)
with open(config_path) as f:
    all_events = yaml.safe_load(f)

if args.event not in all_events:
    available = ', '.join(all_events.keys())
    parser.error(f"Event '{args.event}' not found in {config_path}. Available: {available}")

evt = all_events[args.event]

EVLO = float(evt['lon'])
EVLA = float(evt['lat'])
MW = float(evt.get('mw', 0))
DEPTH_KM = float(evt.get('depth_km', 0))
CASE_DIR = evt['case_dir']

V_OCEAN = 1.48
V_LAND = 3.6

GEOD = Geodesic.WGS84

# 动态路径: 从 step1 输出读台站, 结果写到 step3
STATIONS_CSV = f'/mnt/e/sharedate/case/{CASE_DIR}/process_flow/1_station_distribution/stations.csv'
OUT_DIR = f'/mnt/e/sharedate/case/{CASE_DIR}/process_flow/3_conversion_point'

pygmt.config(GMT_DATA_SERVER='https://oceania.generic-mapping-tools.org')


def find_conversion_point(stlo, stla, grid):
    """
    找大圆路径(震源→台站)与 -1000m 等深线的交点
    返回 (conv_lon, conv_lat, ocean_km, land_km) 或 None
    """
    g = GEOD.Inverse(EVLA, EVLO, stla, stlo)
    azi = g['azi1']
    total_km = g['s12'] / 1000.0
    
    lat_arr = grid.lat.values
    lon_arr = grid.lon.values
    el_arr = grid.values
    
    prev_el = None
    for d_km in np.arange(10, total_km, 2):
        pt = GEOD.Direct(EVLA, EVLO, azi, d_km * 1000)
        i_lat = np.argmin(np.abs(lat_arr - pt['lat2']))
        i_lon = np.argmin(np.abs(lon_arr - pt['lon2']))
        el = float(el_arr[i_lat, i_lon])
        
        if prev_el is not None:
            if (prev_el > -1000 and el <= -1000) or (prev_el <= -1000 and el > -1000):
                conv_lat = (pt['lat2'] + prev_lat) / 2
                conv_lon = (pt['lon2'] + prev_lon) / 2
                ocean_km = (d_km + prev_d_km) / 2
                land_km = total_km - ocean_km
                return conv_lon, conv_lat, ocean_km, land_km
        
        prev_el, prev_lat, prev_lon, prev_d_km = el, pt['lat2'], pt['lon2'], d_km
    return None


def compute_region(stations):
    """根据震源和台站位置动态计算地图范围"""
    all_lons = [EVLO] + [s['lon'] for s in stations]
    all_lats = [EVLA] + [s['lat'] for s in stations]
    lon_min = min(all_lons) - 3
    lon_max = max(all_lons) + 3
    lat_min = min(all_lats) - 3
    lat_max = max(all_lats) + 3
    return [lon_min, lon_max, lat_min, lat_max]


def plot_single_station(sta, stlo, stla, conv_lon, conv_lat, ocean_km, land_km, region):
    """画单站转换点 GMT 图"""
    g = GEOD.Inverse(EVLA, EVLO, stla, stlo)
    azi = g['azi1']
    
    # 大圆路径
    lons, lats = [], []
    for d_km in np.arange(0, g['s12'] / 1000 + 50, 5):
        pt = GEOD.Direct(EVLA, EVLO, azi, d_km * 1000)
        lons.append(pt['lon2']); lats.append(pt['lat2'])
    
    fig = pygmt.Figure()
    
    fig.grdimage(grid='@earth_relief_30s', region=region, projection='M16c',
                 shading='+a45+nt0.5', cmap='gray', transparency=50)
    fig.coast(shorelines='0.3p,gray50', borders=['1/0.3p,gray60'], 
              area_thresh=500, resolution='f')
    fig.basemap(frame=['xa5f2g5', 'ya5f2g5', 
                       f'WSne+t{sta} Conversion Point ({azi:.1f}° Azimuth)'])
    
    # 大圆线
    tf = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False)
    for lo, la in zip(lons, lats):
        tf.write(f'{lo} {la}\n')
    tf.close()
    fig.plot(data=tf.name, pen='1.5p,black,-')
    os.unlink(tf.name)
    
    # 震源
    mw_label = f'Mw {MW}' if MW else ''
    fig.plot(x=EVLO, y=EVLA, style='a0.55c', fill='red', pen='0.8p,black')
    if mw_label:
        fig.text(x=EVLO+1, y=EVLA+0.5, text=mw_label, 
                 font='10p,Helvetica-Bold,red', fill='white@60')
    
    # 台站
    fig.plot(x=stlo, y=stla, style='t0.4c', fill='dodgerblue', pen='0.8p,black')
    fig.text(x=stlo+0.5, y=stla-0.8, text=sta, 
             font='10p,Helvetica-Bold,dodgerblue', fill='white@60')
    
    # 转换点
    fig.plot(x=conv_lon, y=conv_lat, style='d0.5c', fill='yellow', pen='1.2p,black')
    fig.text(x=conv_lon+0.8, y=conv_lat+0.4, 
             text=f'C: {conv_lat:.1f}N, {abs(conv_lon):.1f}W', 
             font='9p,Helvetica-Bold,black', fill='yellow@70')
    
    # 海洋段 / 陆地段
    fig.plot(x=[EVLO, conv_lon], y=[EVLA, conv_lat], pen='2p,blue,..', transparency=60)
    fig.plot(x=[conv_lon, stlo], y=[conv_lat, stla], pen='2p,red,..', transparency=60)
    
    mid_oc = ((EVLO + conv_lon)/2, (EVLA + conv_lat)/2)
    mid_ld = ((conv_lon + stlo)/2, (conv_lat + stla)/2)
    fig.text(x=mid_oc[0], y=mid_oc[1]+1.5, text=f'ocean {ocean_km:.0f}km', 
             font='8p,Helvetica,blue', fill='white@60')
    fig.text(x=mid_ld[0]-2, y=mid_ld[1], text=f'land {land_km:.0f}km (Love v={V_LAND})', 
             font='8p,Helvetica,red', fill='white@60')
    
    out = os.path.join(OUT_DIR, f'{sta}_conversion_point.png')
    fig.savefig(out, dpi=300)
    print(f'  {out}')
    return out


def main():
    # 加载台站 — 从 step1 输出 stations.csv
    # 格式: sta,az,dist,land_km,ocean_km,tt_p,tt_s,tt_tp,tt_ts,snr,lon,lat,baz
    stations = []
    with open(STATIONS_CSV) as f:
        reader = csv.DictReader(
            (line for line in f if not line.startswith('#')),
        )
        for row in reader:
            try:
                stations.append({
                    'sta': row['sta'],
                    'lon': float(row['lon']),
                    'lat': float(row['lat']),
                    'snr': float(row['snr']),
                    'az': float(row['az']),
                    'land_km': float(row['land_km']),
                    'ocean_km': float(row['ocean_km']),
                })
            except (ValueError, KeyError):
                pass
    
    print(f'Event: {args.event} (Mw {MW}, {EVLA}N {EVLO}E)')
    print(f'Stations CSV: {STATIONS_CSV}')
    print(f'台站总数: {len(stations)}')
    print(f'Output dir: {OUT_DIR}')
    
    os.makedirs(OUT_DIR, exist_ok=True)
    
    # 加载 ETOPO — 根据台站分布动态计算区域
    region = compute_region(stations)
    print(f'Map region: {region}')
    grid = pygmt.datasets.load_earth_relief(resolution='02m', region=region)
    
    results = []
    for i, s in enumerate(stations):
        result = find_conversion_point(s['lon'], s['lat'], grid)
        if result is None:
            msg = 'FAILED'
            print(f'  [{i+1}/{len(stations)}] {s["sta"]}: {msg}')
            continue
        
        conv_lon, conv_lat, ocean_km, land_km = result
        ocean_s = ocean_km / V_OCEAN
        land_s_love = land_km / V_LAND  # Love 波陆地走时
        
        g = GEOD.Inverse(EVLA, EVLO, s['lat'], s['lon'])
        azi = g['azi1'] % 360
        
        results.append({
            **s,
            'conv_lon': conv_lon, 'conv_lat': conv_lat,
            'ocean_km': ocean_km, 'land_km': land_km,
            'ocean_s': ocean_s, 'land_s_love': land_s_love,
            't_total_love': ocean_s + land_s_love,
            'azi_geodesic': azi,
        })
        
        if (i + 1) % 15 == 0:
            info = (f'  [{i+1}/{len(stations)}] {s["sta"]}: '
                    f'C=({conv_lat:.2f}N,{conv_lon:.2f}W) '
                    f'ocean={ocean_km:.0f}km land={land_km:.0f}km')
            print(info)
    
    print(f'\n成功: {len(results)}/{len(stations)}')
    
    # 按纬度排序
    results.sort(key=lambda r: r['conv_lat'])
    
    # 输出 CSV
    out_csv = os.path.join(OUT_DIR, 'station_conversion_points.csv')
    with open(out_csv, 'w') as f:
        f.write(f'# Event: {args.event} (Mw {MW}, {EVLA}N {EVLO}E)\n')
        f.write('# Per-station conversion points — great circle intersection with -1000m isobath\n')
        f.write(f'# ocean velocity: {V_OCEAN} km/s (SOFAR), land velocity: {V_LAND} km/s (Love wave)\n')
        f.write('#\n')
        f.write('station,az_deg,conv_lat,conv_lon,ocean_km,land_km,'
                'ocean_s,land_s_love,t_total_love,'
                'station_lon,station_lat,snr\n')
        for r in results:
            f.write(f'{r["sta"]},{r["azi_geodesic"]:.1f},{r["conv_lat"]:.4f},{r["conv_lon"]:.4f},'
                    f'{r["ocean_km"]:.1f},{r["land_km"]:.1f},'
                    f'{r["ocean_s"]:.1f},{r["land_s_love"]:.1f},{r["t_total_love"]:.1f},'
                    f'{r["lon"]:.4f},{r["lat"]:.4f},{r["snr"]:.1f}\n')
    print(f'\n保存: {out_csv}')
    
    # 摘要
    if results:
        lats = [r['conv_lat'] for r in results]
        print(f'\n转换点纬度范围: {min(lats):.2f}°N ~ {max(lats):.2f}°N')
        print(f'转换点经度范围: {min(r["conv_lon"] for r in results):.2f}°W ~ {max(r["conv_lon"] for r in results):.2f}°W')
    
    print('\nDone.')


if __name__ == '__main__':
    main()
