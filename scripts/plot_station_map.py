#!/usr/bin/env python3
"""
通用台站分布图 — 出版级
要素：
- 15s地形图 (geo色标，蓝绿色海洋+棕色陆地)
- SNR着色台站 (viridis)
- 震中距圆弧 + 距离标注
- 两条方位角线（最近/最远海岸）
- 水深1000m等深线（黑色虚线）
- 台站名字标注（SNR>6dB）
- 地震信息标注
- 海洋80%白色覆盖

用法: python3 plot_station_map.py <event>
输出: figure_map.png
"""

import sys, os, csv, yaml, argparse, math, tempfile, warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, "/home/sy/codex")

import numpy as np
import pygmt
from geographiclib.geodesic import Geodesic

CONFIG_PATH = os.path.join(os.path.dirname(__file__), '..', 'config', 'events.yaml')

def plot_azimuth_line(fig, evla, evlo, az, region, pen_style):
    """绘制方位角线并标注"""
    geod = Geodesic.WGS84
    line_lons, line_lats = [], []
    for dkm in range(0, 6000, 50):
        g = geod.Direct(evla, evlo, az, dkm * 1000)
        lo, la = g['lon2'], g['lat2']
        if region[0]-2 <= lo <= region[1]+2 and region[2]-2 <= la <= region[3]+2:
            line_lons.append(lo)
            line_lats.append(la)
        else:
            break
    
    if len(line_lons) > 2:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tf:
            for lo, la in zip(line_lons, line_lats):
                tf.write(f'{lo} {la}\n')
            tn = tf.name
        fig.plot(data=tn, pen=pen_style)
        os.unlink(tn)
        
        # 标注方位角（在圆弧外端）
        label_dist = 4800
        g = geod.Direct(evla, evlo, az, label_dist * 1000)
        if region[0] <= g['lon2'] <= region[1] and region[2] <= g['lat2'] <= region[3]:
            fig.text(
                x=[g['lon2']], y=[g['lat2']], 
                text=f'{az:.0f}°',
                font='8p,Helvetica-Bold,black',
                fill='white@60',
                pen='0.3p'
            )

def main():
    p = argparse.ArgumentParser(description='通用台站分布图 (出版级)')
    p.add_argument('event', help='events.yaml 中的事件 key')
    p.add_argument('--dpi', type=int, default=300, help='输出分辨率 (默认300)')
    p.add_argument('--snr-min-label', type=float, default=6, help='台站名字标注的SNR阈值 (默认6dB)')
    args = p.parse_args()

    # 读取配置
    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f)[args.event]
    case_dir = cfg['case_dir']
    evla = cfg['lat']
    evlo = cfg['lon']
    evdp = cfg['depth_km']
    evmw = cfg['mw']
    evdate = cfg['date']
    
    # 路径
    OUT = f"/mnt/e/sharedate/case/{case_dir}/process_flow/1_station_distribution"
    CSV = f"{OUT}/stations.csv"
    
    # 读取台站数据
    stations = []
    with open(CSV) as f:
        r = csv.DictReader((l for l in f if not l.startswith('#')))
        for row in r:
            stations.append(dict(
                sta=row['sta'], 
                lon=float(row['lon']), 
                lat=float(row['lat']),
                tt_tp=float(row['tt_tp']), 
                snr=float(row['snr']),
                land=float(row['land_km']), 
                az=float(row['az'])
            ))
    
    lons = [s['lon'] for s in stations]
    lats = [s['lat'] for s in stations]
    snr_vals = [s['snr'] for s in stations if np.isfinite(s['snr'])]
    
    snr_min = max(0, min(snr_vals)) if snr_vals else 0
    snr_max = max(snr_vals) if snr_vals else 10
    snr_cpt_min = math.floor(snr_min)
    snr_cpt_max = math.ceil(snr_max) if np.isfinite(snr_max) else 10
    geod = Geodesic.WGS84
    
    # 找到最近和最远海岸距离的台站
    sorted_by_land = sorted(stations, key=lambda s: s['land'])
    min_land_sta = sorted_by_land[0]
    max_land_sta = sorted_by_land[-1]
    
    print(f"事件: {evdate} Mw{evmw}")
    print(f"台站数量: {len(stations)}")
    print(f"SNR范围: {snr_min:.1f}-{snr_max:.1f} dB (色标: {snr_cpt_min}-{snr_cpt_max})")
    print(f"最近海岸: {min_land_sta['sta']} land={min_land_sta['land']:.0f}km az={min_land_sta['az']:.1f}°")
    print(f"最远海岸: {max_land_sta['sta']} land={max_land_sta['land']:.0f}km az={max_land_sta['az']:.1f}°")
    
    # 自动区域 (padding=3度)
    pad = 3
    REGION = [min(min(lons), evlo)-pad, max(max(lons), evlo)+pad,
              min(min(lats), evla)-pad, max(max(lats), evla)+pad]
    w = min(25, max(18, (REGION[1]-REGION[0])*0.4))
    print(f"区域: [{REGION[0]:.1f}, {REGION[1]:.1f}, {REGION[2]:.1f}, {REGION[3]:.1f}]")
    print(f"画布: M{w:.0f}c")
    
    # ── 开始绘图 ──
    fig = pygmt.Figure()
    
    # 1. 地形图 (geo色标)
    print("绘制地形...")
    fig.grdimage(
        grid='@earth_relief_15s', 
        region=REGION, 
        projection=f'M{w:.0f}c',
        shading='+a30+nt1.5',
        cmap='geo'
    )
    
    # 2. 海岸线 + 国界 + 80%白色覆盖
    fig.coast(
        shorelines='0.5p,black',
        borders=['1/0.3p,gray50', '2/0.2p,gray70'],
        area_thresh=500, 
        resolution='h',
        water='white@80'
    )
    
    # 3. 水深1000m等深线（黑色虚线，放在覆盖层上面）
    print("绘制1000m等深线...")
    fig.grdcontour(
        grid='@earth_relief_15s', 
        region=REGION,
        levels=[-1000], 
        pen='0.8p,black,-'
    )
    
    # 4. 震中距圆弧（每500km一条）
    print("绘制震中距圆弧...")
    for rr in [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]:
        radius_deg = rr / 111.12
        fig.plot(
            x=[evlo], y=[evla], 
            style=f'c{radius_deg:.1f}c',
            pen='0.5p,gray40,--'
        )
        # 标注距离（放在右侧）
        g = geod.Direct(evla, evlo, 90, rr * 1000)
        if REGION[0] <= g['lon2'] <= REGION[1] and REGION[2] <= g['lat2'] <= REGION[3]:
            fig.text(
                x=[g['lon2']], y=[g['lat2']], 
                text=f'{rr}',
                font='7p,Helvetica,gray30',
                angle=0,
                justify='LM'
            )
    
    # 5. 方位角参考线（只画两条：最近和最远海岸，深灰色）
    print("绘制方位角参考线...")
    plot_azimuth_line(fig, evla, evlo, min_land_sta['az'], REGION, '1.2p,gray30,-')
    plot_azimuth_line(fig, evla, evlo, max_land_sta['az'], REGION, '1.2p,gray30,-')
    
    # 6. 台站着色 (SNR)
    print("绘制台站...")
    pygmt.makecpt(cmap='viridis', series=[snr_cpt_min, snr_cpt_max, 1])
    
    tf = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False)
    for lo, la, s in zip(lons, lats, snr_vals):
        tf.write(f"{lo} {la} {s}\n")
    tf.close()
    
    fig.plot(
        data=tf.name, 
        style='t0.25c',
        cmap=True, 
        pen='0.3p,black'
    )
    os.unlink(tf.name)
    
    # 6b. 标注台站名字（SNR>阈值，带白色底框）
    print(f"标注台站名字 (SNR>{args.snr_min_label}dB)...")
    for sta in stations:
        if sta['snr'] > args.snr_min_label:
            fig.text(
                x=[sta['lon']], y=[sta['lat']],
                text=sta['sta'],
                font='6p,Helvetica-Bold,black',
                justify='TL',
                offset='0.1c/-0.1c',
                fill='white@70',
                pen='0.2p,gray50'
            )
    
    # 7. 震中（红色星）
    fig.plot(
        x=[evlo], y=[evla], 
        style='a0.4c',
        pen='0.8p,red', 
        fill='red'
    )
    
    # 8. 震源机制解位置预留（留空）
    
    # 9. 地震信息标注（左上角）
    info_text = (
        f"{evdate} UTC\\n"
        f"Mw {evmw}  Depth {evdp:.0f} km\\n"
        f"{abs(evla):.2f}°{'S' if evla<0 else 'N'} {abs(evlo):.2f}°{'E' if evlo>0 else 'W'}"
    )
    fig.text(
        x=[REGION[0]+1], y=[REGION[3]-1.5],
        text=info_text,
        font='10p,Helvetica-Bold,black',
        fill='white@70',
        pen='0.5p,black',
        justify='TL'
    )
    
    # 10. 震中标签
    fig.text(
        x=[evlo+1.5], y=[evla+1.5], 
        text=f'Mw {evmw}',
        font='10p,Helvetica-Bold,red', 
        fill='white@50', 
        pen='0.3p,red'
    )
    
    # 11. 色标
    fig.colorbar(
        frame='xa2f1+l"SNR (dB)"',
        position='jBL+w6c/0.25c+h+o1c/0.6c+ml'
    )
    
    # 12. 坐标轴（带经纬度标注）
    fl = max(5, int((REGION[1]-REGION[0])/8))
    fa = max(5, int((REGION[3]-REGION[2])/8))
    fig.basemap(
        frame=[
            f'xa{fl}f{fl//2}+l"Longitude"',
            f'ya{fa}f{fa//2}+l"Latitude"',
            f'WSne+t{evdate} Mw {evmw}  |  Station Distribution'
        ]
    )
    
    # 13. 保存
    output = f"{OUT}/figure_map.png"
    fig.savefig(output, dpi=args.dpi)
    print(f"\n保存: {output}")

if __name__ == '__main__':
    main()
