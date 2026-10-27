#!/usr/bin/env python3
"""
步骤1: 台站信息提取 + 到时计算 + 台站分布图 (事件无关版)

从 config/events.yaml 读取事件参数，通过 argparse 指定事件名称。

输出:
- stations.xlsx (格式化Excel表格，供统计分析)
- stations.csv  (简洁CSV，供后续脚本处理)
- figure_map_tp.png (台站分布图，T-P理论走时着色)
  地图规范: 15s geo地形, 每10°方位角弧线, 区域自动推算

用法:
    python step1_station_table_and_map.py puerto_rico
    python step1_station_table_and_map.py mindanao_2026 --snr-min 5.0
"""
import sys
import os
import glob
import warnings
import argparse
import yaml

CODE_ROOT = "/home/sy/codex"
for p in [CODE_ROOT, CODE_ROOT + '/t_phase_analysis']:
    if p not in sys.path:
        sys.path.insert(0, p)

import numpy as np
from obspy import read
from obspy.signal.filter import envelope as obs_env
from scipy.ndimage import uniform_filter1d
from lib.arrivals import compute_theoretical_arrivals, dist_az
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
import pygmt

# ── 默认物理常数 (可通过 config 覆盖) ──
DEFAULT_V_OCEAN = 1.48   # km/s (SOFAR通道)
DEFAULT_V_LOVE = 3.6     # km/s (Love波)
DEFAULT_SNR_THRESHOLD = 8.0  # dB

CONFIG_PATH = os.path.join(CODE_ROOT, "config", "events.yaml")


def load_event_config(event_name, config_path=None):
    """从 events.yaml 加载指定事件的配置"""
    if config_path is None:
        config_path = CONFIG_PATH
    with open(config_path, 'r') as f:
        all_events = yaml.safe_load(f)

    if event_name not in all_events:
        available = ', '.join(all_events.keys())
        print(f"错误: 事件 '{event_name}' 不在配置文件中")
        print(f"可用事件: {available}")
        sys.exit(1)

    return all_events[event_name]


def detect_t_wave_peak(data, delta, tp_theory, search_window=120):
    """
    在T-P理论到时附近检测包络峰值，返回实测T波到时(秒)
    """
    env = obs_env(data.astype(np.float64))
    # 平滑包络(2秒窗口)
    smooth_win = max(int(2 / delta), 1)
    env = uniform_filter1d(env, size=smooth_win)
    
    t = np.arange(len(data)) * delta
    t_start = max(0, tp_theory - search_window)
    t_end = min(t[-1], tp_theory + search_window)
    
    mask = (t >= t_start) & (t <= t_end)
    if not np.any(mask):
        return tp_theory
    
    t_win = t[mask]
    env_win = env[mask]
    peak_idx = np.argmax(env_win)
    return t_win[peak_idx]


def compute_snr_welch(data, delta, signal_center, noise_center, signal_width=80, noise_width=300):
    """
    使用Welch PSD计算T波信噪比(dB)
    信号窗: [signal_center, signal_center + signal_width]
    噪声窗: [noise_center - noise_width, noise_center]
    """
    from scipy.signal import welch
    sr = 1.0 / delta
    t = np.arange(len(data)) * delta
    
    # 信号窗
    s_start = int(max(0, signal_center) / delta)
    s_end = int(min(t[-1], signal_center + signal_width) / delta)
    
    # 噪声窗
    n_start = int(max(0, noise_center - noise_width) / delta)
    n_end = int(max(0, noise_center) / delta)
    
    if s_end <= s_start or n_end <= n_start:
        return -999
    
    f_sig, psd_sig = welch(data[s_start:s_end], fs=sr, nperseg=min(256, s_end - s_start))
    f_noise, psd_noise = welch(data[n_start:n_end], fs=sr, nperseg=min(256, n_end - n_start))
    
    # 2-8 Hz频段能量
    e_sig = np.sum(psd_sig[(f_sig >= 2) & (f_sig <= 8)])
    e_noise = np.sum(psd_noise[(f_noise >= 2) & (f_noise <= 8)])
    
    if e_noise > 0 and e_sig > 0:
        return 10 * np.log10(e_sig / e_noise)
    return -999


def main():
    # ── 解析命令行参数 ──
    parser = argparse.ArgumentParser(
        description="步骤1: 台站信息提取 + 到时计算 + 台站分布图")
    parser.add_argument('event', type=str,
                        help="事件名称 (events.yaml 中的 key，如 puerto_rico)")
    parser.add_argument('--snr-min', type=float, default=None,
                        help=f'SNR阈值覆盖 (默认从 config 读取)')
    parser.add_argument('--config', type=str, default=CONFIG_PATH,
                        help=f'配置文件路径 (默认: {CONFIG_PATH})')
    args = parser.parse_args()

    # ── 加载事件配置 ──
    cfg = load_event_config(args.event, config_path=args.config)

    # 从 config 提取事件参数
    EVENT_DATE = cfg['date']
    EVENT_MW = cfg['mw']
    EVLA = cfg['lat']
    EVLO = cfg['lon']
    EVDP = cfg['depth_km']
    case_dir = cfg['case_dir']

    # 构建路径
    PROC_DIR = f"/mnt/f/sharedate/case/{case_dir}/sacprocess"
    OUT_DIR = f"/mnt/f/sharedate/case/{case_dir}/process_flow/1_station_distribution"

    # SNR阈值: CLI > config > 默认值
    if args.snr_min is not None:
        SNR_THRESHOLD = args.snr_min
    else:
        SNR_THRESHOLD = cfg.get('snr_threshold', DEFAULT_SNR_THRESHOLD)

    # 物理常数: config > 默认值
    V_OCEAN = cfg.get('v_ocean', DEFAULT_V_OCEAN)
    V_LOVE = cfg.get('v_love', DEFAULT_V_LOVE)

    os.makedirs(OUT_DIR, exist_ok=True)

    print(f"{'='*70}")
    print(f"步骤1: 台站信息提取 + 到时计算 + 台站分布图")
    print(f"事件: {EVENT_DATE} Mw {EVENT_MW} | {EVLA}°N, {EVLO}°E | 深度 {EVDP} km")
    print(f"配置: {args.event} | SNR≥{SNR_THRESHOLD}dB | V_ocean={V_OCEAN} V_love={V_LOVE}")
    print(f"{'='*70}\n")
    
    # ── 1. 读取所有BHZ SAC文件 ──
    bhz_files = sorted(glob.glob(f"{PROC_DIR}/*.BHZ*_proc.sac"))
    print(f"找到 {len(bhz_files)} 个 BHZ 文件\n")
    
    rows = []
    seen_sta = set()
    for fn in bhz_files:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                tr = read(fn)[0]
            
            sac = tr.stats.sac
            stla = getattr(sac, 'stla', 0)
            stlo = getattr(sac, 'stlo', 0)
            
            if stla == 0 or stlo == 0:
                continue
            
            sta = os.path.basename(fn).split('.')[1]
            if sta in seen_sta:
                continue
            seen_sta.add(sta)
            dist_km, az, baz = dist_az(EVLO, EVLA, stlo, stla)
            
            # 理论到时
            arrs = compute_theoretical_arrivals(EVLO, EVLA, EVDP, stlo, stla)
            tt_p = arrs.get('P', 0)
            tt_s = arrs.get('S', 0)
            tt_tp = arrs.get('T-P', 0)
            tt_ts = arrs.get('T-S', 0)
            ocean_km = arrs.get('ocean_km', 0)
            land_km = arrs.get('land_km', 0)
            
            # 实测T波到时（包络峰值）
            if tt_tp > 0:
                t_obs = detect_t_wave_peak(tr.data, tr.stats.delta, tt_tp)
                
                # SNR计算（以T-S为中心，噪声在P波之前）
                snr = compute_snr_welch(tr.data, tr.stats.delta, 
                                        signal_center=tt_ts, 
                                        noise_center=tt_p - 100)
            else:
                t_obs = 0
                snr = -999
            
            rows.append({
                'sta': sta,
                'az': round(az, 1),
                'dist': round(dist_km, 0),
                'land_km': round(land_km, 0),
                'ocean_km': round(ocean_km, 0),
                'tt_p': round(tt_p, 1),
                'tt_s': round(tt_s, 1),
                'tt_tp': round(tt_tp, 1),
                'tt_ts': round(tt_ts, 1),
                'snr': round(snr, 1),
                'lon': round(stlo, 4),
                'lat': round(stla, 4),
                'baz': round(baz, 1),
            })
            
        except Exception as e:
            continue
    
    # 按SNR降序排列
    rows.sort(key=lambda x: x['snr'], reverse=True)
    
    snr_ok_count = sum(1 for r in rows if r['snr'] >= SNR_THRESHOLD)
    print(f"总台站: {len(rows)} | SNR≥{SNR_THRESHOLD}dB: {snr_ok_count}站\n")
    
    # ── 2. 生成 Excel 表格 ──
    print("生成 Excel 表格...")
    wb = Workbook()
    ws = wb.active
    ws.title = 'Sheet'
    
    # 事件信息行
    ws.merge_cells('A1:M1')
    info_cell = ws.cell(row=1, column=1,
        value=f"事件: {EVENT_DATE} Mw {EVENT_MW}  |  {EVLA}°N, {EVLO}°E  |  "
              f"深度 {EVDP} km  |  总台站: {len(rows)}  |  "
              f"SNR≥{SNR_THRESHOLD}dB: {snr_ok_count}站")
    info_cell.font = Font(bold=True, size=11, color='1f4e79')
    info_cell.alignment = Alignment(horizontal='left')
    
    # 表头（第2行：中文说明）
    hdrs_cn = ['台站', '方位角', '震中距', '陆地距', '海洋距', 
               'P波', 'S波', 'T-P', 'T-S', '信噪比', '经度', '纬度', '反方位角']
    for ci, h in enumerate(hdrs_cn, 1):
        c = ws.cell(row=2, column=ci, value=h)
        c.font = Font(bold=False, size=9, color='808080')
        c.alignment = Alignment(horizontal='center', vertical='center')
    
    # 表头（第3行：英文字段名）
    hdrs_en = ['STA', 'AZ(°)', 'DIST(km)', 'LAND(km)', 'OCEAN(km)',
               'P(s)', 'S(s)', 'T-P(s)', 'T-S(s)', 'SNR(dB)', 'LON', 'LAT', 'BAZ(°)']
    hdr_fill = PatternFill('solid', fgColor='4472C4')
    hdr_font = Font(bold=True, color='FFFFFF', size=10)
    
    for ci, h in enumerate(hdrs_en, 1):
        c = ws.cell(row=3, column=ci, value=h)
        c.font = hdr_font
        c.fill = hdr_fill
        c.alignment = Alignment(horizontal='center', vertical='center')
    
    # 数据行（从第5行开始）
    for ri, r in enumerate(rows, 5):
        vals = [r['sta'], r['az'], r['dist'], r['land_km'], r['ocean_km'],
                r['tt_p'], r['tt_s'], r['tt_tp'], r['tt_ts'], r['snr'],
                r['lon'], r['lat'], r['baz']]
        
        for ci, v in enumerate(vals, 1):
            cell = ws.cell(row=ri, column=ci, value=v)
            cell.alignment = Alignment(horizontal='center', vertical='center')
            cell.font = Font(size=9)
            
            # SNR高亮（>=阈值绿色）
            if ci == 10:  # SNR列
                if isinstance(v, (int, float)) and v >= SNR_THRESHOLD:
                    cell.fill = PatternFill('solid', fgColor='C6EFCE')
                    cell.font = Font(bold=True, size=9, color='006100')
    
    # 列宽
    widths = [12, 10, 12, 12, 12, 10, 10, 12, 12, 10, 14, 14, 12]
    for i, w in enumerate(widths, 1):
        col_letter = chr(64 + i) if i <= 26 else 'A' + chr(64 + i - 26)
        ws.column_dimensions[col_letter].width = w
    
    # 边框
    thin_border = Border(
        left=Side(style='thin', color='D9D9D9'),
        right=Side(style='thin', color='D9D9D9'),
        top=Side(style='thin', color='D9D9D9'),
        bottom=Side(style='thin', color='D9D9D9')
    )
    for row in ws.iter_rows(min_row=3, max_row=len(rows) + 4, max_col=len(hdrs_en)):
        for cell in row:
            cell.border = thin_border
    
    xlsx_path = f"{OUT_DIR}/stations.xlsx"
    try:
        wb.save(xlsx_path)
    except PermissionError:
        xlsx_path = f"{OUT_DIR}/stations_v2.xlsx"
        wb.save(xlsx_path)
    print(f"✓ {xlsx_path}")
    
    # ── 3. 生成 CSV 文件 ──
    print("生成 CSV 文件...")
    csv_path = f"{OUT_DIR}/stations.csv"
    csv_lines = []
    csv_lines.append("# Event: {} Mw {} | {}N, {}E | Depth {}km\n".format(
        EVENT_DATE, EVENT_MW, EVLA, EVLO, EVDP))
    csv_lines.append("sta,az,dist,land_km,ocean_km,tt_p,tt_s,tt_tp,tt_ts,snr,lon,lat,baz\n")
    for r in rows:
        csv_lines.append("{},{:.1f},{:.0f},{:.0f},{:.0f},{:.1f},{:.1f},{:.1f},{:.1f},{:.1f},{:.4f},{:.4f},{:.1f}\n".format(
            r['sta'], r['az'], r['dist'], r['land_km'], r['ocean_km'],
            r['tt_p'], r['tt_s'], r['tt_tp'], r['tt_ts'], r['snr'],
            r['lon'], r['lat'], r['baz']))
    try:
        with open(csv_path, 'w') as f:
            f.writelines(csv_lines)
    except PermissionError:
        csv_path = f"{OUT_DIR}/stations_v2.csv"
        with open(csv_path, 'w') as f:
            f.writelines(csv_lines)
    print(f"✓ {csv_path}")
    
    # ── 4. 生成台站分布图 (通用版) ──
    print("生成台站分布图...")
    
    lons = [r['lon'] for r in rows]
    lats = [r['lat'] for r in rows]
    t_tp_min = [r['tt_tp'] / 60.0 for r in rows]
    
    t_tp_valid = [t for t in t_tp_min if t > 0]
    if t_tp_valid:
        t_min = np.floor(min(t_tp_valid))
        t_max = np.ceil(max(t_tp_valid))
    else:
        t_min, t_max = 0, 30
    print(f"  T-P走时范围: {t_min:.0f}–{t_max:.0f} min")
    
    # 自动推算区域 (台站范围 + 边距 + 包含震中)
    pad = 8  # 度
    reg_lon_min = min(min(lons), EVLO) - pad
    reg_lon_max = max(max(lons), EVLO) + pad
    reg_lat_min = min(min(lats), EVLA) - pad
    reg_lat_max = max(max(lats), EVLA) + pad
    region = [reg_lon_min, reg_lon_max, reg_lat_min, reg_lat_max]
    lon_span = reg_lon_max - reg_lon_min
    
    fig_width = min(30, max(18, lon_span * 0.35))
    proj = f'M{fig_width:.0f}c'
    print(f"  区域: [{reg_lon_min:.0f} {reg_lon_max:.0f} {reg_lat_min:.0f} {reg_lat_max:.0f}]")
    
    fig = pygmt.Figure()
    
    # 地形底图 (使用01d分辨率，避免下载超时)
    try:
        fig.grdimage(grid='@earth_relief_01d', region=region, projection=proj,
                     shading='+a45+nt0.5', cmap='geo')
    except:
        # 如果地形数据下载失败，使用简单的海岸线底图
        fig.coast(region=region, projection=proj, land='lightgray', water='lightblue',
                  shorelines='0.3p,gray30')
    
    fig.coast(shorelines='0.3p,gray30',
              borders=['1/0.3p,gray40', '2/0.2p,gray60'],
              area_thresh=500, resolution='f', water='lightblue@15')
    
    frame_lon = max(5, int(lon_span / 8))
    frame_lat = max(2, int((reg_lat_max - reg_lat_min) / 8))
    fig.basemap(frame=[f'xa{frame_lon}f{frame_lon//2}g{frame_lon}',
                        f'ya{frame_lat}f{frame_lat//2}g{frame_lat}',
                        f'WSne+t{EVENT_DATE} Mw {EVENT_MW}  |  T-P arrival time'])
    
    # 方位角弧线 (每 10°)
    from geographiclib.geodesic import Geodesic
    geod = Geodesic.WGS84
    for azi in range(0, 360, 10):
        clons, clats = [], []
        for d_deg in np.arange(0, 80, 0.2):
            g = geod.Direct(EVLA, EVLO, azi, d_deg * 111.12 * 1000)
            lo, la = g['lon2'], g['lat2']
            if region[0]-2 <= lo <= region[1]+2 and region[2]-2 <= la <= region[3]+2:
                clons.append(lo); clats.append(la)
            else:
                break
        if len(clons) > 3:
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tf:
                for lo, la in zip(clons, clats): tf.write(f"{lo} {la}\n")
                tn = tf.name
            fig.plot(data=tn, pen='0.3p,gray30,--')
            os.unlink(tn)
    
    # 台站 — T-P走时着色
    pygmt.makecpt(cmap='hot', series=[t_min, t_max, 1], reverse=True)
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tf:
        for lo, la, tt in zip(lons, lats, t_tp_min):
            tf.write(f"{lo} {la} {tt}\n")
        tn = tf.name
    fig.plot(data=tn, style='t0.22c', cmap=True, pen='0.15p,black')
    os.unlink(tn)
    
    # 台站名标注 — 只标 SNR 最高的 TOP_N 个，防重叠
    TOP_N = 3
    labeled_positions = []
    for r in sorted(rows, key=lambda x: x.get('snr', 0), reverse=True):
        if len(labeled_positions) >= TOP_N:
            break
        # 防重叠：间距 < 0.8° 只标 SNR 最高的
        too_close = False
        for lx, ly in labeled_positions:
            if abs(r['lon'] - lx) < 0.8 and abs(r['lat'] - ly) < 0.8:
                too_close = True
                break
        if too_close:
            continue
        labeled_positions.append((r['lon'], r['lat']))
        fig.text(x=r['lon']+0.3, y=r['lat']+0.15, text=r['sta'],
                 font='4.5p,Helvetica,black', fill='white@60', pen='0')
    
    # 主要城市标注
    cities = {
        'Perth':(115.86,-31.95), 'Adelaide':(138.60,-34.93),
        'Melbourne':(144.96,-37.81), 'Sydney':(151.21,-33.87),
        'Brisbane':(153.03,-27.47), 'Darwin':(130.84,-12.46),
        'Hobart':(147.33,-42.88),
    }
    for name, (clon, clat) in cities.items():
        if region[0] <= clon <= region[1] and region[2] <= clat <= region[3]:
            fig.plot(x=[clon], y=[clat], style='s0.18c', fill='red', pen='0.3p,black')
            fig.text(x=clon+0.5, y=clat+0.4, text=name, font='7p,Helvetica,black',
                     fill='white@70', pen='0')
    
    # 震源机制 + 标注
    fig.meca(spec={'strike1': 300, 'dip1': 90, 'rake1': 180,
                   'strike2': 30, 'dip2': 90, 'rake2': 0,
                   'mantissa': 3.97, 'exponent': 25},
             longitude=EVLO, latitude=EVLA, depth=EVDP,
             scale='0.45c', compression_fill='red', extension_fill='white',
             outline='0.5p,black')
    fig.text(x=EVLO+1.5, y=EVLA+1.0, text=f'Mw {EVENT_MW} ({EVDP} km)',
             font='9p,Helvetica-Bold,red', fill='white@70', pen='0')
    
    # 色标 + 比例尺
    fig.colorbar(frame='xa2f1+l"T-P arrival (min)"',
                 position='jBL+w8c/0.3c+h+o1c/0.8c+ml')
    fig.basemap(map_scale='jBR+o1c/1c+w500k')
    
    # png_path = f"{OUT_DIR}/figure_map_tp.png"
    # fig.savefig(png_path, dpi=300)
    # print(f"✓ {png_path}")
    
    # 出版级地图 (15s DPI450, 台站名+城市+放大图)
    plot_script = os.path.join(CODE_ROOT, "scripts", "plot_station_map.py")
    if os.path.exists(plot_script):
        import subprocess
        subprocess.run([sys.executable, plot_script, args.event], check=False)
    
    import shutil
    shutil.copy(__file__, f"{OUT_DIR}/step1_station_table_and_map.py")
    
    # ── 总结 ──
    print(f"\n{'='*70}")
    print(f"步骤1 完成 — 事件: {args.event}")
    print(f"{'='*70}")
    print(f"输出文件:")
    print(f"  1. stations.xlsx  — 格式化Excel表格 ({len(rows)}站)")
    print(f"  2. stations.csv   — 简洁CSV ({len(rows)}站)")
    print(f"  3. figure_map.png — 出版级全域台站分布图")
    print(f"\n统计:")
    print(f"  总台站: {len(rows)}")
    print(f"  SNR≥{SNR_THRESHOLD}dB: {snr_ok_count}站")
    print(f"  T-P走时: {t_min:.0f}–{t_max:.0f} min")


if __name__ == '__main__':
    main()
