#!/usr/bin/env python3
"""
步骤1: 台站信息提取 + 到时计算 + 台站分布图 (事件无关版)

从 config/events.yaml 读取事件参数，通过 argparse 指定事件名称。

输出:
- stations.xlsx（格式化台站表）
- stations.csv（后续分析输入）
- 固定5张步1图件：2张台站图 + 3张剖面图

图件必须由事件目录中的 `generate_step1_figures.py` 驱动；旧版内置 T-P 着色地图和单一剖面图逻辑已废弃。

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
from matplotlib.lines import Line2D as MLine
from lib.arrivals import compute_theoretical_arrivals, dist_az
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
import pygmt

# ── 默认物理常数 (可通过 config 覆盖) ──
DEFAULT_V_OCEAN = 1.48   # km/s (SOFAR通道)
DEFAULT_V_LOVE = 3.6     # km/s (Love波)
DEFAULT_SNR_THRESHOLD = 8.0  # dB

CONFIG_PATH = os.path.join(CODE_ROOT, "config", "events.yaml")

# ── 默认排除台站（毛刺/噪声台站）──
DEFAULT_EXCLUDE = {
    # XK 台网 — 规律脉冲毛刺 (2013 SW Indian Ridge)
    'W08KB', 'B04KH', 'W07CR', 'B06OR', 'B15MW', 'B07DX',
    # 数据质量差: TEZI(PAZ灵敏度错误) LONA(DC漂移+饱和)
    'TEZI', 'LONA',
    # 通用噪声测试台站
    'AUJCS', 'AUSMG', 'BN2S', 'LCRK', 'STKA',
}
# 10°S 北部密集噪声簇
EXCLUDE_BOXES = [[32.0, 37.5, -11.5, -7.0]]

def _is_excluded(lon, lat):
    for lo, hi, la, ha in EXCLUDE_BOXES:
        if lo <= lon <= hi and la <= lat <= ha:
            return True
    return False

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


def compute_snr_welch(data, delta, tt_p, tt_tp, tt_ts, noise_width=300):
    """
    使用Welch PSD计算T波信噪比(dB)
    信号窗: [tt_tp - 20, tt_ts + 100]
    噪声窗: 动态取 min(noise_width, 实际可用)（P波之前）
    """
    from scipy.signal import welch
    sr = 1.0 / delta
    t = np.arange(len(data)) * delta

    # 信号窗: T-P−20 s 至 T-S+100 s
    s_start = int(max(0, tt_tp - 20) / delta)
    s_end = int(min(t[-1], tt_ts + 100) / delta)

    # 噪声窗: 动态调整，不少于 20s
    actual_noise = max(20, min(noise_width, tt_p * 0.8))
    n_start = int(max(0, tt_p - actual_noise) / delta)
    n_end = int(max(0, tt_p) / delta)

    if s_end <= s_start or n_end <= n_start:
        return 0.0

    f_sig, psd_sig = welch(data[s_start:s_end].astype(np.float64), sr, nperseg=min(256, s_end - s_start))
    f_noise, psd_noise = welch(data[n_start:n_end].astype(np.float64), sr, nperseg=min(256, n_end - n_start))

    # 2-8 Hz频段能量
    e_sig = np.sum(psd_sig[(f_sig >= 2) & (f_sig <= 8)])
    e_noise = np.sum(psd_noise[(f_noise >= 2) & (f_noise <= 8)])

    if e_noise > 0 and e_sig > 0:
        ratio = e_sig / e_noise
        if ratio > 0 and np.isfinite(ratio):
            return 10 * np.log10(ratio)
    return 0.0  # 回退为 0（非 NaN）


def _gen_waveform_profile(rows, evla, evlo, evdp, evdate, evmw, proc_dir, out_dir, snr_threshold,
                         lon_range=None, profile_label=''):
    """波形剖面图 — 按陆地距离排序，转换点对齐"""
    from scipy.signal import hilbert as sp_hilbert
    from scipy.ndimage import uniform_filter1d as uni_filt

    resp_dir = os.path.join(proc_dir, "remove_response")
    files = sorted(glob.glob(f"{proc_dir}/*.BHZ*_proc.sac"))
    if not files:
        files = sorted(glob.glob(f"{resp_dir}/*BHZ*_resp.sac"))
    
    if not files:
        print("  ⚠ 无波形文件，跳过剖面图")
        return

    sta_data = {r['sta']: r for r in rows if r['snr'] >= snr_threshold}
    # 按经度范围过滤
    if lon_range:
        lo, hi = lon_range
        sta_data = {s: d for s, d in sta_data.items() if lo <= d['lon'] <= hi}
    traces = []
    seen = set()
    
    for fn in files:
        bname = os.path.basename(fn)
        parts = bname.split('.')
        net = parts[0]
        if net == 'GT' and len(parts) >= 3:
            sta = parts[1]
        elif len(parts) >= 2:
            sta = parts[1]
        else:
            continue
        if sta not in sta_data or sta in seen:
            continue
        seen.add(sta)
        
        try:
            tr = read(fn)[0]
            data = tr.data.astype(np.float64)
            delta = tr.stats.delta
            p98 = np.percentile(np.abs(data), 98)
            if p98 > 1e6 or not np.isfinite(p98):
                continue
            info = sta_data[sta]
            sac_b = getattr(tr.stats.sac, 'b', 0)
            t = sac_b + np.arange(len(data)) * delta  # 对齐到发震时刻
            env = obs_env(data)
            env = uni_filt(env, size=max(int(2/delta), 2))
            traces.append({
                'sta': sta, 'land_km': info['land_km'],
                'data': data, 'env': env, 't': t, 'delta': delta, 'sac_b': sac_b,
                'tt_p': info.get('tt_p', 0), 'tt_s': info.get('tt_s', 0),
                'tt_tp': info['tt_tp'], 'tt_ts': info['tt_ts'],
                'snr': info['snr'],
            })
        except:
            continue

    if len(traces) < 3:
        print(f"  ⚠ 波形不足 ({len(traces)}站)，跳过")
        return

    traces.sort(key=lambda x: x['land_km'])
    n = len(traces)
    print(f"  加载 {n} 站波形")

    t_min = 0  # 从发震时刻开始
    t_max = max(t['tt_ts'] for t in traces) + 150

    all_p = np.concatenate([t['data'] for t in traces])
    global_p98 = np.percentile(np.abs(all_p), 98) or 1

    fig_h = min(30, max(10, (traces[-1]['land_km'] - traces[0]['land_km']) / 30))
    land_min = traces[0]['land_km']; land_max = traces[-1]['land_km']
    margin = (land_max - land_min) * 0.05 + 20
    amp = 30

    from matplotlib.ticker import MultipleLocator as MLoc
    from numpy.polynomial import Polynomial
    import matplotlib.pyplot as plt

    for mode, tag, suffix in [('envelope', 'Envelope', ''), ('waveform', 'Bandpass [2-8 Hz]', '_waveform')]:
        fig, ax = plt.subplots(figsize=(16, fig_h * 0.7))
        fig.subplots_adjust(left=0.02)

        for tr in traces:
            y0 = tr['land_km']
            m = (tr['t'] >= t_min) & (tr['t'] <= t_max)
            tw = tr['t'][m]
            if len(tw) == 0: continue

            if mode == 'envelope':
                ey = (tr['env'][m] / global_p98) * amp
                ax.fill_between(tw, y0, y0+ey, color='black', alpha=0.15, lw=0, zorder=2)
                ax.plot(tw, y0+ey, color='black', lw=0.6, alpha=0.9, zorder=4)
            else:
                dy = (tr['data'][m] / global_p98) * amp
                ax.plot(tw, y0+dy, color='black', lw=0.2, alpha=0.7, zorder=4)

            # 到时标记: P(绿) S(橙) T-P(蓝) T-S(紫)
            if mode == 'waveform':
                for tt, color, mrk in [(tr['tt_p'], '#4daf4a', 'P'),
                                        (tr['tt_s'], '#ff9848', 'S'),
                                        (tr['tt_tp'], '#4d9fff', 'T-P'),
                                        (tr['tt_ts'], '#ff4d4d', 'T-S')]:
                    if tt and tt > 0 and t_min <= tt <= t_max:
                        idx = int((tt - tr['sac_b']) / tr['delta'])
                        if 0 <= idx < len(tr['data']):
                            v = (tr['data'][idx] / global_p98) * amp
                            ax.plot(tt, y0+v, 'o', color=color, ms=3, markeredgewidth=0.5,
                                    markeredgecolor='black', zorder=10)

        for pts, color, label in [
            ([(t['land_km'], t['tt_tp']) for t in traces if t['tt_tp'] > 0], 'black', 'T-P'),
            ([(t['land_km'], t['tt_ts']) for t in traces if t['tt_ts'] > 0], 'gray', 'T-S')]:
            if pts:
                pts.sort()
                x, y = [p[1] for p in pts], [p[0] for p in pts]
                if len(x) > 1:
                    p = Polynomial.fit(x, y, 1)
                    xl = np.linspace(max(t_min, 0), min(t_max, 400), 100)
                    yl = p(xl)
                    m_ok = yl <= land_max + margin
                    ax.plot(xl[m_ok], yl[m_ok], color=color, lw=1.5, alpha=0.8,
                            dashes=(10, 6), zorder=27)

        ax.set_xlim(t_min, t_max)
        ax.set_ylim(0, land_max + margin)
        ax.xaxis.set_major_locator(MLoc(200)); ax.xaxis.set_minor_locator(MLoc(50))
        ax.set_xlabel('Time from origin (s)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Land distance (km)', fontsize=14, fontweight='bold')
        ax.yaxis.set_major_locator(MLoc(200)); ax.yaxis.set_minor_locator(MLoc(100))
        ax.tick_params(labelsize=12)
        ax.set_title(f'{evdate} Mw {evmw}  |  {n} stations{profile_label}  |  SNR>{snr_threshold}  |  {tag}',
                     fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.08, linestyle=':')
        leg = [
            MLine([0],[0], color='black', lw=0.8, label=tag),
            MLine([0],[0], marker='o', color='#4daf4a', lw=0, ms=4, label='P'),
            MLine([0],[0], marker='o', color='#ff9848', lw=0, ms=4, label='S'),
            MLine([0],[0], marker='o', color='#4d9fff', lw=0, ms=4, label='T-P'),
            MLine([0],[0], marker='o', color='#ff4d4d', lw=0, ms=4, label='T-S'),
        ]
        ax.legend(handles=leg, fontsize=10, loc='upper left', framealpha=0.9, ncol=3)
        fp = os.path.join(out_dir, f'profile{suffix}{profile_label}.png')
        fig.savefig(fp, dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f"  ✓ {fp}")


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
    if not bhz_files:
        bhz_files = sorted(glob.glob(f"{PROC_DIR}/remove_response/*BHZ*_resp.sac"))
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
            if sta in seen_sta or sta in DEFAULT_EXCLUDE:
                continue
            if _is_excluded(stlo, stla):
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
                
                # SNR计算（信号窗 T-P−20 ~ T-S+100，噪声在P波之前）
                snr = compute_snr_welch(tr.data, tr.stats.delta,
                                        tt_p=tt_p, tt_tp=tt_tp, tt_ts=tt_ts)
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
    info_cell.font = Font(name='微软雅黑', bold=True, size=16)
    info_cell.alignment = Alignment(horizontal='center', vertical='center')
    
    # 表头（第2行：中文）
    hdrs_cn = ['台站', '方位角', '震中距', '陆地距', '海洋距',
               '信噪比', 'P波', 'S波', 'T-P', 'T-S', '经度', '纬度', '反方位角']
    for ci, h in enumerate(hdrs_cn, 1):
        c = ws.cell(row=2, column=ci, value=h)
        c.font = Font(name='微软雅黑', bold=True, size=14, color='000000')
        c.alignment = Alignment(horizontal='center', vertical='center')
    
    # 表头（第3行：英文/单位）
    hdrs_en = ['STA', 'AZ(°)', 'DIST(km)', 'LAND(km)', 'OCEAN(km)',
               'SNR(dB)', 'P(s)', 'S(s)', 'T-P(s)', 'T-S(s)', 'LON', 'LAT', 'BAZ(°)']
    
    for ci, h in enumerate(hdrs_en, 1):
        c = ws.cell(row=3, column=ci, value=h)
        c.font = Font(name='微软雅黑', bold=True, size=14, color='FFFFFF')
        c.fill = PatternFill('solid', fgColor='4472C4')
        c.alignment = Alignment(horizontal='center', vertical='center')
    
    # 数据行（从第4行开始）
    for ri, r in enumerate(rows, 4):
        vals = [r['sta'], r['az'], r['dist'], r['land_km'], r['ocean_km'],
                r['snr'], r['tt_p'], r['tt_s'], r['tt_tp'], r['tt_ts'],
                r['lon'], r['lat'], r['baz']]
        
        for ci, v in enumerate(vals, 1):
            cell = ws.cell(row=ri, column=ci, value=v)
            cell.alignment = Alignment(horizontal='center', vertical='center')
            cell.font = Font(name='微软雅黑', size=14)
            
            # SNR高亮（>=阈值 浅绿色）
            if ci == 6:  # SNR列
                if isinstance(v, (int, float)) and v >= SNR_THRESHOLD:
                    cell.fill = PatternFill('solid', fgColor='C6EFCE')
    
    # 列宽
    widths = [12, 10, 12, 12, 12, 12, 10, 10, 12, 12, 14, 14, 12]
    for i, w in enumerate(widths, 1):
        col_letter = chr(64 + i) if i <= 26 else 'A' + chr(64 + i - 26)
        ws.column_dimensions[col_letter].width = w
    
    # 全框细线（黑色）
    thin_border = Border(
        left=Side(style='thin', color='000000'),
        right=Side(style='thin', color='000000'),
        top=Side(style='thin', color='000000'),
        bottom=Side(style='thin', color='000000')
    )
    for row in ws.iter_rows(min_row=1, max_row=len(rows) + 3, max_col=len(hdrs_en)):
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

    # ── 4. 固定五图输出：只调用事件目录独立驱动 ──
    # 新事件须在 OUT_DIR 放置 generate_step1_figures.py；其中定义两张台站图和三张剖面图。
    import subprocess
    import shutil
    dst_step1 = os.path.join(OUT_DIR, 'step1_station_table_and_map.py')
    if os.path.abspath(__file__) != os.path.abspath(dst_step1):
        shutil.copy2(__file__, dst_step1)
    figure_driver = os.path.join(OUT_DIR, 'generate_step1_figures.py')
    if not os.path.isfile(figure_driver):
        raise FileNotFoundError(f'缺少步1五图驱动脚本: {figure_driver}')
    subprocess.run([sys.executable, figure_driver], check=True)
    print("\n步骤1完成：表格 + CSV + 2张台站图 + 3张剖面图")
    return

    # ── 已废弃：旧版内置 T-P 着色地图/剖面图逻辑（不再执行） ──
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
    
    # 台站名标注 (小字半透底)
    for r in rows:
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
    
    # 期刊级地图 (高精度地形, 大图标)
    plot_script2 = os.path.join(CODE_ROOT, "scripts", "plot_station_map_final.py")
    if os.path.exists(plot_script2):
        subprocess.run([sys.executable, plot_script2, args.event], check=False)
    
    import shutil
    shutil.copy(__file__, f"{OUT_DIR}/step1_station_table_and_map.py")
    
    # ── 波形剖面图 ──
    print("\n生成波形剖面图...")
    _gen_waveform_profile(rows, EVLA, EVLO, EVDP, EVENT_DATE, EVENT_MW,
                          PROC_DIR, OUT_DIR, SNR_THRESHOLD,
                          lon_range=(43, 60), profile_label='_Madagascar')
    _gen_waveform_profile(rows, EVLA, EVLO, EVDP, EVENT_DATE, EVENT_MW,
                          PROC_DIR, OUT_DIR, SNR_THRESHOLD,
                          lon_range=(10, 43), profile_label='_Mainland')
    
    # ── 总结 ──
    print(f"\n{'='*70}")
    print(f"步骤1 完成 — 事件: {args.event}")
    print(f"{'='*70}")
    print(f"输出文件:")
    print(f"  1. stations.xlsx  — 格式化Excel表格 ({len(rows)}站)")
    print(f"  2. stations.csv   — 简洁CSV ({len(rows)}站)")
    print(f"  3. figure_map.png — 台站分布图 (全站名标注)")
    print(f"  4. figure_map_final.png — 期刊级台站图 (高精度地形)")
    print(f"  5. profile.png / profile_waveform.png — 波形剖面图")
    print(f"\n统计:")
    print(f"  总台站: {len(rows)}")
    print(f"  SNR≥{SNR_THRESHOLD}dB: {snr_ok_count}站")
    print(f"  T-P走时: {t_min:.0f}–{t_max:.0f} min")


if __name__ == '__main__':
    main()
