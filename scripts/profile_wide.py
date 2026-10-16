#!/usr/bin/env python3
"""
波形剖面图 — 加宽+变矮+全时段+SNR筛选 + 到时标注优化
支持 --y-land 按陆地距离纵轴排列（允许波形重叠）

期刊出图规范:
  font.size=14 font.weight='bold' axes.labelsize=15 axes.titlesize=16
  xtick/ytick.labelsize=13 legend.fontsize=12  配色黑/灰/白
"""
import sys, os, glob, argparse, yaml
sys.path.insert(0, '/home/sy/codex/t_phase_analysis')

import numpy as np
from obspy.signal.filter import envelope
from obspy.signal.util import smooth
from scipy.ndimage import uniform_filter1d
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D as MLine
from obspy import read
from obspy.taup import TauPyModel
from obspy.signal.filter import envelope as obs_env
from geographiclib.geodesic import Geodesic
from t_phase_analysis.traveltime_profile import get_track, tt_t

# 期刊出图 rcParams
plt.rcParams.update({
    'font.size': 14, 'font.weight': 'bold',
    'axes.labelsize': 15, 'axes.titlesize': 16,
    'xtick.labelsize': 13, 'ytick.labelsize': 13,
    'legend.fontsize': 12,
})

CONFIG_PATH = os.path.join(os.path.dirname(__file__), '..', 'config', 'events.yaml')

def main():
    p = argparse.ArgumentParser(
        description="波形剖面图 — 事件无关，从 events.yaml 读取参数")
    p.add_argument('event', type=str,
                   help='事件名称 (events.yaml 中的 key，如 puerto_rico)')
    p.add_argument('--az-min', type=float, default=None,
                   help='方位角下限 (默认从 config 读取)')
    p.add_argument('--az-max', type=float, default=None,
                   help='方位角上限 (默认从 config 读取)')
    p.add_argument("--snr-min", type=float, default=2.0)
    p.add_argument('--t-min', type=float, default=None,
                   help='时间窗口起点 (默认: min(T-P) - 100s)')
    p.add_argument('--t-max', type=float, default=None,
                   help='时间窗口终点 (默认: max(T-S) + 200s)')
    p.add_argument('--amp', type=float, default=0.10)
    p.add_argument('--spacing', type=float, default=0.35)
    p.add_argument('--sort-by', choices=['land', 'dist'], default='land',
                   help='排序方式: land=陆地距, dist=震中距')
    p.add_argument('--align-tp', action='store_true', default=False,
                   help='海洋走时对齐（t=0 = 波上岸）')
    p.add_argument('--exclude', type=str, default='',
                   help='逗号分隔的排除台站列表')
    p.add_argument('--y-land', action='store_true', default=False,
                   help='纵轴按陆地距离 (km) 排列，允许波形重叠')
    args = p.parse_args()

    # ── 从 events.yaml 加载事件参数 ──
    with open(CONFIG_PATH, 'r') as f:
        events = yaml.safe_load(f)
    if args.event not in events:
        print(f"ERROR: 事件 '{args.event}' 未在 {CONFIG_PATH} 中找到")
        print(f"可用事件: {', '.join(events.keys())}")
        sys.exit(1)
    cfg = events[args.event]

    EVLA = cfg['lat']
    EVLO = cfg['lon']
    EVDP = cfg['depth_km']
    case_dir = cfg['case_dir']
    PROC_DIR = f"/mnt/e/sharedate/case/{case_dir}/sacprocess"
    OUT_DIR  = f"/mnt/e/sharedate/case/{case_dir}/process_flow/2_waveform_profile"

    # 方位角默认从 config 读取；若为空则覆盖全范围
    if args.az_min is None or args.az_max is None:
        az_profiles = cfg.get('azimuth_profiles') or [180]
        az_center = az_profiles[0]
        az_half_width = cfg.get('azimuth_half_width', 15)
        if args.az_min is None:
            args.az_min = max(0, az_center - az_half_width)
        if args.az_max is None:
            args.az_max = min(360, az_center + az_half_width)

    # 使用 config 中的 snr_threshold 作为 --snr-min 的默认值 (仅当未显式指定时)
    if '--snr-min' not in sys.argv:
        args.snr_min = cfg.get('snr_threshold', args.snr_min)

    print(f"事件: {args.event} ({cfg.get('label', args.event)})")
    print(f"  震中: {EVLA:.4f}, {EVLO:.4f}  深度: {EVDP} km")
    print(f"  PROC_DIR: {PROC_DIR}")
    print(f"  OUT_DIR:  {OUT_DIR}")

    taup = TauPyModel(model='iasp91')
    os.makedirs(OUT_DIR, exist_ok=True)

    # ── 收集台站 ──
    files = sorted(glob.glob(f"{PROC_DIR}/*.BHZ*_proc.sac"))
    stations = []; seen = set()
    for fn in files:
        try:
            sta = os.path.basename(fn).split('.')[1]
            if sta in seen: continue; seen.add(sta)
            tr = read(fn, headonly=True)[0]
            stla, stlo = tr.stats.sac.stla, tr.stats.sac.stlo
            g = Geodesic.WGS84.Inverse(EVLA, EVLO, stla, stlo)
            az = g['azi1']; az = az+360 if az<0 else az
            if args.az_min <= az <= args.az_max:
                stations.append({'sta':sta, 'az':az, 'stla':stla, 'stlo':stlo, 'fn':fn})
        except: continue

    print(f"{len(stations)} stations in AZ {args.az_min}°–{args.az_max}°")
    
    if args.exclude:
        exclude_set = set(s.strip() for s in args.exclude.split(',') if s.strip())
        before = len(stations)
        stations = [s for s in stations if s['sta'] not in exclude_set]
        if before != len(stations):
            print(f"  剔除: {before - len(stations)} stations ({', '.join(exclude_set)})")

    # ── 加载波形+计算到时+SNR ──
    traces = []
    seen_traces = set()
    for s in stations:
        # 二次去重保护
        if s['sta'] in seen_traces: continue
        seen_traces.add(s['sta'])
        
        tr = read(s['fn'])[0]; data = tr.data.astype(np.float64)
        delta = tr.stats.delta; t_off = tr.stats.sac.b - tr.stats.sac.o
        t = np.arange(len(data))*delta + t_off
        env = obs_env(data)
        smooth_win = max(1, int(3.0 * sr))  # 3秒窗口 (tpef方法)
        env = smooth(env, smooth_win)

        from obspy.geodetics.base import locations2degrees
        dd = locations2degrees(EVLA, EVLO, s['stla'], s['stlo'])
        arrs = taup.get_travel_times(source_depth_in_km=EVDP, distance_in_degree=dd, phase_list=['P','S'])
        tt_P = next((a.time for a in arrs if a.name=='P'), None)
        tt_S = next((a.time for a in arrs if a.name=='S'), None)

        from obspy.geodetics.base import degrees2kilometers
        dist_km = degrees2kilometers(dd)
        
        land_km = 0; tt_TP = tt_TS = tt_ocean = None
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            try:
                track = get_track(EVLO, EVLA, s['stlo'], s['stla'], 0.02, ref_el=-1000, plot_map=False, plot_elevation=False)
                r = tt_t(track, ref_el=-1000, source_depth_km=EVDP, ignore_small_islands=True)
                tt_ocean = r['tt_ocean']
                tt_TP = r['tt_t']; tt_TS = tt_ocean + r['tt_land_s']
                land_km = r['dist_land'] * 111.0
            except: pass

        # SNR计算 — Welch PSD 频域功率谱法 (2-8 Hz)
        snr = -999
        if tt_TP and tt_TS:
            from scipy.signal import welch
            sr = 1.0 / delta
            # 信号窗: [tt_TP, tt_TS + 60] — 用 t 数组定位索引
            s_mask = (t >= tt_TP) & (t <= tt_TS + 60)
            # 噪声窗: [tt_TP - 300, tt_TP - 60]
            n_mask = (t >= tt_TP - 300) & (t <= tt_TP - 60)
            
            if np.any(s_mask) and np.any(n_mask):
                f_sig, psd_sig = welch(data[s_mask], fs=sr, nperseg=min(256, np.sum(s_mask)))
                f_noise, psd_noise = welch(data[n_mask], fs=sr, nperseg=min(256, np.sum(n_mask)))
                e_sig = np.sum(psd_sig[(f_sig >= 2) & (f_sig <= 8)])
                e_noise = np.sum(psd_noise[(f_noise >= 2) & (f_noise <= 8)])
                if e_noise > 0 and e_sig > 0:
                    snr = round(10 * np.log10(e_sig / e_noise), 2)

        traces.append(dict(sta=s['sta'], az=s['az'], land_km=land_km, dist_km=dist_km,
            data=data, env=env, t=t, delta=delta,
            P=tt_P, S=tt_S, TP=tt_TP, TS=tt_TS, snr=snr, ocean=tt_ocean))

    # ── 海洋走时对齐 ──
    if args.align_tp:
        for tr in traces:
            if tr['ocean'] is not None:
                tr['t'] = tr['t'] - tr['ocean']
                if tr['P'] is not None: tr['P'] -= tr['ocean']
                if tr['S'] is not None: tr['S'] -= tr['ocean']
                if tr['TP'] is not None: tr['TP'] -= tr['ocean']
                if tr['TS'] is not None: tr['TS'] -= tr['ocean']
        print(f"  Aligned to ocean arrival (t=0 = conversion point)")

    # ── SNR筛选+排序 ──
    traces = [t for t in traces if t['snr'] >= args.snr_min]
    traces.sort(key=lambda x: x['land_km'])
    sort_label = 'Land distance'
    n = len(traces)
    print(f"  筛选后: {n} stations (SNR ≥ {args.snr_min})")
    if n == 0: print("无台站"); return
    
    # ── 自动计算时间窗口 ──
    if args.t_min is None or args.t_max is None:
        tp_list = [t['TP'] for t in traces if t['TP'] is not None]
        ts_list = [t['TS'] for t in traces if t['TS'] is not None]
        if tp_list and ts_list:
            if args.t_min is None:
                args.t_min = min(tp_list) - 100
                print(f"  自动时间窗口: T-P(min)-100s = {args.t_min:.0f}s")
            if args.t_max is None:
                args.t_max = max(ts_list) + 200
                print(f"  自动时间窗口: T-S(max)+200s = {args.t_max:.0f}s")
        else:
            if args.t_min is None: args.t_min = 0
            if args.t_max is None: args.t_max = 3600
    
    # ── 振幅裁剪 ──
    if len(traces) > 0:
        p98s = [np.percentile(np.abs(t['data'][(t['t'] >= args.t_min) & (t['t'] <= args.t_max)]), 98)
                for t in traces]
        median_p98 = np.median(p98s)
        clip_ceil = median_p98 * 3.0
        for i, t in enumerate(traces):
            if p98s[i] > clip_ceil:
                scale = clip_ceil / p98s[i]
                traces[i]['data'] = t['data'] * scale
                traces[i]['env'] = t['env'] * scale
                print(f"  Clipped {t['sta']}: {p98s[i]/median_p98:.1f}x → 3x")
    
    # ── 绘图 ──
    # 纵轴 = 陆地距离 (km)
    min_land = traces[0]['land_km']
    max_land = traces[-1]['land_km']
    margin = (max_land - min_land) * 0.05 + 20
    fig_h = min(20, max(6, (max_land - min_land) / 50))
    fig, ax = plt.subplots(figsize=(24, fig_h))
    amp_scale = args.amp * 25  # 放大振幅以适应 km 纵轴
    
    all_plot = np.concatenate([t['data'][(t['t'] >= args.t_min) & (t['t'] <= args.t_max)] 
                               for t in traces])
    global_p98 = np.percentile(np.abs(all_plot), 98)
    print(f"  Global P98: {global_p98:.2e}")

    for tr in traces:
        y0 = tr['land_km']  # 统一用陆地距离定位
        
        m = (tr['t'] >= args.t_min) & (tr['t'] <= args.t_max)
        tw = tr['t'][m]; dw = tr['data'][m]; ew = tr['env'][m]
        if len(tw) == 0: continue

        if global_p98 > 0: dw_n = dw / global_p98; ew_n = ew / global_p98
        else: dw_n = dw; ew_n = ew
        dy = dw_n * amp_scale; ey = ew_n * amp_scale

        ax.plot(tw, y0+dy, 'k-', lw=0.5, alpha=0.9, zorder=3)
        ax.fill_between(tw, y0, y0+ey, color='gray', alpha=0.08, lw=0, zorder=2)
        ax.plot(tw, y0+ey, color='gray', lw=0.5, alpha=0.6, zorder=4)

        y_marker = y0 + 0.75 * amp_scale
        for ph, tt, color, msize in [('T-P', tr['TP'], 'black', 7), ('T-S', tr['TS'], 'gray', 7)]:
            if tt is not None and args.t_min <= tt <= args.t_max:
                ax.plot(tt, y_marker, 'o', color=color, markersize=msize,
                        markeredgecolor='white', markeredgewidth=1.2, zorder=30)

    # ── 图表设置 ──
    ax.set_xlim(args.t_min - 80, args.t_max)  # 左侧留空给台站名标签
    ax.set_ylim(min_land - margin, max_land + margin)
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.xaxis.set_minor_locator(MultipleLocator(20))
    xlabel = 'Time from ocean arrival (s)' if args.align_tp else 'Time from origin (s)'
    ax.set_xlabel(xlabel, fontsize=10, fontweight='bold')
    ax.set_ylabel('Land distance (km)', fontsize=10, fontweight='bold')
    ax.yaxis.set_major_locator(MultipleLocator(100))
    ax.yaxis.set_minor_locator(MultipleLocator(50))
    ax.tick_params(labelsize=7)
    # 显示台站名+陆地距离
    for tr in traces:
        ax.text(args.t_min - 10, tr['land_km'], f"{tr['sta']} ({tr['land_km']:.0f}km)",
                fontsize=6, fontfamily='monospace', ha='right', va='center')
    
    ax.set_title(f'AZ {args.az_min:.0f}° – {args.az_max:.0f}°  |  {n} stations  |  SNR ≥ {args.snr_min}',
                 fontsize=11, fontweight='bold')
    ax.grid(True, alpha=0.08, linestyle=':')
    
    leg = [
        MLine([0],[0], color='black', lw=0.8, label='2-8 Hz waveform'),
        MLine([0],[0], color='gray', lw=0.6, alpha=0.4, label='Envelope'),
        MLine([0],[0], marker='o', color='black', linestyle='None', ms=5, 
              markeredgecolor='white', markeredgewidth=0.8, label='T-P'),
        MLine([0],[0], marker='o', color='gray', linestyle='None', ms=5, 
              markeredgecolor='white', markeredgewidth=0.8, label='T-S'),
    ]
    ax.legend(handles=leg, fontsize=18, loc='upper left', ncol=1, framealpha=0.9)

    plt.subplots_adjust(left=0.08, right=0.97, top=0.94, bottom=0.06)
    suffix = '_yland' if args.y_land else '_optimized'
    tag = f"profile_az{args.az_min:.0f}_{args.az_max:.0f}_snr{args.snr_min}_{args.sort_by}{suffix}"
    out_png = f"{OUT_DIR}/{tag}.png"
    fig.savefig(out_png, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"✅ {out_png}")
    
    import shutil
    shutil.copy(__file__, f"{OUT_DIR}/profile_wide.py")
    print(f"✅ script saved")

if __name__ == '__main__':
    main()
