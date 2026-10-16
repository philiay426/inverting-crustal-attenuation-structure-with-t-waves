#!/usr/bin/env python3
"""
06_profile_section.py — 方位角波形剖面图 (支持 reduced-time / envelope-only)
用法:
  python3 06_profile_section.py <event_key> --az-center 310

模式:
  --time-mode origin             横轴 = 发震后时间 (s) [默认]
  --time-mode conversion_reduced  横轴 = t - L_ocean/1.48

  --plot-mode waveform+envelope  波形 + 包络 [默认]
  --plot-mode envelope_only      只画 RMS 包络

  --envelope-method obspy      Hilbert + 2s轻平滑 (保留峰值波动) [默认]
  --envelope-method hilbert      Hilbert + 10s平滑
  --envelope-method rms          15s滑动RMS + 5s平滑

输出: case/<case_dir>/record_section_<模式>_az<center>_<half>deg.png/pdf
"""
import sys, os, glob, argparse, yaml

sys.path.insert(0, '/mnt/f/sharedate/codex/t_phase_analysis')
for p in ['/home/philiay/shared/codex']:
    if p not in sys.path: sys.path.insert(0, p)

import numpy as np
from scipy.signal import hilbert
from scipy.ndimage import uniform_filter1d
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D as MLine
from obspy import read
from obspy.taup import TauPyModel
from geographiclib.geodesic import Geodesic
from t_phase_analysis.traveltime_profile import get_track, tt_t

CONFIG = yaml.safe_load(open('/home/philiay/shared/codex/config/events.yaml'))

# ── 包络计算 ──
def calc_envelope(data, delta, method='obspy'):
    if method == 'obspy':
        # obspy 自带 Hilbert + 2s 轻平滑（保留峰值波动，去载波振荡）
        env_raw = np.abs(hilbert(data))
        env = uniform_filter1d(env_raw, size=int(2 / delta))
    elif method == 'hilbert':
        env_raw = np.abs(hilbert(data))
        env = uniform_filter1d(env_raw, size=int(10 / delta))
    else:  # rms
        # 15s sliding RMS
        rms_win = int(15 / delta)
        rms = np.sqrt(uniform_filter1d(data**2, size=rms_win))
        # 5s moving average 再平滑
        env = uniform_filter1d(rms, size=int(5 / delta))
    return env


def main():
    p = argparse.ArgumentParser(description='方位角波形剖面图')
    p.add_argument('event', help='事件key')
    p.add_argument('--az-center', type=float, required=True)
    p.add_argument('--az-half', type=float, default=1.0)
    p.add_argument('--t-min', type=float, default=-300)
    p.add_argument('--t-max', type=float, default=3600)
    p.add_argument('--amp', type=float, default=0.15)
    p.add_argument('--spacing', type=float, default=0.55)
    p.add_argument('--time-mode', choices=['origin', 'conversion_reduced'],
                   default='origin')
    p.add_argument('--plot-mode', choices=['waveform+envelope', 'envelope_only'],
                   default='waveform+envelope')
    p.add_argument('--envelope-method', choices=['rms', 'hilbert'],
                   default='rms')
    # reduced-time 横轴范围 (仅在 conversion_reduced 模式下生效)
    p.add_argument('--rt-min', type=float, default=-100)
    p.add_argument('--rt-max', type=float, default=None)
    args = p.parse_args()

    cfg = CONFIG[args.event]
    EVLA, EVLO = cfg['lat'], cfg['lon']
    EVDP = cfg['depth_km']
    case_dir = cfg.get('case_dir', cfg['label'])
    PROC_DIR = f"/home/philiay/shared/case/{case_dir}/sacprocess"
    OUT_DIR  = f"/mnt/f/sharedate/case/{case_dir}/map_{case_dir.replace('-','_')}/data"
    os.makedirs(OUT_DIR, exist_ok=True)

    taup = TauPyModel(model='iasp91')

    # ── 收集台站 ──
    files = sorted(glob.glob(f"{PROC_DIR}/*.BHZ*_proc.sac"))
    stations = []
    seen = set()
    for fn in files:
        try:
            tr = read(fn, headonly=True)[0]
            sta = os.path.basename(fn).split('.')[1]
            if sta in seen: continue
            seen.add(sta)
            stla, stlo = tr.stats.sac.stla, tr.stats.sac.stlo
            g = Geodesic.WGS84.Inverse(EVLA, EVLO, stla, stlo)
            az = g['azi1']
            if az < 0: az += 360
            lo = args.az_center - args.az_half
            hi = args.az_center + args.az_half
            if lo <= az <= hi:
                stations.append({'sta': sta, 'az': az,
                                 'stla': stla, 'stlo': stlo, 'fn': fn})
        except:
            continue

    if len(stations) == 0:
        print("❌ 该方位角范围内无台站"); return
    print(f"{len(stations)} stations in AZ {lo:.1f}° – {hi:.1f}°")

    # ── 加载波形 + 计算到时 + 包络 ──
    traces = []
    for s in stations:
        tr = read(s['fn'])[0]
        data = tr.data.astype(np.float64)
        delta = tr.stats.delta
        t_offset = tr.stats.sac.b - tr.stats.sac.o
        npts = len(data)
        t_origin = np.arange(npts) * delta + t_offset  # 发震后时间

        from obspy.geodetics.base import locations2degrees
        dist_deg = locations2degrees(EVLA, EVLO, s['stla'], s['stlo'])
        arrs = taup.get_travel_times(source_depth_in_km=EVDP,
                                      distance_in_degree=dist_deg,
                                      phase_list=['P', 'S'])
        tt_P = tt_S = None
        for a in arrs:
            if a.name == 'P' and tt_P is None: tt_P = a.time
            elif a.name == 'S' and tt_S is None: tt_S = a.time

        # T-P/T-S + 陆地/海洋距离
        land_km = 0; ocean_km = 0
        tt_TP = tt_TS = None
        try:
            track = get_track(EVLO, EVLA, s['stlo'], s['stla'], 0.02,
                             ref_el=-1000, plot_map=False, plot_elevation=False)
            tr_tt = tt_t(track, ref_el=-1000, source_depth_km=EVDP,
                         ignore_small_islands=True)
            tt_TP = tr_tt['tt_t']
            tt_TS = tr_tt['tt_ocean'] + tr_tt['tt_land_s']
            land_km = tr_tt['dist_land'] * 111.0
            ocean_km = tr_tt['dist_ocean'] * 111.0
        except:
            pass

        # 包络
        env = calc_envelope(data, delta, method=args.envelope_method)

        traces.append({
            'sta': s['sta'], 'az': s['az'],
            'land_km': land_km, 'ocean_km': ocean_km,
            't_origin': t_origin, 'data': data, 'env': env,
            'P': tt_P, 'S': tt_S, 'TP': tt_TP, 'TS': tt_TS,
        })

    # 按陆地距离排序
    traces.sort(key=lambda x: x['land_km'])
    land_info = [(t['sta'], f"{t['land_km']:.0f}km") for t in traces]
    print(f"Sorted by land dist: {land_info}")

    # ── 计算绘图时间坐标 ──
    if args.time_mode == 'origin':
        xlabel = 'Time from origin (s)'
        title_suffix = ''
    else:  # conversion_reduced
        xlabel = 'Reduced time: t − L_ocean/1.48 (s)'
        title_suffix = ' (reduced time)'

    for tr in traces:
        if args.time_mode == 'conversion_reduced' and tr['ocean_km'] > 0:
            t_ocean = tr['ocean_km'] / 1.48  # 海洋走时
            tr['t_plot'] = tr['t_origin'] - t_ocean
            # 理论到时在 reduced-time 坐标中
            if tr['TP'] and tr['TS']:
                tr['TP_plot'] = tr['TP'] - t_ocean
                tr['TS_plot'] = tr['TS'] - t_ocean
            else:
                tr['TP_plot'] = tr['TS_plot'] = None
            # 备用：使用陆地速度直接算
            if tr['land_km'] > 0:
                tr['TP_plot_v2'] = tr['land_km'] / 6.0
                tr['TS_plot_v2'] = tr['land_km'] / 3.6
            else:
                tr['TP_plot_v2'] = tr['TS_plot_v2'] = None
        else:
            tr['t_plot'] = tr['t_origin']
            tr['TP_plot'] = tr['TP']
            tr['TS_plot'] = tr['TS']
            tr['TP_plot_v2'] = tr['TS_plot_v2'] = None

    # ── 自动横轴范围 (reduced-time) ──
    if args.time_mode == 'conversion_reduced':
        if args.rt_max is not None:
            xmin, xmax = args.rt_min, args.rt_max
        else:
            all_ts = [tr['TS_plot'] for tr in traces if tr['TS_plot'] is not None]
            xmin = args.rt_min
            xmax = (max(all_ts) + 300) if all_ts else 700
    else:
        xmin, xmax = args.t_min, args.t_max

    # ── 绘图 ──
    fig, ax = plt.subplots(figsize=(14, args.spacing * len(traces) + 1.5))

    for i, tr in enumerate(traces):
        y0 = i
        m = (tr['t_plot'] >= xmin) & (tr['t_plot'] <= xmax)
        tw = tr['t_plot'][m]; dw = tr['data'][m]; ew = tr['env'][m]
        if len(tw) == 0: continue

        # P98 归一化（始终用原始数据幅值，不依赖窗口裁剪）
        p98 = np.percentile(np.abs(tr['data']), 98)
        if p98 > 0:
            dw_norm = dw / p98
            ew_norm = ew / p98
        else:
            dw_norm = dw; ew_norm = ew

        dy = dw_norm * args.amp
        ey = ew_norm * args.amp

        if args.plot_mode == 'waveform+envelope':
            ax.plot(tw, y0 + dy, 'k-', lw=0.4, alpha=0.8)
            ax.plot(tw, y0 + ey, 'r-', lw=0.6, alpha=0.5)
            ax.plot(tw, y0 - ey, 'r-', lw=0.6, alpha=0.5)
            ax.fill_between(tw, y0 - ey, y0 + ey, color='red', alpha=0.06)
        else:  # envelope_only
            ax.plot(tw, y0 + ey, 'r-', lw=0.7, alpha=0.7)
            ax.plot(tw, y0 - ey, 'r-', lw=0.7, alpha=0.7)
            ax.fill_between(tw, y0 - ey, y0 + ey, color='red', alpha=0.1)

        # 台站标签
        label = f'{tr["sta"]} ({tr["land_km"]:.0f}km)'
        if args.time_mode == 'conversion_reduced' and tr['ocean_km'] > 0:
            label += f' | ocean={tr["ocean_km"]:.0f}km'
        ax.text(xmax + (xmax - xmin) * 0.01, y0, label,
                fontsize=7, va='center', ha='left', color='black')

        # 圆点到时 — reduced-time 模式弱化 P/S
        markers = []
        if args.time_mode == 'origin':
            markers = [('P', tr['P'], '#2ca02c', 3),
                       ('S', tr['S'], '#ff7f0e', 3)]
        # T-P/T-S 始终画
        markers += [('T-P', tr['TP_plot'], '#1f77b4', 3.5),
                    ('T-S', tr['TS_plot'], '#d62728', 3.5)]

        for ph, tt, c, ms in markers:
            if tt is not None and xmin <= tt <= xmax:
                idx = np.argmin(np.abs(tr['t_plot'] - tt))
                amp_at = tr['data'][idx]
                na = amp_at / p98 * args.amp if p98 > 0 else 0
                ax.plot(tt, y0 + na, 'o', color=c, markersize=ms, zorder=5)
                ax.plot(tt, y0 - na, 'o', color=c, markersize=ms, zorder=5)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-1, len(traces))
    ax.xaxis.set_major_locator(MultipleLocator(500 if args.time_mode == 'origin' else 100))
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(f'Station (sorted by land distance{title_suffix})', fontsize=10)

    mode_tag = 'env' if args.plot_mode == 'envelope_only' else 'wf'
    ax.set_title(f'AZ {lo:.1f}°–{hi:.1f}°  {mode_tag}  |  '
                 f'{len(traces)} stations  |  '
                 f'norm=P98  env={args.envelope_method}  |  '
                 f'●T-P ●T-S{"" if args.time_mode=="origin" else "  ●P ●S (hidden)"}',
                 fontsize=9, fontweight='bold')
    ax.grid(True, alpha=0.15, linestyle=':')
    ax.tick_params(labelsize=8)
    ax.set_yticks([])

    # 图例
    leg = []
    if args.plot_mode == 'waveform+envelope':
        leg.append(MLine([0], [0], color='black', lw=0.8, label='Waveform'))
    leg.append(MLine([0], [0], color='red', lw=0.8, alpha=0.6, label='Envelope'))
    leg.append(MLine([0], [0], marker='o', color='#1f77b4', linestyle='None', markersize=4, label='T-P'))
    leg.append(MLine([0], [0], marker='o', color='#d62728', linestyle='None', markersize=4, label='T-S'))
    if args.time_mode == 'origin':
        leg.append(MLine([0], [0], marker='o', color='#2ca02c', linestyle='None', markersize=4, label='P'))
        leg.append(MLine([0], [0], marker='o', color='#ff7f0e', linestyle='None', markersize=4, label='S'))
    ax.legend(handles=leg, fontsize=7, loc='upper right', ncol=2)

    plt.tight_layout()
    tag = f"{args.time_mode}_{mode_tag}_az{args.az_center:.0f}_{args.az_half:.0f}deg"
    fig.savefig(f'{OUT_DIR}/record_section_{tag}.png', dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'✅ {OUT_DIR}/record_section_{tag}.png/pdf')

if __name__ == '__main__':
    main()
