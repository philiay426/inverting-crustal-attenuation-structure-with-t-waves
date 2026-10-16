#!/usr/bin/env python3
"""
全时段波形剖面 — P→S→T, 完整下载窗口 [-10, +60] min
从预处理数据 sacprocess/ 和 stations.csv 读取，参数从 events.yaml
支持 --y-land (陆地距离) / --equal (等间距排列)

期刊出图规范:
  font.size=14 font.weight='bold' axes.labelsize=15 axes.titlesize=16
  xtick/ytick.labelsize=13 legend.fontsize=12  配色黑/灰/白

用法:
  python3 profile_full.py <event> --snr-min 3 --y-land
  python3 profile_full.py <event> --snr-min 3 --equal --land-min 200 --land-max 300
"""
import sys, os, glob, csv, yaml, argparse
sys.path.insert(0, '/home/sy/codex/t_phase_analysis')

import numpy as np
from scipy.signal import hilbert
from obspy.signal.filter import envelope
from obspy.signal.util import smooth
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from obspy import read

plt.rcParams.update({
    'font.size': 14, 'font.weight': 'bold',
    'axes.labelsize': 15, 'axes.titlesize': 16,
    'xtick.labelsize': 13, 'ytick.labelsize': 13,
    'legend.fontsize': 12,
})

CONFIG_PATH = os.path.join(os.path.dirname(__file__), '..', 'config', 'events.yaml')


def main():
    p = argparse.ArgumentParser(description='全时段波形剖面 (P–S–T)')
    p.add_argument('event', help='events.yaml 中的事件 key')
    p.add_argument('--snr-min', type=float, default=3.0, help='SNR 阈值 (dB)')
    p.add_argument('--az-min', type=float, default=0, help='方位角下限')
    p.add_argument('--az-max', type=float, default=360, help='方位角上限')
    p.add_argument('--y-land', action='store_true', help='纵轴按陆地距离排列')
    p.add_argument('--equal', action='store_true', help='等间距排列 (用于窄陆地距离段)')
    p.add_argument('--land-min', type=float, default=0, help='陆地距离下限 (--equal 时生效)')
    p.add_argument('--land-max', type=float, default=1e6, help='陆地距离上限 (--equal 时生效)')
    p.add_argument('--exclude', type=str, default='', help='逗号分隔排除台站')
    args = p.parse_args()

    # 读配置
    with open(CONFIG_PATH) as f:
        cfg = yaml.safe_load(f)[args.event]
    evla, evlo = cfg['lat'], cfg['lon']
    evdate, evmw = cfg['date'], cfg['mw']
    case_dir = cfg['case_dir']
    PROC_DIR = f"/mnt/e/sharedate/case/{case_dir}/sacprocess"
    OUT_DIR  = f"/mnt/e/sharedate/case/{case_dir}/process_flow/2_waveform_profile"
    CSV_PATH = f"/mnt/e/sharedate/case/{case_dir}/process_flow/1_station_distribution/stations.csv"
    os.makedirs(OUT_DIR, exist_ok=True)
    exclude_set = set(s.strip() for s in args.exclude.split(',') if s.strip())

    # 读 CSV 到时
    arrivals = {}
    with open(CSV_PATH) as f:
        reader = csv.DictReader((l for l in f if not l.startswith('#')))
        for row in reader:
            sta = row['sta']
            land = float(row['land_km'])
            snr = float(row['snr'])
            dist = float(row['dist'])
            if snr < args.snr_min: continue
            if land < args.land_min or land > args.land_max: continue
            if sta in exclude_set: continue
            arrivals[sta] = dict(
                tt_p=float(row['tt_p']), tt_s=float(row['tt_s']),
                tt_tp=float(row['tt_tp']), tt_ts=float(row['tt_ts']),
                land_km=land, dist_km=dist, snr=snr, az=float(row['az']),
            )

    # 收集台站文件
    files = sorted(glob.glob(f"{PROC_DIR}/*.BHZ*_proc.sac"))
    stations = []; seen = set()
    for fn in files:
        try:
            sta = os.path.basename(fn).split('.')[1]
            if sta in seen: continue; seen.add(sta)
            if sta not in arrivals: continue
            if not (args.az_min <= arrivals[sta]['az'] <= args.az_max): continue
            tr = read(fn, headonly=True)[0]
            stations.append(dict(sta=sta, fn=fn, o=tr.stats.sac.o, **arrivals[sta]))
        except: continue

    stations.sort(key=lambda x: x['dist_km'])
    print(f"{len(stations)} stations (SNR≥{args.snr_min}, AZ {args.az_min}–{args.az_max})")

    # 时间窗口: 发震时刻到 20 min 后
    t_min, t_max = 0, 1200

    # 加载波形 + 全局 P98
    loaded = []
    for s in stations:
        try: loaded.append((read(s['fn'])[0], s))
        except: continue

    all_amps = []
    for tr, s in loaded:
        t_rel = tr.times() - s['o']
        m = (t_rel >= t_min) & (t_rel <= t_max)
        if m.sum() > 0: all_amps.append(np.percentile(np.abs(tr.data[m]), 98))
    global_p98 = np.median([a for a in all_amps if a > 0]) or 1.0
    print(f"Global P98: {global_p98:.2e}")

    n = len(loaded)
    if args.equal:
        SP = 1.5; AMP_V = 0.04
        fig_h = max(10, n * 0.8)
    else:
        SP = None
        # 振幅自适应: 每条波形占纵轴跨度的 ~40%
        y_min = stations[0]['dist_km']
        y_max = stations[-1]['dist_km']
        y_range = max(1, y_max - y_min)
        AMP_V = y_range / max(n, 1) * 0.40
        fig_h = max(8, n * 0.25)
    print(f"AMP_V={AMP_V:.1f}")

    fig, ax = plt.subplots(figsize=(28, fig_h))
    y_ticks_pos, y_ticks_lab = [], []

    for i, (tr, s) in enumerate(loaded):
        dec = max(1, int(tr.stats.sampling_rate / 5))
        d_norm = tr.data[::dec] / (global_p98 * 3)
        t_rel = tr.times()[::dec] - s['o']

        if args.equal:
            y = i * SP
            if i % 2 == 0:
                y_ticks_pos.append(y)
                y_ticks_lab.append(f"{s['dist_km']:.0f}")
        else:
            y = s['dist_km']

        env = envelope(d_norm)
        ax.fill_between(t_rel, y, y + env*AMP_V, color='gray', alpha=0.30, lw=0)
        ax.plot(t_rel, y + d_norm*AMP_V, 'k-', lw=0.3)

        for tt, color, yf in [
            (s['tt_p'], '#2ca02c', 0.93), (s['tt_s'], '#ff7f0e', 0.79),
            (s['tt_tp'], '#1f77b4', 0.93), (s['tt_ts'], '#d62728', 0.69),
        ]:
            if 0 <= tt <= 1200:
                ax.plot(tt, y + AMP_V*yf, 'o', color=color, ms=5, zorder=5)

        ax.text(5, y + AMP_V*0.30,
                f"{s['sta']} SNR={s['snr']:.1f} az={s['az']:.0f}",
                fontsize=7, va='center', fontweight='bold')

    # 纵轴
    if args.equal:
        ax.set_yticks(y_ticks_pos)
        ax.set_yticklabels(y_ticks_lab, fontsize=11)
        ax.set_ylim(-0.5, (n-1)*SP + 0.5)
        ylabel = 'Epicentral distance (km, equal spacing)'
    else:
        ylabel = 'Epicentral distance (km)'

    ax.grid(True, alpha=0.12)

    # 图例
    handles = [Line2D([0],[0], marker='o', color='w', markerfacecolor=c, markersize=10, label=l)
               for c,l in [('#2ca02c','P'),('#ff7f0e','S'),('#1f77b4','T-P'),('#d62728','T-S')]]
    ax.legend(handles=handles, loc='upper right', fontsize=11, ncol=4, framealpha=0.9)

    # 横轴 (分钟, 每 2min)
    ax.set_xlim(t_min, t_max)
    xt = np.arange(0, 1201, 120)
    ax.set_xticks(xt)
    ax.set_xticklabels([str(int(s/60)) for s in xt], fontsize=12)

    ax.set_xlabel('Time from origin (min)', fontsize=14, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=14, fontweight='bold')

    mode = 'equal' if args.equal else 'land'
    title = f"{evdate} Mw {evmw}  |  Full Waveform [0,+20]min  |  {n} stn (SNR≥{args.snr_min})"
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.tick_params(labelsize=11)

    plt.tight_layout()
    land_tag = f"_land{args.land_min:.0f}_{args.land_max:.0f}" if args.equal else ""
    out = f"{OUT_DIR}/profile_full_snr{args.snr_min:.0f}_{mode}{land_tag}.png"
    fig.savefig(out, dpi=300)
    print(f"OK {out}")
    
    import shutil
    shutil.copy(__file__, f"{OUT_DIR}/profile_full.py")
    print(f"OK script saved")


if __name__ == '__main__':
    main()
