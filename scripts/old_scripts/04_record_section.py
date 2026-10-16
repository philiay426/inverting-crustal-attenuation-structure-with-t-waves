#!/usr/bin/env python3
"""步骤4: 波形剖面视图 — 包络+理论到时圆点"""
import sys, os, glob, yaml, argparse
for p in ['/home/philiay/shared/codex', '/home/philiay/shared/codex/t_phase_analysis']:
    if p not in sys.path: sys.path.insert(0, p)

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D
import numpy as np
from obspy import read
from obspy.signal.filter import envelope as obs_env
from lib.arrivals import compute_theoretical_arrivals, dist_az

CONFIG = yaml.safe_load(open("/home/philiay/shared/codex/config/events.yaml"))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('event')
    p.add_argument('--tmin', type=float, default=1800)
    p.add_argument('--tmax', type=float, default=2500)
    p.add_argument('--az-min', type=float, default=0)
    p.add_argument('--az-max', type=float, default=360)
    args = p.parse_args()

    cfg = CONFIG[args.event]
    case_dir = cfg.get('case_dir', cfg['label'])
    proc_dir = f"/home/philiay/shared/case/{case_dir}/sacprocess"
    out_dir = f"/mnt/f/sharedate/case/{case_dir}"
    os.makedirs(out_dir, exist_ok=True)

    bhz = sorted(glob.glob(f"{proc_dir}/*.BHZ*_proc.sac"))
    stations = []

    for fn in bhz:
        try:
            tr = read(fn)[0]; sac = tr.stats.sac
            stla = getattr(sac,'stla',0); stlo = getattr(sac,'stlo',0)
            if stla == 0: continue
            dist, az, baz = dist_az(cfg['lon'], cfg['lat'], stlo, stla)
            if az < args.az_min or az > args.az_max: continue
            sta = os.path.basename(fn).split('.')[1]
            arrs = compute_theoretical_arrivals(cfg['lon'], cfg['lat'], cfg['depth_km'], stlo, stla)
            # SAC 时间轴偏移
            t_offset = tr.stats.sac.b - tr.stats.sac.o
            env = obs_env(tr.data.astype(np.float64))
            t = np.arange(len(tr.data)) * tr.stats.delta + t_offset
            stations.append(dict(sta=sta, dist=dist, data=tr.data.astype(np.float64),
                                env=env, t=t, delta=tr.stats.delta,
                                tp=arrs.get('T-P',0), ts=arrs.get('T-S',0),
                                tl=arrs.get('T-Love',0)))
        except: continue

    stations.sort(key=lambda x: x['dist'])
    n = len(stations)
    print(f"{n} stations, AZ {args.az_min}-{args.az_max}°, dist {stations[0]['dist']:.0f}-{stations[-1]['dist']:.0f}km")

    dists = [s['dist'] for s in stations]; dmin, dmax = min(dists), max(dists)
    dr = max(dmax-dmin, 100); spacing = dr/max(n-1, 1); amp_s = spacing*0.18

    fig, ax = plt.subplots(figsize=(20, max(8, n*0.25)))
    for s in stations:
        y0 = s['dist']
        m = (s['t'] >= args.tmin) & (s['t'] <= args.tmax)
        tp_ = s['t'][m]; bp = s['data'][m]; en = s['env'][m]
        mx = max(np.max(np.abs(bp)), 1e-30); dn = bp/mx; en_n = en/mx
        ax.fill_between(tp_, y0, y0+en_n*amp_s, color='#e74c3c', alpha=0.12, lw=0)
        ax.plot(tp_, y0+en_n*amp_s, color='#e74c3c', lw=0.3, alpha=0.6)
        ax.plot(tp_, y0+dn*amp_s, 'k-', lw=0.12, alpha=0.4)
        for pk, pc, pm in [('tp','#1f77b4','o'),('ts','#d62728','D'),('tl','#9467bd','s')]:
            tt = s.get(pk,0)
            if args.tmin <= tt <= args.tmax:
                idx = np.argmin(np.abs(tp_ - tt))
                if idx < len(dn):
                    ax.plot(tt, y0+dn[idx]*amp_s, marker=pm, color=pc, markersize=2.5, markeredgewidth=0, zorder=10)
        ax.text(args.tmin+20, y0, s['sta'], fontsize=4.5, fontfamily='monospace', fontweight='bold', va='center', ha='left')

    ax.set_xlim(args.tmin, args.tmax)
    ax.set_ylim(dmin-dr*0.04-amp_s, dmax+dr*0.04+amp_s*1.2)
    ax.set_xlabel('Time (s)', fontsize=12); ax.set_ylabel('Distance (km)', fontsize=12)
    ax.xaxis.set_major_locator(MultipleLocator(200)); ax.grid(True, axis='x', alpha=0.06, linestyle=':')
    leg = [Line2D([0],[0],color='black',lw=0.3,label='2-8 Hz'),
           Line2D([0],[0],color='#e74c3c',lw=0.8,label='Envelope'),
           Line2D([0],[0],marker='o',color='w',markerfacecolor='#1f77b4',markersize=4,label='T-P'),
           Line2D([0],[0],marker='D',color='w',markerfacecolor='#d62728',markersize=4,label='T-S'),
           Line2D([0],[0],marker='s',color='w',markerfacecolor='#9467bd',markersize=4,label='T-Love')]
    ax.legend(handles=leg, loc='upper right', fontsize=8)
    ax.set_title(f'{cfg["label"]} — record section ({n} stations) | 2-8 Hz', fontsize=11, fontweight='bold')
    plt.subplots_adjust(left=0.08, right=0.97, top=0.96, bottom=0.05)
    fig.savefig(f"{out_dir}/record_section.png", dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"\n✅ {out_dir}/record_section.png/pdf")

if __name__ == '__main__':
    main()
