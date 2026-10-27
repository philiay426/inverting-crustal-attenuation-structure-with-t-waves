#!/usr/bin/env python3
"""
波形剖面图 — 双输出（包络 + BP波形）
用法: python3 profile_conversion.py <event> --snr-min 3.0 [--t-min 0 --t-max 2000] [--from-origin]
"""
import sys, os, glob, csv, yaml, argparse
sys.path.insert(0, '/home/sy/codex/t_phase_analysis')

import numpy as np
from scipy.ndimage import uniform_filter1d
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D as MLine
from obspy import read, UTCDateTime
from obspy.signal.filter import envelope as obs_env
from geographiclib.geodesic import Geodesic
from t_phase_analysis.traveltime_profile import get_track, tt_t

CONFIG_PATH = os.path.join(os.path.dirname(__file__), '..', 'config', 'events.yaml')

def main():
    p = argparse.ArgumentParser()
    p.add_argument('event')
    p.add_argument('--snr-min', type=float, default=3.0)
    p.add_argument('--t-min', type=float, default=None)
    p.add_argument('--t-max', type=float, default=None)
    p.add_argument('--amp', type=float, default=0.10)
    p.add_argument('--exclude', type=str, default='')
    p.add_argument('--from-origin', action='store_true', help='横轴=发震时刻')
    args = p.parse_args()

    with open(CONFIG_PATH) as f: cfg = yaml.safe_load(f)[args.event]
    evla, evlo = cfg['lat'], cfg['lon']
    evdp, evdate, evmw = cfg['depth_km'], cfg['date'], cfg['mw']
    case_dir = cfg['case_dir']
    PROC_DIR = f"/mnt/f/sharedate/case/{case_dir}/sacprocess"
    CSV_PATH = f"/mnt/f/sharedate/case/{case_dir}/process_flow/1_station_distribution/stations.csv"
    OUT_DIR = os.path.dirname(__file__)

    # Read CSV
    arrivals = {}
    with open(CSV_PATH) as f:
        for row in csv.DictReader((l for l in f if not l.startswith('#'))):
            sta = row['sta']; snr = float(row['snr'])
            if snr < args.snr_min: continue
            arrivals[sta] = dict(
                tt_p=float(row['tt_p']), tt_s=float(row['tt_s']),
                tt_tp=float(row['tt_tp']), tt_ts=float(row['tt_ts']),
                land_km=float(row['land_km']), ocean_km=float(row['ocean_km']),
                dist_km=float(row['dist']), snr=snr, az=float(row['az']),
                lon=float(row['lon']), lat=float(row['lat']))

    # Compute tt_ocean/tt_land via get_track
    for sta, info in arrivals.items():
        try:
            import signal
            signal.signal(signal.SIGALRM, lambda s,f: (_ for _ in ()).throw(TimeoutError))
            signal.alarm(30)
            track = get_track(evlo, evla, info['lon'], info['lat'], 0.02,
                              ref_el=-1000, plot_map=False, plot_elevation=False)
            tr = tt_t(track, ref_el=-1000, source_depth_km=evdp, ignore_small_islands=True)
            signal.alarm(0)
            info['tt_ocean'] = tr['tt_ocean']
            info['tt_land_p'] = tr['tt_land_p']
            info['tt_land_s'] = tr['tt_land_s']
        except:
            info['tt_ocean'] = info['ocean_km'] / 1.48
            info['tt_land_p'] = info['tt_tp'] - info['tt_ocean']
            info['tt_land_s'] = info['tt_ts'] - info['tt_ocean']

    # Load waveforms
    files = sorted(glob.glob(f"{PROC_DIR}/*.BHZ*_proc.sac"))
    traces, seen = [], set()
    for fn in files:
        try:
            sta = os.path.basename(fn).split('.')[1]
            if sta not in arrivals or sta in seen: continue
            seen.add(sta); info = arrivals[sta]

            tr = read(fn)[0]
            data = tr.data.astype(np.float64)
            delta = tr.stats.delta
            # Skip bad deconvolution
            p98_check = np.percentile(np.abs(data), 98)
            if p98_check > 1e6 or not np.isfinite(p98_check): continue

            # Origin from filename: parts[-4]=year, parts[-3]=julday
            parts = os.path.basename(fn).split('.')
            time_str = parts[-2].replace('_proc','')
            origin = UTCDateTime(year=int(parts[-4]), julday=int(parts[-3]),
                                 hour=int(time_str[0:2]), minute=int(time_str[2:4]),
                                 second=int(time_str[4:6]))
            t_off = tr.stats.starttime - origin
            t = np.arange(len(data)) * delta + t_off

            env = obs_env(data)
            env = uniform_filter1d(env, size=max(int(3/delta), 3))

            tt_ocean = info['tt_ocean']
            t_aligned = t - tt_ocean

            traces.append(dict(
                sta=sta, land_km=info['land_km'], ocean_km=info['ocean_km'],
                data=data, env=env, t=t_aligned, t_origin=t, delta=delta,
                tt_tp=info['tt_land_p'], tt_ts=info['tt_land_s'],
                tt_ocean=tt_ocean, snr=info['snr']))
        except: continue

    traces.sort(key=lambda x: x['land_km'])
    if len(traces) == 0: print("无台站"); return

    # Time window
    if args.t_min is None: args.t_min = min(t['tt_tp'] for t in traces) - 50
    if args.t_max is None: args.t_max = max(t['tt_ts'] for t in traces) + 300

    # Clip extreme
    p98s = [np.percentile(np.abs(t['data'][(t['t']>=args.t_min)&(t['t']<=args.t_max)]),98) for t in traces]
    median_p98 = np.median(p98s)
    for i, t in enumerate(traces):
        if p98s[i] > median_p98 * 10.0:
            scale = median_p98 * 10.0 / p98s[i]
            traces[i]['data'] *= scale; traces[i]['env'] *= scale

    all_plot = np.concatenate([t['data'][(t['t']>=args.t_min)&(t['t']<=args.t_max)] for t in traces])
    global_p98 = np.percentile(np.abs(all_plot), 98)

    # Bin
    bin_width, orig_traces = 20.0, traces[:]
    min_land, max_land = traces[0]['land_km'], traces[-1]['land_km']
    nbins = int((max_land - min_land) / bin_width) + 1
    binned = []
    for ib in range(nbins):
        lo, hi = min_land+ib*bin_width, min_land+(ib+1)*bin_width
        members = [t for t in traces if lo <= t['land_km'] < hi]
        if not members: continue
        ref = members[0]
        env_sum = np.zeros_like(ref['env']); data_sum = np.zeros_like(ref['data'])
        count = np.zeros_like(ref['env'], dtype=int)
        for t in members:
            npts = min(len(ref['env']), len(t['env']))
            env_sum[:npts] += t['env'][:npts]; data_sum[:npts] += t['data'][:npts]
            count[:npts] += 1
        count[count==0] = 1
        binned.append(dict(land_km=(lo+hi)/2, env=env_sum/count, data=data_sum/count,
                           t=ref['t'], t_origin=ref['t_origin'],
                           tt_tp=(lo+hi)/2/6.0, tt_ts=(lo+hi)/2/3.6))
    traces = binned

    # Plotting
    min_land, max_land = traces[0]['land_km'], traces[-1]['land_km']
    margin = (max_land - min_land) * 0.05 + 20
    fig_h = min(30, max(10, (max_land - min_land) / 30))

    def do_plot(mode, outfile):
        fig, ax = plt.subplots(figsize=(16, fig_h * 0.7))
        fig.subplots_adjust(left=0.02)
        amp_scale = args.amp * 25

        for tr in traces:
            y0 = tr['land_km']
            time_axis = tr['t_origin'] if args.from_origin else tr['t']
            m = (time_axis >= args.t_min) & (time_axis <= args.t_max)
            tw = time_axis[m]
            if len(tw) == 0: continue
            if mode == 'envelope':
                ew = tr['env'][m]
                ew_n = ew / global_p98 if global_p98 > 0 else ew
                ey = ew_n * amp_scale
                ax.fill_between(tw, y0, y0+ey, color='black', alpha=0.18, lw=0, zorder=2)
                ax.plot(tw, y0+ey, color='black', lw=0.8, alpha=0.9, zorder=4)
            else:
                dw = tr['data'][m]
                dw_n = dw / global_p98 if global_p98 > 0 else dw
                dy = dw_n * amp_scale
                ax.plot(tw, y0+dy, color='black', lw=0.3, alpha=0.7, zorder=4)

        # T-P/T-S fitted lines (use total travel times for from-origin mode)
        from numpy.polynomial import Polynomial
        # For from-origin: tt_ocean + land_time = total from origin
        tp_pts = [(t['land_km'], t['tt_ocean']+t['tt_tp']) for t in orig_traces if t['tt_tp']]
        ts_pts = [(t['land_km'], t['tt_ocean']+t['tt_ts']) for t in orig_traces if t['tt_ts']]
        for pts, color in [(tp_pts,'black'),(ts_pts,'gray')]:
            pts.sort()
            x, y = [p[1] for p in pts], [p[0] for p in pts]
            if len(x) > 1:
                p = Polynomial.fit(x, y, 1)
                xs = np.linspace(min(x), max(x), 100)
                ax.plot(xs, p(xs), '--', color=color, lw=1.5, alpha=0.8, zorder=27)

        ax.set_xlim(args.t_min, args.t_max)
        ax.set_ylim(0, max_land + margin)
        ax.xaxis.set_major_locator(MultipleLocator(200))
        ax.xaxis.set_minor_locator(MultipleLocator(50))
        xlabel = 'Time from origin (s)' if args.from_origin else 'Time from ocean arrival (s)'
        ax.set_xlabel(xlabel, fontsize=14, fontweight='bold')
        ax.set_ylabel('Land distance (km)', fontsize=14, fontweight='bold')
        ax.yaxis.set_major_locator(MultipleLocator(200))
        ax.yaxis.set_minor_locator(MultipleLocator(100))
        ax.tick_params(labelsize=12)
        tag = 'Envelope (3s smooth)' if mode == 'envelope' else 'Bandpass [2-8 Hz]'
        ax.set_title(f'{evdate} Mw {evmw}  |  {n} stations  |  SNR > {args.snr_min} dB  |  {tag}',
                     fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.08, linestyle=':')
        ax.legend(handles=[
            MLine([0],[0], color='black', lw=0.8, label=tag),
            MLine([0],[0], color='black', linestyle='--', lw=1.5, label='T-P'),
            MLine([0],[0], color='gray', linestyle='--', lw=1.5, label='T-S'),
        ], fontsize=14, loc='upper left', framealpha=0.9)
        fig.savefig(outfile, dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved {outfile}")

    suffix = '_from_origin' if args.from_origin else ''
    do_plot('envelope', os.path.join(OUT_DIR, f'profile_envelope{suffix}.png'))
    do_plot('waveform', os.path.join(OUT_DIR, f'profile_waveform{suffix}.png'))

if __name__ == '__main__':
    main()
