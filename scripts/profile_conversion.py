#!/home/sy/tools/tphase_env/bin/python3
"""
波形剖面图 — 按海岸转换点对齐，陆地距离排序 [通用版 v2]
横轴: t - tt_ocean (转换点对齐, t=0 = 波上岸)
纵轴: 陆地距离 (km) | 从0起算
包络: 1s滑动RMS + 1s平滑, 单边红色 | T-P/T-S虚线从y=0起画

用法: python3 profile_conversion.py <event> --snr-min 3.0
     python3 profile_conversion.py sio2022
     python3 profile_conversion.py sio2022 --snr-min 6.0 --amp 0.15
"""
import sys, os, glob, csv, yaml, argparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 't_phase_analysis'))

import numpy as np
from scipy.ndimage import uniform_filter1d
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D as MLine
from obspy import read

# 期刊出图 rcParams
plt.rcParams.update({
    'font.size': 14, 'font.weight': 'bold',
    'axes.labelsize': 15, 'axes.titlesize': 16,
    'xtick.labelsize': 13, 'ytick.labelsize': 13,
    'legend.fontsize': 12,
})

CODE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_PATH = os.path.join(CODE_ROOT, 'config', 'events.yaml')
DATA_BASE = '/mnt/f/sharedate/case'
SOFAR_V = 1.48


def main():
    p = argparse.ArgumentParser(description='波形剖面图 — 转换点对齐 [通用版]')
    p.add_argument('event', help='events.yaml 中的事件 key')
    p.add_argument('--snr-min', type=float, default=6.0)
    p.add_argument('--t-min', type=float, default=None)
    p.add_argument('--t-max', type=float, default=None)
    p.add_argument('--amp', type=float, default=0.10)
    p.add_argument('--exclude', type=str, default='')
    p.add_argument('--rms-win', type=float, default=1.0, help='RMS窗口(秒), 默认1.0')
    p.add_argument('--smooth-win', type=float, default=1.0, help='平滑窗口(秒), 默认1.0')
    p.add_argument('--spacing', type=float, default=20.0, help='陆距分箱间隔(km), 每箱只保留SNR最高的台站')
    args = p.parse_args()

    # 读配置
    if not os.path.exists(CONFIG_PATH):
        # 回退到事件目录下的 config
        cfg_path = os.path.join(os.path.dirname(__file__), '..', 'config', 'events.yaml')
    else:
        cfg_path = CONFIG_PATH

    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)[args.event]
    evla, evlo = cfg['lat'], cfg['lon']
    evdp = cfg['depth_km']
    evdate, evmw = cfg['date'], cfg['mw']
    case_dir = cfg['case_dir']
    PROC_DIR = f"{DATA_BASE}/{case_dir}/sacprocess"
    CSV_PATH = f"{DATA_BASE}/{case_dir}/process_flow/1_station_distribution/stations.csv"
    OUT_DIR = f"{DATA_BASE}/{case_dir}/process_flow/2_waveform_profile"
    os.makedirs(OUT_DIR, exist_ok=True)

    exclude_set = set(s.strip() for s in args.exclude.split(',') if s.strip())

    # 读 CSV (SNR + 到时)
    print(f"事件: {evdate} Mw {evmw} | {evla}°, {evlo}° | 深度 {evdp} km")
    arrivals = {}
    with open(CSV_PATH) as f:
        reader = csv.DictReader((l for l in f if not l.startswith('#')))
        for row in reader:
            sta = row['sta']
            snr = float(row['snr'])
            if snr < args.snr_min: continue
            if sta in exclude_set: continue
            arrivals[sta] = dict(
                tt_p=float(row['tt_p']), tt_s=float(row['tt_s']),
                tt_tp=float(row['tt_tp']), tt_ts=float(row['tt_ts']),
                land_km=float(row['land_km']), ocean_km=float(row['ocean_km']),
                dist_km=float(row['dist']), snr=snr, az=float(row['az']),
            )
    print(f"  CSV SNR≥{args.snr_min}: {len(arrivals)} 站")

    for sta, info in arrivals.items():
        info['tt_ocean'] = info['ocean_km'] / SOFAR_V

    # 收集波形
    print("  加载波形...")
    files = sorted(glob.glob(f"{PROC_DIR}/*.BHZ*_proc.sac"))
    traces = []
    seen = set()
    for fn in files:
        try:
            sta = os.path.basename(fn).split('.')[1]
            if sta not in arrivals or sta in seen: continue
            seen.add(sta)
            info = arrivals[sta]

            tr = read(fn)[0]
            data = tr.data.astype(np.float64)
            delta = tr.stats.delta
            t_off = tr.stats.sac.b - tr.stats.sac.o
            t = np.arange(len(data)) * delta + t_off

            # 包络: RMS窗口 + 移动平均平滑
            rms_win = max(int(args.rms_win / delta), 1)
            rms = np.sqrt(uniform_filter1d(data**2, size=rms_win))
            env = uniform_filter1d(rms, size=max(int(args.smooth_win / delta), 1))

            # 转换点对齐: t' = t - tt_ocean
            tt_ocean = info['tt_ocean']
            t_aligned = t - tt_ocean

            traces.append(dict(
                sta=sta, land_km=info['land_km'], ocean_km=info['ocean_km'],
                data=data, env=env, t=t_aligned, delta=delta,
                tt_p=info['tt_p'] - tt_ocean,
                tt_s=info['tt_s'] - tt_ocean,
                tt_tp=info['tt_tp'] - tt_ocean,
                tt_ts=info['tt_ts'] - tt_ocean,
                snr=info['snr'], tt_ocean=tt_ocean
            ))
        except: continue

    print(f"  加载 {len(traces)} 站波形")

    # 20km分箱: 每箱只保留SNR最高的台站
    if args.spacing > 0 and len(traces) > 0:
        binned = {}
        for tr in traces:
            bin_key = int(tr['land_km'] // args.spacing)
            if bin_key not in binned or tr['snr'] > binned[bin_key]['snr']:
                binned[bin_key] = tr
        traces = list(binned.values())
        print(f"  分箱(spacing={args.spacing}km)后: {len(traces)} 站")

    traces.sort(key=lambda x: x['land_km'])
    n = len(traces)
    if n == 0:
        print("无台站"); return

    # 时间窗口
    if args.t_min is None:
        args.t_min = min(t['tt_tp'] for t in traces) - 200
    if args.t_max is None:
        args.t_max = max(t['tt_ts'] for t in traces) + 300

    # 振幅裁剪 (3x 中位数 P98)
    p98s = [np.percentile(np.abs(t['env'][(t['t'] >= args.t_min) & (t['t'] <= args.t_max)]), 98)
            for t in traces]
    median_p98 = np.median(p98s)
    clip_ceil = median_p98 * 3.0
    for i, t in enumerate(traces):
        if p98s[i] > clip_ceil:
            scale = clip_ceil / p98s[i]
            traces[i]['data'] = t['data'] * scale
            traces[i]['env'] = t['env'] * scale
            print(f"  Clipped {t['sta']}: {p98s[i]/median_p98:.1f}x → 3x")

    # 全局 P98 (基于包络)
    all_env = np.concatenate([t['env'][(t['t'] >= args.t_min) & (t['t'] <= args.t_max)]
                              for t in traces])
    global_p98 = np.percentile(np.abs(all_env), 98)
    print(f"  Global P98: {global_p98:.2e}")

    # 绘图
    min_land = traces[0]['land_km']
    max_land = traces[-1]['land_km']
    fig_h = min(30, max(10, (max_land - min_land) / 30))
    fig, ax = plt.subplots(figsize=(24, fig_h))
    amp_scale = args.amp * 25

    for tr in traces:
        y0 = tr['land_km']
        m = (tr['t'] >= args.t_min) & (tr['t'] <= args.t_max)
        tw = tr['t'][m]; dw = tr['data'][m]; ew = tr['env'][m]
        if len(tw) == 0: continue

        if global_p98 > 0: dw_n = dw / global_p98; ew_n = ew / global_p98
        else: dw_n = dw; ew_n = ew
        dy = dw_n * amp_scale; ey = ew_n * amp_scale

        ax.plot(tw, y0+dy, 'k-', lw=0.5, alpha=0.9, zorder=3)
        ax.plot(tw, y0+ey, 'r-', lw=0.6, alpha=0.7, zorder=4)

    # T-P / T-S 理论走时线（从 y=0 起画的虚线）
    v_p, v_s = 6.0, 3.6
    y_ref = np.linspace(0, max_land + 100, 200)
    ax.plot(y_ref / v_p, y_ref, 'k--', lw=2.0, alpha=0.6, zorder=25)
    ax.plot(y_ref / v_s, y_ref, 'k:', lw=2.0, alpha=0.6, zorder=25)

    # 图表设置
    ax.set_xlim(args.t_min, args.t_max)
    ax.set_ylim(0, max_land + 100)
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.xaxis.set_minor_locator(MultipleLocator(20))
    ax.set_xlabel('Time from ocean arrival (s)', fontsize=10, fontweight='bold')
    ax.set_ylabel('Land distance (km)', fontsize=10, fontweight='bold')
    ax.yaxis.set_major_locator(MultipleLocator(100))
    ax.yaxis.set_minor_locator(MultipleLocator(50))
    ax.tick_params(labelsize=7)

    ax.set_title(f'{evdate} Mw {evmw}  |  {n} stations  |  SNR ≥ {args.snr_min} dB  |  Aligned by conversion point',
                 fontsize=11, fontweight='bold')
    ax.grid(True, alpha=0.08, linestyle=':')

    # 图例
    leg = [
        MLine([0],[0], color='black', lw=0.8, label='2-8 Hz waveform'),
        MLine([0],[0], color='red', lw=1.5, alpha=0.8, label='Envelope'),
        MLine([0],[0], color='black', lw=2.0, linestyle='--', alpha=0.6, label=f'T-P ({v_p} km/s)'),
        MLine([0],[0], color='black', lw=2.0, linestyle=':', alpha=0.6, label=f'T-S ({v_s} km/s)'),
    ]
    ax.legend(handles=leg, fontsize=18, loc='upper left', ncol=1, framealpha=0.9)

    out = os.path.join(OUT_DIR, f'profile_conversion_snr{args.snr_min}.png')
    fig.savefig(out, dpi=200, bbox_inches='tight')
    print(f"✅ {out}")


if __name__ == '__main__':
    main()
