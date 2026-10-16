#!/usr/bin/env python3
"""步骤7: F-K分析 — PMCC方位角/视速度反演 (T_Phase_Analysis)"""
import sys, os, glob, argparse, yaml

for p in ['/home/philiay/shared/codex', '/home/philiay/shared/codex/t_phase_analysis']:
    if p not in sys.path: sys.path.insert(0, p)

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from obspy import read, Stream

from t_phase_analysis import pmcc, plot_pmcc_diagnostics
from lib.arrivals import compute_theoretical_arrivals, dist_az

CONFIG = yaml.safe_load(open(f"{CODE_ROOT}/config/events.yaml"))
OUT_ROOT = "/mnt/f/sharedate/oi"  # removed, using CODE_ROOT instead


def main():
    p = argparse.ArgumentParser(description='PMCC F-K 方位角/视速度分析')
    p.add_argument('event')
    p.add_argument('--az-min', type=float, default=305)
    p.add_argument('--az-max', type=float, default=315)
    p.add_argument('--window-center', type=float, default=None,
                   help='分析窗口中心 (s from origin), 默认理论T-P中值')
    p.add_argument('--window-half', type=float, default=60,
                   help='分析窗口半宽 (s), 默认120s窗口')
    p.add_argument('--fmin', type=float, default=2.0)
    p.add_argument('--fmax', type=float, default=8.0)
    p.add_argument('--max-stations', type=int, default=20,
                   help='最多取多少个台站（PMCC需要3个，取N个选最好的3个）')
    args = p.parse_args()

    cfg = CONFIG[args.event]
    case_dir = cfg.get('case_dir', cfg['label'])
    proc_dir = f"/home/philiay/shared/case/{case_dir}/sacprocess"
    out_dir = f"/mnt/f/sharedate/case/{case_dir}"
    os.makedirs(out_dir, exist_ok=True)

    # ── 找到方位角范围内的台站 ──
    bhz = sorted(glob.glob(f"{proc_dir}/*.BHZ*_proc.sac"))
    stations = []
    for fn in bhz:
        try:
            tr = read(fn, headonly=True)[0]; sac = tr.stats.sac
            stla, stlo = getattr(sac,'stla',0), getattr(sac,'stlo',0)
            if stla == 0: continue
            dist, az, baz = dist_az(cfg['lon'], cfg['lat'], stlo, stla)
            if az < args.az_min or az > args.az_max: continue
            sta = os.path.basename(fn).split('.')[1]
            arrs = compute_theoretical_arrivals(cfg['lon'], cfg['lat'], cfg['depth_km'], stlo, stla)
            tp = arrs.get('T-P', 0)
            stations.append(dict(sta=sta, stla=stla, stlo=stlo, az=az, dist=dist, tp=tp, fn=fn))
        except: continue

    stations.sort(key=lambda x: x['az'])
    stations = stations[:args.max_stations]
    print(f"{len(stations)} stations, AZ {args.az_min}-{args.az_max}°")

    if len(stations) < 3:
        print("❌ PMCC需要至少3个台站")
        return

    # ── 确定分析窗口 ──
    if args.window_center is None:
        tps = [s['tp'] for s in stations if s['tp'] > 0]
        window_center = np.median(tps) if tps else 2000
    else:
        window_center = args.window_center

    time_start = window_center - args.window_half
    time_length = args.window_half * 2

    print(f"PMCC窗口: {time_start:.0f}-{time_start+time_length:.0f}s "
          f"(中心={window_center:.0f}s)")
    print(f"台站群: {', '.join(s['sta'] for s in stations[:10])}")

    # ── 组3分量Stream ──
    st = Stream()
    for s in stations:
        fn_bhz = s['fn']
        # 加载BHZ数据用于PMCC（proc.sac 已做2-8 Hz滤波）
        tr_z = read(fn_bhz)[0]
        # 确保有台站坐标头段
        tr_z.stats.sac.stla = s['stla']
        tr_z.stats.sac.stlo = s['stlo']
        tr_z.stats.sac.t0 = 0  # PMCC参考时刻
        st.append(tr_z)

    # ── PMCC ──
    print("Running PMCC...")
    try:
        result = pmcc(st, time_start=time_start,
                      time_one_slice=20.0, time_slide=5.0,
                      time_length=time_length)
    except Exception as e:
        print(f"❌ PMCC失败: {e}")
        # 尝试只用前3个台站
        print("Retrying with first 3 stations...")
        st3 = Stream([st[0], st[1], st[2]])
        try:
            result = pmcc(st3, time_start=time_start,
                          time_one_slice=20.0, time_slide=5.0,
                          time_length=time_length)
        except Exception as e2:
            print(f"❌ Retry failed: {e2}")
            return

    # ── 诊断图 ──
    fig, axs = plot_pmcc_diagnostics(result)
    fig.suptitle(f'{cfg["label"]} — PMCC F-K {args.az_min}-{args.az_max}° '
                 f'({window_center:.0f}s) | {len(st)} stas | {args.fmin:.0f}-{args.fmax:.0f} Hz',
                 fontsize=11, fontweight='bold')
    plt.subplots_adjust(top=0.92)
    fig.savefig(f"{out_dir}/pmcc_az{args.az_min:.0f}_{args.az_max:.0f}.png",
                dpi=300, bbox_inches='tight')
    plt.close(fig)

    # ── 统计 ──
    print(f"\n✅ {out_dir}/pmcc_az{args.az_min:.0f}_{args.az_max:.0f}.png")
    if 'azimuth_mode1' in result:
        print(f"  方位角(mode1): {result['azimuth_mode1']:.1f}°")
        print(f"  视速度(mode1): {result['speed_mode1']:.1f} m/s")
    if 'azimuth_mode2' in result:
        print(f"  方位角(mode2): {result['azimuth_mode2']:.1f}°")
        print(f"  视速度(mode2): {result['speed_mode2']:.1f} m/s")

if __name__ == '__main__':
    main()
