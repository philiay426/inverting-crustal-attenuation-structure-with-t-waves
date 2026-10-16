#!/usr/bin/env python3
"""
步骤8: 单台ZRT三分量包络+时频图
左=Z/R/T波形(黑色)+包络(红色)，右=各分量的时频图
用法:
  python3 scripts/08_station_zrt_spec.py <event> <station>
示例:
  python3 scripts/08_station_zrt_spec.py puerto_rico K62A
  python3 scripts/08_station_zrt_spec.py mid_atlantic_2014 K62A
"""
import sys, os, yaml, argparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 't_phase_analysis'))
import numpy as np
from scipy.signal import hilbert
from scipy.ndimage import uniform_filter1d
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D
from obspy import read
from obspy.taup import TauPyModel
from obspy.geodetics.base import locations2degrees
from geographiclib.geodesic import Geodesic
from t_phase_analysis.traveltime_profile import get_track, tt_t
from t_phase_analysis.waveform import compute_spectrogram, plot_spectrogram_panel


def main():
    p = argparse.ArgumentParser(description='单台ZRT三分量包络+时频图')
    p.add_argument('event', help='事件名称（对应 events.yaml 中的 key）')
    p.add_argument('station', help='台站代码，如 K62A')
    p.add_argument('--network', default='TA', help='台网代码（默认 TA）')
    p.add_argument('--tmin', type=float, default=-500, help='起始时间 (s, 默认 -500)')
    p.add_argument('--tmax', type=float, default=3500, help='结束时间 (s, 默认 3500)')
    p.add_argument('--freq-max', type=float, default=10, help='时频图频率上限 Hz (默认 10)')
    p.add_argument('--ylim', type=float, default=0.5, help='y轴上下限 (默认 0.5)')
    p.add_argument('--amp', type=float, default=0.08, help='振幅缩放因子 (默认 0.08)')
    p.add_argument('--no-rotation', action='store_true', help='不旋转ZNE→ZRT，直接显示ZNE')
    args = p.parse_args()

    # ── 加载事件配置 ──
    cfg_path = os.path.join(os.path.dirname(__file__), '..', 'config', 'events.yaml')
    with open(cfg_path) as f:
        config = yaml.safe_load(f)
    if args.event not in config:
        print(f"错误: 事件 '{args.event}' 不在 events.yaml 中")
        print(f"可用事件: {', '.join(config.keys())}")
        sys.exit(1)
    ev = config[args.event]
    evla, evlo = ev['lat'], ev['lon']
    evdp = ev.get('depth_km', 10.0)
    sac_date = ev.get('sac_date', '')
    mw = ev.get('mw', '')
    case_dir = ev.get('case_dir', ev['label'])

    proc_dir = f"/mnt/f/sharedate/case/{case_dir}/sacprocess"
    out_dir  = f"/mnt/f/sharedate/case/{case_dir}/waveform_contrast"
    os.makedirs(out_dir, exist_ok=True)

    # ── 加载三分量 SAC ──
    sta, net = args.station.upper(), args.network.upper()
    trs = {}
    found_chs = []
    # 尝试 BHZ/BHN/BHE, 也尝试 BH1/BH2/BHZ
    for ch in ['BHZ', 'BHN', 'BHE', 'BH1', 'BH2']:
        pattern = f"{proc_dir}/{net}.{sta}..{ch}"
        import glob
        matches = sorted(glob.glob(pattern + "*_proc.sac"))
        if not matches:
            matches = sorted(glob.glob(pattern + "*.sac"))
        if matches:
            trs[ch] = read(matches[0])[0]
            found_chs.append(ch)

    for req_ch in ['BHZ', 'BHN', 'BHE']:
        if req_ch not in trs and req_ch not in ['BH1', 'BH2']:
            # 尝试 BH1/BH2 替代 BHN/BHE
            alt = {'BHN': 'BH1', 'BHE': 'BH2'}
            if req_ch in alt and alt[req_ch] in trs:
                trs[req_ch] = trs[alt[req_ch]]
                del trs[alt[req_ch]]

    if 'BHZ' not in trs or 'BHN' not in trs or 'BHE' not in trs:
        print(f"错误: 找不到 {net}.{sta} 的三分量 SAC 文件 ({proc_dir})")
        print(f"找到: {list(trs.keys())}")
        sys.exit(1)

    # 对齐时间轴
    t0 = max(tr.stats.starttime for tr in trs.values())
    t1 = min(tr.stats.endtime for tr in trs.values())
    for ch in trs: trs[ch].trim(t0, t1)

    delta = trs['BHZ'].stats.delta
    sac = trs['BHZ'].stats.sac
    try: stla = sac.stla
    except AttributeError: stla = evla
    try: stlo = sac.stlo
    except AttributeError: stlo = evlo
    t_off = (sac.b - sac.o) if hasattr(sac, 'b') and hasattr(sac, 'o') else 0
    npts = len(trs['BHZ'].data)
    t = np.arange(npts) * delta + t_off

    # ── 方位角 ──
    g = Geodesic.WGS84.Inverse(evla, evlo, stla, stlo)
    az = g['azi1']; az = az+360 if az<0 else az
    baz = (az+180)%360; br = np.radians(baz)

    # ── ZNE → ZRT 旋转 ──
    z = trs['BHZ'].data.astype(np.float64)
    n = trs['BHN'].data.astype(np.float64)
    e = trs['BHE'].data.astype(np.float64)
    if args.no_rotation:
        comps = [('Z', z), ('N', n), ('E', e)]
    else:
        r = n*np.cos(br) + e*np.sin(br)
        tc = -n*np.sin(br) + e*np.cos(br)
        comps = [('Z', z), ('R', r), ('T', tc)]

    # ── 理论到时 ──
    taup = TauPyModel(model='iasp91')
    dd = locations2degrees(evla, evlo, stla, stlo)
    arrs = taup.get_travel_times(source_depth_in_km=evdp, distance_in_degree=dd,
                                 phase_list=['P', 'S'])
    tt_P = next((a.time for a in arrs if a.name == 'P'), None)
    tt_S = next((a.time for a in arrs if a.name == 'S'), None)
    track = get_track(evlo, evla, stlo, stla, 0.02, ref_el=-1000,
                      plot_map=False, plot_elevation=False)
    tt_res = tt_t(track, ref_el=-1000, source_depth_km=evdp,
                  ignore_small_islands=True)
    tt_TP = tt_res.get('tt_t', None)
    tt_TS = tt_res['tt_ocean'] + tt_res['tt_land_s'] if 'tt_land_s' in tt_res else None
    land_km = tt_res.get('dist_land', 0) * 111.0
    ocean_km = tt_res.get('dist_ocean', 0) * 111.0

    # ── 包络 ──
    envs = {}
    for name, data in comps:
        envs[name] = uniform_filter1d(
            np.abs(hilbert(data)), size=max(int(2/delta), 1))

    # ── 时窗 ──
    tmin_f, tmax_f = args.tmin, args.tmax
    m = (t >= tmin_f) & (t <= tmax_f)
    tw = t[m]

    arrivals = [('P', tt_P, '#2ca02c'), ('S', tt_S, '#ff7f0e'),
                ('T-P', tt_TP, '#1f77b4'), ('T-S', tt_TS, '#d62728')]

    # ── 绘图: 3行×2列 ──
    fig, axes = plt.subplots(3, 2, figsize=(20, 10),
                             gridspec_kw={'width_ratios': [1, 1.3],
                                          'hspace': 0.28, 'wspace': 0.30})
    plt.subplots_adjust(left=0.06, right=0.98, top=0.93, bottom=0.06)

    # ── 公共归一化因子（取三个分量 P98 的最大值，保持振幅可比）──
    common_p98 = max(np.percentile(np.abs(data), 98) for _, data in comps)

    for i, (name, data) in enumerate(comps):
        ax_w = axes[i, 0]
        env = envs[name]

        dw = data[m]; ew = env[m]
        p98 = common_p98
        dn = dw/p98 if p98 > 0 else dw
        en = ew/p98 if p98 > 0 else ew
        amp = args.amp

        # 波形 + 包络
        ax_w.plot(tw, dn*amp, 'k-', lw=0.35, alpha=0.6)
        ax_w.fill_between(tw, 0, en*amp, color='red', alpha=0.12, lw=0)
        ax_w.plot(tw, en*amp, 'r-', lw=0.7, alpha=0.6)
        ax_w.axhline(0, color='gray', lw=0.2)

        # 到时标注
        for ph, tt, clr in arrivals:
            if tt and tmin_f <= tt <= tmax_f:
                ax_w.axvline(x=tt, color=clr, lw=0.8, ls='--', alpha=0.5)
                y_pos = {'P': 0.92, 'S': 0.78, 'T-P': 0.92, 'T-S': 0.64}
                ha_side = {'P': 'left', 'S': 'left', 'T-P': 'right', 'T-S': 'right'}
                ax_w.annotate(f'{ph}={tt:.0f}s', xy=(tt, 0.5),
                            xytext=(tt, y_pos.get(ph, 0.8)),
                            fontsize=6.5, color=clr, fontweight='bold',
                            ha=ha_side.get(ph, 'center'), va='bottom',
                            xycoords=('data', 'axes fraction'),
                            textcoords=('data', 'axes fraction'),
                            arrowprops=dict(arrowstyle='-', color=clr, lw=0.3, alpha=0.5))

        ax_w.set_ylabel(name, fontsize=12, fontweight='bold', color='black')
        ax_w.set_ylim(-args.ylim, args.ylim)
        ax_w.set_xlim(tmin_f, tmax_f)
        ax_w.xaxis.set_major_locator(MultipleLocator(500))
        ax_w.xaxis.set_minor_locator(MultipleLocator(100))
        ax_w.tick_params(labelsize=7)
        if i < 2:
            ax_w.tick_params(labelbottom=False)
        else:
            ax_w.set_xlabel('Time from origin (s)', fontsize=9)
        ax_w.grid(True, alpha=0.06, linestyle=':')

        # 时频图
        ax_s = axes[i, 1]
        f_spec, t_spec, Sxx_db = compute_spectrogram(data, 1/delta)
        t_spec_rel = t_spec + t_off
        plot_spectrogram_panel(ax_s, t_spec_rel, f_spec, Sxx_db,
                               tt_p=tt_P, tt_s=tt_S, tt_t=tt_TP,
                               has_twave=True, freq_max=args.freq_max,
                               db_range=50, db_percentile=98, cmap='inferno')
        ax_s.set_xlim(tmin_f, tmax_f)
        ax_s.tick_params(labelsize=7)
        if i < 2:
            ax_s.tick_params(labelbottom=False)
        else:
            ax_s.set_xlabel('Time from origin (s)', fontsize=9)
        if i == 1:
            ax_s.set_ylabel('Frequency (Hz)', fontsize=9)
        ax_s.xaxis.set_major_locator(MultipleLocator(500))
        ax_s.xaxis.set_minor_locator(MultipleLocator(100))

    # 图例
    leg = [Line2D([0],[0], color='k', lw=0.5, label='Waveform'),
           Line2D([0],[0], color='red', lw=0.6, alpha=0.6, label='Envelope'),
           Line2D([0],[0], marker='o', color='#2ca02c', linestyle='None', ms=4, label='P'),
           Line2D([0],[0], marker='o', color='#ff7f0e', linestyle='None', ms=4, label='S'),
           Line2D([0],[0], marker='o', color='#1f77b4', linestyle='None', ms=4, label='T-P'),
           Line2D([0],[0], marker='o', color='#d62728', linestyle='None', ms=4, label='T-S')]
    axes[0, 0].legend(handles=leg, fontsize=6.5, loc='upper right',
                      ncol=6, framealpha=0.85)

    fig.suptitle(f'{net}.{sta}  |  AZ={az:.1f}°  LAND={land_km:.0f}km  '
                 f'OCEAN={ocean_km:.0f}km  |  2-8 Hz',
                 fontsize=11, fontweight='bold')

    out = f"{out_dir}/{sta}_ZRT_spectrogram.png"
    fig.savefig(out, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"✅ {out}")


if __name__ == '__main__':
    main()
