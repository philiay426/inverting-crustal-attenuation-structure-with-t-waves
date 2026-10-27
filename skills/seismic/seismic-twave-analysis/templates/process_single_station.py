#!/usr/bin/env python3
"""
process_single_station.py — 单台三分量 SAC 预处理 + 并排图 (2026-08-01)

核心: tr.simulate(paz_remove=paz_dict) 直接出速度 (BOSA 验证方案)
用 --plot 输出左=带通包络 + 右=频谱图的 3×2 并排图

用法:
  python3 process_single_station.py <STATION> --network GT --raw-dir <DIR>
  python3 process_single_station.py <STATION> --network GT --raw-dir <DIR> --plot

输出路径: process_flow/2_single_station_analysis/{STA}/

边缘情况:
  - 依赖 codex/lib/processing.py (parse_sacpz, remove_response)
  - 震中默认 (-42.78, 42.19) 是 Prince Edward 区域，其他事件需 --evlo --evla
  - 走时仅计算 P/S（TauP iasp91），T-P/T-S 需从 stations.csv 读取再传入
"""
import sys, os, glob, argparse
sys.path.insert(0, '/home/sy/codex/lib')
import numpy as np

from obspy import read, UTCDateTime
from obspy.geodetics.base import locations2degrees, gps2dist_azimuth
from obspy.io.sac import SACTrace
from obspy.taup import TauPyModel

from scipy.signal import butter, filtfilt, hilbert
from processing import parse_sacpz, remove_response, find_sacpz

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ===== 默认配置 =====
FMIN, FMAX = 2.0, 8.0
PRE_FILT = (0.001, 0.005, 45.0, 50.0)
WATER_LEVEL = 60.0
TAUP_MODEL = 'iasp91'
EVENT_DEPTH = 10.0  # km
TAPER = 0.05


def bandpass_filter(tr, freqmin=FMIN, freqmax=FMAX, corners=4):
    """零相位带通滤波 (scipy butter + filtfilt)"""
    b, a = butter(corners, [freqmin, freqmax], btype='band', fs=tr.stats.sampling_rate)
    data = tr.data.astype(np.float64)
    filtered = filtfilt(b, a, data)
    tr_out = tr.copy()
    tr_out.data = filtered
    return tr_out


def write_sac(tr, filepath):
    """保存 SAC, 保留头信息"""
    sac = SACTrace.from_obspy_trace(tr)
    sac.write(filepath)


def compute_arrivals(stla, stlo, evlo, evla, depth=10.0):
    """TauP 计算 P/S 走时"""
    arrivals = {}
    dist_deg = locations2degrees(evla, evlo, stla, stlo)
    taup = TauPyModel(model=TAUP_MODEL)
    try:
        arrs = taup.get_travel_times(source_depth_in_km=depth,
                                      distance_in_degree=dist_deg,
                                      phase_list=['P', 'S'])
        for a in arrs:
            if a.name == 'P' and 'P' not in arrivals:
                arrivals['P'] = a.time
            elif a.name == 'S' and 'S' not in arrivals:
                arrivals['S'] = a.time
    except Exception as e:
        print(f"  [WARN] TauP failed: {e}")
    return arrivals


def main():
    p = argparse.ArgumentParser(
        description='单台三分量 SAC 预处理: 去响应 + 带通滤波 + 波形预览')
    p.add_argument('station', help='台站码 (如 BOSA)')
    p.add_argument('--network', default='GT', help='台网码')
    p.add_argument('--raw-dir', required=True,
                   help='原始 SAC 目录 (含 SACPZ/ 子目录)')
    p.add_argument('--out-dir', default=None,
                   help='输出目录 (默认: <raw-dir>/../../process_flow/2_single_station_analysis/<STA>/)')
    p.add_argument('--evlo', type=float, default=42.19,
                   help='震中经度')
    p.add_argument('--evla', type=float, default=-42.78,
                   help='震中纬度')
    p.add_argument('--evdp', type=float, default=10.0,
                   help='震源深度 (km)')
    p.add_argument('--fmin', type=float, default=FMIN)
    p.add_argument('--fmax', type=float, default=FMAX)
    p.add_argument('--plot', action='store_true',
                   help='生成 ZRT 三分量波形预览图')
    args = p.parse_args()

    evlo, evla = args.evlo, args.evla
    evdp = args.evdp
    raw_dir = os.path.abspath(args.raw_dir)
    sacpz_dir = os.path.join(raw_dir, 'SACPZ')
    out_dir = args.out_dir or os.path.join(os.path.dirname(raw_dir), '..',
                                            'process_flow', '2_single_station_analysis', args.station)
    resp_dir = os.path.join(out_dir, 'remove_response')
    os.makedirs(resp_dir, exist_ok=True)

    print(f"处理 {args.network}.{args.station}")
    print(f"  震中: ({evla:.2f}°N, {evlo:.2f}°E), 深度={evdp}km")
    print(f"  原始 SAC: {raw_dir}")
    print(f"  输出: {out_dir}")
    print(f"  带通: {args.fmin}-{args.fmax} Hz")

    # 查找三分量原始文件
    raw_files = {}
    for comp in ['BHZ', 'BHN', 'BHE']:
        for pat in [f'{args.network}.{args.station}..{comp}.*.SAC',
                     f'{args.network}.{args.station}.{comp}.*.SAC',
                     f'{args.network}.{args.station}.*.{comp}.*.SAC']:
            m = glob.glob(os.path.join(raw_dir, pat))
            if m:
                raw_files[comp] = m[0]
                break

    if 'BHZ' not in raw_files:
        print("错误: 找不到 BHZ 原始文件")
        sys.exit(1)

    # 获取台站坐标
    tr_bhz = read(raw_files['BHZ'])[0]
    try:
        stla = tr_bhz.stats.sac.stla
        stlo = tr_bhz.stats.sac.stlo
    except AttributeError:
        print("错误: BHZ SAC 头缺少 stla/stlo")
        sys.exit(1)

    dist_deg = locations2degrees(evla, evlo, stla, stlo)
    _, az, baz = gps2dist_azimuth(evla, evlo, stla, stlo)
    arrivals = compute_arrivals(stla, stlo, evlo, evla, evdp)
    print(f"  台站: {stla:.3f}°N, {stlo:.3f}°E")
    print(f"  距离: {dist_deg:.2f}°  AZ={az:.1f}°  BAZ={baz:.1f}°")
    print(f"  到时: {arrivals}")

    # 逐分量处理
    processed = {}
    for comp in ['BHZ', 'BHN', 'BHE']:
        if comp not in raw_files:
            continue

        raw_path = raw_files[comp]
        print(f"\n--- {comp} ---")

        # 查找 SACPZ
        sacpz_path = None
        for sp in [f'SACPZ.{args.network}.{args.station}.*.{comp}',
                    f'SACPZ.*.{args.station}.*.{comp}',
                    f'SACPZ..{args.station}..{comp}']:
            sm = glob.glob(os.path.join(sacpz_dir, sp))
            if sm:
                sacpz_path = sm[0]
                break
        if not sacpz_path:
            print(f"  跳过: 缺少 SACPZ")
            continue

        tr = read(raw_path)[0]
        tr_resp = remove_response(tr, sacpz_path)
        resp_path = os.path.join(resp_dir,
            os.path.basename(raw_path).replace('.SAC', '_resp.sac'))
        write_sac(tr_resp, resp_path)
        print(f"  → 去响应: {os.path.basename(resp_path)}")

        tr_proc = bandpass_filter(tr_resp, args.fmin, args.fmax)
        proc_path = os.path.join(out_dir,
            os.path.basename(raw_path).replace('.SAC', '_proc.sac'))
        write_sac(tr_proc, proc_path)
        print(f"  → 带通: {os.path.basename(proc_path)}")

        processed[comp] = {'raw': tr, 'resp': tr_resp, 'proc': tr_proc}

    if not processed:
        print("错误: 无分量处理成功")
        sys.exit(1)

    # 并排图：左=带通包络(ZRT) + 右=频谱图(0-10Hz) — 单台分析标准图
    if args.plot and len(processed) == 3:
        from scipy.signal import spectrogram as _sp_spectrogram

        br = np.radians(baz)
        comps_zrt = {}
        for comp_name, src_comp in [('Z', 'BHZ'), ('N', 'BHN'), ('E', 'BHE')]:
            if src_comp in processed:
                comps_zrt[src_comp] = bandpass_filter(
                    processed[src_comp]['resp'], args.fmin, args.fmax)

        ml = min(len(comps_zrt[c].data) for c in ['BHZ', 'BHN', 'BHE'])
        raw_z, raw_n, raw_e = (comps_zrt['BHZ'].data[:ml].astype(np.float64),
                                comps_zrt['BHN'].data[:ml].astype(np.float64),
                                comps_zrt['BHE'].data[:ml].astype(np.float64))
        rz = raw_n * np.cos(br) + raw_e * np.sin(br)
        tc = -raw_n * np.sin(br) + raw_e * np.cos(br)
        comps = [('Z', raw_z), ('R', rz), ('T', tc)]

        sr = processed['BHZ']['resp'].stats.sampling_rate
        t_axis = np.arange(ml) / sr

        # ——— 左侧: 带通包络 ———
        p98_bp = np.percentile(np.abs(np.concatenate([raw_z, rz, tc])), 98)
        NP = 3 if p98_bp > 0 else 1

        nperseg = min(int(sr * 12.8), max(256, ml // 8))
        noverlap = int(nperseg * 0.9375)
        nfft = nperseg * 2

        spec_cache = {}
        all_db = []
        for lbl, d in comps:
            f_s, t_s, Sxx = _sp_spectrogram(d, fs=sr, window='hann',
                                             nperseg=nperseg, noverlap=noverlap,
                                             nfft=nfft, scaling='spectrum', mode='psd')
            Sxx_db = 10 * np.log10(np.maximum(Sxx, 1e-30))
            f_mask = f_s <= 10
            spec_cache[lbl] = (f_s, t_s, Sxx_db, f_mask)
            all_db.append(Sxx_db[f_mask, :].ravel())
        all_db = np.concatenate(all_db)
        vmin, vmax = np.percentile(all_db, 2), np.percentile(all_db, 95)
        if vmax <= vmin: vmax = vmin + 10

        from matplotlib.gridspec import GridSpec
        fig_w = 19.0 / 2.54
        fig_h = fig_w * 0.60
        fig = plt.figure(figsize=(fig_w, fig_h), dpi=150)
        gs = GridSpec(3, 2, fig, hspace=0.10, wspace=0.06,
                       top=0.90, bottom=0.12, left=0.06, right=0.92)

        for i, (lbl, d) in enumerate(comps):
            # 左列: 带通包络
            ax_l = fig.add_subplot(gs[i, 0])
            env = np.abs(hilbert(d))
            env_s = np.convolve(env, np.ones(60)/60, mode='same')
            ax_l.fill_between(t_axis, 0, env_s / NP, color='#d62728', alpha=0.18, lw=0)
            ax_l.plot(t_axis, env_s / NP, '#d62728', lw=0.6)
            ax_l.plot(t_axis, d / NP, '#555555', lw=0.25, alpha=0.5)
            ax_l.text(0.985, 0.88, lbl, transform=ax_l.transAxes,
                       fontsize=8, fontweight='bold', va='top', ha='right')
            ax_l.set_ylim(-2.5, 2.5); ax_l.axhline(0, color='gray', lw=0.3)
            for ph, clr in [('P', '#4daf4a'), ('S', '#ff9848'),
                            ('T-P', '#4d9fff'), ('T-S', '#ff4d4d')]:
                if ph in arrivals:
                    ax_l.axvline(arrivals[ph], color=clr, lw=0.6, ls='--', alpha=0.7)
            ax_l.set_xlim(0, t_axis[-1])
            if i < 2: ax_l.tick_params(labelbottom=False)
            if i == 1: ax_l.set_ylabel('BP 2-8 Hz', fontsize=8)
            for sp in ('top', 'right'): ax_l.spines[sp].set_visible(False)

            # 右列: 频谱图
            ax_r = fig.add_subplot(gs[i, 1])
            f_s, ts_s, Sxx_db, fm = spec_cache[lbl]
            im = ax_r.pcolormesh(ts_s, f_s[fm], Sxx_db[fm, :],
                                  cmap='viridis', vmin=vmin, vmax=vmax, shading='auto')
            ax_r.text(0.985, 0.90, lbl, transform=ax_r.transAxes,
                       fontsize=8, fontweight='bold', color='white', va='top', ha='right')
            ax_r.set_xlim(ax_l.get_xlim()); ax_r.set_ylim(0, 10)
            if i < 2: ax_r.tick_params(labelbottom=False)

        axes = fig.get_axes()
        axes[-2].set_xlabel('Time from origin (s)', fontsize=8)
        axes[-1].set_xlabel('Time from origin (s)', fontsize=8)

        title = (f'{args.network}.{args.station}  |  '
                 f'dist={dist_deg:.1f}°  AZ={az:.1f}°  BAZ={baz:.1f}°')
        fig.suptitle(title, fontsize=9, fontweight='bold')

        cbar_ax = fig.add_axes([0.94, 0.12, 0.012, 0.78])
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.set_label('dB', fontsize=7); cbar.ax.tick_params(labelsize=6)

        img_path = os.path.join(out_dir, f'{args.station}_bandpass_spectrogram.png')
        fig.savefig(img_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"\n✅ 并排图: {img_path}")

    # 台站信息
    with open(os.path.join(out_dir, 'station_info.txt'), 'w') as f:
        f.write(f"Station: {args.network}.{args.station}\n")
        f.write(f"Location: {stla:.4f}°N, {stlo:.4f}°E\n")
        f.write(f"Event: ({evla:.2f}°N, {evlo:.2f}°E), depth={evdp}km\n")
        f.write(f"Dist: {dist_deg:.2f}°  AZ={az:.1f}°  BAZ={baz:.1f}°\n")
        f.write(f"Bandpass: {args.fmin}-{args.fmax} Hz\n")
        for ph, t_arr in sorted(arrivals.items()):
            f.write(f"  {ph}: {t_arr:.1f}s\n")
    print(f"✅ 台站信息: {os.path.join(out_dir, 'station_info.txt')}")
    print(f"\n完成. 输出目录: {out_dir}/")


if __name__ == '__main__':
    main()
