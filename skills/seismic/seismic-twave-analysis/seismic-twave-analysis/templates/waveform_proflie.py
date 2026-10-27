#!/usr/bin/env python3
"""
波形剖面图（独立脚本，从 step1_station_table_and_map.py 拆分出来）

读取 step1 生成的 stations.csv + sacprocess 目录下的波形，
按震中距(dist)排序，每 bin_km 一档，档内多台站波形插值到公共时间轴后叠加(平均)。

输出:
- profile{_waveform}{profile_label}.png

用法:
    python plot_waveform_profile.py puerto_rico
    python plot_waveform_profile.py mindanao_2026 --snr-min 5.0 --bin-km 30
    python plot_waveform_profile.py sw_indian_ridge --regions "10:43:Mainland,43:60:Madagascar"
"""
import os
import sys
import glob
import csv
import argparse
import yaml

import numpy as np
from obspy import read, UTCDateTime
from obspy.signal.filter import envelope as obs_env
from scipy.ndimage import uniform_filter1d as uni_filt
from matplotlib.lines import Line2D as MLine
from matplotlib.ticker import MultipleLocator as MLoc
from numpy.polynomial import Polynomial
import matplotlib.pyplot as plt

CODE_ROOT = "/home/sy/codex"
CONFIG_PATH = os.path.join(CODE_ROOT, "config", "events.yaml")
DEFAULT_SNR_THRESHOLD = 8.0
DEFAULT_BIN_KM = 10
DEFAULT_REGIONS = "43:60:Madagascar"

# 剖面图排除台站（异常噪声/毛刺）
PROFILE_EXCLUDE = {'BITY', 'Z01TG', 'W12MB'}


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


def load_stations_csv(csv_path):
    """读取 step1 生成的 stations.csv"""
    rows = []
    with open(csv_path) as f:
        for row in csv.DictReader((l for l in f if not l.startswith('#'))):
            rows.append({
                'sta': row['sta'],
                'dist': float(row['dist']),
                'land_km': float(row['land_km']),
                'ocean_km': float(row['ocean_km']),
                'tt_p': float(row['tt_p']),
                'tt_s': float(row['tt_s']),
                'tt_tp': float(row['tt_tp']),
                'tt_ts': float(row['tt_ts']),
                'snr': float(row['snr']),
                'lon': float(row['lon']),
                'lat': float(row['lat']),
            })
    return rows


def _parse_origin_from_filename(bname):
    """从 SAC 文件名解析发震时刻 (如 NET.STA..BHZ.D.2013.203.070100.SAC → 2013-07-22T07:01:00)。
    解析失败返回 None。"""
    import re
    from obspy import UTCDateTime
    m = re.search(r'\.(\d{4})\.(\d{3})\.(\d{6})', bname)
    if not m:
        return None
    y, jd, hm = int(m.group(1)), int(m.group(2)), m.group(3)
    return UTCDateTime(year=y, julday=jd,
                       hour=int(hm[0:2]), minute=int(hm[2:4]), second=int(hm[4:6]))


def gen_waveform_profile(rows, evdate, evmw, proc_dir, out_dir, snr_threshold,
                          lon_range=None, profile_label='', bin_km=DEFAULT_BIN_KM,
                          ykey='dist', band_fl=0.5, band_fh=2.0):
    """波形剖面图 — 按 ykey(dist/land_km) 排序，每 bin_km 一档，档内取 SNR 最高台站"""
    resp_dir = os.path.join(proc_dir, "remove_response")
    envelope_dir = os.path.join(proc_dir, "envelope")
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
        if sta not in sta_data or sta in seen or sta in PROFILE_EXCLUDE:
            continue
        seen.add(sta)

        try:
            tr = read(fn)[0]
            delta = tr.stats.delta
            info = sta_data[sta]
            env_npz = None
            if abs(band_fl - 2.0) < 1e-6 and abs(band_fh - 8.0) < 1e-6:
                env_npz = os.path.join(envelope_dir, f'{sta}_envelope_2-8Hz.npz')
            if env_npz is not None and os.path.exists(env_npz):
                # 优先使用预处理阶段保存的包络，保证与四分量图完全一致
                ep = np.load(env_npz)
                data = ep['waveform_z'].astype(np.float64)
                env = ep['envelope_z'].astype(np.float64)
                t = ep['time'].astype(np.float64)
                delta = float(ep['delta'])
                print(f'  {sta}: envelope cache')
            else:
                # 兼容旧事件：没有缓存时才重新计算
                data = tr.data.astype(np.float64)
                p98 = np.percentile(np.abs(data), 98)
                if p98 > 1e6 or not np.isfinite(p98):
                    continue
                sac_b = getattr(tr.stats.sac, 'b', 0)
                sac_o = getattr(tr.stats.sac, 'o', None)
                if sac_o is not None and 0 < abs(sac_o) < 1e6:
                    t = sac_b - sac_o + np.arange(len(data)) * delta
                else:
                    _origin = _parse_origin_from_filename(bname)
                    t = ((tr.stats.starttime - _origin) if _origin is not None else 0.0) + np.arange(len(data)) * delta
                from scipy.signal import butter as _butter, filtfilt as _filtfilt
                _b, _a = _butter(4, [band_fl, band_fh], btype='band', fs=1/delta)
                data = _filtfilt(_b, _a, data)
                env = obs_env(data)
                env = uni_filt(env, size=max(int(3 / delta), 2))
            # 陆地走时对齐: t_land = t - tt_ocean, t=0 = 波上岸
            tt_ocean = info['ocean_km'] / 1.48
            t_land = t - tt_ocean
            traces.append({
                'sta': sta, 'dist': info['dist'], 'land_km': info['land_km'],
                'data': data, 'env': env, 't': t_land, 't_abs': t,
                'delta': delta, 'tt_ocean': tt_ocean,
                # 到时全部转到陆地走时坐标
                'tt_p': (info.get('tt_p', 0) - tt_ocean) if info.get('tt_p', 0) else 0,
                'tt_s': (info.get('tt_s', 0) - tt_ocean) if info.get('tt_s', 0) else 0,
                'tt_tp': info['tt_tp'] - tt_ocean, 'tt_ts': info['tt_ts'] - tt_ocean,
                'snr': info['snr'],
            })
        except Exception:
            continue

    if len(traces) < 3:
        print(f"  ⚠ 波形不足 ({len(traces)}站)，跳过")
        return

    n_raw = len(traces)

    # ── 按 ykey 分箱：每 bin_km 一档，档内只取 SNR 最高的台站 ──
    ylabel = {'dist': 'Epicentral distance', 'land_km': 'Land distance'}[ykey]
    ysfx = {'dist': '', 'land_km': '_land'}[ykey]
    for tr in traces:
        tr['bin'] = round(tr[ykey] / bin_km) * bin_km
    bins = {}
    for tr in traces:
        bins.setdefault(tr['bin'], []).append(tr)

    # 每档选 SNR 最高者
    best = []
    for b, group in sorted(bins.items()):
        g = max(group, key=lambda x: x['snr'])
        g[ykey] = b  # 对齐到档位中心
        best.append(g)
    traces = best
    n = len(traces)
    print(f"  {n_raw} 站 → {n} 个{ylabel}档 (每{bin_km}km, 取SNR最高)")

    t_min = min(t['tt_tp'] for t in traces if t['tt_tp'] > 0) - 30  # T 波时段起点(收紧)
    t_max = max(t['tt_ts'] for t in traces if t['tt_ts'] > 0) + 60   # T 波时段终点(收紧)

    all_p = np.concatenate([t['data'] for t in traces])
    global_p98 = np.percentile(np.abs(all_p), 98) or 1

    # 包络单独归一化：每条 trace 先扣信号窗口前的噪声基线（t < t_min 中位数），
    # 再基于扣基线后的包络算统一 P98 基准。不能复用波形 global_p98：
    # 包络是恒正平滑量，与振荡波形的 P98 尺度不同，混用会把包络压成直线。
    for _tr in traces:
        _pre = _tr['env'][_tr['t'] < t_min]
        _tr['env_base'] = float(np.median(_pre)) if len(_pre) > 0 else 0.0
        _tr['env_c'] = np.maximum(_tr['env'] - _tr['env_base'], 0.0)
    all_env = np.concatenate([t['env_c'] for t in traces])
    global_env_p98 = np.percentile(all_env, 98) or 1
    for _tr in traces:
        _m = (_tr['t'] >= t_min) & (_tr['t'] <= t_max)
        _w = _tr['env_c'][_m]
        _tr['p98_self'] = float(np.percentile(_w, 98)) if len(_w) > 0 else 1.0
    print(f"  Global P98 (waveform) = {global_p98:.2e} | Global P98 (envelope) = {global_env_p98:.2e}")

    fig_h = min(30, max(10, (traces[-1][ykey] - traces[0][ykey]) / 30))
    dist_min = traces[0][ykey]; dist_max = traces[-1][ykey]
    margin = (dist_max - dist_min) * 0.05 + 20

    # 分箱后相邻档间距已固定为 bin_km，振幅按此留出安全余量
    amp = min(30, max(5, bin_km * 0.45))
    amp_env = amp * 3.0          # 包络主图放大：突出 T 波隆起
    amp_env_pb = amp * 1.6       # 包络 per-bin 辅助图
    print(f"  振幅缩放: amp={amp:.1f}km (分箱宽度 {bin_km}km) | envelope amp={amp_env:.1f}km | per-bin amp={amp_env_pb:.1f}km")

    # 仅输出一张纯包络主图（全局 P98，保持台站间相对振幅）
    plot_specs = [
        ('envelope', f'Envelope (global env P98) {band_fl:.1f}-{band_fh:.0f}Hz', '_env', 'env_global', amp_env, False),
    ]

    for mode, tag, suffix, norm_mode, amp_use, overlay in plot_specs:
        fig, ax = plt.subplots(figsize=(16, fig_h * 0.7))
        fig.subplots_adjust(left=0.02)

        for tr in traces:
            y0 = tr[ykey]
            m = (tr['t'] >= t_min) & (tr['t'] <= t_max)
            tw = tr['t'][m]
            if len(tw) == 0:
                continue

            if mode == 'envelope':
                # 包络: 已扣基线 env_c, 归一化基准可选全局 env P98 或每档自身 P98
                if norm_mode == 'env_perbin':
                    ey = (tr['env_c'][m] / tr['p98_self']) * amp_use
                else:
                    ey = (tr['env_c'][m] / global_env_p98) * amp_use
                if overlay:
                    # 主图: 灰色填充包络 + 黑色细波形叠加
                    dw = (tr['data'][m] / global_p98) * amp
                    ax.fill_between(tw, y0, y0 + ey, color='#888888', alpha=0.45, lw=0, zorder=3)
                    ax.plot(tw, y0 + dw, color='black', lw=0.25, alpha=0.7, zorder=4)
                else:
                    ax.plot(tw, y0 + ey, color='black', lw=0.6, alpha=1.0, zorder=4)
            else:
                dy = (tr['data'][m] / global_p98) * amp_use
                ax.plot(tw, y0 + dy, color='black', lw=0.25, alpha=0.8, zorder=4)

        # 理论到时虚线: T-P(蓝) T-S(红)
        for pts, color, lbl in [
            ([(t[ykey], t['tt_tp']) for t in traces if t['tt_tp'] > 0], '#4d9fff', 'T-P'),
            ([(t[ykey], t['tt_ts']) for t in traces if t['tt_ts'] > 0], '#ff4d4d', 'T-S')]:
            if pts:
                pts.sort()
                x, y = [p[1] for p in pts], [p[0] for p in pts]
                if len(x) > 1:
                    p = Polynomial.fit(x, y, 1)
                    xl = np.linspace(min(x), max(x), 200)
                    yl = p(xl)
                    m_ok = (yl >= 0) & (yl <= dist_max + margin)
                    ax.plot(xl[m_ok], yl[m_ok], color=color, lw=1.2, alpha=0.8,
                            dashes=(8, 4), zorder=27)

        ax.set_xlim(t_min, t_max)
        ax.set_ylim(0, dist_max + margin)
        ax.xaxis.set_major_locator(MLoc(200)); ax.xaxis.set_minor_locator(MLoc(50))
        ax.set_xlabel('Time from shore arrival (s)', fontsize=12, fontweight='bold')
        ax.set_ylabel(f'{ylabel} (km)', fontsize=14, fontweight='bold')
        ax.yaxis.set_major_locator(MLoc(200)); ax.yaxis.set_minor_locator(MLoc(100))
        ax.tick_params(labelsize=12)
        ax.set_title(f'{evdate} Mw {evmw}  |  {n} bins ({bin_km}km){profile_label}  |  SNR>{snr_threshold}  |  {tag}',
                     fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.08, linestyle=':')
        leg = [
            MLine([0], [0], color='black', lw=0.8, label=tag),
            MLine([0], [0], color='#4d9fff', lw=1.2, ls='--', label='T-P'),
            MLine([0], [0], color='#ff4d4d', lw=1.2, ls='--', label='T-S'),
        ]
        ax.legend(handles=leg, fontsize=10, loc='upper left', framealpha=0.9, ncol=3)
        fp = os.path.join(out_dir, f'profile{suffix}{profile_label}{ysfx}.png')
        fig.savefig(fp, dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f"  ✓ {fp}")


def main():
    parser = argparse.ArgumentParser(description="波形剖面图（按震中距分箱叠加）")
    parser.add_argument('event', type=str, help="事件名称 (events.yaml 中的 key)")
    parser.add_argument('--snr-min', type=float, default=None,
                        help=f'SNR阈值覆盖 (默认从 config 读取，缺省 {DEFAULT_SNR_THRESHOLD}dB)')
    parser.add_argument('--bin-km', type=float, default=DEFAULT_BIN_KM,
                        help=f'震中距分箱宽度 (km, 默认 {DEFAULT_BIN_KM})')
    parser.add_argument('--regions', type=str, default=DEFAULT_REGIONS,
                        help='lon_min:lon_max:label，逗号分隔多个区域，每个区域出一份剖面图 '
                             f'(默认 "{DEFAULT_REGIONS}")')
    parser.add_argument('--band', type=str, default='0.5-2.0',
                        help='带通频段 (Hz), 默认 0.5-2.0 低频(T波低频衰减慢,远处可见); '
                             '近处可选 2.0-8.0')
    parser.add_argument('--config', type=str, default=CONFIG_PATH,
                        help=f'配置文件路径 (默认: {CONFIG_PATH})')
    args = parser.parse_args()

    try:
        band_fl, band_fh = [float(x) for x in args.band.split('-')]
    except ValueError:
        print(f"错误: --band 格式应为 'fl-fh'，如 0.5-2.0，收到: {args.band}")
        sys.exit(1)

    cfg = load_event_config(args.event, config_path=args.config)
    EVENT_DATE = cfg['date']
    EVENT_MW = cfg['mw']
    case_dir = cfg['case_dir']

    PROC_DIR = f"/mnt/f/sharedate/case/{case_dir}/sacprocess"
    OUT_DIR = f"/mnt/f/sharedate/case/{case_dir}/process_flow/1_station_distribution"

    csv_path = f"{OUT_DIR}/stations.csv"
    if not os.path.exists(csv_path):
        csv_path = f"{OUT_DIR}/stations_v2.csv"
    if not os.path.exists(csv_path):
        print(f"错误: 找不到 stations.csv，请先运行 step1_station_table_and_map.py")
        sys.exit(1)

    rows = load_stations_csv(csv_path)
    SNR_THRESHOLD = args.snr_min if args.snr_min is not None else cfg.get('snr_threshold', DEFAULT_SNR_THRESHOLD)

    print(f"{'='*70}")
    print(f"波形剖面图 | 事件: {args.event} | {len(rows)}站 | SNR≥{SNR_THRESHOLD}dB | bin={args.bin_km}km")
    print(f"{'='*70}\n")

    for spec in args.regions.split(','):
        lo, hi, label = spec.split(':')
        print(f"[{label}] lon {lo}~{hi}°")
        gen_waveform_profile(rows, EVENT_DATE, EVENT_MW, PROC_DIR, OUT_DIR, SNR_THRESHOLD,
                              lon_range=(float(lo), float(hi)), profile_label=f'_{label}',
                              bin_km=args.bin_km, ykey='land_km',
                              band_fl=band_fl, band_fh=band_fh)
        print()


if __name__ == '__main__':
    main()