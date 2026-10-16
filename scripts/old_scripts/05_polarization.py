#!/usr/bin/env python3
"""
05_polarization.py — T波偏振分析 (T-P/T-S/全T波三窗口)
用法:
  python3 05_polarization.py mid_atlantic_2014 --sta P61A
"""
import sys, os, glob, argparse, yaml
import numpy as np

CODE_ROOT = "/home/philiay/shared/codex"
sys.path.insert(0, CODE_ROOT)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator

from obspy import read
from obspy.geodetics.base import gps2dist_azimuth, locations2degrees
from lib.arrivals import compute_theoretical_arrivals, dist_az, compute_snr
from lib.rotation import rotate_zne_to_zrt, energy_ratios, amplitude_rose

CONFIG = yaml.safe_load(open(f"{CODE_ROOT}/config/events.yaml"))
FMIN, FMAX = 2.0, 8.0
OUT_DIR = os.path.join(CODE_ROOT, "case")
PAD_AFTER  = 40
PAD_BEFORE = 60  # 截窗前 padding (s)
WIN_DEFAULT = 25  # ±25s
WIN_NARROW  = 12  # ±12s 当重叠时


def load_component(station, network, comp, loc='*'):
    """加载已处理的 proc 文件，返回 (data, delta, t)
    用 glob 匹配 location（支持 00/10/-- 等）"""
    import glob as _glob
    pattern = os.path.join(PROC_DIR, f'{network}.{station}.{loc}.{comp}.*_proc.sac')
    matches = _glob.glob(pattern)
    if not matches:
        # 再试双点格式
        pattern = os.path.join(PROC_DIR, f'{network}.{station}..{comp}.*_proc.sac')
        matches = _glob.glob(pattern)
    if matches:
        path = matches[0]
        tr = read(path)[0]
        data = tr.data.astype(np.float64)
        delta = tr.stats.delta
        t = np.arange(len(data)) * delta
        return data, delta, t
    return None, None, None


def find_raw_coords(station, network, loc='*'):
    """从原始 SAC 头段获取台站坐标，用 glob 匹配 location"""
    import glob as _glob
    pattern = os.path.join(RAW_DIR, f'{network}.{station}.{loc}.BHZ.*.SAC')
    matches = _glob.glob(pattern)
    if not matches:
        pattern = os.path.join(RAW_DIR, f'{network}.{station}..BHZ.*.SAC')
        matches = _glob.glob(pattern)
    if matches:
        try:
            tr = read(matches[0], headonly=True)[0]
            return tr.stats.sac.stla, tr.stats.sac.stlo
        except:
            pass
    return None, None


def rotate_zne_to_zrt(z, n, e, baz):
    """ZNE → ZRT"""
    br = np.radians(baz)
    r = n * np.cos(br) + e * np.sin(br)
    t = -n * np.sin(br) + e * np.cos(br)
    return z, r, t


def compute_energy_ratios(zw, rw, tw):
    """未归一化能量比"""
    ez = np.sum(zw**2); er = np.sum(rw**2); et = np.sum(tw**2)
    eh = er + et
    return ez, er, et, eh


def rose_full(theta, weights, nbins=36):
    """全圆玫瑰图直方图 (θ ∈ [-π, π])"""
    bins = np.linspace(-np.pi, np.pi, nbins + 1)
    centers = (bins[:-1] + bins[1:]) / 2
    hist = np.zeros(nbins)
    for i in range(nbins):
        m = (theta >= bins[i]) & (theta < bins[i+1])
        hist[i] = np.sum(weights[m])
    if np.max(hist) > 0:
        hist = hist / np.max(hist)
    return centers, hist


def pca_3d(z, r, t):
    """3分量去均值 PCA，返回 (linearity, planarity, λ1, λ2, λ3)"""
    # 去均值
    z0 = z - np.mean(z)
    r0 = r - np.mean(r)
    t0 = t - np.mean(t)
    # 3×N 矩阵
    X = np.vstack([z0, r0, t0])
    # 协方差矩阵 (3×3)
    C = X @ X.T / X.shape[1]
    # 特征值分解
    eigvals = np.linalg.eigvalsh(C)  # 升序
    λ3, λ2, λ1 = eigvals  # λ1 ≥ λ2 ≥ λ3
    # 指标
    denom = λ1 if λ1 > 0 else 1e-30
    linearity = 1 - λ2 / denom
    planarity = 1 - 2 * λ3 / (λ1 + λ2) if (λ1 + λ2) > 0 else 0
    return linearity, planarity, λ1, λ2, λ3


def detect_channel_components(station, network):
    """自动检测可用的水平分量 (优先 BHN/BHE，回退 BH1/BH2)"""
    # 检查 BHN/BHE
    z, _, _ = load_component(station, network, 'BHZ')
    n, _, _ = load_component(station, network, 'BHN')
    e, _, _ = load_component(station, network, 'BHE')
    if z is not None and n is not None and e is not None:
        return 'BHZ', 'BHN', 'BHE', None
    # 检查 BH1/BH2
    h1, _, _ = load_component(station, network, 'BH1')
    h2, _, _ = load_component(station, network, 'BH2')
    if z is not None and h1 is not None and h2 is not None:
        return 'BHZ', 'BH1', 'BH2', None
    print(f"  ⚠ 未找到完整三分量: {station}")
    return None, None, None, None


def detect_orientation(comp_h1, comp_h2):
    """判断水平分量名称 → (cmpaz_H1, cmpaz_H2)"""
    if comp_h1 == 'BHN':
        return 0.0, 90.0  # BHN=N, BHE=E
    elif comp_h1 == 'BH1':
        # BH1/BH2 的方位角需从 SACPZ 或头段读取，默认 N/E
        return 0.0, 90.0
    return 0.0, 90.0


def run(station, network='TA'):
    """主分析函数"""
    
    # ── 0. 检测分量 ──
    print(f"\n{'='*60}")
    print(f"  {network}.{station} 偏振分析")
    print(f"{'='*60}")
    
    comp_z, comp_n, comp_e, _ = detect_channel_components(station, network)
    if comp_z is None:
        return
    
    # ── 1. 加载三分量 ──
    z, delta, t = load_component(station, network, comp_z)
    n, _, _ = load_component(station, network, comp_n)
    e, _, _ = load_component(station, network, comp_e)
    
    min_len = min(len(z), len(n), len(e))
    z, n, e = z[:min_len], n[:min_len], e[:min_len]
    t = t[:min_len]
    delta = delta
    
    # ── 2. 计算到时 ──
    stla, stlo = find_raw_coords(station, network)
    if stla is None:
        print(f"  ❌ 无法获取台站坐标")
        return
    
    dist_deg = locations2degrees(EVLA, EVLO, stla, stlo)
    arrivals = _compute_arrivals(stla, stlo, dist_deg)
    tTP = arrivals.get('T-P', 0)
    tTS = arrivals.get('T-S', 0)
    tP  = arrivals.get('P', 0)
    tS  = arrivals.get('S', 0)
    
    if tTP == 0 or tTS == 0:
        print(f"  ❌ 无法计算 T-P/T-S 到时")
        return
    
    # 方位角
    _, az, baz = gps2dist_azimuth(EVLA, EVLO, stla, stlo)
    baz_deg = baz
    
    print(f"  AZ={az:.1f}°  BAZ={baz:.1f}°")
    print(f"  P={tP:.0f}s  S={tS:.0f}s  T-P={tTP:.0f}s  T-S={tTS:.0f}s")
    
    # ── 3. 数据已是 2-8 Hz BPF（来自 proc.sac），直接使用 ──
    # (已在 process_obs.py 中去响应+滤波，不再二次滤波)
    
    # ── 4. ZNE → ZRT ──
    z_rt, r, t_comp = rotate_zne_to_zrt(z, n, e, baz_deg)
    
    # ── 5. 动态窗口：噪声水平 → SNR<1.5 截断 ──
    # 噪声窗: P 波之前
    noise_start = 100
    noise_end = min(tP - 50, 250) if tP > 100 else 100
    if noise_end <= noise_start:
        noise_start, noise_end = 50, 150
    idx_noise = (t >= noise_start) & (t <= noise_end)
    noise_rms = np.sqrt(np.mean(r[idx_noise]**2 + t_comp[idx_noise]**2))
    print(f"  噪声窗: {noise_start:.0f}-{noise_end:.0f}s  RMS_horiz={noise_rms:.3e}")
    
    env_total = np.sqrt(z_rt**2 + r**2 + t_comp**2)
    # 平滑总能量包络 (20s 滑动 RMS)
    win_len = int(20 / delta)
    env_smooth = np.full_like(env_total, np.nan)
    half = win_len // 2
    for i in range(len(env_total)):
        i0 = max(0, i - half)
        i1 = min(len(env_total), i + half + 1)
        env_smooth[i] = np.sqrt(np.mean(env_total[i0:i1]**2))
    
    # 找 SNR<1.5 的时间，最长不超过 120s
    snr_thresh = 1.5
    max_after = 120.0
    below_thresh = env_smooth < snr_thresh * noise_rms
    
    # 台站特定窗口覆盖
    station_overrides = {
        'HRV': {'tp_end': 1800, 'ts_end': 1800},
    }
    
    # T-P 窗口: [tTP-50, min(SNR<1.5, tTP+120, 台站覆盖)]
    tp_start_snr = tTP - 50
    after_tp = t > tTP
    idx_tp_end = np.where(below_thresh & after_tp)[0]
    tp_end_snr = t[idx_tp_end[0]] if len(idx_tp_end) > 0 else t[-1]
    tp_end_snr = min(tp_end_snr, tTP + max_after)
    if station in station_overrides and 'tp_end' in station_overrides[station]:
        tp_end_snr = station_overrides[station]['tp_end']
    
    # T-S 窗口: [tTS-50, min(SNR<1.5, tTS+120, 台站覆盖)]
    ts_start_snr = tTS - 50
    after_ts = t > tTS
    idx_ts_end = np.where(below_thresh & after_ts)[0]
    ts_end_snr = t[idx_ts_end[0]] if len(idx_ts_end) > 0 else t[-1]
    ts_end_snr = min(ts_end_snr, tTS + max_after)
    if station in station_overrides and 'ts_end' in station_overrides[station]:
        ts_end_snr = station_overrides[station]['ts_end']
    
    print(f"  T-P窗口: {tp_start_snr:.0f}–{tp_end_snr:.0f}s  ({tp_end_snr-tTP:.0f}s after T-P)")
    print(f"  T-S窗口: {ts_start_snr:.0f}–{ts_end_snr:.0f}s  ({ts_end_snr-tTS:.0f}s after T-S)")
    
    windows = [
        ('T-P', tp_start_snr, tp_end_snr, tTP),
        ('T-S', ts_start_snr, ts_end_snr, tTS),
    ]
    
    # ── 6. 各窗口计算 ──
    results = {}
    for wname, ws, we, wc in windows:
        idx = (t >= ws) & (t <= we)
        zw = z_rt[idx]; rw = r[idx]; tw = t_comp[idx]
        
        # 能量比
        ez, er, et, eh = compute_energy_ratios(zw, rw, tw)
        er_et = er/et if et > 0 else 0
        ez_eh = ez/eh if eh > 0 else 0
        
        # 玫瑰图
        theta = np.arctan2(tw, rw)
        weights = rw**2 + tw**2
        rose_c, rose_h = rose_full(theta, weights)
        
        # PCA
        lin, pla, l1, l2, l3 = pca_3d(zw, rw, tw)
        
        results[wname] = {
            'ws': ws, 'we': we, 'wc': wc,
            'er_et': er_et, 'ez_eh': ez_eh,
            'log_er_et': np.log10(er_et) if er_et > 0 else -999,
            'rose_c': rose_c, 'rose_h': rose_h,
            'linearity': lin, 'planarity': pla,
            'λ1': l1, 'λ2': l2, 'λ3': l3,
        }
        
        print(f"\n  [{wname}] {ws:.0f}-{we:.0f}s  (center={wc:.0f}s)")
        print(f"    ER/ET={er_et:.2f}  EZ/EH={ez_eh:.3f}  log10(ER/ET)={np.log10(er_et):+.2f}" if er_et>0 else "")
        print(f"    PCA: linearity={lin:.3f}  planarity={pla:.3f}  λ1={l1:.3e}  λ2={l2:.3e}  λ3={l3:.3e}")
    
    # ── 7. 判断标注 ──
    tp_er = results['T-P']['er_et']
    ts_er = results['T-S']['er_et']
    
    if tp_er > 1:
        tp_judge = '偏 radial (R>T)'
    elif tp_er < 0.8:
        tp_judge = '偏 tangential (T>R)'
    else:
        tp_judge = 'ambiguous (R≈T)'
    
    if ts_er < 1:
        ts_judge = '偏 tangential (T>R)'
    elif ts_er > 1.2:
        ts_judge = '偏 radial (R>T)'
    else:
        ts_judge = 'ambiguous (R≈T)'
    
    print(f"\n  T-P判断: {tp_judge}")
    print(f"  T-S判断: {ts_judge}")
    
    # ── 8. 绘图 ──
    x_min = tTP - PAD_BEFORE
    x_max = tTS + PAD_AFTER
    
    fig = plt.figure(figsize=(14, 10))
    gs = gridspec.GridSpec(2, 2, figure=fig,
                           hspace=0.30, wspace=0.30,
                           height_ratios=[2, 1.5],
                           left=0.07, right=0.96, top=0.93, bottom=0.08)
    
    # ── Row 0: Z/R/T Waveforms — each window in its column ──
    for iw, (wname, wdata) in enumerate(results.items()):
        ax_wf = fig.add_subplot(gs[0, iw])
        
        ws, we = wdata['ws'], wdata['we']
        idx_win = (t >= ws) & (t <= we)
        tw_win = t[idx_win]
        zw_win = z_rt[idx_win]
        rw_win = r[idx_win]
        tw_win2 = t_comp[idx_win]
        
        # 窗口内归一化
        max_win = max(np.max(np.abs(zw_win)), np.max(np.abs(rw_win)),
                     np.max(np.abs(tw_win2)))
        if max_win > 0:
            z_plot = zw_win / max_win
            r_plot = rw_win / max_win
            t_plot = tw_win2 / max_win
        else:
            z_plot, r_plot, t_plot = zw_win, rw_win, tw_win2
        
        offset = 3
        ax_wf.plot(tw_win, z_plot + offset, 'k-', lw=0.5, label='Z')
        ax_wf.plot(tw_win, r_plot, 'g-', lw=0.5, label='R')
        ax_wf.plot(tw_win, t_plot - offset, 'r-', lw=0.5, label='T')
        for yv in [offset, 0, -offset]:
            ax_wf.axhline(y=yv, color='gray', lw=0.2, alpha=0.3)
        
        # 窗口背景色
        ax_wf.axvspan(ws, we, color='#1f77b4' if wname=='T-P' else '#d62728', alpha=0.04)
        
        # ── T-P / T-S 到时虚线标注 ──
        dot_phases = []
        if wname == 'T-P':
            dot_phases = [('T-P', tTP, '#1f77b4')]
        else:
            dot_phases = [('T-S', tTS, '#d62728')]
        
        for ph_name, ph_time, ph_color in dot_phases:
            if ws <= ph_time <= we:
                # 竖虚线贯穿三个分量
                ax_wf.axvline(x=ph_time, color=ph_color, lw=0.8, ls='--', alpha=0.6, zorder=8)
                ax_wf.text(ph_time, ax_wf.get_ylim()[1] * 0.92,
                          f'{ph_name}={ph_time:.0f}s',
                          fontsize=7, fontweight='bold', color=ph_color,
                          ha='center', va='top',
                          bbox=dict(boxstyle='round,pad=0.2', fc='white', alpha=0.7))
        
        ax_wf.set_xlim(ws, we)
        ax_wf.set_ylabel('Amplitude (norm)', fontsize=9)
        ax_wf.set_title(f'{wname} ({ws:.0f}-{we:.0f}s)', fontsize=10, fontweight='bold')
        ax_wf.legend(fontsize=7, ncol=3, loc='upper right')
        ax_wf.tick_params(labelsize=7)
        ax_wf.grid(True, alpha=0.08, linestyle=':')
    
    # ── Row 1: Rose diagrams ──
    for iw, (wname, wdata) in enumerate(results.items()):
        ax_rd = fig.add_subplot(gs[1, iw], polar=True)
        
        centers = wdata['rose_c']
        hist = wdata['rose_h']
        width = np.pi / 36
        
        ax_rd.bar(centers, hist, width=width, bottom=0.1,
                  color='#1f77b4' if wname == 'T-P' else '#d62728',
                  alpha=0.75, edgecolor='w', lw=0.3)
        
        # Reference lines (full circle)
        for ang, lbl, clr in [(0, 'R', 'green'), (np.pi/2, 'T', 'red'),
                              (np.pi, 'R', 'green'), (-np.pi/2, 'T', 'red')]:
            ax_rd.plot([ang, ang], [0, 1.2], color=clr, lw=0.8, alpha=0.5)
            # 标签只标一次（避免重叠）
            if ang <= 0:
                ax_rd.text(ang, 1.35, lbl, fontsize=8, fontweight='bold',
                          ha='center', color=clr)
        
        ax_rd.set_title(
            f'{wname} ({wdata["ws"]:.0f}-{wdata["we"]:.0f}s)\n'
            f'ER/ET={wdata["er_et"]:.2f}  EZ/EH={wdata["ez_eh"]:.3f}',
            fontsize=8, fontweight='bold', pad=20)
        ax_rd.tick_params(labelsize=6)
        ax_rd.set_theta_zero_location('E')
        ax_rd.set_theta_direction(-1)
        ax_rd.set_ylim(0, 1.5)
    
    # ── 图标题 ──
    tp_half = (tp_end_snr - tp_start_snr) / 2
    ts_half = (ts_end_snr - ts_start_snr) / 2
    fig.suptitle(
        f'{network}.{station}  |  T-P: {tp_judge}  |  T-S: {ts_judge}  |  '
        f'{FMIN:.0f}-{FMAX:.0f} Hz\n'
        f'T-P={tTP:.0f}s → {tp_end_snr:.0f}s (±{tp_half:.0f}s)  |  '
        f'T-S={tTS:.0f}s → {ts_end_snr:.0f}s (±{ts_half:.0f}s)',
        fontsize=10, fontweight='bold', y=0.97
    )
    
    os.makedirs(OUT_DIR, exist_ok=True)
    out_png = os.path.join(OUT_DIR, f'polar_analysis_{station}.png')
    fig.savefig(out_png, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"\n✅ {out_png}")
    
    # ── 返回结果用于后续分析 ──
    return results


def _compute_arrivals(stla, stlo, dist_deg):
    """计算 P/S (TauP) 和 T-P/T-S 理论到时"""
    from lib.arrivals import compute_theoretical_arrivals as cta
    from t_phase_analysis.traveltime_profile import get_track, tt_t
    import sys as _sys
    _sys.path.insert(0, '/mnt/f/sharedate/codex/t_phase_analysis')
    return cta(EVLO, EVLA, EVDP, stlo, stla)


def plot_batch_roses(stations, network, out_dir):
    """
    批量: 每个台站一个全T波窗偏振玫瑰图，合并为一张图。
    SNR 前 N 个台站的简化偏振分析。
    """
    import sys as _sys
    _sys.path.insert(0, '/mnt/f/sharedate/codex/t_phase_analysis')
    from t_phase_analysis.traveltime_profile import get_track, tt_t

    fig, axes = plt.subplots(1, len(stations), figsize=(5*len(stations), 5),
                             subplot_kw={'polar': True})
    if len(stations) == 1:
        axes = [axes]

    for iax, sta in enumerate(stations):
        ax = axes[iax]
        print(f"\n{'='*50}")
        print(f"  {network}.{sta}")
        print(f"{'='*50}")

        # 加载三分量
        comp_z, comp_n, comp_e, _ = detect_channel_components(sta, network)
        if comp_z is None:
            ax.text(0, 0, f'{sta}: no 3C', ha='center', fontsize=8)
            continue

        z, delta, t = load_component(sta, network, comp_z)
        n, _, _ = load_component(sta, network, comp_n)
        e, _, _ = load_component(sta, network, comp_e)
        if z is None:
            ax.text(0, 0, f'{sta}: no data', ha='center', fontsize=8)
            continue

        min_len = min(len(z), len(n), len(e))
        z, n, e = z[:min_len], n[:min_len], e[:min_len]
        t = t[:min_len]

        # SAC 时间偏移
        tr = read(f"{PROC_DIR}/{network}.{sta}..BHZ.*_proc.sac")[0]
        t_offset = tr.stats.sac.b - tr.stats.sac.o
        t = t + t_offset  # fix previously loaded t (was from 0)

        # 到时
        stla, stlo = find_raw_coords(sta, network)
        dist_deg = locations2degrees(EVLA, EVLO, stla, stlo)
        arrs = _compute_arrivals(stla, stlo, dist_deg)
        tTP = arrs.get('T-P', 0)
        tTS = arrs.get('T-S', 0)
        if tTP == 0 or tTS == 0:
            ax.text(0, 0, f'{sta}: no T arr', ha='center', fontsize=8); continue

        # BAZ
        _, _, baz = gps2dist_azimuth(EVLA, EVLO, stla, stlo)

        # ZNE → ZRT
        z_rt, r, tc = rotate_zne_to_zrt(z, n, e, baz)

        # 全T波窗: T-P - 30s 到 T-S + 120s
        ws = max(tTP - 30, t[0])
        we = min(tTS + 120, t[-1])
        idx = (t >= ws) & (t <= we)
        rw = r[idx]; tcw = tc[idx]

        if len(rw) < 10:
            ax.text(0, 0, f'{sta}: short win', ha='center', fontsize=8); continue

        # 玫瑰图
        theta = np.arctan2(tcw, rw)
        weights = rw**2 + tcw**2
        centers, hist = rose_full(theta, weights)

        width = np.pi / 36
        ax.bar(centers, hist, width=width, bottom=0.1,
               color='#1f77b4', alpha=0.75, edgecolor='w', lw=0.3)

        # R/T 参考线
        for ang, lbl, clr in [(0,'R','green'),(np.pi/2,'T','red'),
                              (np.pi,'R','green'),(-np.pi/2,'T','red')]:
            ax.plot([ang,ang], [0,1.2], color=clr, lw=0.8, alpha=0.5)
            if ang <= 0:
                ax.text(ang, 1.35, lbl, fontsize=8, fontweight='bold',
                        ha='center', color=clr)

        # BAZ 虚线
        baz_rad = 0
        ax.plot([baz_rad, baz_rad], [0, 1.4], 'k--', lw=1.0, alpha=0.5)
        ax.plot([baz_rad+np.pi, baz_rad+np.pi], [0, 1.4], 'k--', lw=1.0, alpha=0.5)

        # 能量比
        er_et = np.sum(rw**2)/np.sum(tcw**2) if np.sum(tcw**2) > 0 else 0
        ez_eh = np.sum(z_rt[idx]**2)/(np.sum(rw**2)+np.sum(tcw**2)) if (np.sum(rw**2)+np.sum(tcw**2)) > 0 else 0

        ax.set_title(f'{sta}\\nfull T-wave ({ws:.0f}-{we:.0f}s)\\n'
                     f'ER/ET={er_et:.2f}  EZ/EH={ez_eh:.3f}',
                     fontsize=9, fontweight='bold', pad=20)
        ax.tick_params(labelsize=6)
        ax.set_theta_zero_location('E')
        ax.set_theta_direction(-1)
        ax.set_ylim(0, 1.5)

        print(f"  全T波窗: {ws:.0f}-{we:.0f}s  ER/ET={er_et:.2f}  EZ/EH={ez_eh:.3f}")

    fig.suptitle(f'Top {len(stations)} SNR stations — Polarization (full T-wave window)',
                 fontsize=11, fontweight='bold', y=1.02)
    plt.tight_layout()
    out_png = os.path.join(out_dir, f'polar_top{len(stations)}_fulltwave.png')
    fig.savefig(out_png, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"\n✅ {out_png}")


def main():
    p = argparse.ArgumentParser(description='T波偏振分析')
    p.add_argument('event', help='事件key (见 config/events.yaml)')
    p.add_argument('--sta', help='单个台站名')
    p.add_argument('--stations', help='多个台站，逗号分隔 (如 P61A,J59A,K60A)')
    p.add_argument('--top-n', type=int, default=0,
                   help='自动选 SNR 前 N 个台站做偏振玫瑰图（全T波窗）')
    p.add_argument('--network', default='TA', help='台网')
    args = p.parse_args()

    cfg = CONFIG[args.event]
    global EVLA, EVLO, EVDP
    EVLA, EVLO = cfg['lat'], cfg['lon']
    EVDP = cfg['depth_km']
    global PROC_DIR, RAW_DIR
    case_dir = cfg.get('case_dir', cfg['label'])
    PROC_DIR = f"/home/philiay/shared/case/{case_dir}/sacprocess"
    RAW_DIR  = f"/mnt/f/sharedate/case/{case_dir}/raw"
    out_dir  = f"/home/philiay/shared/case/{case_dir}/polarization"
    os.makedirs(out_dir, exist_ok=True)

    # 确定待处理台站列表
    stations = []
    if args.stations:
        stations = [s.strip() for s in args.stations.split(',') if s.strip()]
    elif args.top_n > 0:
        # 从 proc 目录扫描所有台站，计算 SNR，取前 N
        import glob
        bhz_files = sorted(glob.glob(f"{PROC_DIR}/*.BHZ*_proc.sac"))
        scored = []
        for fn in bhz_files:
            try:
                sta = os.path.basename(fn).split('.')[1]
                if sta in [s['sta'] for s in scored]:
                    continue
                tr = read(fn)[0]
                data = tr.data.astype(np.float64)
                delta = tr.stats.delta
                t_offset = tr.stats.sac.b - tr.stats.sac.o
                t = np.arange(len(data)) * delta + t_offset

                # 简单 SNR: T-P 窗口 / P 波前噪声
                stla, stlo = find_raw_coords(sta, args.network)
                if stla is None: continue
                dist_deg = locations2degrees(EVLA, EVLO, stla, stlo)
                arrs = _compute_arrivals(stla, stlo, dist_deg)
                tTP = arrs.get('T-P', 0)
                tP = arrs.get('P', 0)
                if tTP == 0 or tP == 0: continue

                # 噪声窗: P波前
                n_start = 50; n_end = max(int(tP - 30), n_start + 10)
                # 信号窗: T-P 前后
                sig_start = int(tTP - 30)
                sig_end   = int(tTP + 90)
                if sig_start < 0 or sig_end > len(data): continue

                noise_rms = np.std(data[n_start:n_end])
                sig_rms = np.std(data[sig_start:sig_end])
                snr = sig_rms / noise_rms if noise_rms > 0 else 0

                _, az, baz = gps2dist_azimuth(EVLA, EVLO, stla, stlo)
                scored.append({'sta': sta, 'snr': snr, 'az': az, 'baz': baz})
                print(f"  {sta}: SNR={snr:.1f} AZ={az:.0f}°")
            except:
                continue

        scored.sort(key=lambda x: x['snr'], reverse=True)
        stations = [s['sta'] for s in scored[:args.top_n]]
        print(f"\nSNR top {args.top_n}: {stations}")
        for s in scored[:args.top_n]:
            print(f"  {s['sta']}: SNR={s['snr']:.1f} AZ={s['az']:.0f}° BAZ={s['baz']:.0f}°")
    elif args.sta:
        stations = [args.sta]

    if not stations:
        print("请指定 --sta, --stations, 或 --top-n")
        return

    # 逐个处理或批量出图
    if len(stations) == 1 and not args.top_n:
        # 单台: 原有详细模式 (2×2 布局)
        run(args.sta, args.network)
    else:
        # 批量: 简化的全T波窗玫瑰图模式
        plot_batch_roses(stations, args.network, out_dir)

    print(f"\n输出目录: {out_dir}")


if __name__ == '__main__':
    main()
