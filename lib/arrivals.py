#!/usr/bin/env python3
"""到时计算 — TauP P/S + tt_t T-P/T-S + 包络检测"""
import sys, numpy as np
for p in ['/home/sy/codex', '/home/sy/codex/t_phase_analysis']:
    if p not in sys.path: sys.path.insert(0, p)
from obspy.taup import TauPyModel
from obspy.signal.filter import envelope as obs_env
from obspy.signal.trigger import classic_sta_lta
from scipy.ndimage import uniform_filter1d
from geographiclib.geodesic import Geodesic
from t_phase_analysis.traveltime_profile import get_track, tt_t

TAUP = TauPyModel(model='iasp91')
GEOD = Geodesic.WGS84


def compute_theoretical_arrivals(evlo, evla, evdp, stlo, stla):
    """返回 {P, S, T-P, T-S, T-Love, land_km, ocean_km} 秒"""
    from obspy.geodetics.base import locations2degrees, degrees2kilometers
    dist_deg = locations2degrees(evla, evlo, stla, stlo)
    results = {}
    
    # P/S via TauP
    try:
        arrs = TAUP.get_travel_times(source_depth_in_km=evdp, distance_in_degree=dist_deg,
                                      phase_list=['P','S'])
        for a in arrs:
            if a.name == 'P' and 'P' not in results: results['P'] = a.time
            elif a.name == 'S' and 'S' not in results: results['S'] = a.time
    except: pass
    
    # T-P/T-S via tt_t
    try:
        track = get_track(evlo, evla, stlo, stla, 0.02, ref_el=-1000,
                         plot_map=False, plot_elevation=False)
        tr = tt_t(track, ref_el=-1000, source_depth_km=evdp, ignore_small_islands=True)
        if tr:
            results['T-P'] = tr['tt_t']
            results['T-S'] = tr['tt_ocean'] + tr['tt_land_s']
            land_deg = tr.get('dist_land', 0)
            ocean_deg = tr.get('dist_ocean', 0)
            results['land_km'] = degrees2kilometers(land_deg)
            results['ocean_km'] = degrees2kilometers(ocean_deg)
            results['T-Love'] = tr['tt_ocean'] + results['land_km'] / 3.6
    except: pass
    
    return results


def detect_arrival_envelope(data, delta, theory_s, search_half=120, celerity_max=1.8,
                             dist_km=None):
    """
    包络峰值检测T波到时。
    在 theory_s ± search_half 内找包络峰值, 用celerity检查。
    """
    env = obs_env(data.astype(np.float64))
    t = np.arange(len(data)) * delta
    ws = max(0, theory_s - search_half); we = min(t[-1], theory_s + search_half)
    win = (t >= ws) & (t <= we)
    if not np.any(win): return theory_s
    
    tw, ew = t[win], env[win]
    peak_idx = np.argmax(ew); t_peak = tw[peak_idx]
    
    # Celerity check
    if dist_km is not None:
        celerity = dist_km / t_peak if t_peak > 0 else 999
        if celerity >= celerity_max: return theory_s  # likely P/S, not T
    
    return t_peak


def compute_snr(data, delta, signal_center, noise_center=None, signal_start=None,
                noise_pre=300, noise_post=60, t_offset=0.0):
    """
    T波 SNR (dB) — 频域功率谱法 (Welch PSD, 2-8Hz积分)
    信号窗: [signal_start, signal_center+80]
    噪声窗: [noise_center-noise_pre, noise_center-noise_post] (P波到达前)
    SNR = 10*log10(integrated_power_sig / integrated_power_noise)
    """
    from scipy.signal import welch
    sr = 1.0 / delta
    t = np.arange(len(data)) * delta + t_offset
    
    ss = signal_start if signal_start is not None else noise_center
    s_start = max(0, int((ss - t_offset) / delta))
    s_end = min(len(data), int((signal_center + 80 - t_offset) / delta))
    
    nc = noise_center if noise_center is not None else signal_center
    n_start = max(0, int((nc - noise_pre - t_offset) / delta))
    n_end = max(0, int((nc - noise_post - t_offset) / delta))
    
    if s_end <= s_start or n_end <= n_start:
        return -999
    
    f_sig, psd_sig = welch(data[s_start:s_end], fs=sr, nperseg=min(256, s_end-s_start))
    f_noise, psd_noise = welch(data[n_start:n_end], fs=sr, nperseg=min(256, n_end-n_start))
    
    e_sig = np.sum(psd_sig[(f_sig >= 2) & (f_sig <= 8)])
    e_noise = np.sum(psd_noise[(f_noise >= 2) & (f_noise <= 8)])
    
    if e_noise > 0 and e_sig > 0:
        return 10 * np.log10(e_sig / e_noise)
    return -999


def dist_az(evlo, evla, stlo, stla):
    """返回 (dist_km, az, baz)。

    az  = 震源处指向台站的方位角（geographiclib azi1）—— 从震源出发投影/画射线用这个。
    baz = **台站处指向震源**的方位角（= 到达端 azi2 + 180）—— 旋转 N/E→R/T、标注方位用这个。

    ⚠️ 2026-09-16 修 bug：旧实现 baz = azi1 + 180（震源端方位角取反），不是台站端方位角，
       两者差一个"大圆弧方位角沿路径的旋转量" ≈ Δλ·sinφ_mid（南非 XA 案例实测 6.0–13.0°，
       均值 +9.0°）。用旧值旋转 R/T 会整体偏 ~11°，且与用 gps2dist_azimuth(台,震)[1] 的
       F-K 链不一致（R/T 波形差 15–28%）。
    """
    g = GEOD.Inverse(evla, evlo, stla, stlo)
    az = g['azi1']
    if az < 0: az += 360
    baz = (g['azi2'] + 180) % 360      # 台站端反向方位角（= 台站→震源）
    return g['s12'] / 1000.0, az, baz
