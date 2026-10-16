#!/usr/bin/env python3
"""
T波衰减分析模板
用法:
  python3 attenuation_analysis.py EVLA EVLO EVDP BASE [--az-min 310 --az-max 350 --snr-min 2]

功能:
  1. 读取 stations.xlsx 中 SNR>=阈值 的台站
  2. 多台站子频段(2-3, 3-4, 4-5, 5-6, 6-8 Hz) T波振幅测量
  3. ln(A) vs land_km 线性拟合 → 衰减系数 α, 估计 Q
  4. ln(A) vs frequency → t* 测量
  5. 输出: 衰减曲线图 + 数据表

理论基础 (Shearer 2019 Ch.6.6):
  A(f, x) = A₀ · exp(-π·f·x / (Q·V))
  ln(A) = ln(A₀) - (π·f / Q·V) · x
  → 斜率 = -π·f / (Q·V) = -α(f)
  → Q = π·f / (α·V)

结合 codex 已有模块:
  - tpef.tpef(): T波窗口自动检测
  - waveform.compute_spectrogram(): 时频谱
  - snr.snr(): 信噪比
"""
import sys, os, glob, argparse
import numpy as np
from scipy.signal import hilbert
from scipy.ndimage import uniform_filter1d
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from obspy import read
from obspy.signal.filter import bandpass as bp_filter
from obspy.taup import TauPyModel
from obspy.geodetics.base import locations2degrees
from geographiclib.geodesic import Geodesic

_codex = '/home/philiay/shared/codex/t_phase_analysis'
if os.path.isdir(_codex): sys.path.insert(0, _codex)
from t_phase_analysis.traveltime_profile import get_track, tt_t

# ─── CLI ───
parser = argparse.ArgumentParser(description='T-wave attenuation analysis')
parser.add_argument('evla', type=float)
parser.add_argument('evlo', type=float)
parser.add_argument('evdp', type=float)
parser.add_argument('base', help='event root directory')
parser.add_argument('--az-min', type=float, default=0)
parser.add_argument('--az-max', type=float, default=360)
parser.add_argument('--snr-min', type=float, default=2.0)
parser.add_argument('--v-land', type=float, default=6.0, help='assumed land P velocity (km/s)')
parser.add_argument('--sofar-v', type=float, default=1.48, help='SOFAR velocity (km/s)')
parser.add_argument('--out', default=None, help='output directory')
args = parser.parse_args()

EVLA, EVLO, EVDP = args.evla, args.evlo, args.evdp
BASE = args.base
PROC_DIR = os.path.join(BASE, 'sacprocess')
OUT_DIR = args.out or os.path.join(BASE, 'process_flow', '4_attenuation')
os.makedirs(OUT_DIR, exist_ok=True)

# ─── 子频段定义 ───
SUBBANDS = [(2, 3), (3, 4), (4, 5), (5, 6), (6, 8)]
SUBBAND_COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

# ─── 台站信息 ───
from openpyxl import load_workbook
wb = load_workbook(os.path.join(BASE, 'stations.xlsx'))
ws = wb.active
stations = []
for row in ws.iter_rows(min_row=4, values_only=True):
    if row[0] is None: continue
    # 新列序: STA(0), AZ(1), DIST(2), LAND(3), OCEAN(4), P(5), S(6), T-P(7), T-S(8), SNR(9), LON(10), LAT(11), BAZ(12)
    sta = str(row[0]).strip()
    az = row[1]; dist = row[2]; land_km = row[3]; ocean_km = row[4]
    lon = row[10]; lat = row[11]
    snr_val = row[9]
    if snr_val is None or snr_val < args.snr_min: continue
    if az is None or not (args.az_min <= az <= args.az_max): continue
    stations.append({
        'sta': sta, 'land_km': land_km, 'ocean_km': ocean_km,
        'dist_km': dist, 'az': az, 'snr': snr_val
    })
stations.sort(key=lambda x: x['land_km'])
print(f'Selected {len(stations)} stations (SNR>={args.snr_min}, az={args.az_min}-{args.az_max})')

# ─── 振幅测量 ───

def measure_twave_amplitude(fn, fmin, fmax, tt_TP, tt_TS):
    """测量T波窗口内子频段包络峰值振幅"""
    try:
        tr = read(fn)[0]
    except Exception:
        return None
    data = tr.data.astype(np.float64)
    delta = tr.stats.delta; sr = 1/delta
    t_off = tr.stats.sac.b - tr.stats.sac.o
    t = np.arange(len(data))*delta + t_off

    # 带通滤波
    data_bp = bp_filter(data, fmin, fmax, sr, corners=4, zerophase=True)

    # T波窗: T-P to T-S+60s
    sig_start, sig_end = tt_TP, tt_TS + 60
    m_sig = (t >= sig_start) & (t <= sig_end)
    if np.sum(m_sig) < 10: return None

    # 包络峰值振幅
    env = uniform_filter1d(np.abs(hilbert(data_bp[m_sig])), size=max(int(2/delta), 1))
    amp_peak = np.max(env)

    # 噪声窗: T-P前300s到T-P前60s
    noise_start, noise_end = tt_TP - 300, tt_TP - 60
    m_noise = (t >= noise_start) & (t <= noise_end)
    if np.sum(m_noise) < 10: return None
    env_noise = uniform_filter1d(np.abs(hilbert(data_bp[m_noise])), size=max(int(2/delta), 1))
    noise_rms = np.std(env_noise)

    return {
        'amp_peak': amp_peak,
        'amp_rms': np.sqrt(np.mean(data_bp[m_sig]**2)),
        'noise_rms': noise_rms,
        'snr_band': amp_peak / noise_rms if noise_rms > 0 else 0,
    }

# ─── 主循环: 测量所有台站×所有子频段 ───
taup = TauPyModel('iasp91')
results = []

for s in stations:
    fn_pat = f"{PROC_DIR}/TA.{s['sta']}.*.BHZ*_proc.sac"
    fns = glob.glob(fn_pat)
    if not fns: continue
    fn = fns[0]

    # 计算到时
    tr = read(fn, headonly=True)[0]
    stla, stlo = tr.stats.sac.stla, tr.stats.sac.stlo
    dd = locations2degrees(EVLA, EVLO, stla, stlo)
    try:
        track = get_track(EVLO, EVLA, stlo, stla, 0.02, ref_el=-1000, plot_map=False, plot_elevation=False)
        tt_r = tt_t(track, ref_el=-1000, source_depth_km=EVDP, ignore_small_islands=True)
        tt_TP = tt_r['tt_t']; tt_TS = tt_r['tt_ocean'] + tt_r['tt_land_s']
    except Exception:
        continue

    row = {'sta': s['sta'], 'land_km': s['land_km'], 'ocean_km': s['ocean_km']}
    for fmin, fmax in SUBBANDS:
        res = measure_twave_amplitude(fn, fmin, fmax, tt_TP, tt_TS)
        if res and res['snr_band'] >= 1.5:
            row[f'amp_{fmin}_{fmax}'] = res['amp_peak']
            row[f'snr_{fmin}_{fmax}'] = res['snr_band']
        else:
            row[f'amp_{fmin}_{fmax}'] = np.nan
            row[f'snr_{fmin}_{fmax}'] = np.nan
    results.append(row)

print(f'Measured {len(results)} stations')

# ─── 衰减拟合: ln(A) vs land_km ───

fig, axes = plt.subplots(len(SUBBANDS), 1, figsize=(10, 3*len(SUBBANDS)), sharex=True)
fig.subplots_adjust(hspace=0.25)

q_results = []
for i, (fmin, fmax) in enumerate(SUBBANDS):
    ax = axes[i]
    key = f'amp_{fmin}_{fmax}'
    land = np.array([r['land_km'] for r in results if not np.isnan(r[key])])
    amp = np.array([r[key] for r in results if not np.isnan(r[key])])

    if len(land) < 3:
        ax.set_title(f'{fmin}-{fmax} Hz: insufficient data ({len(land)} pts)', fontsize=10)
        q_results.append({'band': f'{fmin}-{fmax}', 'Q': np.nan, 'alpha': np.nan, 'n': len(land)})
        continue

    ln_amp = np.log(amp)

    # 线性拟合: ln(A) = intercept + slope * land_km
    A = np.vstack([land, np.ones_like(land)]).T
    slope, intercept = np.linalg.lstsq(A, ln_amp, rcond=None)[0]

    # Q = π·f_center / (|slope| · V_land)
    f_center = (fmin + fmax) / 2
    alpha = abs(slope)  # 1/km
    Q_est = np.pi * f_center / (alpha * args.v_land) if alpha > 1e-10 else np.inf

    ax.scatter(land, ln_amp, s=20, c=SUBBAND_COLORS[i], alpha=0.7, edgecolors='k', linewidths=0.3)
    x_fit = np.linspace(land.min(), land.max(), 100)
    ax.plot(x_fit, intercept + slope*x_fit, '--', c=SUBBAND_COLORS[i], lw=1.5,
            label=f'slope={slope:.4f}/km, Q={Q_est:.0f}')
    ax.set_ylabel('ln(A)', fontsize=9)
    ax.set_title(f'{fmin}-{fmax} Hz  (n={len(land)}, Q≈{Q_est:.0f})', fontsize=10, fontweight='bold')
    ax.legend(fontsize=8); ax.grid(True, alpha=0.15)
    q_results.append({'band': f'{fmin}-{fmax}', 'Q': Q_est, 'alpha': alpha, 'n': len(land), 'slope': slope})

axes[-1].set_xlabel('Land distance (km)', fontsize=10)
fig.suptitle(f'T-wave attenuation: ln(A) vs land distance\n'
             f'V_land={args.v_land} km/s, {len(results)} stations', fontsize=11, fontweight='bold')
out_fn = os.path.join(OUT_DIR, 'attenuation_lnA_vs_land.png')
fig.savefig(out_fn, dpi=200, bbox_inches='tight'); plt.close(fig)
print(f'OK {out_fn}')

# ─── 输出 Q 汇总 ───
print('\n=== Q estimates ===')
print(f'{"Band":>8s} {"N":>4s} {"slope":>10s} {"alpha":>10s} {"Q":>8s}')
for q in q_results:
    print(f'{q["band"]:>8s} {q.get("n",0):>4d} {q.get("slope",0):>10.5f} {q.get("alpha",0):>10.5f} {q["Q"]:>8.0f}')

# ─── t* 测量图 (单台站频谱斜率) ───
# 选取中间距离的一个台站做 t* 示例
if results:
    mid_idx = len(results) // 2
    mid_sta = results[mid_idx]
    fn_pat = f"{PROC_DIR}/TA.{mid_sta['sta']}.*.BHZ*_proc.sac"
    fns = glob.glob(fn_pat)
    if fns:
        tr = read(fns[0])[0]
        data = tr.data.astype(np.float64)
        delta = tr.stats.delta; sr = 1/delta
        t_off = tr.stats.sac.b - tr.stats.sac.o
        t_arr = np.arange(len(data))*delta + t_off
        stla, stlo = tr.stats.sac.stla, tr.stats.sac.stlo
        track = get_track(EVLO, EVLA, stlo, stla, 0.02, ref_el=-1000, plot_map=False, plot_elevation=False)
        tt_r2 = tt_t(track, ref_el=-1000, source_depth_km=EVDP, ignore_small_islands=True)
        tt_TP2 = tt_r2['tt_t']; tt_TS2 = tt_r2['tt_ocean'] + tt_r2['tt_land_s']

        # T波窗内FFT
        m_sig = (t_arr >= tt_TP2) & (t_arr <= tt_TS2+60)
        sig_data = data[m_sig]
        if len(sig_data) > 0:
            from scipy.fft import rfft, rfftfreq
            N = len(sig_data)
            freqs = rfftfreq(N, delta)
            spec = np.abs(rfft(sig_data))
            # 2-8 Hz 范围内拟合 ln|A| vs f
            f_mask = (freqs >= 2) & (freqs <= 8)
            f_fit = freqs[f_mask]
            ln_spec = np.log(spec[f_mask] + 1e-30)

            A2 = np.vstack([f_fit, np.ones_like(f_fit)]).T
            slope_tstar, intercept_tstar = np.linalg.lstsq(A2, ln_spec, rcond=None)[0]
            tstar = -slope_tstar / np.pi

            fig2, ax2 = plt.subplots(figsize=(8, 4))
            ax2.plot(f_fit, ln_spec, 'k-', lw=0.8, alpha=0.6)
            ax2.plot(f_fit, intercept_tstar + slope_tstar*f_fit, 'r--', lw=1.5,
                     label=f't*={tstar:.2f}s (slope={slope_tstar:.3f})')
            ax2.set_xlabel('Frequency (Hz)', fontsize=10)
            ax2.set_ylabel('ln|A(f)|', fontsize=10)
            ax2.set_title(f't* measurement: {mid_sta["sta"]} (land={mid_sta["land_km"]:.0f}km)',
                         fontsize=10, fontweight='bold')
            ax2.legend(fontsize=9); ax2.grid(True, alpha=0.15)
            out_fn2 = os.path.join(OUT_DIR, f'tstar_{mid_sta["sta"]}.png')
            fig2.savefig(out_fn2, dpi=200, bbox_inches='tight'); plt.close(fig2)
            print(f'OK {out_fn2}')
            print(f't* = {tstar:.3f} s for {mid_sta["sta"]}')
