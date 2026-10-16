#!/usr/bin/env python3
"""
偏振分析验证：对 K62A 的 P 波和 S 波做短时窗偏振分析
理论预期：
  P波: Z-R 平面线性极化，偏振方向指向震源 (BAZ)
  S波: 水平面径向极化 (R分量为主)
"""
import sys, os
sys.path.insert(0, '/mnt/f/sharedate/codex/t_phase_analysis')
sys.path.insert(0, '/mnt/f/sharedate/T_Phase_Analysis-main')
sys.path.insert(0, '/mnt/f/sharedate/case/code')
import numpy as np
from obspy import read
from obspy.signal.filter import bandpass
from obspy.geodetics.base import locations2degrees, gps2dist_azimuth
from obspy.taup import TauPyModel
from geographiclib.geodesic import Geodesic
from t_phase_analysis.traveltime_profile import get_track, tt_t
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D
from scipy.signal import hilbert
from scipy.ndimage import uniform_filter1d

# ── 参数 ──
STA, NET = 'K62A', 'TA'
EVLA, EVLO = 18.99, -66.87
EVDP = 10.0
PROC_DIR = "/mnt/f/sharedate/case/2014-01-13-mw63-process/sacprocess"
OUT_DIR = "/mnt/f/sharedate/case/2014-01-13-mw63-process/waveform_contrast"
FMIN, FMAX = 2.0, 8.0

# ── 加载数据（去仪器响应后，未滤波）──
RESP_DIR = "/mnt/f/sharedate/case/2014-01-13-mw63-process/sacprocess/remove_response"
trs = {}
for ch in ['BHZ','BHN','BHE']:
    fn = f"{RESP_DIR}/{NET}.{STA}..{ch}.M.2014.013.035104_resp.sac"
    trs[ch] = read(fn)[0]
t0 = max(tr.stats.starttime for tr in trs.values())
t1 = min(tr.stats.endtime for tr in trs.values())
for ch in trs: trs[ch].trim(t0, t1)
delta = trs['BHZ'].stats.delta
sac = trs['BHZ'].stats.sac
stla, stlo = sac.stla, sac.stlo
t_off = sac.b - sac.o
t = np.arange(len(trs['BHZ'].data)) * delta + t_off

# ── 理论到时 ──
taup = TauPyModel(model='iasp91')
dd = locations2degrees(EVLA, EVLO, stla, stlo)
arrs = taup.get_travel_times(source_depth_in_km=EVDP, distance_in_degree=dd, phase_list=['P','S'])
tt_P = next((a.time for a in arrs if a.name=='P'), None)
tt_S = next((a.time for a in arrs if a.name=='S'), None)
_, az, baz = gps2dist_azimuth(EVLA, EVLO, stla, stlo)
print(f"K62A: P={tt_P:.0f}s  S={tt_S:.0f}s  AZ={az:.1f}°  BAZ={baz:.1f}°")

# ── 无滤波（直接用原始速度波形）──
z = trs['BHZ'].data.astype(np.float64)
n = trs['BHN'].data.astype(np.float64)
e = trs['BHE'].data.astype(np.float64)

br = np.radians(baz)
r = n*np.cos(br) + e*np.sin(br)
tc = -n*np.sin(br) + e*np.cos(br)

# ── 定义时窗 ──
windows = [
    ('P波', tt_P, 15.0),   # ±15s around P
    ('S波', tt_S, 15.0),   # ±15s around S
]

# ── PCA函数 ──
def pca_3d(z, r, t):
    z0 = z - np.mean(z); r0 = r - np.mean(r); t0 = t - np.mean(t)
    X = np.vstack([z0, r0, t0])
    C = X @ X.T / X.shape[1]
    eigvals = np.linalg.eigvalsh(C)
    λ3, λ2, λ1 = eigvals
    denom = λ1 if λ1 > 0 else 1e-30
    linearity = 1 - λ2 / denom
    planarity = 1 - 2 * λ3 / (λ1 + λ2) if (λ1 + λ2) > 0 else 0
    return linearity, planarity, λ1, λ2, λ3

# ── 分析 ──
results = {}
print(f"\n{'='*60}")
print(f"偏振验证分析: {NET}.{STA}")
print(f"{'='*60}")

for label, t_center, win_half in windows:
    ws, we = t_center - win_half, t_center + win_half
    m = (t >= ws) & (t <= we)
    zw, rw, tw = z[m], r[m], tc[m]
    t_win = t[m]

    # 能量
    ez = np.sum(zw**2); er = np.sum(rw**2); et = np.sum(tw**2); eh = er + et
    er_et = er/et if et > 0 else 0
    ez_eh = ez/eh if eh > 0 else 0
    log_er_et = np.log10(er_et) if er_et > 0 else -999

    # PCA
    lin, pla, l1, l2, l3 = pca_3d(zw, rw, tw)

    # 玫瑰图(水平面：arctan(T/R))
    theta = np.arctan2(tw, rw)
    weights = rw**2 + tw**2
    bins = np.linspace(-np.pi, np.pi, 37)
    hist, _ = np.histogram(theta, bins=bins, weights=weights)
    if np.max(hist) > 0: hist = hist / np.max(hist)

    # Z-R 平面玫瑰图 (P波偏振方向)
    theta_zr = np.arctan2(rw, zw)  # 相对垂直方向的夹角
    weights_zr = zw**2 + rw**2
    hist_zr, _ = np.histogram(theta_zr, bins=bins, weights=weights_zr)
    if np.max(hist_zr) > 0: hist_zr = hist_zr / np.max(hist_zr)

    results[label] = {
        'ws': ws, 'we': we,
        'er_et': er_et, 'ez_eh': ez_eh, 'log_er_et': log_er_et,
        'linearity': lin, 'planarity': pla,
        'rose_theta': theta, 'rose_w': weights,
        'hist': hist, 'bins': bins,
        'hist_zr': hist_zr,
        'zw': zw, 'rw': rw, 'tw': tw, 't_win': t_win,
    }

    # 判断
    pol_judge = 'radial (R>T)' if er_et > 1 else ('tangential (T>R)' if er_et < 0.8 else 'ambiguous (R≈T)')
    pca_judge = 'linear' if lin > 0.7 else 'elliptical'

    print(f"\n  【{label}】 {ws:.0f}-{we:.0f}s")
    print(f"    ER/ET={er_et:.2f}  → {pol_judge}")
    print(f"    PCA: linearity={lin:.3f} ({pca_judge})  planarity={pla:.3f}")
    print(f"    EZ/EH={ez_eh:.3f}")


# ── 绘图：左P右S，每列上半=波形+包络，下半=玫瑰图 ──
fig = plt.figure(figsize=(14, 10))
# 使用 GridSpec 精确控制，玫瑰图用 polar projection
gs = fig.add_gridspec(3, 2, height_ratios=[1.5, 1.2, 1.2],
                      hspace=0.35, wspace=0.35,
                      left=0.07, right=0.97, top=0.92, bottom=0.08)

for col, (label, res) in enumerate(results.items()):
    # ── Row 0: Z/R/T 波形+包络 ──
    ax_w = fig.add_subplot(gs[0, col])
    tw_w = res['t_win']
    # 先包络（Hilbert），然后各自归一化用于图面美观
    env_z_all = uniform_filter1d(np.abs(hilbert(res['zw'])), size=max(int(2/delta),1))
    env_r_all = uniform_filter1d(np.abs(hilbert(res['rw'])), size=max(int(2/delta),1))
    env_t_all = uniform_filter1d(np.abs(hilbert(res['tw'])), size=max(int(2/delta),1))
    # 每个窗口自己归一化（显示波形形态用）
    max_all = max(np.max(np.abs(res['zw'])), np.max(np.abs(res['rw'])), np.max(np.abs(res['tw'])))
    if max_all > 0:
        off = 3
        ax_w.plot(tw_w, res['zw']/max_all + off, 'k-', lw=0.6, label='Z')
        ax_w.plot(tw_w, res['rw']/max_all, 'g-', lw=0.6, label='R')
        ax_w.plot(tw_w, res['tw']/max_all - off, 'r-', lw=0.6, label='T')
    for yv in [off, 0, -off]:
        ax_w.axhline(y=yv, color='gray', lw=0.3, alpha=0.3)
    ax_w.set_xlim(res['ws'], res['we'])
    ax_w.set_title(f'{label}  ({res["ws"]:.0f}-{res["we"]:.0f}s)', fontsize=11, fontweight='bold')
    ax_w.legend(fontsize=7, ncol=3, loc='upper right')
    ax_w.set_ylabel('Norm. Amp.', fontsize=9)
    ax_w.grid(True, alpha=0.08, linestyle=':')
    ax_w.tick_params(labelsize=7)
    # 到时虚线
    ax_w.axvline(x=tt_P if 'P' in label else tt_S, color='gray', lw=0.8, ls='--', alpha=0.4)

    # ── Row 1: 水平面玫瑰图 (R-T) ──
    ax_rt = fig.add_subplot(gs[1, col], projection='polar')
    centers = (res['bins'][:-1] + res['bins'][1:]) / 2
    ax_rt.bar(centers, res['hist'], width=np.pi/18, bottom=0.1,
              color='gray', alpha=0.75, edgecolor='w', lw=0.3)
    for ang, lbl, clr in [(0, 'R', 'green'), (np.pi/2, 'T', 'red'),
                          (np.pi, 'R', 'green'), (-np.pi/2, 'T', 'red')]:
        ax_rt.plot([ang, ang], [0, 1.2], color=clr, lw=0.8, alpha=0.5)
        if ang <= 0:
            ax_rt.text(ang, 1.35, lbl, fontsize=8, fontweight='bold', ha='center', color=clr)
    ax_rt.set_title(f'R-T玫瑰图  |  ER/ET={res["er_et"]:.2f}', fontsize=9, fontweight='bold', pad=15)
    ax_rt.tick_params(labelsize=6)
    ax_rt.set_theta_zero_location('E')
    ax_rt.set_theta_direction(-1)
    ax_rt.set_ylim(0, 1.5)

    # ── Row 2: Z-R平面玫瑰图 ──
    ax_zr = fig.add_subplot(gs[2, col], projection='polar')
    # 重新定义Z-R玫瑰图（角度=arctan2(R, Z)，R=0°为径向水平，90°为垂直）
    bins_90 = np.linspace(-np.pi/2, np.pi/2, 19)  # -90° ~ +90°
    theta_zr = np.arctan2(res['rw'], res['zw'])  # R/Z的反正切
    weights_zr = res['zw']**2 + res['rw']**2
    hist_zr, _ = np.histogram(theta_zr, bins=bins_90, weights=weights_zr)
    if np.max(hist_zr) > 0: hist_zr = hist_zr / np.max(hist_zr)
    centers_zr = (bins_90[:-1] + bins_90[1:]) / 2
    ax_zr.bar(centers_zr, hist_zr, width=np.pi/18, bottom=0.1,
              color='purple', alpha=0.75, edgecolor='w', lw=0.3)
    # 标注方向
    ax_zr.plot([0, 0], [0, 1.2], color='green', lw=0.8, alpha=0.5)  # R=0°
    ax_zr.plot([np.pi/2, np.pi/2], [0, 1.2], color='black', lw=0.8, alpha=0.5)  # Z=90°
    ax_zr.text(0, 1.35, 'R', fontsize=8, fontweight='bold', ha='center', color='green')
    ax_zr.text(np.pi/2, 1.35, 'Z', fontsize=8, fontweight='bold', ha='center', color='black')
    ax_zr.set_title(f'Z-R玫瑰图  |  lin={res["linearity"]:.2f}', fontsize=9, fontweight='bold', pad=15)
    ax_zr.tick_params(labelsize=6)
    ax_zr.set_theta_zero_location('E')
    ax_zr.set_theta_direction(-1)
    ax_zr.set_ylim(0, 1.5)

# 总标题
fig.suptitle(f'{NET}.{STA}  偏振分析验证 (未滤波)  |  AZ={az:.1f}°  BAZ={baz:.1f}°',
             fontsize=12, fontweight='bold', y=0.98)

out = f"{OUT_DIR}/polar_validation_{STA}.png"
fig.savefig(out, dpi=200, bbox_inches='tight')
plt.close(fig)
print(f"\n✅ {out}")

# ── 结论 ──
print(f"\n{'='*60}")
print("理论预期验证：")
p_res = results['P波']
s_res = results['S波']
print(f"  P波: ER/ET={p_res['er_et']:.2f}, linearity={p_res['linearity']:.3f}  → "
      f"{'✓ P波径向偏振' if p_res['er_et'] > 1 else '✗ 不符合预期'}")
print(f"  S波: ER/ET={s_res['er_et']:.2f}, linearity={s_res['linearity']:.3f}  → "
      f"{'✓ S波径向偏振' if s_res['er_et'] > 1 else '✗ 不符合预期'}")
print(f"{'='*60}")
