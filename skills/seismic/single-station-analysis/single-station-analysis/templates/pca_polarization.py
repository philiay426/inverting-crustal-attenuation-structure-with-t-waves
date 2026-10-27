#!/usr/bin/env python3
"""
用 t_phase_analysis 包跑单台偏振 — 观测到时版 (2026-08-26 偏振引导定稿)
- 拾取: 理论到时 ±15s 内扫描 1s 子窗(75%重叠)，T-P 选最接近径向 R(0°/±180°)、
        T-S 选最接近横向 T(±90°) 的时窗中心；最终窗口 = 中心 ±2s（2026-09-02 定稿）
- 子窗: 1.0s, 75% 重叠
- 旋转: ZRT (back_azimuth), 仅 R-T 平面输出
- 协方差: analytic (Hilbert, 窄带2-8Hz更稳)
输出: SOE_pca_TP_TS.png
"""
import sys, os, glob, csv
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from obspy import read, UTCDateTime
from obspy.signal.filter import bandpass

sys.path.insert(0, '/home/sy/codex/t_phase_analysis')
from t_phase_analysis.pca import pca, rotate_stream_to_zrt

EV_ORIGIN = UTCDateTime(2019, 1, 22, 19, 1, 43)  # 发震时刻 (D.2019.022.190143)
BASE = "/mnt/f/sharedate/case/Southern_Africa/2019-01-22-mww67-prince-edward-islands-station_low"
RESP = f'{BASE}/sacprocess/remove_response'
OUT = f'{BASE}/process_flow/2_single_station_analysis/SOE'
CSV = f'{BASE}/process_flow/1_station_distribution/stations.csv'

STA, NET = 'SOE', 'AF'

# 读台站参数
with open(CSV) as f:
    for row in csv.DictReader((l for l in f if not l.startswith('#'))):
        if row['sta'] == STA:
            land = float(row['land_km']); ttp = float(row['tt_tp']); tts = float(row['tt_ts'])
            tp = float(row['tt_p']); bt = float(row['baz'])
            ocean_km = float(row['ocean_km'])
            break

# 陆地走时基准: t_land = t - tt_ocean, t=0 = 波上岸
TT_OCEAN = ocean_km / 1.48
TTP_LAND = ttp - TT_OCEAN   # T-P 陆地走时
TTS_LAND = tts - TT_OCEAN   # T-S 陆地走时
print(f'陆地走时: T-P={TTP_LAND:.0f}s T-S={TTS_LAND:.0f}s (ocean={TT_OCEAN:.0f}s)')

# 加载三分量
trs = {}
for ch in ['BHZ','BHN','BHE']:
    m = glob.glob(f'{RESP}/{NET}.{STA}.*.{ch}*_resp.sac')
    if m: trs[ch] = read(m[0])[0]
t0 = max(tr.stats.starttime for tr in trs.values())
t1 = min(tr.stats.endtime for tr in trs.values())
for ch in trs: trs[ch].trim(t0, t1)
delta = trs['BHZ'].stats.delta

# 时间轴: SAC o 字段可靠时用 b-o (数据常从发震前开始记录)
sac = trs['BHZ'].stats.sac
o_val = getattr(sac, 'o', None)
if o_val is not None and 0 < abs(o_val) < 1e6:
    t_off = sac.b - o_val
else:
    t_off = trs['BHZ'].stats.starttime - EV_ORIGIN

# 时间轴(相对发震)
t = np.arange(trs['BHZ'].stats.npts)*delta + t_off

# 2-8 Hz 带通 + 包络(平滑8s) 仅用于显示
sr = 1/delta
from scipy.signal import hilbert
from obspy.signal.util import smooth

def env(d):
    return smooth(np.abs(hilbert(d)), max(1, int(8*sr)))

zf = bandpass(trs['BHZ'].data.astype(float), 2, 8, df=sr, corners=4)
envz = env(zf)
t_land = t - TT_OCEAN  # 陆地走时坐标

# 组装 Stream（pca 内部自动 ZNE→ZRT 旋转）
from obspy import Stream
st = Stream([trs['BHZ'], trs['BHN'], trs['BHE']])

def run_pca(w0_l, w1_l):
    # 传入陆地走时坐标，内部 +TT_OCEAN 转绝对坐标（⚠️ 不要传绝对坐标，会双重偏移）
    w0_abs = w0_l + TT_OCEAN; w1_abs = w1_l + TT_OCEAN
    st_win = st.slice(starttime=trs['BHZ'].stats.starttime + (w0_abs - t_off),
                      endtime=trs['BHZ'].stats.starttime + (w1_abs - t_off))
    # 1s 窗 75% 重叠 → step 0.25s
    res = pca(st_win, window_length_s=1.0, window_step_s=0.25,
              back_azimuth=bt, covariance_method='analytic')
    return res

def scan_rt_angle(t0_l, t1_l):
    """在陆地走时 [t0_l, t1_l] 内扫描 1s 子窗，返回 (子窗中心陆地时刻, R-T 角, 线性度)。"""
    res = run_pca(t0_l, t1_l)
    tw = t0_l + np.asarray(res['time_windows_s']) + 0.5  # 子窗中心（陆地走时）
    v = np.asarray(res['principal_vector'])
    ang = np.degrees(np.arctan2(v[:, 2], v[:, 1]))
    lin = np.asarray(res['linearity'])
    valid = np.asarray(res['valid_mask']) & np.isfinite(ang)
    return tw[valid], ang[valid], lin[valid]

def pick_by_polarization(t0_l, t1_l, target='radial', half=2.0):
    """偏振引导拾取：在搜索范围内扫描候选时窗，选整体偏振最接近目标的时窗中心。

    每个候选时窗 = 中心 ±half（默认 ±2s，2026-09-02 前十站敏感性定稿；旧值 ±3s），用窗内所有 1s 子窗的 R-T 角中值代表其偏振；
    target='radial' → 最接近径向 R (0°/±180°)，用于 T-P；
    target='transverse' → 最接近横向 T (±90°)，用于 T-S。
    返回 (中心陆地时刻, 窗口 R-T 角中值, 窗口线性度中值)。
    """
    tw, ang, lin = scan_rt_angle(t0_l, t1_l)
    if len(tw) == 0:
        raise RuntimeError('扫描窗口无有效子窗')
    if target == 'radial':
        # 角度距离折叠到 [0,90]：径向 0/180 → 0，横向 90/-90 → 90
        d_ang = np.abs(ang)
        d_ang = np.minimum(d_ang, np.abs(d_ang - 180))
        d_ang = np.minimum(d_ang, np.abs(d_ang + 180))
        d_ang = np.minimum(d_ang, 180 - d_ang)
    else:
        d_ang = np.abs(ang)
        d_ang = np.minimum(d_ang, np.abs(d_ang - 180))
        d_ang = np.minimum(d_ang, np.abs(d_ang + 180))
        d_ang = 90 - np.minimum(d_ang, 180 - d_ang)
    step = np.median(np.diff(tw)) if len(tw) > 1 else 0.25
    half_n = max(1, int(round(half / step)))
    # 每个候选中心：取 ±half 邻域子窗的中值偏振/线性度
    best = None
    for i in range(len(tw)):
        i0, i1 = max(0, i - half_n), min(len(tw), i + half_n + 1)
        med_ang = np.median(d_ang[i0:i1])
        med_lin = np.median(lin[i0:i1])
        # 与目标距离 = 中值角度距离；线性度作为次要判据
        score = (med_ang, -med_lin)
        if best is None or score < best[0]:
            best = (score, i)
    _, i = best
    i0, i1 = max(0, i - half_n), min(len(tw), i + half_n + 1)
    return tw[i], float(np.median(ang[i0:i1])), float(np.median(lin[i0:i1]))

# 偏振引导观测拾取: 搜索范围 = 理论到时 ±15s (避免选到远处偶然偏振相似的时段)
SEARCH_HALF = 15.0
t_obs_tp_l, ang_tp_pick, lin_tp_pick = pick_by_polarization(TTP_LAND-SEARCH_HALF, TTP_LAND+SEARCH_HALF, 'radial')
t_obs_ts_l, ang_ts_pick, lin_ts_pick = pick_by_polarization(
    max(TTS_LAND-SEARCH_HALF, t_obs_tp_l+5), TTS_LAND+SEARCH_HALF, 'transverse')
# 最终窗口: 观测时窗中心 ±2s (与理论版口径一致, 2026-09-02 定稿)
w_tp0_l, w_tp1_l = t_obs_tp_l-2, t_obs_tp_l+2
w_ts0_l, w_ts1_l = t_obs_ts_l-2, t_obs_ts_l+2
if w_tp1_l > w_ts0_l:
    shift = w_tp1_l - w_ts0_l + 0.5
    w_ts0_l, w_ts1_l = w_ts0_l+shift, w_ts1_l+shift
# 陆地坐标 → 绝对坐标 (用于画图)
t_obs_tp = t_obs_tp_l + TT_OCEAN; w_tp0 = w_tp0_l + TT_OCEAN; w_tp1 = w_tp1_l + TT_OCEAN
t_obs_ts = t_obs_ts_l + TT_OCEAN; w_ts0 = w_ts0_l + TT_OCEAN; w_ts1 = w_ts1_l + TT_OCEAN
print(f'T-P: th={ttp:.0f}(陆地{TTP_LAND:.0f}) obs={t_obs_tp:.0f} 窗=[{w_tp0:.1f},{w_tp1:.1f}] '
      f'拾取角={ang_tp_pick:.0f}° lin={lin_tp_pick:.2f}')
print(f'T-S: th={tts:.0f}(陆地{TTS_LAND:.0f}) obs={t_obs_ts:.0f} 窗=[{w_ts0:.1f},{w_ts1:.1f}] '
      f'拾取角={ang_ts_pick:.0f}° lin={lin_ts_pick:.2f}')

res_tp = run_pca(w_tp0_l, w_tp1_l)
res_ts = run_pca(w_ts0_l, w_ts1_l)

for lbl, res in [('T-P', res_tp), ('T-S', res_ts)]:
    v = np.asarray(res['principal_vector'])
    # R-T 平面角: 在 ZRT 中 principal_vector 的行 = [Z,R,T]
    ang_rt = np.degrees(np.arctan2(v[:,2], v[:,1]))
    med = np.median(ang_rt)
    lin = np.median(np.asarray(res['linearity']))
    pla = np.median(np.asarray(res['planarity']))
    print(f'{lbl}: R-T角中值={med:.1f}° lin={lin:.3f} pla={pla:.3f} 窗数={len(ang_rt)}')

# ── 图: 2×2 ──
fig = plt.figure(figsize=(14, 10))
fig.suptitle(f'{NET}.{STA} PCA (t_phase_analysis pkg)  T-P={t_obs_tp:.0f}s T-S={t_obs_ts:.0f}s  '
             f'1s window/75% overlap/analytic', fontsize=13, fontweight='bold')

def rose_rt(ax, res, title):
    v = np.asarray(res['principal_vector'])
    ang = np.degrees(np.arctan2(v[:,2], v[:,1]))
    w = np.asarray(res['eigenvalues'])[:,0]
    bins = np.linspace(-180, 180, 37)
    h, _ = np.histogram(ang, bins=bins, weights=w)
    th = (bins[:-1]+bins[1:])/2
    ax.bar(np.radians(th), h, width=np.radians(10), alpha=0.7, color='#1f77b4')
    mx = max(h)*1.3 if len(h) else 1
    ax.plot([0,0],[0,mx],'g--',lw=1); ax.text(0,mx*1.15,'R',c='green',ha='center')
    ax.plot([np.pi/2,np.pi/2],[0,mx],'r--',lw=1); ax.text(np.pi/2,mx*1.15,'T',c='red',ha='center')
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_theta_zero_location('N'); ax.set_theta_direction(-1)

ax1 = fig.add_subplot(2,2,1, polar=True)
rose_rt(ax1, res_tp, f'T-P R-T  med={np.median(np.degrees(np.arctan2(np.asarray(res_tp["principal_vector"])[:,2], np.asarray(res_tp["principal_vector"])[:,1]))):.0f}°')
ax2 = fig.add_subplot(2,2,2, polar=True)
rose_rt(ax2, res_ts, f'T-S R-T  med={np.median(np.degrees(np.arctan2(np.asarray(res_ts["principal_vector"])[:,2], np.asarray(res_ts["principal_vector"])[:,1]))):.0f}°')

# 包络+窗口
ax3 = fig.add_subplot(2,2,3)
w0, w1 = ttp-100, tts+200
mm = (t>=w0)&(t<=w1)
ax3.plot(t[mm], envz[mm]/max(envz[mm]), 'r-', lw=1.2)
ax3.plot(t[mm], zf[mm]/max(abs(zf[mm])), 'k-', lw=0.3, alpha=0.5)
ax3.axvspan(w_tp0, w_tp1, color='#1f77b4', alpha=0.15)
ax3.axvspan(w_ts0, w_ts1, color='#9467bd', alpha=0.15)
ax3.axvline(ttp, c='#1f77b4', ls='--'); ax3.axvline(tts, c='#9467bd', ls='--')
ax3.set_xlim(w0,w1); ax3.set_title('Envelope + pick windows', fontsize=12, fontweight='bold')
ax3.set_xlabel('Time from origin (s)'); ax3.grid(alpha=0.08)

# 线性度
ax4 = fig.add_subplot(2,2,4)
for lbl, res, clr in [('T-P', res_tp, '#1f77b4'), ('T-S', res_ts, '#9467bd')]:
    tw = np.asarray(res['time_windows_s'])
    tw = tw - tw[0]
    ax4.plot(tw, np.asarray(res['linearity']), color=clr, lw=1.2, label=f'{lbl} L')
ax4.axhline(0.5, c='gray', ls=':'); ax4.set_ylim(0,1.05)
ax4.legend(fontsize=9); ax4.set_title('Linearity within windows', fontsize=12, fontweight='bold')
ax4.set_xlabel('Time in window (s)'); ax4.grid(alpha=0.08)

plt.tight_layout(rect=[0,0,1,0.96])
fn = f'{OUT}/{STA}_pca_TP_TS.png'
fig.savefig(fn, dpi=200, bbox_inches='tight'); plt.close(fig)
print(f'=> {fn}')
