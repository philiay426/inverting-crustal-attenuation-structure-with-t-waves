#!/usr/bin/env python3
"""
用 t_phase_analysis 包跑单台偏振 — 理论到时版 (T-P / T-S)
- 窗口: 理论 T-P/T-S ±2s (不重叠, 用 w_tp[1] > w_ts[0] 检查)
- 方法开关 --method: cov | amp | count, 默认 count(采样点等权/投票);
  优先级: --method > 环境变量 POLAR_METHOD > 默认 count。方向+权重见 polar_methods.py
  · cov   = t_phase_analysis pca() 子窗(1s/75%重叠, analytic)主轴方位 + λ1 权重 (数值与旧版一致)
  · amp   = 2-8Hz R/T 逐采样点 atan2(T,R) + (R²+T²) 权重
  · count = 2-8Hz R/T 逐采样点 atan2(T,R), 每点 1 票
- 输出文件名保持 {STA}_pca_theory_TP_TS.png (不改名, 避免断链);
  实际方法名与 default 标记写在图标题里
- 下方 linearity 面板始终是 PCA 诊断量 (与 polarization_attributes.py 同定义)
- 旋转: ZNE→ZRT (back_azimuth), 仅 R-T 平面
输出: {STA}_pca_theory_TP_TS.png
"""
import sys, os, glob, csv, argparse
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from obspy import read, UTCDateTime, Stream
from obspy.signal.filter import bandpass
from obspy.signal.util import smooth
from scipy.signal import hilbert

for _codex in ('/home/sy/codex/t_phase_analysis', '/home/sy/codex'):
    if os.path.isdir(_codex): sys.path.insert(0, _codex)
from t_phase_analysis.pca import pca, rotate_stream_to_zrt

# 共享模块 polar_methods.py 与本脚本同目录(驱动脚本会一起拷到台站目录)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import polar_methods

EV_ORIGIN = UTCDateTime(2019, 1, 22, 19, 1, 43)  # 发震时刻 (D.2019.022.190143)
BASE = "/mnt/f/sharedate/case/Southern_Africa/2019-01-22-mww67-prince-edward-islands-station_low"
RESP = f'{BASE}/sacprocess/remove_response'
OUT = f'{BASE}/process_flow/2_single_station_analysis/SOE'
CSV = f'{BASE}/process_flow/1_station_distribution/stations.csv'

STA, NET = 'SOE', 'AF'

# ─── CLI: --method (cov|amp|count, 默认 count) / --out-dir ───
def _parse_args(argv=None):
    p = argparse.ArgumentParser(description='T-P/T-S theory-window R-T polarization')
    polar_methods.add_common_args(p)
    return p.parse_args(argv)


ARGS = _parse_args()
METHOD = polar_methods.resolve_method(ARGS.method)
METHOD_LBL = polar_methods.method_label(METHOD)
OUT = polar_methods.resolve_out_dir(ARGS.out_dir, OUT)
os.makedirs(OUT, exist_ok=True)

# 读台站参数
with open(CSV) as f:
    for row in csv.DictReader((l for l in f if not l.startswith('#'))):
        if row['sta'] == STA:
            land = float(row['land_km']); ttp = float(row['tt_tp']); tts = float(row['tt_ts'])
            bt = float(row['baz'])
            ocean_km = float(row['ocean_km'])
            break

# 陆地走时基准: 窗口在陆地走时坐标 (t=0 波上岸), 避免绝对到时偏差
TT_OCEAN = ocean_km / 1.48
TTP_LAND = ttp - TT_OCEAN; TTS_LAND = tts - TT_OCEAN
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
t = np.arange(trs['BHZ'].stats.npts)*delta + t_off
sr = 1/delta

# 2-8 Hz 带通 + 包络(平滑8s) 仅用于显示
zf = bandpass(trs['BHZ'].data.astype(float), 2, 8, df=sr, corners=4)
envz = smooth(np.abs(hilbert(zf)), max(1, int(8*sr)))

# 理论到时窗口: T-P/T-S 均按计算到时 ±2s（2026-09-02 前十站窗口敏感性测试定稿：±2s 为最短稳定窗；T-S 不延后，与 T-P 同口径）
w_tp = (TTP_LAND-2, TTP_LAND+2)
w_ts = (TTS_LAND-2, TTS_LAND+2)
if w_tp[1] > w_ts[0]:
    shift = w_tp[1] - w_ts[0] + 0.5
    w_ts = (w_ts[0]+shift, w_ts[1]+shift)
print(f'T-P: th={ttp:.0f}s(陆地{TTP_LAND:.0f}) 窗=[{w_tp[0]:.1f},{w_tp[1]:.1f}]')
print(f'T-S: th={tts:.0f}s(陆地{TTS_LAND:.0f}) 窗=[{w_ts[0]:.1f},{w_ts[1]:.1f}]')

st = Stream([trs['BHZ'], trs['BHN'], trs['BHE']])

def run_pca(w0, w1):
    # 陆地走时坐标 → 绝对坐标 (相对发震)
    w0_abs = w0 + TT_OCEAN; w1_abs = w1 + TT_OCEAN
    st_win = st.slice(starttime=trs['BHZ'].stats.starttime + (w0_abs - t_off),
                      endtime=trs['BHZ'].stats.starttime + (w1_abs - t_off))
    return pca(st_win, window_length_s=1.0, window_step_s=0.25,
               back_azimuth=bt, covariance_method='analytic')

res_tp = run_pca(*w_tp)
res_ts = run_pca(*w_ts)

for lbl, res in [('T-P', res_tp), ('T-S', res_ts)]:
    v = np.asarray(res['principal_vector'])
    ang_rt = np.degrees(np.arctan2(v[:,2], v[:,1]))
    lin = np.median(np.asarray(res['linearity']))
    pla = np.median(np.asarray(res['planarity']))
    print(f'{lbl}: R-T角中值={np.median(ang_rt):.1f}° lin={lin:.3f} pla={pla:.3f} 窗数={len(ang_rt)}')

# ── 玫瑰图方向+权重 (按 --method) ──
# 窗口在陆地走时坐标定义; 取窗/画图用绝对坐标(加回海洋走时)
w_tp_abs = (w_tp[0] + TT_OCEAN, w_tp[1] + TT_OCEAN)
w_ts_abs = (w_ts[0] + TT_OCEAN, w_ts[1] + TT_OCEAN)


def window_rose_data(res, w_abs):
    """返回 (angles_deg, weights) —— R-T 系角度, 0=径向 R, ±90=横向 T"""
    if METHOD == polar_methods.METHOD_COV:
        # cov: 沿用 pca() 子窗主轴 + λ1 权重 (与旧版逐位一致)
        v = np.asarray(res['principal_vector'])
        return np.degrees(np.arctan2(v[:, 2], v[:, 1])), np.asarray(res['eigenvalues'])[:, 0]
    # amp / count: 2-8Hz R/T 逐采样点 atan2(T,R) (= atan2(E,N) - BAZ)
    st_zrt_ = rotate_stream_to_zrt(st.copy(), bt)
    rt_bp = {}
    for ch in ('R', 'T'):
        trc = next(x for x in st_zrt_ if x.stats.channel.endswith(ch))
        rt_bp[ch] = bandpass(trc.data.astype(float), 2, 8, df=sr, corners=4)
    mw = (t >= w_abs[0]) & (t <= w_abs[1])
    ang, w = polar_methods.sample_direction_weights(METHOD, rt_bp['R'][mw], rt_bp['T'][mw])
    return np.degrees(ang), w


ang_tp_rose, wt_tp_rose = window_rose_data(res_tp, w_tp_abs)
ang_ts_rose, wt_ts_rose = window_rose_data(res_ts, w_ts_abs)
for lbl, ang, w in [('T-P', ang_tp_rose, wt_tp_rose), ('T-S', ang_ts_rose, wt_ts_rose)]:
    sec = polar_methods.sector_fractions(np.radians(ang), w)
    print(f'{lbl} [{METHOD}]: 玫瑰中值={np.median(ang):.1f}° '
          f'径向={sec["frac_radial"]*100:.0f}% 横向={sec["frac_transverse"]*100:.0f}% 样本={len(ang)}')

# ── 图: 2×2 ──
fig = plt.figure(figsize=(14, 10))
# 图标题只保留台站、SNR与理论到时；子窗/平滑/窗口策略保留在本代码注释中。
fig.suptitle(f'{STA}  |  SNR {float(row["snr"]):.1f} dB  |  T-P {ttp:.0f} s  |  T-S {tts:.0f} s'
             f'  |  method: {METHOD_LBL}',
             fontsize=13, fontweight='bold')

def rose_rt(ax, ang, w, title):
    bins = np.linspace(-180, 180, 37)
    h, _ = np.histogram(ang, bins=bins, weights=w)
    th = (bins[:-1]+bins[1:])/2
    ax.bar(np.radians(th), h, width=np.radians(10), alpha=0.7, color='#1f77b4')
    mx = max(h)*1.3 if len(h) else 1
    # R/T 参考虚线与图内加粗标签；标签接近图外缘，不与方位角刻度重叠。
    ax.plot([0,0],[0,mx],'g--',lw=1)
    ax.text(np.radians(5), mx*0.92, 'R', color='green', ha='center', va='center',
            fontweight='bold', fontsize=12)
    ax.plot([np.pi/2,np.pi/2],[0,mx],'r--',lw=1)
    ax.text(np.radians(85), mx*0.92, 'T', color='red', ha='center', va='center',
            fontweight='bold', fontsize=12)
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_theta_zero_location('N'); ax.set_theta_direction(-1)

ax1 = fig.add_subplot(2,2,1, polar=True)
rose_rt(ax1, ang_tp_rose, wt_tp_rose,
        f'T-P R-T  med={np.median(ang_tp_rose):.0f}°  [{METHOD_LBL}]')
ax2 = fig.add_subplot(2,2,2, polar=True)
rose_rt(ax2, ang_ts_rose, wt_ts_rose,
        f'T-S R-T  med={np.median(ang_ts_rose):.0f}°  [{METHOD_LBL}]')

ax3 = fig.add_subplot(2,2,3)
w0, w1 = ttp-100, tts+200
mm = (t>=w0) & (t<=w1)
# 理论窗口按陆地走时定义；波形图横轴为相对发震时间，必须加回海洋走时(w_tp_abs/w_ts_abs 已在上方定义)。
# Z/R/T 共用 P98，保留真实分量相对振幅；灰色波形+彩色4 s平滑包络。
st_zrt = rotate_stream_to_zrt(st.copy(), bt)
colors = {'Z':'#333333', 'R':'#1f77b4', 'T':'#d62728'}
plot_data = {}
for ch in ('Z', 'R', 'T'):
    trc = next(x for x in st_zrt if x.stats.channel.endswith(ch))
    plot_data[ch] = bandpass(trc.data.astype(float), 2, 8, df=sr, corners=4)
global_scale = np.percentile(np.abs(np.concatenate([plot_data[ch][mm] for ch in ('Z','R','T')])), 98) or 1.0
for ch, offset in [('Z', 1.0), ('R', 0.0), ('T', -1.0)]:
    dat = plot_data[ch]
    env_ch = smooth(np.abs(hilbert(dat)), max(1, int(4*sr))) / global_scale * 0.32
    ax3.plot(t[mm], dat[mm]/global_scale*0.32 + offset, color='0.70', lw=0.35, alpha=0.65)
    ax3.plot(t[mm], env_ch[mm] + offset, color=colors[ch], lw=0.85, label=ch)
    ax3.fill_between(t[mm], offset, env_ch[mm] + offset, color=colors[ch], alpha=0.12, lw=0)
ax3.axvspan(w_tp_abs[0], w_tp_abs[1], color='#4d9fff', alpha=0.22, label='Theory T-P ±2 s')
ax3.axvspan(w_ts_abs[0], w_ts_abs[1], color='#9467bd', alpha=0.18, label='Theory T-S window')
ax3.axvline(ttp, c='#4d9fff', ls='--', lw=0.8)
ax3.axvline(tts, c='#9467bd', ls='--', lw=0.8)
ax3.set_yticks([-1,0,1], ['T','R','Z'])
ax3.set_xlim(w0,w1); ax3.set_title('Z/R/T waveforms + envelopes + PCA windows', fontsize=12, fontweight='bold')
ax3.set_xlabel('Time from origin (s)'); ax3.grid(alpha=0.08)
ax3.legend(loc='upper left', fontsize=7, ncol=2, framealpha=0.9)

ax4 = fig.add_subplot(2,2,4)
# 该面板是 PCA 诊断量(定义同 polarization_attributes.py), 与 --method 无关
for lbl, res, clr in [('T-P', res_tp, '#1f77b4'), ('T-S', res_ts, '#9467bd')]:
    tw = np.asarray(res['time_windows_s']); tw = tw - tw[0]
    ax4.plot(tw, np.asarray(res['linearity']), color=clr, lw=1.2, label=f'{lbl} L')
ax4.axhline(0.5, c='gray', ls=':'); ax4.set_ylim(0,1.05)
ax4.legend(fontsize=9); ax4.set_title('Linearity within windows (PCA diagnostic)', fontsize=12, fontweight='bold')
ax4.set_xlabel('Time in window (s)'); ax4.grid(alpha=0.08)

plt.tight_layout(rect=[0,0,1,0.96])
fn = f'{OUT}/{STA}_pca_theory_TP_TS.png'
fig.savefig(fn, dpi=200, bbox_inches='tight'); plt.close(fig)
print(f'=> {fn}')
