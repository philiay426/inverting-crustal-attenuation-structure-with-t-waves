#!/usr/bin/env python3
"""
偏振属性: T波时段 线性度/平面度 (仅T波窗口, v3)
- 窗口: T-P-50s ~ T-S+150s, 滑动 1s 子窗 75% 重叠 (analytic 协方差)
- 输出: {STA}_linearity_planarity_3D.png
- 面板: (a) 线性度 vs 时间  (b) 平面度 vs 时间  (c) 3D 散点(L, P, 时间着色)
"""
import sys, os, glob, csv
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from obspy import read, UTCDateTime, Stream
from obspy.signal.filter import bandpass

sys.path.insert(0, '/home/sy/codex/t_phase_analysis')
from t_phase_analysis.pca import pca, rotate_stream_to_zrt

EV_ORIGIN = UTCDateTime(2019, 1, 22, 19, 1, 43)  # ← 必改: 事件发震时刻 (SAC 缺 o 字段时时间轴靠它)
BASE = "/mnt/f/sharedate/case/Southern_Africa/2019-01-22-mww67-prince-edward-islands-station_low"  # ← 必改
RESP = f'{BASE}/sacprocess/remove_response'
OUT = f'{BASE}/process_flow/2_single_station_analysis/SOE'  # ← 必改: OUT 硬编码台站名, 批量跑只改 STA,NET 不够
CSV = f'{BASE}/process_flow/1_station_distribution/stations.csv'

STA, NET = 'SOE', 'AF'  # ← 必改

# 读台站参数
with open(CSV) as f:
    for row in csv.DictReader((l for l in f if not l.startswith('#'))):
        if row['sta'] == STA:
            ttp = float(row['tt_tp']); tts = float(row['tt_ts'])
            bt = float(row['baz'])
            break

trs = {}
for ch in ['BHZ','BHN','BHE']:
    m = glob.glob(f'{RESP}/{NET}.{STA}.*.{ch}*_resp.sac')
    if m: trs[ch] = read(m[0])[0]
t0 = max(tr.stats.starttime for tr in trs.values())
t1 = min(tr.stats.endtime for tr in trs.values())
for ch in trs: trs[ch].trim(t0, t1)
delta = trs['BHZ'].stats.delta
sr = 1/delta
# 时间轴: SAC o 字段可靠时用 b-o (数据常从发震前开始记录)
sac = trs['BHZ'].stats.sac
o_val = getattr(sac, 'o', None)
if o_val is not None and 0 < abs(o_val) < 1e6:
    t_off = sac.b - o_val
else:
    t_off = trs['BHZ'].stats.starttime - EV_ORIGIN

# ZNE 流 — pca() 内部用 back_azimuth 自动旋转 ZRT, 不要手动 rotate 后传入(会触发 Non-ZNE 报错)
st = Stream([trs['BHZ'], trs['BHN'], trs['BHE']])

# T 波时段: T-P-50 ~ T-S+150
w0 = ttp - 50
w1 = tts + 150

st_win = st.slice(starttime=trs['BHZ'].stats.starttime + (w0 - t_off),
                  endtime=trs['BHZ'].stats.starttime + (w1 - t_off))
res = pca(st_win, window_length_s=1.0, window_step_s=0.25,
          back_azimuth=bt, covariance_method='analytic')

tw = np.asarray(res['time_windows_s'])  # 相对窗口起点
tw_abs = w0 + tw  # 相对发震
lin = np.asarray(res['linearity'])
pla = np.asarray(res['planarity'])
v = np.asarray(res['principal_vector'])
ang_rt = np.degrees(np.arctan2(v[:,2], v[:,1]))

# 分段统计: T-P 段 / T-S 段
m_tp = (tw_abs >= ttp-20) & (tw_abs <= ttp+40)
m_ts = (tw_abs >= tts-10) & (tw_abs <= tts+60)
print(f'T波时段 [{w0:.0f},{w1:.0f}]s, 窗数={len(tw)}')
print(f'T-P 段: lin中值={np.median(lin[m_tp]):.3f} pla中值={np.median(pla[m_tp]):.3f} R-T角={np.median(ang_rt[m_tp]):.1f}°')
print(f'T-S 段: lin中值={np.median(lin[m_ts]):.3f} pla中值={np.median(pla[m_ts]):.3f} R-T角={np.median(ang_rt[m_ts]):.1f}°')

# ── 图 ──
fig = plt.figure(figsize=(15, 5))
fig.suptitle(f'{NET}.{STA} Polarization attributes (T-wave window {w0:.0f}-{w1:.0f}s)',
             fontsize=13, fontweight='bold')

# (a) 线性度
ax1 = fig.add_subplot(1,3,1)
ax1.plot(tw_abs, lin, 'k-', lw=1.0)
ax1.axvline(ttp, c='#4d9fff', ls='--', lw=1.2, label='T-P')
ax1.axvline(tts, c='#ff4d4d', ls='--', lw=1.2, label='T-S')
ax1.axvspan(ttp, tts, color='gray', alpha=0.08)
ax1.set_xlim(w0, w1); ax1.set_ylim(0, 1.05)
ax1.set_xlabel('Time from origin (s)', fontsize=11, fontweight='bold')
ax1.set_ylabel('Linearity', fontsize=11, fontweight='bold')
ax1.legend(fontsize=8); ax1.grid(alpha=0.15)
ax1.set_title('(a) Linearity (rectilinearity)', fontsize=11, fontweight='bold')

# (b) 平面度
ax2 = fig.add_subplot(1,3,2)
ax2.plot(tw_abs, pla, 'b-', lw=1.0)
ax2.axvline(ttp, c='#4d9fff', ls='--', lw=1.2, label='T-P')
ax2.axvline(tts, c='#ff4d4d', ls='--', lw=1.2, label='T-S')
ax2.axvspan(ttp, tts, color='gray', alpha=0.08)
ax2.set_xlim(w0, w1); ax2.set_ylim(0, 1.05)
ax2.set_xlabel('Time from origin (s)', fontsize=11, fontweight='bold')
ax2.set_ylabel('Planarity', fontsize=11, fontweight='bold')
ax2.legend(fontsize=8); ax2.grid(alpha=0.15)
ax2.set_title('(b) Planarity', fontsize=11, fontweight='bold')

# (c) 3D 散点: R-T角 vs 时间 vs 线性度
ax3 = fig.add_subplot(1,3,3, projection='3d')
sc = ax3.scatter(ang_rt, tw_abs, lin, c=tw_abs, cmap='viridis', s=8, alpha=0.8)
ax3.set_xlabel('R-T angle (°)', fontsize=10)
ax3.set_ylabel('Time (s)', fontsize=10)
ax3.set_zlabel('Linearity', fontsize=10)
ax3.set_title('(c) 3D: R-T / time / linearity', fontsize=11, fontweight='bold')
fig.colorbar(sc, ax=ax3, shrink=0.6, label='Time (s)')

plt.tight_layout(rect=[0,0,1,0.94])
fn = f'{OUT}/{STA}_linearity_planarity_3D.png'
fig.savefig(fn, dpi=200, bbox_inches='tight'); plt.close(fig)
print(f'=> {fn}')
