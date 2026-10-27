#!/usr/bin/env python3
"""
T波三分量窗口 v2 (2026-08-04 定稿): 仅 T 波三分量窗口
- 不画 P/S 窗口
- 时间范围: T-P-100 ~ T-S+200
- 三分量: Z / R / T (2-8 Hz bandpass) 三行并排
- 样式: 波形灰色 (color='0.55'), 包络红色平滑8s
用法: python3 plot_waveform_windows.py
"""
import sys, os, glob
import numpy as np
from scipy.signal import hilbert
from scipy.ndimage import uniform_filter1d
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from obspy import read, UTCDateTime
from obspy.signal.filter import bandpass as bp_filter
from obspy.taup import TauPyModel
from obspy.geodetics.base import locations2degrees
from geographiclib.geodesic import Geodesic

sys.path.insert(0, '/home/sy/codex/t_phase_analysis')
from t_phase_analysis.traveltime_profile import get_track, tt_t

# ====== 事件参数 (改这里) ======
EVLA, EVLO, EVDP = -43.1219, 42.3568, 13.0
EV_ORIGIN = UTCDateTime(2019, 1, 22, 19, 1, 43)  # ← 发震时刻（SAC 缺 o 字段时必需）
STA, NET = 'SOE', 'AF'
BASE = "/mnt/f/sharedate/case/Southern_Africa/2019-01-22-mww67-prince-edward-islands-station_low"
# ================================
RESP_DIR = os.path.join(BASE, "sacprocess", "remove_response")
OUT = os.path.join(BASE, "process_flow", "2_single_station_analysis", STA)
os.makedirs(OUT, exist_ok=True)

SMOOTH_S = 8.0   # 包络平滑 8s

trs = {}
for ch in ['BHZ','BHN','BHE','BH1','BH2']:
    for pat in [f"{RESP_DIR}/{NET}.{STA}.*.{ch}*_resp.sac",
                f"{RESP_DIR}/{NET}.{STA}..{ch}*_resp.sac",
                f"{RESP_DIR}/{NET}.{STA}.{ch}*_resp.sac"]:
        m = glob.glob(pat)
        if m: trs[ch] = read(m[0])[0]; break
if len(trs) < 3: print("缺分量"); exit()

# 映射 BH1/BH2 到 BHN/BHE
if "BHN" not in trs and "BH1" in trs:
    trs["BHN"] = trs["BH1"]
if "BHE" not in trs and "BH2" in trs:
    trs["BHE"] = trs["BH2"]

t0 = max(tr.stats.starttime for tr in trs.values())
t1 = min(tr.stats.endtime for tr in trs.values())
for ch in trs: trs[ch].trim(t0, t1)

delta = trs['BHZ'].stats.delta; sr = 1/delta
sac = trs['BHZ'].stats.sac
# 时间轴: SAC o 字段可靠时用 b-o (数据常从发震前开始记录)
# 2013 南非: o=600 = 发震在数据起点后600s, b-o=-600 → t=0=发震 ✓
# 缺 o 时 fallback: starttime - EV_ORIGIN
o_val = getattr(sac, 'o', None)
if o_val is not None and 0 < abs(o_val) < 1e6:
    t_off = sac.b - o_val
else:
    t_off = trs['BHZ'].stats.starttime - EV_ORIGIN
g = Geodesic.WGS84.Inverse(EVLA, EVLO, sac.stla, sac.stlo)
az = g['azi1']; az = az+360 if az<0 else az
baz = (az+180)%360; br = np.radians(baz)

z = trs['BHZ'].data.astype(np.float64)
n = trs['BHN'].data.astype(np.float64)
e = trs['BHE'].data.astype(np.float64)
r = n*np.cos(br) + e*np.sin(br)
tc = -n*np.sin(br) + e*np.cos(br)
z_bp = bp_filter(z.copy(), 2, 8, sr, corners=4, zerophase=True)
r_bp = bp_filter(r.copy(), 2, 8, sr, corners=4, zerophase=True)
t_bp = bp_filter(tc.copy(), 2, 8, sr, corners=4, zerophase=True)
t = np.arange(len(z))*delta + t_off

taup = TauPyModel('iasp91')
dd = locations2degrees(EVLA, EVLO, sac.stla, sac.stlo)
arrs = taup.get_travel_times(source_depth_in_km=EVDP, distance_in_degree=dd, phase_list=['P','S'])
tt_P = next((a.time for a in arrs if a.name=='P'), None)
tt_S = next((a.time for a in arrs if a.name=='S'), None)
track = get_track(EVLO, EVLA, sac.stlo, sac.stla, 0.02, ref_el=-1000, plot_map=False, plot_elevation=False)
tt_r = tt_t(track, ref_el=-1000, source_depth_km=EVDP, ignore_small_islands=True)
tt_TP = tt_r["tt_t"]; tt_TS = tt_r["tt_ocean"] + tt_r["tt_land_s"]
print(f"{STA}: P={tt_P:.0f}s S={tt_S:.0f}s T-P={tt_TP:.0f}s T-S={tt_TS:.0f}s")

def calc_env(d):
    """包络, 平滑 SMOOTH_S 秒"""
    return uniform_filter1d(np.abs(hilbert(d.astype(np.float64))), size=max(int(SMOOTH_S/delta),1))

# ── Figure: T波三分量窗口 (T-P-100 ~ T-S+200) ──
tw_min, tw_max = int(tt_TP-100), int(tt_TS+200)
comps = [('Z', z_bp), ('R', r_bp), ('T', t_bp)]
fig, axes = plt.subplots(3, 1, figsize=(14, 9), sharex=True)
for k, (cname, cwf) in enumerate(comps):
    ax = axes[k]
    m = (t>=tw_min) & (t<=tw_max)
    p98 = np.percentile(np.abs(cwf[m]), 98) or 1
    env = calc_env(cwf)[m]/p98
    ax.plot(t[m], cwf[m]/p98, color='0.55', lw=0.7)     # 波形灰色
    ax.plot(t[m], env, 'r-', lw=1.2, alpha=0.8)          # 包络红色
    ax.fill_between(t[m], -env, env, color='red', alpha=0.08)
    ax.set_xlim(tw_min, tw_max)
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.set_ylabel('Norm. Amp.', fontsize=10)
    ax.set_title(f'{NET}.{STA}  {cname}  (BP 2-8Hz, env smooth {SMOOTH_S:.0f}s)', fontsize=11, fontweight='bold', loc='left')
    ax.grid(True, alpha=0.15, linestyle=':')
    if k == 2:
        ax.set_xlabel('Time from origin (s)', fontsize=10)
fig.suptitle(f'T-wave 3-component window  {tw_min}-{tw_max}s  (T-P-100 ~ T-S+200)',
             fontsize=13, fontweight='bold', y=0.995)
plt.tight_layout(rect=[0,0,1,0.98])
fn = f'{OUT}/{STA}_Twave_3comp_{tw_min}-{tw_max}s.png'
fig.savefig(fn, dpi=200, bbox_inches='tight')
plt.close(fig)
print(f"✓ {STA}_Twave_3comp_{tw_min}-{tw_max}s.png")
