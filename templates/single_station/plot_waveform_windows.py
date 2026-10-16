#!/usr/bin/env python3
"""
单台波形窗: PS窗 + T波窗 (2-8 Hz)
用法:
  python3 plot_waveform_windows.py STA NET EVLA EVLO EVDP BASE [--resp-dir DIR]

生成: {STA}_PS_window_{t1}-{t2}s.png, {STA}_Twave_window_{t1}-{t2}s.png
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

parser = argparse.ArgumentParser(description='Single-station waveform windows')
parser.add_argument('sta')
parser.add_argument('net')
parser.add_argument('evla', type=float)
parser.add_argument('evlo', type=float)
parser.add_argument('evdp', type=float)
parser.add_argument('base')
parser.add_argument('--resp-dir', default=None)
parser.add_argument('--filter-fmin', type=float, default=2.0)
parser.add_argument('--filter-fmax', type=float, default=8.0)
args = parser.parse_args()

STA, NET = args.sta, args.net
EVLA, EVLO, EVDP = args.evla, args.evlo, args.evdp
BASE = args.base
RESP_DIR = args.resp_dir or os.path.join(BASE, 'sacprocess', 'remove_response')
OUT_DIR = os.path.join(BASE, 'process_flow', '3_single_station_analysis', STA)
os.makedirs(OUT_DIR, exist_ok=True)

# ─── Load ───
chs = ['BHZ','BHN','BHE']
trs = {}
for ch in chs:
    for pat in [f"{RESP_DIR}/{NET}.{STA}.*.{ch}*_resp.sac",
                f"{RESP_DIR}/{NET}.{STA}..{ch}*_resp.sac",
                f"{RESP_DIR}/{NET}.{STA}.{ch}*_resp.sac"]:
        m = glob.glob(pat)
        if m: trs[ch] = read(m[0])[0]; break
if len(trs) < 3: sys.exit('ERROR: missing components')

t0 = max(tr.stats.starttime for tr in trs.values())
t1 = min(tr.stats.endtime for tr in trs.values())
for ch in trs: trs[ch].trim(t0, t1)

delta = trs['BHZ'].stats.delta; sr = 1/delta
sac = trs['BHZ'].stats.sac; t_off = sac.b - sac.o
g = Geodesic.WGS84.Inverse(EVLA, EVLO, sac.stla, sac.stlo)
az = g['azi1']; az = az+360 if az<0 else az
baz = (az+180)%360; br = np.radians(baz)

z = trs['BHZ'].data.astype(np.float64)
n = trs['BHN'].data.astype(np.float64)
e = trs['BHE'].data.astype(np.float64)
# 2026-09-02: 统一 obspy 约定（R 背离震源）
from obspy.signal.rotate import rotate_ne_rt
r, tc = rotate_ne_rt(n, e, float(baz))
z_bp = bp_filter(z.copy(), args.filter_fmin, args.filter_fmax, sr, corners=4, zerophase=True)
t = np.arange(len(z))*delta + t_off

# ─── Arrivals ───
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
    return uniform_filter1d(np.abs(hilbert(d.astype(np.float64))), size=max(int(2/delta),1))

# ── PS window ──
ps_min, ps_max = int(tt_P-90), int(tt_S+100)
fig, ax = plt.subplots(1, 1, figsize=(12, 3))
m = (t>=ps_min)&(t<=ps_max); p98 = np.percentile(np.abs(z_bp[m]),98) or 1
amp = 0.25
ax.plot(t[m], z_bp[m]/p98*amp, 'k-', lw=0.6)
env = calc_env(z_bp)[m]/p98*amp
ax.fill_between(t[m], -env, env, color='red', alpha=0.10)
ax.plot(t[m], env, 'r-', lw=0.8, alpha=0.5)
ax.plot(t[m], -env, 'r-', lw=0.8, alpha=0.5)
ax.set_xlim(ps_min, ps_max); ax.set_ylim(-0.35, 0.35)
ax.xaxis.set_major_locator(MultipleLocator(50)); ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.set_xlabel('Time from origin (s)', fontsize=10)
ax.set_ylabel('Norm. Amp.', fontsize=10)
for ph, tt_v, clr in [('P',tt_P,'#2ca02c'),('S',tt_S,'#ff7f0e'),
                       ('T-P',tt_TP,'#1f77b4'),('T-S',tt_TS,'#d62728')]:
    if tt_v and ps_min<=tt_v<=ps_max:
        ax.axvline(tt_v, color=clr, ls='--', lw=0.8, alpha=0.5)
        ax.text(tt_v, 0.92, f'{ph}={tt_v:.0f}s', color=clr, fontsize=8, fontweight='bold',
                ha='center', va='top', transform=ax.get_xaxis_transform())
ax.set_title(f'{NET}.{STA}  PS window  {ps_min}-{ps_max}s', fontsize=10, fontweight='bold')
ax.grid(True, alpha=0.15, linestyle=':')
fn_ps = f'{OUT_DIR}/{STA}_PS_window_{ps_min}-{ps_max}s.png'
fig.savefig(fn_ps, dpi=200, bbox_inches='tight'); plt.close(fig)
print(f'OK {fn_ps}')

# ── T-wave window ──
tw_min, tw_max = int(tt_TP-80), int(tt_TS+80)
fig, ax = plt.subplots(1, 1, figsize=(12, 3))
m = (t>=tw_min)&(t<=tw_max); p98 = np.percentile(np.abs(z_bp[m]),98) or 1
ax.plot(t[m], z_bp[m]/p98*amp, 'k-', lw=0.6)
env = calc_env(z_bp)[m]/p98*amp
ax.fill_between(t[m], -env, env, color='red', alpha=0.10)
ax.plot(t[m], env, 'r-', lw=0.8, alpha=0.5)
ax.plot(t[m], -env, 'r-', lw=0.8, alpha=0.5)
ax.set_xlim(tw_min, tw_max); ax.set_ylim(-0.35, 0.35)
ax.xaxis.set_major_locator(MultipleLocator(50)); ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.set_xlabel('Time from origin (s)', fontsize=10)
ax.set_ylabel('Norm. Amp.', fontsize=10)
for ph, tt_v, clr in [('T-P',tt_TP,'#1f77b4'),('T-S',tt_TS,'#d62728')]:
    if tt_v and tw_min<=tt_v<=tw_max:
        ax.axvline(tt_v, color=clr, ls='--', lw=0.8, alpha=0.5)
        ax.text(tt_v, 0.92, f'{ph}={tt_v:.0f}s', color=clr, fontsize=8, fontweight='bold',
                ha='center', va='top', transform=ax.get_xaxis_transform())
ax.set_title(f'{NET}.{STA}  T-wave window  {tw_min}-{tw_max}s', fontsize=10, fontweight='bold')
ax.grid(True, alpha=0.15, linestyle=':')
fn_tw = f'{OUT_DIR}/{STA}_Twave_window_{tw_min}-{tw_max}s.png'
fig.savefig(fn_tw, dpi=200, bbox_inches='tight'); plt.close(fig)
print(f'OK {fn_tw}')
