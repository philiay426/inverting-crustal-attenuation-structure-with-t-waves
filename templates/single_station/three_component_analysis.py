#!/usr/bin/env python3
"""
单台三分量波形+时频谱 (2-8 Hz)
用法:
  python3 three_component_analysis.py STA NET EVLA EVLO EVDP BASE [--resp-dir DIR]

示例:
  # Barbados 2014
  python3 three_component_analysis.py O51A TA 14.6682 -58.9272 15.0 /mnt/e/sharedate/case/2014-02-18-mw65-north-atlantic-ocean
"""
import sys, os, glob, argparse
import numpy as np
from scipy.signal import hilbert
from scipy.ndimage import uniform_filter1d
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.lines import Line2D
from obspy import read
from obspy.signal.filter import bandpass as bp_filter
from obspy.taup import TauPyModel
from obspy.geodetics.base import locations2degrees
from geographiclib.geodesic import Geodesic

_codex = '/home/philiay/shared/codex/t_phase_analysis'
if os.path.isdir(_codex): sys.path.insert(0, _codex)
from t_phase_analysis.traveltime_profile import get_track, tt_t
from t_phase_analysis.waveform import compute_spectrogram

# ─── CLI ───
parser = argparse.ArgumentParser(description='Single-station 3-component spectrogram')
parser.add_argument('sta')
parser.add_argument('net')
parser.add_argument('evla', type=float)
parser.add_argument('evlo', type=float)
parser.add_argument('evdp', type=float)
parser.add_argument('base', help='event root directory')
parser.add_argument('--resp-dir', default=None, help='override remove_response dir (default: BASE/sacprocess/remove_response)')
parser.add_argument('--filter-fmin', type=float, default=2.0)
parser.add_argument('--filter-fmax', type=float, default=8.0)
parser.add_argument('--out', default=None, help='output filename (default: auto)')
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
if len(trs) < 3:
    sys.exit(f'ERROR: missing components, found {list(trs.keys())}')

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
t = np.arange(len(z))*delta + t_off

z = bp_filter(z, args.filter_fmin, args.filter_fmax, sr, corners=4, zerophase=True)
r = bp_filter(r, args.filter_fmin, args.filter_fmax, sr, corners=4, zerophase=True)
tc = bp_filter(tc, args.filter_fmin, args.filter_fmax, sr, corners=4, zerophase=True)

# ─── Arrivals ───
taup = TauPyModel('iasp91')
dd = locations2degrees(EVLA, EVLO, sac.stla, sac.stlo)
arrs = taup.get_travel_times(source_depth_in_km=EVDP, distance_in_degree=dd, phase_list=['P','S'])
tt_P = next((a.time for a in arrs if a.name=='P'), None)
tt_S = next((a.time for a in arrs if a.name=='S'), None)
track = get_track(EVLO, EVLA, sac.stlo, sac.stla, 0.02, ref_el=-1000, plot_map=False, plot_elevation=False)
tt_r = tt_t(track, ref_el=-1000, source_depth_km=EVDP, ignore_small_islands=True)
tt_TP = tt_r["tt_t"]; tt_TS = tt_r["tt_ocean"] + tt_r["tt_land_s"]
print(f"{STA}: P={tt_P:.0f}s S={tt_S:.0f}s T-P={tt_TP:.0f}s T-S={tt_TS:.0f}s AZ={az:.1f}")

# ─── Figure ───
def calc_env(d):
    return uniform_filter1d(np.abs(hilbert(d)), size=max(int(2/delta),1))
env_z, env_r, env_tc = calc_env(z), calc_env(r), calc_env(tc)

tmin_f, tmax_f = -300, 3600
m = (t>=tmin_f)&(t<=tmax_f); tw = t[m]

fig, axes = plt.subplots(3, 2, figsize=(20, 10),
    gridspec_kw={'width_ratios':[1,1.3],'hspace':0.28,'wspace':0.30})
plt.subplots_adjust(left=0.06, right=0.98, top=0.93, bottom=0.06)

comps = [('Z',z,env_z),('R',r,env_r),('T',tc,env_tc)]
arrivals = [('P',tt_P,'#2ca02c'),('S',tt_S,'#ff7f0e'),('T-P',tt_TP,'#1f77b4'),('T-S',tt_TS,'#d62728')]
amp = 0.12
p98_global = max(np.percentile(np.abs(z), 98), np.percentile(np.abs(r), 98), np.percentile(np.abs(tc), 98)) or 1

for i,(name,data,env) in enumerate(comps):
    ax_w = axes[i,0]
    dw=data[m]; ew=env[m]
    dn=dw/p98_global*amp; en=ew/p98_global*amp
    ax_w.plot(tw,dn,'k-',lw=0.6,alpha=0.6)
    ax_w.fill_between(tw,0,en,color='red',alpha=0.10,lw=0)
    ax_w.plot(tw,en,'r-',lw=0.9,alpha=0.5)
    ax_w.axhline(0,color='gray',lw=0.2)
    ax_w.set_ylabel(name,fontsize=12,fontweight='bold')
    ax_w.set_ylim(-0.5,0.5); ax_w.set_xlim(tmin_f,tmax_f)
    ax_w.xaxis.set_major_locator(MultipleLocator(500))
    ax_w.tick_params(labelsize=7)
    if i<2: ax_w.tick_params(labelbottom=False)
    else: ax_w.set_xlabel('Time from origin (s)',fontsize=9)
    ax_w.grid(alpha=0.06,linestyle=':')
    for ph,tt_val,clr in arrivals:
        if tt_val:
            ax_w.axvline(x=tt_val,color=clr,lw=0.6,ls='--',alpha=0.4)
            ax_w.text(tt_val, ax_w.get_ylim()[1]*0.92, f'{ph}={tt_val:.0f}s', fontsize=7, fontweight='bold', color=clr, ha='center', va='top')

    ax_s = axes[i,1]
    f_spec, t_spec, Sxx_db = compute_spectrogram(data, sr)
    t_spec_rel = t_spec + t_off
    f_mask = f_spec <= 10
    vmin=np.percentile(Sxx_db[f_mask,:],5); vmax=np.percentile(Sxx_db[f_mask,:],95)
    im=ax_s.pcolormesh(t_spec_rel,f_spec[f_mask],Sxx_db[f_mask,:],cmap='inferno',vmin=vmin,vmax=vmax,shading='auto')
    ax_s.set_xlim(tmin_f,tmax_f); ax_s.set_ylim(0,10)
    ax_s.tick_params(labelsize=7)
    if i<2: ax_s.tick_params(labelbottom=False)
    else: ax_s.set_xlabel('Time from origin (s)',fontsize=9)
    if i==1: ax_s.set_ylabel('Frequency (Hz)',fontsize=9)
    ax_s.xaxis.set_major_locator(MultipleLocator(500))
    for ph,tt_val,clr in arrivals:
        if tt_val: ax_s.axvline(x=tt_val,color=clr,lw=0.6,ls='--',alpha=0.4)

cbar_ax = fig.add_axes([0.985,0.08,0.008,0.84])
cbar = fig.colorbar(im,cax=cbar_ax,orientation='vertical')
cbar.set_label('dB',fontsize=8); cbar.ax.tick_params(labelsize=6)

leg=[Line2D([0],[0],color='k',lw=0.5,label='Waveform'),
     Line2D([0],[0],color='red',lw=0.6,alpha=0.6,label='Envelope'),
     Line2D([0],[0],marker='o',color='#2ca02c',linestyle='None',ms=4,label='P'),
     Line2D([0],[0],marker='o',color='#ff7f0e',linestyle='None',ms=4,label='S'),
     Line2D([0],[0],marker='o',color='#1f77b4',linestyle='None',ms=4,label='T-P'),
     Line2D([0],[0],marker='o',color='#d62728',linestyle='None',ms=4,label='T-S')]
axes[0,0].legend(handles=leg,fontsize=6.5,loc='upper right',ncol=6,framealpha=0.85)

fig.suptitle(f'{NET}.{STA}  AZ={az:.1f}°  {args.filter_fmin}-{args.filter_fmax} Hz', fontsize=11, fontweight='bold')

out_fn = args.out or f"{OUT_DIR}/{STA}_three_filter.png"
fig.savefig(out_fn, dpi=200, bbox_inches='tight')
plt.close(fig)
print(f'OK {out_fn}')
