#!/usr/bin/env python3
"""
单台偏振分析 (Flinn 1965 协方差法)
用法:
  python3 polarization_analysis.py STA NET EVLA EVLO EVDP BASE [--resp-dir DIR]

输出4张图:
  1. P/S ±5s 偏振 (0.01-50Hz)  — P: Z-R+R-T, S: R-T
  2. P/S 全窗偏振 (0.01-50Hz)  — P: Z-R+R-T, S: R-T
  3. T-P/T-S 偏振 (2-8Hz)       — T-P: Z-R+R-T, T-S: R-T
  4. T波全窗偏振 (2-8Hz)         — T-wave: R-T

窗口:
  P: tt_P±5s, S: tt_S±5s (窄)  /  P: tt_P-30~+40, S: tt_S-30~+80 (全)
  T-P: tt_TP±10s, T-S: tt_TS-5~+15s
  T-wave: tt_TP-30~tt_TS+80
"""
import sys, os, glob, argparse
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from obspy import read
from obspy.signal.filter import bandpass as bp_filter
from obspy.taup import TauPyModel
from obspy.geodetics.base import locations2degrees
from geographiclib.geodesic import Geodesic

_codex = '/home/philiay/shared/codex/t_phase_analysis'
if os.path.isdir(_codex): sys.path.insert(0, _codex)
from t_phase_analysis.traveltime_profile import get_track, tt_t

# ─── CLI ───
parser = argparse.ArgumentParser(description='Single-station polarization (Flinn 1965)')
parser.add_argument('sta')
parser.add_argument('net')
parser.add_argument('evla', type=float)
parser.add_argument('evlo', type=float)
parser.add_argument('evdp', type=float)
parser.add_argument('base')
parser.add_argument('--resp-dir', default=None)
parser.add_argument('--sub-win', type=float, default=5.0, help='sliding sub-window (s)')
parser.add_argument('--overlap', type=float, default=0.75, help='sub-window overlap fraction')
args = parser.parse_args()

STA, NET = args.sta, args.net
EVLA, EVLO, EVDP = args.evla, args.evlo, args.evdp
BASE = args.base
RESP_DIR = args.resp_dir or os.path.join(BASE, 'sacprocess', 'remove_response')
OUT_DIR = os.path.join(BASE, 'process_flow', '3_single_station_analysis', STA)
os.makedirs(OUT_DIR, exist_ok=True)

# ─── Flinn (1965) 协方差极化 ───

def covariance_polarization(t, z, r, tc, ws, we):
    """滑动子窗口 → 3×3协方差 → 特征分解 → 方向+线性度"""
    sub_win = args.sub_win
    step = sub_win * (1 - args.overlap)
    n_sub = max(1, int((we - ws - sub_win) / step) + 1)
    if n_sub == 1:
        sub_win = we - ws

    angles_zr, weights_zr = [], []
    angles_rt, weights_rt = [], []

    for i in range(n_sub):
        t0 = ws + i * step if n_sub > 1 else ws
        t1 = min(t0 + sub_win, we)
        m = (t >= t0) & (t <= t1)
        zw, rw, tw_c = z[m], r[m], tc[m]
        if len(zw) < 3: continue

        X = np.vstack([zw-np.mean(zw), rw-np.mean(rw), tw_c-np.mean(tw_c)])
        C = X @ X.T / X.shape[1]
        ev, evec = np.linalg.eigh(C)
        l1, l2, l3 = ev[2], ev[1], ev[0]
        v1 = evec[:, 2]

        lin = 1.0 - (l2 + l3) / (2.0 * l1) if l1 > 1e-30 else 0.0
        angles_zr.append(np.arctan2(v1[1], v1[0]))
        weights_zr.append(lin)
        angles_rt.append(np.arctan2(v1[2], v1[1]))
        weights_rt.append(lin)

    return {
        'angles_zr': np.array(angles_zr), 'weights_zr': np.array(weights_zr),
        'angles_rt': np.array(angles_rt), 'weights_rt': np.array(weights_rt),
    }

def rose_from_covar(angles, weights, nbins=36):
    """协方差子窗方向 → 玫瑰图直方图"""
    bins = np.linspace(-np.pi, np.pi, nbins+1)
    centers = (bins[:-1]+bins[1:])/2
    hist = np.zeros(nbins)
    for a, w in zip(angles, weights):
        i = np.digitize(a, bins) - 1
        if i >= nbins: i = 0
        hist[i] += w
    if np.max(hist) > 0: hist /= np.max(hist)
    return centers, hist

# ─── 数据加载 ───

def load_station(sta, net, resp_dir, baz_deg, fmin=0.01, fmax=50.0):
    chs = ['BHZ','BHN','BHE']
    trs = {}
    for ch in chs:
        for pat in [f"{resp_dir}/{net}.{sta}.*.{ch}*_resp.sac",
                    f"{resp_dir}/{net}.{sta}..{ch}*_resp.sac",
                    f"{resp_dir}/{net}.{sta}.{ch}*_resp.sac"]:
            m = glob.glob(pat)
            if m: trs[ch] = read(m[0])[0]; break
    if len(trs) < 3: return None
    t0 = max(tr.stats.starttime for tr in trs.values())
    t1 = min(tr.stats.endtime for tr in trs.values())
    for ch in trs: trs[ch].trim(t0, t1)
    delta = trs['BHZ'].stats.delta; sr = 1/delta
    sac = trs['BHZ'].stats.sac; t_off = sac.b - sac.o
    z = trs['BHZ'].data.astype(np.float64)
    n = trs['BHN'].data.astype(np.float64)
    e = trs['BHE'].data.astype(np.float64)
    br = np.radians(baz_deg)
    # 2026-09-02: 统一 obspy 约定（R 背离震源）
    from obspy.signal.rotate import rotate_ne_rt
    r, tc = rotate_ne_rt(n, e, float(baz_deg))
    t = np.arange(len(z))*delta + t_off
    z_bp = bp_filter(z, fmin, fmax, sr, corners=4, zerophase=True)
    r_bp = bp_filter(r, fmin, fmax, sr, corners=4, zerophase=True)
    tc_bp = bp_filter(tc, fmin, fmax, sr, corners=4, zerophase=True)
    return {'t': t, 'z': z_bp, 'r': r_bp, 'tc': tc_bp, 'delta': delta}

# ─── 主绘图 ───

def plot_polarization(sta, net, resp_dir, windows, baz_deg, out_fn, tt_P, tt_S, tt_TP, tt_TS, fmin=0.01, fmax=50.0):
    d = load_station(sta, net, resp_dir, baz_deg, fmin=fmin, fmax=fmax)
    if d is None: return
    t, z, r, tc = d['t'], d['z'], d['r'], d['tc']

    rose_list = []
    for wname, ws, we, _ in windows:
        cov = covariance_polarization(t, z, r, tc, ws, we)
        m_win = (t>=ws) & (t<=we)
        ez = np.sum(z[m_win]**2); er = np.sum(r[m_win]**2); et_val = np.sum(tc[m_win]**2)
        er_et = er/et_val if et_val>0 else 0
        ez_eh = ez/(er+et_val) if (er+et_val)>0 else 0

        # R-T rose (always)
        rc, rh = rose_from_covar(cov['angles_rt'], cov['weights_rt'])
        rose_list.append({'name': wname, 'plane': 'R-T', 'ws': ws, 'we': we,
                          'er_et': er_et, 'ez_eh': ez_eh, 'rose_c': rc, 'rose_h': rh,
                          'refs': [(0,'R','green'),(np.pi/2,'T','red'),(np.pi,'R','green'),(-np.pi/2,'T','red')]})
        # Z-R rose (only P, T-P)
        if wname in ('P', 'T-P'):
            rc, rh = rose_from_covar(cov['angles_zr'], cov['weights_zr'])
            rose_list.append({'name': wname, 'plane': 'Z-R', 'ws': ws, 'we': we,
                              'er_et': er_et, 'ez_eh': ez_eh, 'rose_c': rc, 'rose_h': rh,
                              'refs': [(0,'Z','blue'),(np.pi/2,'R','green'),(np.pi,'Z','blue'),(-np.pi/2,'R','green')]})

    n_roses = len(rose_list); n_waveforms = len(windows); ncols = max(n_roses, n_waveforms)
    fig = plt.figure(figsize=(6*ncols, 8.5))
    gs = gridspec.GridSpec(2, ncols, figure=fig, hspace=0.30, wspace=0.25,
                           height_ratios=[2, 1.2], left=0.04, right=0.98, top=0.92, bottom=0.06)

    for i_r, rd in enumerate(rose_list):
        ax = fig.add_subplot(gs[0, i_r], polar=True)
        wname = rd['name']; color = '#1f77b4' if ('T-P' in wname or wname=='P') else '#d62728'
        ax.bar(rd['rose_c'], rd['rose_h'], width=np.pi/36, bottom=0.1, color=color, alpha=0.75, edgecolor='w', lw=0.3)
        for ang, lbl, cl in rd['refs']:
            ax.plot([ang,ang],[0,1.2],color=cl,lw=0.8,alpha=0.5)
            if ang<=0: ax.text(ang,1.35,lbl,fontsize=8,fontweight='bold',ha='center',color=cl)
        ax.set_title(f'{wname} ({rd["plane"]})\n{rd["ws"]:.0f}-{rd["we"]:.0f}s\nER/ET={rd["er_et"]:.2f}  EZ/EH={rd["ez_eh"]:.3f}',
                     fontsize=9, fontweight='bold', pad=18)
        ax.tick_params(labelsize=6); ax.set_theta_zero_location('N'); ax.set_theta_direction(-1); ax.set_ylim(0,1.5)

    arrivals_map = {'P':tt_P, 'S':tt_S, 'T-P':tt_TP, 'T-S':tt_TS}
    for i_w, (wname, ws, we, _) in enumerate(windows):
        ax_wf = fig.add_subplot(gs[1, i_w])
        m_wf = (t>=ws)&(t<=we)
        tw_w, zw_w, rw_w, tw_w2 = t[m_wf], z[m_wf], r[m_wf], tc[m_wf]
        win_p98 = max(np.percentile(np.abs(zw_w),98), np.percentile(np.abs(rw_w),98), np.percentile(np.abs(tw_w2),98)) or 1
        amp=0.15; off=0.8
        ax_wf.plot(tw_w, zw_w/win_p98*amp+off,'k-',lw=0.5)
        ax_wf.plot(tw_w, rw_w/win_p98*amp,'k-',lw=0.5)
        ax_wf.plot(tw_w, tw_w2/win_p98*amp-off,'k-',lw=0.5)
        dx = (we-ws)*0.03
        ax_wf.text(ws-dx,off,'Z',fontsize=8,fontweight='bold',ha='right',va='center')
        ax_wf.text(ws-dx,0,'R',fontsize=8,fontweight='bold',ha='right',va='center')
        ax_wf.text(ws-dx,-off,'T',fontsize=8,fontweight='bold',ha='right',va='center')
        for yv in[off,0,-off]: ax_wf.axhline(y=yv,c='gray',lw=0.2,alpha=0.3)
        for pname,pt in arrivals_map.items():
            if ws<=pt<=we:
                clr = {'P':'#2ca02c','S':'#ff7f0e','T-P':'#1f77b4','T-S':'#d62728'}[pname]
                ax_wf.axvline(pt,c=clr,lw=0.8,ls='--',alpha=0.6)
                ax_wf.text(pt,ax_wf.get_ylim()[1]*0.92,f'{pname}={pt:.0f}s',fontsize=7,fontweight='bold',
                          c=clr,ha='center',va='top',bbox=dict(boxstyle='round,pad=0.2',fc='white',alpha=0.7))
        ax_wf.set_xlim(ws,we); ax_wf.set_ylabel('Amp(norm)',fontsize=9)
        ax_wf.set_title(f'{wname} ({ws:.0f}-{we:.0f}s)',fontsize=10,fontweight='bold')
        ax_wf.tick_params(labelsize=7); ax_wf.grid(True,alpha=0.08,linestyle=':')

    for i_e in range(n_waveforms, ncols):
        fig.add_subplot(gs[1,i_e]).axis('off')
    fig.savefig(out_fn, dpi=200, bbox_inches='tight'); plt.close(fig)
    print(f'OK {out_fn}')

# ─── Main ───
tr = read(glob.glob(f"{RESP_DIR}/{NET}.{STA}.*.BHZ*_resp.sac")[0])[0]
sac = tr.stats.sac; stla, stlo = sac.stla, sac.stlo
dd = locations2degrees(EVLA, EVLO, stla, stlo)
g = Geodesic.WGS84.Inverse(EVLA, EVLO, stla, stlo)
az = g['azi1']; az = az+360 if az<0 else az; baz = (az+180)%360
taup = TauPyModel('iasp91')
arrs = taup.get_travel_times(source_depth_in_km=EVDP, distance_in_degree=dd, phase_list=['P','S'])
tt_P = next((a.time for a in arrs if a.name=='P'), 0)
tt_S = next((a.time for a in arrs if a.name=='S'), 0)
track = get_track(EVLO, EVLA, stlo, stla, 0.02, ref_el=-1000, plot_map=False, plot_elevation=False)
tt_r = tt_t(track, ref_el=-1000, source_depth_km=EVDP, ignore_small_islands=True)
tt_TP = tt_r["tt_t"]; tt_TS = tt_r["tt_ocean"] + tt_r["tt_land_s"]
print(f"{STA}: P={tt_P:.0f}s S={tt_S:.0f}s T-P={tt_TP:.0f}s T-S={tt_TS:.0f}s AZ={az:.1f} BAZ={baz:.1f}")

# 1. P/S ±5s (0.01-50 Hz)
wps = [('P', tt_P-5, tt_P+5, 'Z-R'), ('S', tt_S-5, tt_S+5, 'R-T')]
plot_polarization(STA, NET, RESP_DIR, wps, baz, f"{OUT_DIR}/{STA}_p_s_polar_P{tt_P-5:.0f}-{tt_P+5:.0f}_S{tt_S-5:.0f}-{tt_S+5:.0f}.png", tt_P, tt_S, tt_TP, tt_TS)

# 2. P/S full (0.01-50 Hz)
wpf = [('P', tt_P-30, tt_P+40, 'Z-R'), ('S', tt_S-30, tt_S+80, 'R-T')]
plot_polarization(STA, NET, RESP_DIR, wpf, baz, f"{OUT_DIR}/{STA}_p_s_fullwindow.png", tt_P, tt_S, tt_TP, tt_TS)

# 3. T-P/T-S (2-8 Hz): TP±10s, TS(-5~+15s)
wts = [('T-P', tt_TP-10, tt_TP+10, 'Z-R'), ('T-S', tt_TS-5, tt_TS+15, 'R-T')]
plot_polarization(STA, NET, RESP_DIR, wts, baz, f"{OUT_DIR}/{STA}_tp_ts_polar_TP{tt_TP-10:.0f}-{tt_TP+10:.0f}_TS{tt_TS-5:.0f}-{tt_TS+15:.0f}.png", tt_P, tt_S, tt_TP, tt_TS, fmin=2.0, fmax=8.0)

# 4. T波全窗 (2-8 Hz)
wtf = [('T-wave', tt_TP-30, tt_TS+80, 'R-T')]
plot_polarization(STA, NET, RESP_DIR, wtf, baz, f"{OUT_DIR}/{STA}_tp_ts_fullwindow.png", tt_P, tt_S, tt_TP, tt_TS, fmin=2.0, fmax=8.0)
