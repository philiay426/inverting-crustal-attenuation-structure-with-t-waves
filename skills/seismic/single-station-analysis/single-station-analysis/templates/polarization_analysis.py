#!/usr/bin/env python3
"""
单台偏振分析: P/S 理论到时窗 (±2 s) 的 R-T 平面偏振
- 输入: 去响应后 ZNE SAC 文件
- 默认输出: 只出 {STA}_p_s_polar_P*_S*.png 一张(驱动重命名为 {STA}_PS_polar.png)
  其余旧图 (p_s_fullwindow / tp_ts_polar / p_s_onset10s) 需加 --all
- 方法开关 --method: cov | amp | count, 默认 count(采样点等权/投票);
  优先级: --method > 环境变量 POLAR_METHOD > 默认 count。方向+权重见 polar_methods.py
- P波窗: P±2s, R-T 平面 (0.01-50 Hz)
- S波窗: S±2s, R-T 平面 (0.01-50 Hz)
- T-P窗: T-P±10s, R-T 平面 (2-8 Hz)      (--all)
- T-S窗: T-S-5/+15s, R-T 平面 (2-8 Hz)   (--all)
- cov 路径算式在 polar_methods.cov_subwindow_direction_weights() 中逐位保留,
  数值与旧版完全一致(向后兼容)
- 画图标准: 见 seismic/plotting-standards skill
- v3: 移除 Z-R/Z-T 平面，仅保留 R-T 平面
"""



import sys, os, glob, argparse
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from obspy import read, UTCDateTime
from obspy.signal.filter import bandpass as bp_filter
from obspy.taup import TauPyModel
from obspy.geodetics.base import locations2degrees
from geographiclib.geodesic import Geodesic
from scipy.signal import hilbert, welch

for _codex in ('/home/sy/shared/codex/t_phase_analysis', '/home/sy/codex/t_phase_analysis',
               '/home/sy/codex'):
    if os.path.isdir(_codex): sys.path.insert(0, _codex)
from t_phase_analysis.traveltime_profile import get_track, tt_t

# 共享模块 polar_methods.py 与本脚本同目录(驱动脚本会一起拷到台站目录)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import polar_methods
from polar_methods import METHOD_COV

EVLA, EVLO, EVDP = 19.0428, -66.8096, 20.0  # ← 修改为事件坐标
EV_ORIGIN = UTCDateTime(2019, 1, 22, 19, 1, 43)  # ← 修改为事件发震时刻（SAC 缺 o 字段时必需）

# ─── Flinn (1965) 协方差极化 (cov 路径) ───
# 算式已移入 polar_methods.cov_subwindow_direction_weights(), 逐位保留、数值不变；
# 这里只保留同名薄封装, 供外部沿用旧调用方式。

def covariance_polarization(t, z, r, tc, ws, we, sub_win=0.5, overlap=0.75):
    """Flinn 1965: 滑动协方差 PCA, 仅 R-T 平面
    返回 {'angles_rt': rad, 'weights_rt': lin*sqrt(lambda1)}"""
    angles_rt, weights_rt = polar_methods.cov_subwindow_direction_weights(
        t, z, r, tc, ws, we, sub_win=sub_win, overlap=overlap)
    return {'angles_rt': angles_rt, 'weights_rt': weights_rt}

def rose_from_covar(angles, weights, nbins=36):
    return polar_methods.rose_hist(angles, weights, nbins=nbins)

def compute_snr_db(z_data, t, delta, tt_P, tt_TP, tt_TS):
    """Welch PSD SNR: 10*log10(E_sig/E_noise) dB over 2-8 Hz
    Signal: T-P ~ T-S+80s, Noise: P-300 ~ P-60s"""
    sr = 1/delta
    m_sig = (t >= tt_TP) & (t <= tt_TS + 80)
    m_noi = (t >= tt_P - 300) & (t <= tt_P - 60)
    sig_data = z_data[m_sig]; noi_data = z_data[m_noi]
    if len(sig_data) < 64 or len(noi_data) < 64: return 0
    nperseg = min(256, min(len(sig_data), len(noi_data)) // 2)
    f_sig, p_sig = welch(sig_data, sr, nperseg=nperseg)
    f_noi, p_noi = welch(noi_data, sr, nperseg=nperseg)
    fm = (f_sig >= 2) & (f_sig <= 8)
    e_sig = np.trapezoid(p_sig[fm], f_sig[fm])
    e_noi = np.trapezoid(p_noi[fm], f_noi[fm])
    return 10 * np.log10(e_sig / e_noi) if e_noi > 0 else 0

def load_station(sta, net, resp_dir, baz_deg, fmin=0.01, fmax=50.0):
    chs = ['BHZ','BHN','BHE']; trs = {}
    for ch in chs:
        for pat in [f"{resp_dir}/{net}.{sta}.*.{ch}*_resp.sac",
                    f"{resp_dir}/{net}.{sta}..{ch}*_resp.sac",
                    f"{resp_dir}/{net}.{sta}.{ch}*_resp.sac"]:
            m = glob.glob(pat)
            if m: trs[ch] = read(m[0])[0]; break
    # Fallback: IU/US networks use BH1/BH2 instead of BHN/BHE
    if 'BHN' not in trs or 'BHE' not in trs:
        for ch_src, ch_dst in [('BH1','BHN'), ('BH2','BHE')]:
            if ch_dst not in trs:
                for pat in [f"{resp_dir}/{net}.{sta}.*.{ch_src}*_resp.sac",
                            f"{resp_dir}/{net}.{sta}..{ch_src}*_resp.sac",
                            f"{resp_dir}/{net}.{sta}.{ch_src}*_resp.sac"]:
                    m = glob.glob(pat)
                    if m: trs[ch_dst] = read(m[0])[0]; break
    if len(trs) < 3: return None
    t0 = max(tr.stats.starttime for tr in trs.values())
    t1 = min(tr.stats.endtime for tr in trs.values())
    for ch in trs: trs[ch].trim(t0, t1)
    delta = trs['BHZ'].stats.delta; sr = 1/delta
    sac = trs['BHZ'].stats.sac
    # 时间轴: SAC o 字段可靠时用 b-o (数据常从发震前开始记录)
    o_val = getattr(sac, 'o', None)
    if o_val is not None and 0 < abs(o_val) < 1e6:
        t_off = sac.b - o_val
    else:
        t_off = trs['BHZ'].stats.starttime - EV_ORIGIN
    z = trs['BHZ'].data.astype(np.float64)
    n = trs['BHN'].data.astype(np.float64); e = trs['BHE'].data.astype(np.float64)
    br = np.radians(baz_deg)
    r = n*np.cos(br) + e*np.sin(br); tc = -n*np.sin(br) + e*np.cos(br)
    t = np.arange(len(z))*delta + t_off
    z_bp = bp_filter(z, fmin, fmax, sr, corners=4, zerophase=True)
    r_bp = bp_filter(r, fmin, fmax, sr, corners=4, zerophase=True)
    tc_bp = bp_filter(tc, fmin, fmax, sr, corners=4, zerophase=True)
    return {'t': t, 'z': z_bp, 'r': r_bp, 'tc': tc_bp, 'delta': delta}

# ─── 标准偏振绘图 ───

def plot_polarization(sta, net, resp_dir, windows, baz_deg, out_fn, tt_P, tt_S, tt_TP, tt_TS,
                      fmin=0.01, fmax=50.0, station_info=None, method=polar_methods.DEFAULT_METHOD):
    d = load_station(sta, net, resp_dir, baz_deg, fmin=fmin, fmax=fmax)
    if d is None: return
    t, z, r, tc = d['t'], d['z'], d['r'], d['tc']
    sub_win = max(1.25, min(5.0, min(w[2]-w[1] for w in windows)/3))
    method_lbl = polar_methods.method_label(method)
    rose_list = []
    for wname, ws, we, _ in windows:
        if method == METHOD_COV:
            # cov: 子窗 3x3 协方差最大特征向量方位, 权重 lin*sqrt(lambda1) (旧口径)
            ang_win, w_win = polar_methods.cov_subwindow_direction_weights(
                t, z, r, tc, ws, we, sub_win=sub_win)
        else:
            # amp/count: 逐采样点 atan2(E,N)(=atan2(T,R), R-T 系), 权重见 polar_methods
            m_win0 = (t>=ws)&(t<=we)
            ang_win, w_win = polar_methods.sample_direction_weights(method, r[m_win0], tc[m_win0])
        m_win = (t>=ws)&(t<=we)
        ez = np.sum(z[m_win]**2); er = np.sum(r[m_win]**2); et_val = np.sum(tc[m_win]**2)
        er_et = er/et_val if et_val>0 else 0; ez_eh = ez/(er+et_val) if (er+et_val)>0 else 0
        rc, rh = rose_from_covar(ang_win, w_win)
        sec = polar_methods.sector_fractions(ang_win, w_win)
        rose_list.append({'name':wname,'plane':'R-T','ws':ws,'we':we,'er_et':er_et,'ez_eh':ez_eh,
                          'rose_c':rc,'rose_h':rh,'sec':sec,
                          'refs':[(0,'R','green'),(np.pi/2,'T','red'),(np.pi,'R','green'),(-np.pi/2,'T','red')]})
    n_roses = len(rose_list); n_waveforms = len(windows); ncols = max(n_roses, n_waveforms)
    fig = plt.figure(figsize=(6*ncols, 8.5))
    gs = gridspec.GridSpec(2, ncols, figure=fig, hspace=0.30, wspace=0.25,
                           height_ratios=[2,1.2], left=0.04, right=0.98, top=0.92, bottom=0.06)
    for i_r, rd in enumerate(rose_list):
        ax = fig.add_subplot(gs[0,i_r], polar=True)
        color = '#1f77b4' if ('T-P' in rd['name'] or rd['name']=='P') else '#d62728'
        ax.bar(rd['rose_c'], rd['rose_h'], width=np.pi/36, bottom=0.1, color=color, alpha=0.75, edgecolor='w', lw=0.3)
        for ang, lbl, cl in rd['refs']:
            ax.plot([ang,ang],[0,1.2],color=cl,lw=0.8,alpha=0.5)
            if ang<=0: ax.text(ang,1.35,lbl,fontsize=12,fontweight='bold',ha='center',color=cl)
        ax.set_title(f'{rd["name"]} ({rd["plane"]}) [{method_lbl}]\n{rd["ws"]:.0f}-{rd["we"]:.0f}s\nER/ET={rd["er_et"]:.2f}  EZ/EH={rd["ez_eh"]:.3f}',
                     fontsize=12, fontweight='bold', pad=18)
        ax.tick_params(labelsize=9); ax.set_theta_zero_location('N'); ax.set_theta_direction(-1); ax.set_ylim(0,1.5)
    arrivals_map = {'P':tt_P,'S':tt_S,'T-P':tt_TP,'T-S':tt_TS}
    for i_w, (wname, ws, we, _) in enumerate(windows):
        ax_wf = fig.add_subplot(gs[1,i_w])
        m_wf = (t>=ws)&(t<=we); tw_w, zw_w, rw_w, tw_w2 = t[m_wf], z[m_wf], r[m_wf], tc[m_wf]
        win_p98 = max(np.percentile(np.abs(zw_w),98), np.percentile(np.abs(rw_w),98), np.percentile(np.abs(tw_w2),98)) or 1
        amp=0.15; off=0.8
        ax_wf.plot(tw_w, zw_w/win_p98*amp+off,'k-',lw=0.6)
        ax_wf.plot(tw_w, rw_w/win_p98*amp,'k-',lw=0.6)
        ax_wf.plot(tw_w, tw_w2/win_p98*amp-off,'k-',lw=0.6)
        dx = (we-ws)*0.03
        for yv,lab in[(off,'Z'),(0,'R'),(-off,'T')]:
            ax_wf.text(ws-dx,yv,lab,fontsize=12,fontweight='bold',ha='right',va='center')
            ax_wf.axhline(y=yv,c='gray',lw=0.2,alpha=0.3)
        for pname,pt in arrivals_map.items():
            if ws<=pt<=we:
                clr = {'P':'#2ca02c','S':'#ff7f0e','T-P':'#1f77b4','T-S':'#d62728'}[pname]
                ax_wf.axvline(pt,c=clr,lw=0.8,ls='--',alpha=0.6)
                ax_wf.text(pt,ax_wf.get_ylim()[1]*0.92,f'{pname}={pt:.0f}s',fontsize=11,fontweight='bold',
                          c=clr,ha='center',va='top',bbox=dict(boxstyle='round,pad=0.2',fc='white',alpha=0.7))
        ax_wf.set_xlim(ws,we); ax_wf.set_ylabel('Amp (norm)',fontsize=12)
        ax_wf.set_title(f'{wname} ({ws:.0f}-{we:.0f}s)',fontsize=14,fontweight='bold')
        ax_wf.tick_params(labelsize=10); ax_wf.grid(True,alpha=0.08,linestyle=':')
    for i_e in range(n_waveforms, ncols):
        if i_e == ncols - 1 and station_info is not None:
            ax_info = fig.add_subplot(gs[1, i_e]); ax_info.axis('off')
            band_str = f'{fmin:.2f}-{fmax:.0f} Hz' if fmax < 50 else '0.01-50 Hz'
            info_lines = [f'{sta} ({net})',
                          f'Land: {station_info["land_km"]:.0f} km',
                          f'Ocean: {station_info["ocean_km"]:.0f} km',
                          f'AZ: {station_info["az"]:.1f}°',
                          f'SNR: {station_info["snr_db"]:.1f} dB',
                          f'Band: {band_str}',
                          f'Method: {method_lbl}']
            ax_info.text(0.05, 0.95, '\n'.join(info_lines),
                        transform=ax_info.transAxes, fontsize=16, fontweight='bold',
                        verticalalignment='top', linespacing=1.4,
                        bbox=dict(boxstyle='round,pad=0.6', fc='whitesmoke', alpha=0.8))
        else:
            fig.add_subplot(gs[1, i_e]).axis('off')
    fig.savefig(out_fn, dpi=200, bbox_inches='tight'); plt.close(fig)
    print(f'OK {out_fn}')

# ─── CLI ───────────────────────────────────────────────────────────────────
# 优先级: --method > 环境变量 POLAR_METHOD > 默认 count; --out-dir > OUT_DIR > 脚本内 OUT
def _parse_args(argv=None):
    p = argparse.ArgumentParser(description='P/S (and legacy T-P/T-S) R-T plane polarization')
    polar_methods.add_common_args(p)
    p.add_argument('--all', action='store_true',
                   help='also write legacy extra figures: p_s_fullwindow, tp_ts_polar, p_s_onset10s')
    return p.parse_args(argv)


ARGS = _parse_args()
METHOD = polar_methods.resolve_method(ARGS.method)

# ─── Main ──
STA, NET = 'K62A', 'TA'  # ← 修改台站
BASE = "/mnt/f/sharedate/case/2014-01-13-mw63-process"
RESP_DIR = os.path.join(BASE, "sacprocess", "remove_response")
OUT = os.path.join(BASE, "process_flow", "2_single_station_analysis", STA)
OUT = polar_methods.resolve_out_dir(ARGS.out_dir, OUT)
os.makedirs(OUT, exist_ok=True)

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
print(f"{STA}: P={tt_P:.0f}s S={tt_S:.0f}s T-P={tt_TP:.0f}s T-S={tt_TS:.0f}s AZ={az:.1f} BAZ={baz:.1f} method={METHOD}")

# ── SNR (Welch PSD dB) ──
land_km = tt_r.get('dist_land', 0) * 111.0
ocean_km = tt_r.get('dist_ocean', 0) * 111.0
dist_km = dd * 111.19
d_snr = load_station(STA, NET, RESP_DIR, baz, fmin=2.0, fmax=8.0)
if d_snr:
    snr_db = compute_snr_db(d_snr['z'], d_snr['t'], d_snr['delta'], tt_P, tt_TP, tt_TS)
else:
    snr_db = 0
station_info = {'land_km': land_km, 'ocean_km': ocean_km, 'dist_km': dist_km,
                'az': az, 'baz': baz, 'snr_db': snr_db}

# 1. P/S ±2s 0.01-50Hz —— **默认唯一输出**(驱动重命名为 {STA}_PS_polar.png)
wps = [('P', tt_P-2, tt_P+2, 'R-T'), ('S', tt_S-2, tt_S+2, 'R-T')]
plot_polarization(STA, NET, RESP_DIR, wps, baz,
    f"{OUT}/{STA}_p_s_polar_P{tt_P-2:.0f}-{tt_P+2:.0f}_S{tt_S-2:.0f}-{tt_S+2:.0f}.png",
    tt_P, tt_S, tt_TP, tt_TS, station_info=station_info, method=METHOD)

if ARGS.all:
    # 2. P/S full 0.01-50Hz
    wpf = [('P', tt_P, tt_P+35, 'R-T'), ('S', tt_S, tt_S+80, 'R-T')]
    plot_polarization(STA, NET, RESP_DIR, wpf, baz, f"{OUT}/{STA}_p_s_fullwindow.png",
        tt_P, tt_S, tt_TP, tt_TS, station_info=station_info, method=METHOD)

    # 3. T-P/T-S 2-8Hz
    wts = [('T-P', tt_TP-10, tt_TP+10, 'R-T'), ('T-S', tt_TS-5, tt_TS+15, 'R-T')]
    plot_polarization(STA, NET, RESP_DIR, wts, baz,
        f"{OUT}/{STA}_tp_ts_polar_TP{tt_TP-10:.0f}-{tt_TP+10:.0f}_TS{tt_TS-5:.0f}-{tt_TS+15:.0f}.png",
        tt_P, tt_S, tt_TP, tt_TS, fmin=2.0, fmax=8.0, station_info=station_info, method=METHOD)

    # 4. P/S onset+10s 0.01-50Hz
    wps10 = [('P', tt_P+15, tt_P+25, 'R-T'), ('S', tt_S, tt_S+10, 'R-T')]
    plot_polarization(STA, NET, RESP_DIR, wps10, baz,
        f"{OUT}/{STA}_p_s_onset10s_P{tt_P+15:.0f}-{tt_P+25:.0f}_S{tt_S:.0f}-{tt_S+10:.0f}.png",
        tt_P, tt_S, tt_TP, tt_TS, station_info=station_info, method=METHOD)
