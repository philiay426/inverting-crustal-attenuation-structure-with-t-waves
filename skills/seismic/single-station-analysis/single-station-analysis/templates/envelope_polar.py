#!/usr/bin/env python3
"""
T-P/T-S 包络拾取偏振分析 v2 (2026-08-04 定稿)
- 包络检测: 峰值法 + 30%阈值, 偏振窗口 ±5s, T-P/T-S 不重叠
- 偏振平面: 仅 R-T (去掉 Z-R)
- 波形: T波三分量 (Z/R/T) 窗口 T-P-100 ~ T-S+200
- 样式: 波形灰色, 包络红色平滑8s, 包络范围用彩色方块标注
- 不画 P/S 窗口
用法: python3 envelope_polar.py <STA> <NET>
"""
import sys, os, glob, csv
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from obspy import read, UTCDateTime
from obspy.signal.filter import bandpass, envelope
from obspy.signal.util import smooth
from scipy.signal import welch

# ====== 事件参数 (改这里) ======
EVLA, EVLO, EVDP = -43.1219, 42.3568, 13.0
EV_ORIGIN = UTCDateTime(2019, 1, 22, 19, 1, 43)  # ← 发震时刻（SAC 缺 o 字段时必需）
BASE = "/mnt/f/sharedate/case/Southern_Africa/2019-01-22-mww67-prince-edward-islands-station_low"
# ================================
RESP = f'{BASE}/sacprocess/remove_response'
OUT = f'{BASE}/process_flow/2_single_station_analysis'
CSV = f'{BASE}/process_flow/1_station_distribution/stations.csv'

SMOOTH_S = 8.0    # 包络平滑 8s
THRESH    = 0.30  # 包络阈值 30%
HALF      = 5.0   # 偏振窗口 ±5s
SUB_WIN   = 1.0   # 子窗 1s
OVERLAP   = 0.75  # 子窗重叠 75%

def polarize_rt(t, z, r, tc, ws, we, sub_win=SUB_WIN, overlap=OVERLAP):
    """仅 R-T 平面滑动 PCA (Flinn 1965)"""
    step = sub_win*(1-overlap)
    n = max(1, int((we-ws-sub_win)/step)+1)
    if n == 1: sub_win = we-ws
    arts, weights, lins, plas = [], [], [], []
    for i in range(n):
        t0=ws+i*step; t1=min(t0+sub_win,we); m=(t>=t0)&(t<=t1)
        zz,rw,tw=z[m],r[m],tc[m]
        if len(zz)<3: continue
        X=np.vstack([zz-np.mean(zz),rw-np.mean(rw),tw-np.mean(tw)])
        C=X@X.T/X.shape[1]; ev,evc=np.linalg.eigh(C); l3,l2,l1=ev; v1=evc[:,2]
        lin=1-(l2+l3)/(2*l1) if l1>1e-30 else 0
        pla=1-(2*l3)/(l1+l2) if l1+l2>1e-30 else 0
        w=np.sqrt(l1)
        arts.append(np.arctan2(v1[2],v1[1])); weights.append(w)
        lins.append(lin); plas.append(pla)
    return arts, weights, (np.mean(lins) if lins else 0), (np.mean(plas) if plas else 0), n

def rose(ax, angles, weights, title, refs):
    bins=np.linspace(-np.pi,np.pi,37); h,_=np.histogram(angles,bins=bins,weights=weights)
    th=(bins[:-1]+bins[1:])/2
    ax.bar(th, h, width=np.pi/18, alpha=0.7, color='#1f77b4')
    mx=max(h)*1.3 if len(h)>0 else 1
    for a,l,c in refs: ax.plot([a,a],[0,mx],c,ls='--',lw=1); ax.text(a,mx*1.15,l,c=c,fontsize=10,ha='center')
    ax.set_title(title,fontsize=12,fontweight='bold')
    ax.set_theta_zero_location('N'); ax.set_theta_direction(-1)

def pick_peak(t, env, t0, t1, half=HALF, thr=THRESH):
    """峰值法拾取 + 30%阈值包络范围, 限制 ±half, 返回 (t_obs, w0, w1)"""
    m=(t>=t0)&(t<=t1)
    tseg, eseg = t[m], env[m]
    if len(tseg)==0: return None, None, None
    i0=np.argmax(eseg); t_obs=tseg[i0]; emax=eseg[i0]
    mw=(tseg>=t_obs-half)&(tseg<=t_obs+half)
    e2=eseg[mw]; t2=tseg[mw]
    above=e2 > thr*emax
    if above.any():
        idx=np.where(above)[0]
        w0, w1 = t2[idx[0]], t2[idx[-1]]
    else:
        w0, w1 = t_obs-half, t_obs+half
    return t_obs, w0, w1

def process_sta(sta, net):
    print(f'\n{net}.{sta}')
    with open(CSV) as f:
        for row in csv.DictReader((l for l in f if not l.startswith('#'))):
            if row['sta']==sta:
                land=float(row['land_km']); ttp=float(row['tt_tp']); tts=float(row['tt_ts'])
                tp=float(row['tt_p']); ts=float(row['tt_s']); bt=float(row['baz'])
                break

    trs={}
    for ch in ['BHZ','BHN','BHE','BH1','BH2']:
        for pat in [f'{RESP}/{net}.{sta}.*.{ch}*_resp.sac',f'{RESP}/{net}.{sta}..{ch}*_resp.sac']:
            m=glob.glob(pat)
            if m: trs[ch]=read(m[0])[0]; break
    if 'BHN' not in trs and 'BH1' in trs: trs['BHN']=trs['BH1']
    if 'BHE' not in trs and 'BH2' in trs: trs['BHE']=trs['BH2']

    t0=max(tr.stats.starttime for tr in trs.values())
    t1=min(tr.stats.endtime for tr in trs.values())
    for ch in trs: trs[ch].trim(t0,t1)

    delta=trs['BHZ'].stats.delta; sr=1/delta
    sac=trs['BHZ'].stats.sac
    o_val = getattr(sac, 'o', None)
    if o_val is not None:
        t_off = sac.b - o_val
    else:
        # ⚠️ 缺 o 字段：时间轴必须相对发震时刻（stations.csv 到时即相对发震）
        t_off = trs['BHZ'].stats.starttime - EV_ORIGIN
    t=np.arange(trs['BHZ'].stats.npts)*delta+t_off
    z=trs['BHZ'].data.astype(np.float64)
    n=trs['BHN'].data.astype(np.float64)
    e=trs['BHE'].data.astype(np.float64)

    br=np.radians(bt); r=n*np.cos(br)+e*np.sin(br); tc=-n*np.sin(br)+e*np.cos(br)
    zf=bandpass(z,2,8,df=sr,corners=4)
    rf=bandpass(r,2,8,df=sr,corners=4)
    tf=bandpass(tc,2,8,df=sr,corners=4)

    # 包络: 平滑 8s
    env_z = smooth(envelope(zf), max(1,int(SMOOTH_S*sr)))
    env_r = smooth(envelope(rf), max(1,int(SMOOTH_S*sr)))
    env_t = smooth(envelope(tf), max(1,int(SMOOTH_S*sr)))

    # 峰值法拾取 (T-S 从 T-P 观测+5s 后搜索, 避免同一峰)
    t_obs_tp, w_tp0, w_tp1 = pick_peak(t, env_z, ttp, ttp+40)
    ts_start = max(tts, (t_obs_tp or ttp) + 5)
    t_obs_ts, w_ts0, w_ts1 = pick_peak(t, env_z, ts_start, ts_start+30)

    # 窗口不重叠检查
    w_tp=(w_tp0, w_tp1); w_ts=(w_ts0, w_ts1)
    if w_tp[1] > w_ts[0]:
        shift = w_tp[1] - w_ts[0] + 0.5
        w_ts = (w_ts[0]+shift, w_ts[1]+shift)
        w_ts0, w_ts1 = w_ts
    print(f'  T-P: theory={ttp:.0f}s obs={t_obs_tp:.0f}s 窗口=[{w_tp[0]:.1f},{w_tp[1]:.1f}]s')
    print(f'  T-S: theory={tts:.0f}s obs={t_obs_ts:.0f}s 窗口=[{w_ts[0]:.1f},{w_ts[1]:.1f}]s')

    m_sig=(t>=ttp)&(t<=tts+80); m_noi=(t>=tp-300)&(t<=tp-60)
    fs,ps=welch(zf[m_sig],fs=sr,nperseg=256); fn,pn=welch(zf[m_noi],fs=sr,nperseg=256)
    snr=10*np.log10(np.sum(ps[(fs>=2)&(fs<=8)])/np.sum(pn[(fn>=2)&(fn<=8)]))

    # 偏振 (仅 R-T)
    res_tp=polarize_rt(t,zf,rf,tf,*w_tp)
    res_ts=polarize_rt(t,zf,rf,tf,*w_ts)
    print(f'  T-P: {res_tp[4]}子窗, R-T角中值={np.degrees(np.median(res_tp[0])):.1f}°, lin={res_tp[2]:.3f}')
    print(f'  T-S: {res_ts[4]}子窗, R-T角中值={np.degrees(np.median(res_ts[0])):.1f}°, lin={res_ts[2]:.3f}')

    # ── 图: 2行×3列 ──
    fig=plt.figure(figsize=(18, 10))
    fig.suptitle(f'{net}.{sta}  |  Land={land:.0f}km  |  SNR={snr:.1f}dB  |  '
                 f'T-P th={ttp:.0f}s obs={t_obs_tp:.0f}s  T-S th={tts:.0f}s obs={t_obs_ts:.0f}s  |  '
                 f'pick: peak+{int(THRESH*100)}%, ±{HALF:.0f}s, sub={SUB_WIN:.0f}s/{OVERLAP*100:.0f}%',
                 fontsize=13, fontweight='bold', y=0.99)

    # Row1: T-P R-T玫瑰 | T-S R-T玫瑰 | Z包络拾取(带彩色方块)
    ax1=fig.add_subplot(2,3,1,polar=True)
    rose(ax1,res_tp[0],res_tp[1],f'T-P R-T ({res_tp[4]}sub)',[(0,'R','green'),(np.pi/2,'T','red')])
    ax2=fig.add_subplot(2,3,2,polar=True)
    rose(ax2,res_ts[0],res_ts[1],f'T-S R-T ({res_ts[4]}sub)',[(0,'R','green'),(np.pi/2,'T','red')])

    ax3=fig.add_subplot(2,3,3)
    w0=ttp-100; w1=tts+200
    mm=(t>=w0)&(t<=w1)
    # 彩色方块: T-P 窗口蓝, T-S 窗口紫
    ax3.axvspan(w_tp[0], w_tp[1], color='#1f77b4', alpha=0.15)
    ax3.axvspan(w_ts[0], w_ts[1], color='#9467bd', alpha=0.15)
    ax3.plot(t[mm], env_z[mm]/max(env_z[mm]),'r-',lw=1.2,alpha=0.8)
    ax3.plot(t[mm], zf[mm]/max(abs(zf[mm])),'k-',lw=0.3,alpha=0.5)
    for tt_v, clr, lbl in [(ttp,'#1f77b4','T-P th'),(t_obs_tp,'#1f77b4','T-P obs'),
                           (tts,'#9467bd','T-S th'),(t_obs_ts,'#9467bd','T-S obs')]:
        ax3.axvline(tt_v,c=clr,ls='--' if 'th' in lbl else '-',lw=1.2,alpha=0.7)
        ax3.text(tt_v, 1.03, f'{lbl}={tt_v:.0f}s', c=clr, fontsize=8, fontweight='bold',
                 ha='center', va='bottom', rotation=0)
    ax3.set_xlim(w0,w1); ax3.set_ylim(-1.2,1.25)
    ax3.set_title('Z Envelope + picks (gray: BP waveform)',fontsize=12,fontweight='bold')
    ax3.set_xlabel('Time from origin (s)',fontsize=10); ax3.grid(alpha=0.08)

    # Row2: Z/R/T 三分量 T波窗口 (T-P-100 ~ T-S+200), 灰波形+红包络8s+彩色方块
    comps=[('Z', zf, env_z), ('R', rf, env_r), ('T', tf, env_t)]
    for k,(cname, cwf, cenv) in enumerate(comps):
        ax=fig.add_subplot(2,3,4+k)
        mm=(t>=w0)&(t<=w1)
        p98=np.percentile(np.abs(cwf[mm]),98) or 1
        ax.plot(t[mm], cwf[mm]/p98, color='0.55', lw=0.7)          # 波形灰色
        ax.plot(t[mm], cenv[mm]/p98, 'r-', lw=1.2, alpha=0.8)      # 包络红色
        ax.fill_between(t[mm], -cenv[mm]/p98, cenv[mm]/p98, color='red', alpha=0.08)
        ax.axvspan(w_tp[0], w_tp[1], color='#1f77b4', alpha=0.15)
        ax.axvspan(w_ts[0], w_ts[1], color='#9467bd', alpha=0.15)
        for tt_v, clr in [(ttp,'#1f77b4'),(tts,'#9467bd')]:
            ax.axvline(tt_v,c=clr,ls='--',lw=0.9,alpha=0.6)
        ax.set_xlim(w0,w1)
        ax.set_title(f'{cname} comp  (BP 2-8Hz, env smooth {SMOOTH_S:.0f}s)',fontsize=12,fontweight='bold')
        ax.set_xlabel('Time from origin (s)',fontsize=10)
        ax.set_ylabel('Norm. Amp.',fontsize=10)
        ax.grid(alpha=0.08)

    plt.tight_layout(rect=[0,0,1,0.96])
    fn=f'{OUT}/{sta}/{sta}_envelope_polar.png'
    fig.savefig(fn,dpi=200,bbox_inches='tight'); plt.close(fig)
    print(f'  => {fn}')

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("用法: python3 envelope_polar.py <STA> <NET>")
        sys.exit(1)
    process_sta(sys.argv[1], sys.argv[2])
