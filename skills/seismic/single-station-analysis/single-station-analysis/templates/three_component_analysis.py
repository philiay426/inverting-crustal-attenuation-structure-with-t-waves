#!/usr/bin/env python3
"""
单台三分量频谱图 — 直接调用预处理脚本的 save_three_component_spectrogram
用法: 修改事件参数和台站名后运行
"""

import sys, os, glob

sys.path.insert(0, '/home/sy/codex/t_phase_analysis')

# ====== 事件参数 ======
EVLA, EVLO = -41.6794, 42.4486
EVDP = 10.0
STA = 'BOSA'
BASE = "/mnt/f/sharedate/case/Southern_Africa/2020-08-06-mww63-prince-edward-islands-region"
# ======================

# 加载预处理脚本中的函数
proc_obs = os.path.join(BASE, 'process_obs.py')
if not os.path.exists(proc_obs):
    print(f"process_obs.py not found at {proc_obs}"); sys.exit(1)

with open(proc_obs) as f:
    code = f.read()
code = code.split("if __name__ == \"__main__\":")[0]
ns = {}
exec(code, ns)

save_three_component_spectrogram = ns['save_three_component_spectrogram']

from obspy import read, UTCDateTime
from obspy.geodetics.base import locations2degrees, degrees2kilometers
from obspy.taup import TauPyModel
from t_phase_analysis.traveltime_profile import get_track, tt_t

RESP_DIR = os.path.join(BASE, "sacprocess", "remove_response")
OUT = os.path.join(BASE, "process_flow", "2_single_station_analysis", STA)
os.makedirs(OUT, exist_ok=True)

# 加载三分量
trs = {}
for ch in ['BHZ', 'BHN', 'BHE']:
    m = glob.glob(f"{RESP_DIR}/*{STA}*{ch}*_resp.sac")
    if m: trs[ch] = read(m[0])[0]
if len(trs) < 3:
    print(f"缺分量: {list(trs.keys())}"); sys.exit(1)

bhz = trs['BHZ']
stla = bhz.stats.sac.stla; stlo = bhz.stats.sac.stlo
origin = UTCDateTime(year=bhz.stats.sac.nzyear, julday=bhz.stats.sac.nzjday,
                     hour=bhz.stats.sac.nzhour, minute=bhz.stats.sac.nzmin,
                     second=bhz.stats.sac.nzsec,
                     microsecond=(getattr(bhz.stats.sac,'nzmsec',0) or 0)*1000)
baz = bhz.stats.sac.baz or 0

dist_deg = locations2degrees(EVLA, EVLO, stla, stlo)
dist_km = degrees2kilometers(dist_deg)

model = TauPyModel('iasp91')
arrivals = {}
for a in model.get_travel_times(source_depth_in_km=EVDP, distance_in_degree=dist_deg,
                                 phase_list=['P','S']):
    arrivals[a.name] = a.time

tt = get_track(EVLO, EVLA, stlo, stla, 0.02, ref_el=-1000, plot_map=False, plot_elevation=False)
tr_tt = tt_t(tt, ref_el=-1000, source_depth_km=EVDP, ignore_small_islands=True)
if tr_tt:
    arrivals['T-P'] = tr_tt['tt_t']
    arrivals['T-S'] = tr_tt['tt_ocean'] + tr_tt['tt_land_s']
    arrivals['land_km'] = tr_tt.get('dist_land', 0) * 111

# 生成频谱图
fp = os.path.join(OUT, f"{STA}_three_filter.png")
save_three_component_spectrogram(STA, trs, origin, arrivals, baz, fp, dist_km, dist_deg)
print(f"✓ {fp}")
