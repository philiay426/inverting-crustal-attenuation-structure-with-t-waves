#!/usr/bin/env python3
"""
Love 波验证分析模板 — 三个独立方法
(1) R-T 偏振方向 vs 横向预期（最直接）
(2) 内陆传播速度估计（沿方位角剖面）
(3) H/V 谱比

用法：复制到事件目录，修改 EVLA/EVLO/EVDP/BASE/RESP_DIR，运行。
"""
import os, sys, glob
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from obspy import read
from obspy.signal.filter import bandpass as bp_filter
from obspy.taup import TauPyModel
from obspy.geodetics.base import locations2degrees
from geographiclib.geodesic import Geodesic

sys.path.insert(0, '/home/sy/shared/codex/t_phase_analysis')
from t_phase_analysis.traveltime_profile import get_track, tt_t

# ═══ 配置 ═══
EVLA, EVLO, EVDP = 19.0428, -66.8096, 20.0
BASE = "/mnt/f/sharedate/case/2014-01-13-mw63-process"
RESP_DIR = os.path.join(BASE, "sacprocess", "remove_response")
TMP = os.path.join(BASE, "tmp")
os.makedirs(TMP, exist_ok=True)

# ═══ Flinn 协方差 ═══
def flinn_polarization(t, z, r, tc, ws, we, sub_win=5.0, overlap=0.75):
    step = sub_win * (1 - overlap)
    n_sub = max(1, int((we - ws - sub_win) / step) + 1)
    if n_sub == 1: sub_win = we - ws
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
    return {'angles_zr': np.array(angles_zr), 'weights_zr': np.array(weights_zr),
            'angles_rt': np.array(angles_rt), 'weights_rt': np.array(weights_rt)}

def weighted_mean_dir(angles, weights):
    if len(angles) == 0: return np.nan
    xs = np.sum(weights * np.cos(angles)); ys = np.sum(weights * np.sin(angles))
    return np.arctan2(ys, xs)

# ═══ 方法1: R-T 偏振方向 vs 横向 ═══
def verify_rt_direction():
    """Love波: R-T偏振方向应接近±90° (纯横向)"""
    df = pd.read_excel(os.path.join(BASE, "stations.xlsx"), header=2)
    # ... (full implementation in tmp/verify_love_wave.py)
    pass

# ═══ 方法2: 内陆传播速度 ═══
def verify_propagation_velocity():
    """Love波: V < 4.5 km/s; P波: V > 5 km/s"""
    pass

# ═══ 方法3: H/V 谱比 ═══
def verify_hv_ratio():
    """Love波: H/V >> 1; Rayleigh: H/V ≈ 1"""
    pass
