#!/usr/bin/env python3
"""
核心处理工具库 — 基于 BOSA test 验证方案 (2026-08-01 更新)

核心变更:
  - 去响应: tr.simulate(paz_remove=paz_dict, remove_sensitivity=True) — 直接出速度 m/s
  - 带通滤波: scipy butter + filtfilt (零相位)
  - SACPZ 解析: 统一 parse_sacpz() 返回 (zeros, poles, a0, sensitivity)

旧方案 (已废弃): simulate_seismometer() + np.gradient() 手动转速度
"""
import os, numpy as np
from obspy import read
from scipy.signal import butter, filtfilt

TAPER = 0.05
PRE_FILT = (0.001, 0.005, 45.0, 50.0)
WATER_LEVEL = 60.0  # dB
BP_FMIN, BP_FMAX = 2.0, 8.0
BP_CORNERS = 4


def parse_sacpz(filepath):
    """解析SACPZ, 返回 (zeros, poles, a0, sensitivity)"""
    zeros, poles = [], []
    constant = 1.0
    a0 = None
    sensitivity = None
    found_first = False
    in_zeros = False
    in_poles = False
    done = False

    with open(filepath) as f:
        for line in f:
            if done:
                break
            s = line.strip()

            # 注释行: 提取 A0, SENSITIVITY
            if s.startswith('*'):
                if 'SENSITIVITY' in s.upper():
                    parts = s.split()
                    for i, p in enumerate(parts):
                        if 'SENSITIVITY' in p.upper() and i + 2 < len(parts):
                            try:
                                sensitivity = float(parts[i + 2])
                            except ValueError:
                                pass
                if s.upper().lstrip('* ').startswith('A0'):
                    parts = s.split()
                    for i, p in enumerate(parts):
                        if p.upper() == 'A0' and i + 2 < len(parts):
                            try:
                                a0 = float(parts[i + 2])
                            except ValueError:
                                pass
                continue

            if s.startswith('ZEROS'):
                if not found_first:
                    found_first = True
                    in_zeros = True
                    in_poles = False
                else:
                    done = True
                continue
            if s.startswith('POLES') and found_first:
                in_zeros = False
                in_poles = True
                continue
            if s.startswith('CONSTANT') and found_first:
                parts = s.split()
                if len(parts) >= 2:
                    constant = float(parts[1])
                continue

            parts = s.split()
            if len(parts) >= 2 and found_first and not done:
                try:
                    c = complex(float(parts[0]), float(parts[1]))
                    if in_zeros:
                        zeros.append(c)
                    elif in_poles:
                        poles.append(c)
                except ValueError:
                    continue

    if a0 is None and sensitivity is not None:
        a0 = constant / sensitivity
    elif sensitivity is None and a0 is not None:
        sensitivity = constant / a0
    elif a0 is None and sensitivity is None:
        a0 = constant
        sensitivity = 1.0

    return zeros, poles, a0, sensitivity


def remove_response(tr, paz_file):
    """
    去仪器响应，直接输出速度 (m/s).
    使用 tr.simulate(paz_remove=...) — BOSA 验证方案.
    """
    zeros, poles, a0, sensitivity = parse_sacpz(paz_file)
    paz = {'zeros': zeros, 'poles': poles, 'gain': a0, 'sensitivity': sensitivity}
    tr_out = tr.copy()
    tr_out.data = tr_out.data.astype(np.float64)
    tr_out.detrend('demean')
    tr_out.detrend('linear')
    tr_out.taper(max_percentage=TAPER)
    tr_out.simulate(paz_remove=paz, remove_sensitivity=True,
                    water_level=WATER_LEVEL, pre_filt=PRE_FILT)
    return tr_out


def apply_bandpass(tr, fmin=BP_FMIN, fmax=BP_FMAX, corners=BP_CORNERS, fs=None):
    """
    零相位带通滤波 (scipy butter + filtfilt).
    fs 可从 tr.stats.sampling_rate 自动获取.
    """
    if fs is None:
        fs = tr.stats.sampling_rate
    b, a = butter(corners, [fmin, fmax], btype='band', fs=fs)
    data = tr.data.astype(np.float64)
    filtered = filtfilt(b, a, data)
    tr_out = tr.copy()
    tr_out.data = filtered
    return tr_out


def find_sacpz(sta, ch, sacpz_dir):
    """
    查找 SACPZ 文件，依次尝试不同 location code.
    返回文件路径，或 None.
    """
    for pat in [f"SACPZ.*.{sta}..{ch}",
                f"SACPZ..{sta}..{ch}",
                f"SACPZ.*.{sta}.*.{ch}",
                f"SACPZ..{sta}.*.{ch}"]:
        import glob
        m = glob.glob(os.path.join(sacpz_dir, pat))
        if m:
            return m[0]
    return None


def process_single_file(raw_path, sacpz_dir, out_dir,
                        do_bandpass=True, save_resp=True):
    """
    处理单个 SAC 文件: detrend→taper→去响应→(可选带通滤波)→保存.

    参数:
      raw_path: 原始 SAC 路径
      sacpz_dir: SACPZ 目录
      out_dir: 输出目录
      do_bandpass: 是否做带通滤波
      save_resp: 是否保存去响应后的 SAC

    返回: (resp_trace, proc_trace) — proc 可能是 resp 的副本(若 do_bandpass=False)
    """
    from obspy.io.sac import SACTrace

    tr_raw = read(raw_path)[0]

    # 从 SAC 头获取台站信息
    sac = tr_raw.stats.sac
    net = sac.knetwk.strip() if sac.knetwk else ''
    sta = sac.kstnm.strip() if sac.kstnm else ''
    ch = sac.kcmpnm.strip() if sac.kcmpnm else ''

    paz_file = find_sacpz(sta, ch, sacpz_dir)
    if paz_file is None:
        # fallback: 从文件名解析
        fname = os.path.basename(raw_path).replace('.SAC', '')
        parts = fname.split('.')
        net2, sta2, ch2 = parts[0], parts[1], parts[3]
        paz_file = find_sacpz(sta2, ch2, sacpz_dir)
    if paz_file is None:
        print(f"  ❌ SACPZ 未找到: {os.path.basename(raw_path)}")
        return None, None

    # 预处理
    tr_resp = remove_response(tr_raw, paz_file)

    base = os.path.basename(raw_path).replace('.SAC', '').replace('.sac', '')

    # 保存去响应
    if save_resp:
        resp_dir = os.path.join(out_dir, 'remove_response')
        os.makedirs(resp_dir, exist_ok=True)
        resp_out = os.path.join(resp_dir, f"{base}_resp.sac")
        sac_out = SACTrace.from_obspy_trace(tr_resp)
        sac_out.write(resp_out)

    # 带通滤波
    if do_bandpass:
        tr_proc = apply_bandpass(tr_resp)
        proc_out = os.path.join(out_dir, f"{base}_proc.sac")
        sac_p = SACTrace.from_obspy_trace(tr_proc)
        sac_p.write(proc_out)
        return tr_resp, tr_proc

    return tr_resp, tr_resp
