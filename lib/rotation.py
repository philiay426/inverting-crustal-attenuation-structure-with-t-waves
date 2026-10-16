#!/usr/bin/env python3
"""ZNE→ZRT 分量旋转（统一 obspy 约定，2026-09-02）"""
import numpy as np


def rotate_zne_to_zrt(z, n, e, baz):
    """ZNE → ZRT。R 背离震源（obspy rotate_ne_rt 约定，r = -n·cos - e·sin）。

    2026-09-02: 原手写公式 r=n·cos+e·sin（R 指向震源）与 obspy 差 180° 反相，
    已统一到 obspy 标准（调用 rotate_ne_rt），保证与单台偏振分析一致。
    """
    from obspy.signal.rotate import rotate_ne_rt
    r, t = rotate_ne_rt(np.asarray(n, dtype=float), np.asarray(e, dtype=float), float(baz))
    return np.asarray(z, dtype=float), r, t


def energy_ratios(zw, rw, tw):
    """窗口内未归一化能量比"""
    ez = np.sum(zw**2)
    er = np.sum(rw**2)
    et = np.sum(tw**2)
    eh = er + et
    er_et = er / et if et > 0 else 0
    ez_eh = ez / eh if eh > 0 else 0
    return ez, er, et, eh, er_et, ez_eh


def amplitude_rose(theta, weights, nbins=36):
    """全圆振幅加权玫瑰图直方图"""
    bins = np.linspace(-np.pi, np.pi, nbins+1)
    centers = (bins[:-1] + bins[1:]) / 2
    hist = np.zeros(nbins)
    for i in range(nbins):
        m = (theta >= bins[i]) & (theta < bins[i+1])
        hist[i] = np.sum(weights[m])
    if np.max(hist) > 0: hist /= np.max(hist)
    return centers, hist
