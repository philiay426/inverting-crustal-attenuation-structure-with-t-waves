#!/usr/bin/env python3
"""
偏振方法共享模块 (cov / amp / count) — 方向 + 权重
=================================================
三种口径 (与 SKILL.md 三方法表一致):

  method | 每个样本的方向                                        | 权重
  -------|------------------------------------------------------|--------------------
  cov    | 子窗 3x3 协方差最大特征向量的水平(R-T)主轴方位         | lin*sqrt(lambda1)
  amp    | 逐采样点 atan2(E, N)                                  | N^2 + E^2
  count  | 逐采样点 atan2(E, N)                                  | 1 (每点 1 票) ← 默认

角度约定: 地理 N/E, 0 deg = N、顺时针 = atan2(E, N)。
统一返回到 **R-T 坐标系** 时减去 BAZ: ang_rt = wrap180(atan2(E, N) - BAZ),
R-T 系里 0 deg = 径向 R(P/SV), +-90 deg = 横向 T(SH/Love)。
注意 arctan2(T, R) 恒等于 atan2(E, N) - BAZ(旋转保角), 所以只拿到 R/T 两道的
脚本可直接用 angles_rt_from_rt(), 不必再取 BHN/BHE。

默认方法: count(采样点等权 / 每点 1 票)。优先级: 命令行 > 环境变量 POLAR_METHOD > 默认。
为什么默认 count: 不看幅度、无偏, 背景噪声与波包同权; 用于"方向分布形状"比较时
是最保守的口径(见 SKILL.md 三方法对比小节)。amp 会被强尖峰主导, cov 假定子窗内单一偏振。

线性度/平面度定义(与 templates/polarization_attributes.py 及 t_phase_analysis.pca()
完全一致, lam1 >= lam2 >= lam3):
    linearity = clip((lam1 - lam2) / lam1, 0, 1)
    planarity = clip((lam2 - lam3) / lam1, 0, 1)

用法:
    from polar_methods import resolve_method, sample_direction_weights, method_label
    METHOD = resolve_method(args.method)                  # CLI > POLAR_METHOD > count
    ang, w = sample_direction_weights(METHOD, r, t_comp)  # ang 弧度, R-T 系
"""
import os

import numpy as np

METHODS = ('cov', 'amp', 'count')
METHOD_COV = 'cov'
DEFAULT_METHOD = 'count'
ENV_METHOD = 'POLAR_METHOD'
ENV_OUT_DIR = 'OUT_DIR'

SECTOR_HALF_DEG = 22.5   # 径向/横向扇区半宽 (SKILL.md 三方法口径)


# ─── 方法选择 / 输出目录 ───

def resolve_method(cli_method=None, default=DEFAULT_METHOD, env=ENV_METHOD):
    """方法解析: 命令行 > 环境变量 > 默认 count。非法值直接报错退出(不静默回退)。"""
    for src, val in (('--method', cli_method), (env, os.environ.get(env))):
        if val is None or str(val).strip() == '':
            continue
        m = str(val).strip().lower()
        if m not in METHODS:
            raise SystemExit(f'[polar_methods] unknown method {val!r} from {src}; '
                             f'choose one of {METHODS}')
        return m
    return default


def method_label(method, default=DEFAULT_METHOD):
    """图/日志用标签: 默认方法带 default 标记 -> 'count (default)'。"""
    return f'{method} (default)' if method == default else str(method)


def resolve_out_dir(cli_out_dir=None, fallback=None, env=ENV_OUT_DIR):
    """输出目录: --out-dir > OUT_DIR 环境变量 > 脚本内 OUT。"""
    return cli_out_dir or os.environ.get(env) or fallback


def add_common_args(parser):
    """给各脚本 argparse 加通用开关(--method / --out-dir)。"""
    parser.add_argument('--method', choices=list(METHODS), default=None,
                        help='polarization method: cov (covariance) | amp (amplitude-weighted) '
                             '| count (one vote per sample, DEFAULT); env POLAR_METHOD also works')
    parser.add_argument('--out-dir', default=None,
                        help='output directory override (--out-dir > env OUT_DIR > in-script OUT)')
    return parser


# ─── 方向 + 权重 ───

def wrap180(deg):
    """角度折叠到 (-180, 180]。"""
    return (np.asarray(deg, dtype=float) + 180.0) % 360.0 - 180.0


def angles_rt_from_rt(r, t_comp):
    """R/T 两道 → R-T 系角(弧度)。恒等于 atan2(E, N) - BAZ。"""
    return np.arctan2(np.asarray(t_comp, dtype=float), np.asarray(r, dtype=float))


def angles_rt_from_ne(n, e, baz_deg):
    """地理 N/E → R-T 系角(弧度)。0 = 径向 R。"""
    ang = np.degrees(np.arctan2(np.asarray(e, dtype=float), np.asarray(n, dtype=float)))
    return np.radians(wrap180(ang - float(baz_deg)))


def sample_direction_weights(method, r, t_comp):
    """amp / count: 逐采样点方向 + 权重。返回 (angles_rad, weights)。

    amp   -> 权重 N^2 + E^2 (旋转保范数, 用 R^2 + T^2 等价)
    count -> 权重恒为 1 (每点 1 票)
    """
    if method == METHOD_COV:
        raise ValueError('sample_direction_weights() is for amp/count; use '
                         'cov_subwindow_direction_weights() for cov')
    ang = angles_rt_from_rt(r, t_comp)
    if method == 'amp':
        a = np.asarray(r, dtype=float)
        b = np.asarray(t_comp, dtype=float)
        return ang, a * a + b * b
    return ang, np.ones_like(np.asarray(r, dtype=float))


def cov_subwindow_direction_weights(t, z, r, tc, ws, we, sub_win=0.5, overlap=0.75):
    """cov: Flinn 1965 滑动协方差, 子窗 3x3 协方差最大特征向量的 R-T 方位。

    权重 = lin * sqrt(lambda1)(逐位复用 2026-08 之前的 covariance_polarization(),
    数值与旧版完全一致, 不得改动本函数里的算式)。
    返回 (angles_rad, weights)。
    """
    win_len = we - ws
    if win_len < sub_win * 3:
        sub_win = max(0.25, win_len / 8)
    step = sub_win * (1 - overlap)
    n_sub = max(1, int((we - ws - sub_win) / step) + 1)
    if n_sub == 1:
        sub_win = we - ws
    angles_rt, weights_rt = [], []
    for i in range(n_sub):
        t0 = ws + i * step if n_sub > 1 else ws
        t1 = min(t0 + sub_win, we)
        m = (t >= t0) & (t <= t1)
        zw, rw, tw_c = z[m], r[m], tc[m]
        if len(zw) < 3:
            continue
        X = np.vstack([zw - np.mean(zw), rw - np.mean(rw), tw_c - np.mean(tw_c)])
        C = X @ X.T / X.shape[1]
        ev, evec = np.linalg.eigh(C)
        l1, l2, l3 = ev[2], ev[1], ev[0]
        v1 = evec[:, 2]
        lin = 1.0 - (l2 + l3) / (2.0 * l1) if l1 > 1e-30 else 0.0
        weight = lin * np.sqrt(l1)
        angles_rt.append(np.arctan2(v1[2], v1[1]))
        weights_rt.append(weight)
    return np.array(angles_rt), np.array(weights_rt)


def rose_hist(angles_rad, weights, nbins=36):
    """36 bin 玫瑰直方(归一化到峰值 1)。逐位复用旧 rose_from_covar()。"""
    bins = np.linspace(-np.pi, np.pi, nbins + 1)
    centers = (bins[:-1] + bins[1:]) / 2
    hist = np.zeros(nbins)
    for a, w in zip(angles_rad, weights):
        i = np.digitize(a, bins) - 1
        if i >= nbins:
            i = 0
        hist[i] += w
    if np.max(hist) > 0:
        hist /= np.max(hist)
    return centers, hist


def sector_fractions(angles_rad, weights, half_deg=SECTOR_HALF_DEG):
    """径向/横向扇区权重占比: 径向 = |ang| <= half, 横向 = ||ang| - 180| <= half。

    ⚠️ 只有占比可跨方法比较; 柱高/权重之和量纲各不相同(见 SKILL.md)。
    随机基线: 径向 12.5%、横向 25%。
    """
    ang = np.degrees(np.asarray(angles_rad, dtype=float))
    w = np.asarray(weights, dtype=float)
    tot = float(np.nansum(w))
    if not np.isfinite(tot) or tot <= 0:
        return {'frac_radial': float('nan'), 'frac_transverse': float('nan')}
    a = np.abs(wrap180(ang))
    m_r = a <= half_deg
    m_t = np.abs(a - 180.0) <= half_deg
    return {'frac_radial': float(np.nansum(w[m_r]) / tot),
            'frac_transverse': float(np.nansum(w[m_t]) / tot)}


def polarization_attributes(eigenvalues):
    """线性度/平面度(与 templates/polarization_attributes.py 定义一致)。

    eigenvalues: (..., 3), 升序或降序均可 —— 内部按 lam1 >= lam2 >= lam3 处理。
    返回 (linearity, planarity)。仅供自检/交叉核对; 主流程直接用 pca() 的同名列。
    """
    ev = np.sort(np.asarray(eigenvalues, dtype=float), axis=-1)
    lam3, lam2, lam1 = ev[..., 0], ev[..., 1], ev[..., 2]
    lam1_d = lam1 + np.finfo(float).eps
    lin = np.clip((lam1 - lam2) / lam1_d, 0.0, 1.0)
    pla = np.clip((lam2 - lam3) / lam1_d, 0.0, 1.0)
    return lin, pla
