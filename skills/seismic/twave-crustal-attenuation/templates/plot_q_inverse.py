#!/usr/bin/env python3
# =============================================================================
# Q^-1 随频率图（定案版式）：1/Q_S = (1/Q0) * f^(-alpha)，对数-对数，黑点黑线
# -----------------------------------------------------------------------------
# 输入 = twave_q_inversion.py 产出的 q_<tag>_detail.csv（逐台逐频点振幅表）
# 输出 = 单张 PNG（无标题、无误差棒、拟合线旁标 a = alpha +- sigma）
# 版式要点（2026-10-20 用户定稿）：
#   · 纵轴 Q^-1，统一科学计数法刻度（0.2~3.0 x10^-3 取整数倍，不出现裸 10^n）
#   · 横轴 log，主刻度 2.0 / 2.5 / 3.0 / 4.0，标签写 "log frequency (Hz)"
#   · 数据点黑色实心圆 ms 5.5；拟合线黑色实线 lw 1.7
#   · a 的标注沿拟合线旋转（角度按显示坐标算），自动避开数据点、不出框
# 环境变量：
#   Q_DETAIL=<detail.csv>   Q_OUTPNG=<png 路径>
#   Q_CODA=2S               尾波口径列名后缀（A_coda_<口径>）
#   Q_FREQ_EXCL=2.25,2.75   剔除频点（须有 |k|/sigma_k 判据支撑）
#   Q_FLO=2.0 Q_FHI=4.0     参与拟合的频段
#   Q_V=3.6 Q_HM=80         陆段速度 / 几何扩散特征长度
# 用法：
#   Q_DETAIL=.../q_s13_detail.csv Q_FREQ_EXCL=2.25,2.75 \
#     /home/sy/miniforge3/envs/obspy/bin/python plot_q_inverse.py
# =============================================================================
import csv, math, os
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9,
    'axes.labelsize': 10.5,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'axes.linewidth': 0.6,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

# ============================== ★ 参数区
HERE   = os.path.dirname(os.path.abspath(__file__))
DET    = os.environ.get('Q_DETAIL') or os.path.join(HERE, 'q_s13_detail.csv')
OUTPNG = os.environ.get('Q_OUTPNG') or os.path.join(HERE, os.pardir, 'Q倒数_拟合线.png')
CODA   = os.environ.get('Q_CODA', '2S')      # ★ 尾波口径（列名后缀 A_coda_<CODA>）
CODA_OF_STA = {}                             # 逐台例外口径（只在有判据时用）
FREQ_EXCL = tuple(float(v) for v in os.environ.get('Q_FREQ_EXCL', '').split(',') if v.strip())
FLO    = float(os.environ.get('Q_FLO', 2.0))
FHI    = float(os.environ.get('Q_FHI', 4.0))
V      = float(os.environ.get('Q_V', 3.6))   # 陆段视速度 (km/s)
H_M    = float(os.environ.get('Q_HM', 80.0)) # 几何扩散 h_m (km)；常数因子，不影响斜率
SNR_MIN = float(os.environ.get('Q_SNR_MIN', 1.0))   # 逐台入选门槛 A_d/A_c >= 1
COL     = 'k'                                # 数据点/拟合线颜色
DPI     = 300

# ============================== 读表
rows = list(csv.DictReader(open(DET)))
freqs = sorted({float(r['f_hz']) for r in rows})
freqs = [f for f in freqs if FLO - 1e-9 <= f <= FHI + 1e-9
         and not any(abs(f - fe) < 1e-9 for fe in FREQ_EXCL)]

fs, Qs, sQs = [], [], []
for f in freqs:
    x, y = [], []
    for r in rows:
        if float(r['f_hz']) != f:
            continue
        ac = float(r['A_coda_' + CODA_OF_STA.get(r['sta'], CODA)]); ad = float(r['A_theory'])
        r_km = float(r['r_km']); land = float(r['land_km'])
        if ac <= 0 or ad <= 0 or ad / ac < SNR_MIN:
            continue
        # 尾波归一化后的振幅比，乘回几何扩散 sqrt(r*h_m)
        x.append(land); y.append(math.log(ad * math.sqrt(r_km * H_M) / ac))
    if len(x) < 4:
        continue
    x = np.asarray(x); y = np.asarray(y); n = len(x)
    p = np.polyfit(x, y, 1)                       # y = k*x + b
    res = y - np.polyval(p, x)
    s2 = float(res @ res) / (n - 2)
    Sxx = float(((x - x.mean()) ** 2).sum())
    sk = math.sqrt(s2 / Sxx)                      # sigma_k
    k = abs(p[0])
    Q = math.pi * f / (V * k)                     # Q = pi f / (v |k|)
    fs.append(f); Qs.append(Q); sQs.append(Q * sk / k)
    print('f=%.2f Hz  n=%d  |k|=%.3e +- %.1e /km  Q=%.0f +- %.0f' % (f, n, k, sk, Q, Q * sk / k))

fs = np.asarray(fs); Qs = np.asarray(Qs); sQs = np.asarray(sQs)
inv = 1.0 / Qs; sinv = sQs / Qs ** 2                # sigma(1/Q) = sigma_Q / Q^2

# 在 1/Q 上做不加权 log-log 回归：ln(1/Q) = ln(1/Q0) - alpha*ln f
lx = np.log(fs); ly = np.log(inv)
p = np.polyfit(lx, ly, 1)                           # 斜率 = -alpha
alpha = -p[0]; Q0 = math.exp(-p[1])
res = ly - np.polyval(p, lx)
s2 = float(res @ res) / (len(fs) - 2)
Sxx = float(((lx - lx.mean()) ** 2).sum())
sea = math.sqrt(s2 / Sxx)
# a 的误差按论文/导师口径：95% 置信区间的一半再除以 1.96
sig_a = sea * stats.t.ppf(0.975, len(fs) - 2) / 1.96
# Q0 的 1 sigma：截距标准误按同一口径（t(0.975,n-2)*se_intercept/1.96）再指数化
n_f = len(fs)
s2_i = float(res @ res) / (n_f - 2)
se_b = math.sqrt(s2_i * (1.0 / n_f + lx.mean() ** 2 / Sxx))
sig_Q0 = Q0 * se_b * stats.t.ppf(0.975, n_f - 2) / 1.96
print('\nQ^-1 = (1/%.0f) f^(%+.2f)   即 Q = (%.0f +- %.0f) f^(%+.2f +- %.2f)   (n=%d)'
      % (Q0, -alpha, Q0, sig_Q0, alpha, sig_a, n_f))

# ============================== 画图
fig, ax = plt.subplots(figsize=(5.0, 3.8))
ax.plot(fs, inv, 'o', ms=5.5, mfc='k', mec='k', mew=0.8, zorder=3)   # 黑色实心圆点
_fl = np.linspace(FLO - 0.05, FHI + 0.05, 60)
ax.plot(_fl, (1.0 / Q0) * _fl ** (-alpha), '-', lw=1.7, color=COL, zorder=2)

ax.set_xscale('log'); ax.set_yscale('log')
_tk = [2.0, 2.5, 3.0, 4.0]
ax.xaxis.set_major_locator(FixedLocator(_tk))
ax.xaxis.set_major_formatter(FixedFormatter(['%g' % t for t in _tk]))
ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_xlim(FLO - 0.08, FHI + 0.10)

# 纵轴：统一科学计数法（带尾数 m x 10^n，不出现裸 10^n）
_ymin = min(inv.min(), (1.0 / Q0) * FLO ** (-alpha))
_ymax = max(inv.max(), (1.0 / Q0) * FHI ** (-alpha))
_lo, _hi = _ymin / 1.45, _ymax * 1.55          # ★ 纵轴留足余量（2026-10-20 用户要求画大一些）
ax.set_ylim(_lo, _hi)
_yt = [v for v in np.arange(0.2e-3, 3.01e-3, 0.1e-3) if _lo < v < _hi]
ax.yaxis.set_major_locator(FixedLocator(_yt))
ax.yaxis.set_major_formatter(FixedFormatter([r'$%.1f{\times}10^{-3}$' % (v / 1e-3) for v in _yt]))
ax.yaxis.set_minor_formatter(NullFormatter())

ax.set_xlabel('log frequency (Hz)', fontsize=10.5)
ax.set_ylabel(r'$Q_S^{-1}$', fontsize=11.5)
ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=9)
ax.grid(alpha=0.22, lw=0.4, which='major')

# ---- 沿拟合线标注 a 值：自动挑一个不压数据点、不出框的位置 ----
# 角度按显示坐标算（对数轴下斜率会变，必须用 transform 后的坐标）
fig.canvas.draw()
_p1 = ax.transData.transform((FLO + 0.03, (1.0 / Q0) * (FLO + 0.03) ** (-alpha)))
_p2 = ax.transData.transform((FHI - 0.03, (1.0 / Q0) * (FHI - 0.03) ** (-alpha)))
_ang = math.degrees(math.atan2(_p2[1] - _p1[1], _p2[0] - _p1[0]))
_pts = ax.transData.transform(np.column_stack([fs, inv]))      # 数据点的显示坐标
_axbb = ax.get_window_extent()
_rd = fig.canvas.get_renderer()
_lab = 'a = %.2f \u00b1 %.2f' % (alpha, sig_a)
_found, _tries = [], []
for _fx in np.arange(FLO + 0.12, FHI - 0.06, 0.06):
    _yl = (1.0 / Q0) * _fx ** (-alpha)
    for _dy in (10, -10, 14, -14, 18, -18, 24, -24, 30, -30, 38, -38, 48, -48):
        _t = ax.annotate(_lab, xy=(_fx, _yl), xycoords='data',
                         xytext=(0, _dy), textcoords='offset points',
                         ha='center', va='center', fontsize=8.5, color='k',
                         rotation=_ang, rotation_mode='anchor')
        _bb = _t.get_window_extent(renderer=_rd)
        _inside = (_bb.x0 > _axbb.x0 + 10 and _bb.x1 < _axbb.x1 - 10 and
                   _bb.y0 > _axbb.y0 + 10 and _bb.y1 < _axbb.y1 - 10)
        if not _inside:
            _t.remove(); continue
        _dmin = 1e9
        for _px, _py in _pts:
            _dx = max(_bb.x0 - _px, 0.0, _px - _bb.x1)
            _dyy = max(_bb.y0 - _py, 0.0, _py - _bb.y1)
            _dmin = min(_dmin, math.hypot(_dx, _dyy))
        _tries.append((_dmin, abs(_dy), _fx, _dy, _t))
# 先满足"离数据点 >= 9 pt"，再取偏移最小（即最贴线）的那个
_ok = [c for c in _tries if c[0] >= 9.0]
_pick = (min(_ok, key=lambda c: c[1]) if _ok else max(_tries, key=lambda c: c[0]))
for c in _tries:
    if c[4] is not _pick[4]:
        c[4].remove()
print('标注位置: f=%.2f Hz, 竖直偏移 %+d pt, 与最近数据点间距 %.1f pt' % (_pick[2], _pick[3], _pick[0]))

fig.savefig(OUTPNG, dpi=DPI, bbox_inches='tight')
print('图 -> %s' % os.path.normpath(OUTPNG))
