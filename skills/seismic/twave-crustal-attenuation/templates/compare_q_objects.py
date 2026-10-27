#!/usr/bin/env python3
# =============================================================================
# Q_S 的倒数（Q_S^-1）随频率 —— 多个对象对比（本文结果 / 对照事件 / 文献值）
# -----------------------------------------------------------------------------
# 每条曲线都是 Q_S = Q0 * f^a 的等价写法  Q_S^-1 = (1/Q0) * f^(-a)，只靠颜色 + 线型区分。
# ★ 曲线参数由使用者填进 CASES（Q0、a、sigma_a 从 Q 反演输出或文献抄）；
#   本模板**不含任何具体事件或文献的数值**。
# 误差：a 的误差按论文口径 = 95% 置信区间 / 1.96
# 输出：Q_OUTPNG（默认 ../Q倒数_对比.png），无标题，黑点黑线为主
# 用法：Q_OUTPNG=<out.png> python compare_q_objects.py
# =============================================================================
import math, os
import numpy as np
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
OUTPNG = os.environ.get('Q_OUTPNG') or os.path.join(HERE, os.pardir, 'Q倒数_对比.png')
DPI    = 300
FLO, FHI = 2.00, 4.00        # 横轴范围 (Hz)

# (名称, Q0, a, sigma_a, 颜色, 线型, 拟合频段右端)
CASES = [
    # (图例名称, Q0, a, sigma_a, 颜色, 线型, 拟合频段右端 Hz)
    # ★ 按自己的数据/文献填；示例占位值不可直接当结论用
    ('This study',        1000.0, 0.50, 0.20, 'k',       '-', 4.00),
    ('Reference (20XX)',   800.0, 0.60, 0.15, '#d62728', '-', 4.00),
]

# ============================== 画图
fig, ax = plt.subplots(figsize=(5.0, 3.8))
fl = np.linspace(FLO, FHI, 200)

for name, Q0, a, sa, col, ls, fend in CASES:
    m = fl <= fend + 1e-9
    ax.plot(fl[m], (1.0 / Q0) * fl[m] ** (-a), ls, color=col, lw=1.7, zorder=2, label=name)

ax.set_xscale('log'); ax.set_yscale('log')
_tk = [2.0, 2.5, 3.0, 4.0]
ax.xaxis.set_major_locator(FixedLocator(_tk))
ax.xaxis.set_major_formatter(FixedFormatter(['%g' % t for t in _tk]))
ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_xlim(FLO - 0.08, FHI + 0.10)

# 纵轴范围：数据上下各留 ~45%/55% 余量
_yv = np.concatenate([(1.0 / Q0) * np.array([FLO, fend]) ** (-a) for _, Q0, a, _, _, _, fend in CASES])
_lo, _hi = _yv.min() / 1.45, _yv.max() * 1.55
ax.set_ylim(_lo, _hi)
_yt = [v for v in (0.4e-3, 0.6e-3, 0.8e-3, 1.0e-3, 1.2e-3, 1.4e-3, 1.6e-3, 1.8e-3, 2.0e-3) if _lo < v < _hi]
ax.yaxis.set_major_locator(FixedLocator(_yt))
ax.yaxis.set_major_formatter(FixedFormatter([r'$%.1f{\times}10^{-3}$' % (v / 1e-3) for v in _yt]))
ax.yaxis.set_minor_formatter(NullFormatter())

ax.set_xlabel('log frequency (Hz)', fontsize=10.5)
ax.set_ylabel(r'$Q_S^{-1}$', fontsize=11.5)
ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=9)
ax.grid(alpha=0.22, lw=0.4, which='major')
ax.legend(loc='lower left', fontsize=8, frameon=False, handlelength=2.2, borderaxespad=0.5)

# ---- 逐条线标注 a 值：自动挑不压其它线、不出框的位置 ----
fig.canvas.draw()
_rd = fig.canvas.get_renderer()
_axbb = ax.get_window_extent()
# 三条线的显示坐标（用于间距检查）
_lines = []
for name, Q0, a, sa, col, ls, fend in CASES:
    ff = np.linspace(FLO, fend, 240)
    _lines.append(ax.transData.transform(np.column_stack([ff, (1.0 / Q0) * ff ** (-a)])))

# 三个标注分别放在不同横位置，避免互压（按曲线由低到高分配）
_XTARGET = {'Mw 6.2 (this study)': 2.35, 'Zhou et al. (2021)': 3.15, 'Mw 5.8': 3.60}
_placed = []
for k, (name, Q0, a, sa, col, ls, fend) in enumerate(CASES):
    _p1 = ax.transData.transform((FLO + 0.03, (1.0 / Q0) * (FLO + 0.03) ** (-a)))
    _p2 = ax.transData.transform((fend - 0.03, (1.0 / Q0) * (fend - 0.03) ** (-a)))
    _ang = math.degrees(math.atan2(_p2[1] - _p1[1], _p2[0] - _p1[0]))
    _lab = 'a = %.2f \u00b1 %.2f' % (a, sa)
    _tries = []
    _xc = _XTARGET.get(name, 0.5 * (FLO + fend))
    for _fx in np.arange(_xc - 0.22, _xc + 0.23, 0.04):
        _fx = min(max(_fx, FLO + 0.08), fend - 0.06)
        _yl = (1.0 / Q0) * _fx ** (-a)
        for _dy in (8, -8, 12, -12, 16, -16, 20, -20, 26, -26, 32, -32, 40, -40):
            _t = ax.annotate(_lab, xy=(_fx, _yl), xycoords='data',
                             xytext=(0, _dy), textcoords='offset points',
                             ha='center', va='center', fontsize=8.5, color=col,
                             rotation=_ang, rotation_mode='anchor')
            _bb = _t.get_window_extent(renderer=_rd)
            _inside = (_bb.x0 > _axbb.x0 + 4 and _bb.x1 < _axbb.x1 - 4 and
                       _bb.y0 > _axbb.y0 + 4 and _bb.y1 < _axbb.y1 - 4)
            if not _inside:
                _t.remove(); continue
            # 与其它两条线的最小间距（pt）
            _d = 1e9
            for j, _pl in enumerate(_lines):
                if j == k:
                    continue
                for _px, _py in _pl:
                    if _bb.x0 - 40 < _px < _bb.x1 + 40 and _bb.y0 - 40 < _py < _bb.y1 + 40:
                        _dx = max(_bb.x0 - _px, 0.0, _px - _bb.x1)
                        _dyy = max(_bb.y0 - _py, 0.0, _py - _bb.y1)
                        _d = min(_d, math.hypot(_dx, _dyy))
            # 与已放好的标注的最小间距（pt）
            for _ob in _placed:
                _dx = max(_bb.x0 - _ob.x1, 0.0, _ob.x0 - _bb.x1)
                _dyy = max(_bb.y0 - _ob.y1, 0.0, _ob.y0 - _bb.y1)
                _d = min(_d, math.hypot(_dx, _dyy))
            _tries.append((_d, abs(_dy) + 3 * abs(_fx - _xc), _fx, _dy, _t, _bb))
    _ok = [c for c in _tries if c[0] >= 7.0]
    _pick = (min(_ok, key=lambda c: c[1]) if _ok else max(_tries, key=lambda c: c[0]))
    for c in _tries:
        if c[4] is not _pick[4]:
            c[4].remove()
    _placed.append(_pick[5])
    print('%s: f=%.2f Hz, 偏移 %+d pt, 最小间距 %.1f pt' % (name, _pick[2], _pick[3], _pick[0]))

fig.savefig(OUTPNG, dpi=DPI, bbox_inches='tight')
print('图 -> %s' % os.path.normpath(OUTPNG))
