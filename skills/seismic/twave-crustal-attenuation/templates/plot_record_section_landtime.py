# -*- coding: utf-8 -*-
"""
T 分量剖面图（13 台，定案台组）：陆地计时 + 标出信号窗 / 尾波窗 / 理论到时
--------------------------------------------------------------------
横轴：陆地计时 t - L_ocean/1.48（0 = 该台 T 波上岸），海洋段传播时间已扣除，
      这样 T 波随陆地距离单调排列。
标注：
  · 信号窗（直达窗）= [T-S, T-S + 20 s]；陆地计时下 = [L_land/3.6, +20]
  · 尾波窗 = 全台同一绝对窗 [1200, 1220] s；ZAKA 例外用 T-S+400 = [1643, 1663] s
             （2026-10-20 定案）。绝对时间换算到陆地计时要逐台减 L_ocean/1.48，
             所以在图上是一条逐台错开的短段，不是竖带。
  · 理论到时：T-P 蓝短虚线 (陆地计时 = L_land/6.0)、T-S 紫短虚线 (= L_land/3.6)
数据：0_data_processing/envelope/<sta>_envelope_2-8Hz.npz 的 waveform_t
      （2-8 Hz 带通 T 分量，速度 m/s）；逐台按自身峰值归一化（只看形态，不比较振幅）。
版式：纵轴用台站序号等间距、刻度标真实陆距（台站两两陆距很接近时——例如同一子阵里
      几对台相距 <20 km——真实距离轴下波形必然压在一起，所以用等间距纵轴）。
"""
import os
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# ============================ ★ 参数区 ============================
HERE = os.path.dirname(os.path.abspath(__file__))
_b = os.environ.get('Q_BASE')
if not _b:
    raise SystemExit('请用环境变量 Q_BASE 指定事件的 process_flow 目录（其下有 1_station_distribution/stations.csv）')
BASE = os.path.abspath(_b)
CSV_PATH = os.path.join(BASE, '1_station_distribution', 'stations.csv')
# 包络目录：默认 BASE/0_data_processing/envelope；新布局/旧布局不同时用 Q_ENV_DIR 指过去
ENV_DIR = os.environ.get('Q_ENV_DIR') or os.path.join(BASE, '0_data_processing', 'envelope')

STAS = [x.strip() for x in os.environ.get('Q_STAS', '').split(',') if x.strip()]
if not STAS:
    raise SystemExit('请用环境变量 Q_STAS=A,B,C 指定要画的台站（顺序即纵轴顺序）')
EXPECTED_N = len(STAS)

V_OCEAN = 1.48           # 海洋段 T 波速度 (km/s)：扣海洋段、算陆地计时用
V_P     = 6.0            # 陆段 P 波速度 (km/s)：T-P 在陆地计时下 = L_land/V_P
V_S     = 3.6            # 陆段 S 波速度 (km/s)：T-S 在陆地计时下 = L_land/V_S

SIG_LEN = 20.0                                   # 信号窗长度 (s)
CODA_T0 = 1200.0                                 # 公共尾波窗起点（绝对时间, s）
CODA_LEN = 20.0                                  # 尾波窗长度 (s)
# ★ 尾波口径（定案 = '2S'）：'2S' | 'abs1200' | 'T-S+400'
CODA_MODE = os.environ.get('Q_CODA_MODE', '2S')
CODA_DT = float(os.environ.get('Q_CODA_DT', 400.0))
CODA_T0_OF_STA = {}        # 逐台例外，如 {'ZAKA': 1643.1}（只在有判据时用）

# ★ 两套时间轴各出一版：'land' = 陆地计时（扣海洋段），'abs' = 绝对计时（相对发震）
MODES = [('land', (-1400.0, 2400.0), (-100.0, 900.0)),   # (a)全时段, (b)放大
         ('abs',  (-600.0, 2200.0),  (800.0, 1800.0))]   # 绝对计时: 放大到直达窗~尾波窗

AMP_FRAC = 0.45                 # 单道显示幅度 = AMP_FRAC x 道间距
COL = 'k'                       # 波形
LW = 0.25
C_SIG = '#ff7f0e'               # 信号窗（橙）
C_CODA = '#d62728'              # 尾波窗（红）
C_TP = '#1f77b4'                # T-P（蓝）
C_TS = '#9467bd'                # T-S（紫）
BAR_HALF = 0.035                # 窗条半高（道间距的比例）
TICK_HALF = 0.33                # 到时短虚线半高

FIG_W, FIG_H = 13.0, 9.0
OUT = os.path.join(HERE, os.pardir, '台集S13_T分量剖面_陆地计时_带时窗.png')
# ==================================================================

plt.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
                     'font.size': 8, 'axes.labelsize': 9,
                     'xtick.labelsize': 8, 'ytick.labelsize': 8,
                     'axes.linewidth': 0.6, 'xtick.direction': 'in', 'ytick.direction': 'in',
                     'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                     'pdf.fonttype': 42, 'ps.fonttype': 42})


def load_rows():
    """读 stations.csv（按列名取；表头前有注释行，先定位以 sta, 开头那一行）"""
    raw = open(CSV_PATH).read().splitlines()
    i = [k for k, l in enumerate(raw) if l.startswith('sta,')][0]
    rows = [r for r in csv.DictReader(raw[i:]) if r.get('sta')]
    rows = [r for r in rows if r['sta'] in STAS]
    for r in rows:
        r['land'] = float(r['land_km'])
        r['ocean'] = float(r['ocean_km'])
        r['tt_s'] = float(r['tt_s'])
    rows.sort(key=lambda r: r['land'])              # ★ 按陆地距离升序
    return rows


def coda_t0(r):
    """尾波窗起点（绝对时间，相对发震）：'2S' = 2 x 该台 S 走时（lapse 从发震算）"""
    if CODA_MODE == '2S':
        return 2.0 * r['tt_s']
    if CODA_MODE.startswith('T-S+'):
        return r['ocean'] / V_OCEAN + r['land'] / V_S + CODA_DT
    return CODA_T0_OF_STA.get(r['sta'], CODA_T0)


def coda_label():
    if CODA_MODE == '2S':
        return 'Coda window [2 x S travel time, +%g s]' % CODA_LEN
    if CODA_MODE.startswith('T-S+'):
        return 'Coda window [T-S + %g s, +%g s]' % (CODA_DT, CODA_LEN)
    return 'Coda window [%g, %g] s' % (CODA_T0, CODA_T0 + CODA_LEN)


def draw(ax, rows, xlim, tag, mode, show_labels=True):
    n = len(rows)
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.9, n - 1 + 0.9)
    ax.set_yticks(np.arange(n))
    ax.set_yticklabels(['%s (%g km)' % (r['sta'], r['land']) for r in rows] if show_labels else [])
    ax.set_xlabel('Land time: t - ocean travel (s)   [0 = T-wave landfall]' if mode == 'land'
                  else 'Time from event origin (s)')
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.text(0.004, 1.004, '(%s)' % tag, transform=ax.transAxes, ha='left', va='bottom', fontsize=9)

    for k, r in enumerate(rows):
        f = os.path.join(ENV_DIR, '%s_envelope_2-8Hz.npz' % r['sta'])
        z = np.load(f)
        t = z['time'].astype(float)
        w = z['waveform_t'].astype(float)
        # ★ 陆地计时 = 绝对时间 - 海洋段传播时间；绝对计时模式不减
        x = t - (r['ocean'] / V_OCEAN if mode == 'land' else 0.0)
        m = (x >= xlim[0]) & (x <= xlim[1])
        if not np.any(m):
            continue
        xx, ww = x[m], w[m]
        peak = float(np.abs(ww).max())
        if peak > 0:
            ax.plot(xx, k + ww / peak * AMP_FRAC, color=COL, lw=LW, zorder=2)

        # ---- 理论到时：T-P（蓝）/ T-S（紫），陆地计时下分别 = L/V_P、L/V_S
        _off = r['ocean'] / V_OCEAN if mode == 'land' else 0.0   # 到时同样按时间轴换算
        for tv, cv, ls in ((r['ocean'] / V_OCEAN + r['land'] / V_P - _off, C_TP, (0, (4, 3))),
                           (r['ocean'] / V_OCEAN + r['land'] / V_S - _off, C_TS, (0, (1.5, 1.5)))):
            if xlim[0] <= tv <= xlim[1]:
                ax.plot([tv, tv], [k - TICK_HALF, k + TICK_HALF], color=cv, lw=0.9,
                        ls=ls, zorder=4)
        # ---- 信号窗 [T-S, T-S+20]（橙条，画在基线上）
        s0 = r['ocean'] / V_OCEAN + r['land'] / V_S - _off
        ax.plot([s0, s0 + SIG_LEN], [k, k], color=C_SIG, lw=4.2, solid_capstyle='butt', zorder=5)
        # ---- 尾波窗（红条）：绝对窗 [t0, t0+20] 换算到陆地计时
        c0 = coda_t0(r) - _off
        ax.plot([c0, c0 + CODA_LEN], [k, k], color=C_CODA, lw=4.2, solid_capstyle='butt', zorder=5)


rows = load_rows()
assert len(rows) == EXPECTED_N, '台站数 %d != %d' % (len(rows), EXPECTED_N)

handles = [Line2D([], [], color=COL, lw=0.9, label='T component (2-8 Hz)'),
           Line2D([], [], color=C_SIG, lw=4.2, label='Signal window [T-S, T-S+20 s]'),
           Line2D([], [], color=C_CODA, lw=4.2, label=coda_label()),
           Line2D([], [], color=C_TP, lw=0.9, ls=(0, (4, 3)), label='T-P theory'),
           Line2D([], [], color=C_TS, lw=0.9, ls=(0, (1.5, 1.5)), label='T-S theory')]

for mode, xf, xz in MODES:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(FIG_W, FIG_H), sharey=True,
                                   gridspec_kw=dict(height_ratios=[1.0, 1.35], hspace=0.16))
    draw(ax1, rows, xf, 'a', mode)
    draw(ax2, rows, xz, 'b', mode)
    ax1.set_ylabel('Station (land distance)')
    ax1.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.0, 1.03), ncol=3,
               frameon=False, fontsize=7.5, handlelength=2.0, columnspacing=1.4, borderaxespad=0)
    fig.subplots_adjust(left=0.105, right=0.985, top=0.895, bottom=0.065)
    FIGDIR = os.environ.get('Q_FIG_DIR') or os.path.join(HERE, os.pardir)
    out = os.path.join(FIGDIR,
                       '台集_T分量剖面_陆地计时_带时窗.png' if mode == 'land'
                       else '台集_T分量剖面_绝对计时_带时窗.png')
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print('saved: %s' % os.path.abspath(out))

# ---- 数值自检（不看图）：两套时间轴下窗的位置
for mode, _xf, _xz in MODES:
    print('\n[%s] %-5s %8s %10s %11s' % (mode, 'sta', 'T-P', 'T-S', '信号窗', ))
    for r in rows:
        _off = r['ocean'] / V_OCEAN if mode == 'land' else 0.0
        c0 = coda_t0(r) - _off
        print('   %-5s %8.1f %10.1f %6.0f-%-3.0f 尾波 %6.1f-%-5.1f' % (
            r['sta'], r['ocean'] / V_OCEAN + r['land'] / V_P - _off,
            r['ocean'] / V_OCEAN + r['land'] / V_S - _off,
            r['ocean'] / V_OCEAN + r['land'] / V_S - _off,
            r['ocean'] / V_OCEAN + r['land'] / V_S - _off + SIG_LEN, c0, c0 + CODA_LEN))
