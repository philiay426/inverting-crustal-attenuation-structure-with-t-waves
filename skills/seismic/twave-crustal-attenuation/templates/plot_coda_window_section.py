# -*- coding: utf-8 -*-
"""
2S 尾波窗检查图（13 台）：把 2S 窗画在波形上，看它落在什么位置
--------------------------------------------------------------------
CODA_MODE 控制画哪个尾波窗（★ 参数区）：
  '2S'      : 起点 = 2 x 该台 S 波到时 tt_s（stations.csv 的 tt_s 列），窗长 20 s
  'T-S+400' : 起点 = 该台 T-S 到时 + 400 s（= L_ocean/1.48 + L_land/3.6 + 400）
  'abs1200' : 全台同一绝对窗 [1200, 1220] s
图上同时画出：
  · 信号窗（橙）= [T-S, T-S+20] s，T-S = L_ocean/1.48 + L_land/3.6
  · 2S 尾波窗（红）= [2*tt_s, 2*tt_s + 20] s
  · T-P 蓝短虚线、T-S 紫短虚线
横轴：绝对计时（相对发震），因为 2S 本身是按绝对到时定义的。
数据：0_data_processing/envelope/<sta>_envelope_2-8Hz.npz 的 waveform_t（2-8 Hz T 分量）。
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
BASE = os.environ.get('Q_BASE') or os.path.abspath(
    os.path.join(HERE, '..', '..', '..'))      # process_flow/（可用 Q_BASE 覆盖）
CSV_PATH = os.path.join(BASE, '1_station_distribution', 'stations.csv')
# 包络目录：默认 BASE/0_data_processing/envelope；新布局/旧布局不同时用 Q_ENV_DIR 指过去
ENV_DIR = os.environ.get('Q_ENV_DIR') or os.path.join(BASE, '0_data_processing', 'envelope')

STAS = [x.strip() for x in os.environ.get('Q_STAS', '').split(',') if x.strip()]
if not STAS:
    raise SystemExit('请用环境变量 Q_STAS=A,B,C 指定要画的台站（顺序即纵轴顺序）')
EXPECTED_N = len(STAS)
V_OCEAN, V_P, V_S = 1.48, 6.0, 3.6
SIG_LEN = 20.0
CODA_LEN = 20.0
CODA_MODE = '2S'             # ★ '2S' / 'T-S+400' / 'abs1200'（定案 = 2S）
CODA_DT = 600.0              # 'T-S+<dt>' 里的 dt (s)
CODA_T0 = 1200.0             # 'abs<t0>' 里的 t0 (s)
# ★ 纵向排序依据：
#   'tt_ts' = 按该台理论 T-S 到时（= 传播顺序，T 波在图上依次排开）  <- 现用
#   'land'  = 按陆地距离（陆段衰减分析用，但 T 波不会连成斜线）
#   'ocean' = 按海洋路径长度（与 tt_ts 几乎同序）
#   'tt_s'  = 按 2S 参考窗所在时刻（S 到时）
SORT_BY = 'tt_ts'

ONE_PANEL = True                 # ★ True = 单面板（放 PPT 左半栏用）
XLIM = (0.0, 1500.0)             # ★ 显示时段 (s)：0-1500 覆盖 2S 参考窗与信号窗

BAR_PAD = 0.0                # 时窗色块横向不加宽（长度保持实际 20 s）
BAR_LW  = 10.0               # ★ 时窗色块的线宽 (pt)：加粗 = 方块在竖直方向更高/更醒目
AMP_FRAC = 0.45
COL, LW = 'k', 0.25
C_SIG, C_CODA, C_TP, C_TS = '#ff7f0e', '#d62728', '#1f77b4', '#9467bd'
BAR_HALF, TICK_HALF = 0.035, 0.33
FIG_W, FIG_H = 5.4, 7.6          # ★ 竖版：够高，适合 PPT 左侧一半栏
_SUF = {'2S': '2S', 'T-S+400': 'TS%g' % CODA_DT, 'abs1200': 'abs1200'}[CODA_MODE]
FIGDIR = os.environ.get('Q_FIG_DIR') or os.path.join(HERE, os.pardir)
OUT = os.path.join(FIGDIR, '台集_T分量剖面_%s时窗_按%s.png' % (_SUF, SORT_BY))
# ==================================================================

plt.rcParams.update({'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
                     'font.size': 8, 'axes.labelsize': 9,
                     'xtick.labelsize': 8, 'ytick.labelsize': 8,
                     'axes.linewidth': 0.6, 'xtick.direction': 'in', 'ytick.direction': 'in',
                     'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                     'pdf.fonttype': 42, 'ps.fonttype': 42})


def load_rows():
    raw = open(CSV_PATH).read().splitlines()
    i = [k for k, l in enumerate(raw) if l.startswith('sta,')][0]
    rows = [r for r in csv.DictReader(raw[i:]) if r.get('sta')]
    rows = [r for r in rows if r['sta'] in STAS]
    for r in rows:
        r['land'] = float(r['land_km'])
        r['ocean'] = float(r['ocean_km'])
        r['tt_s'] = float(r['tt_s'])                 # ★ S 波到时（2S 窗的锚点）
        r['tt_ts'] = r['ocean'] / V_OCEAN + r['land'] / V_S   # T-S 理论到时
        r['tt_tp'] = r['ocean'] / V_OCEAN + r['land'] / V_P   # T-P 理论到时
    _KEY = {'tt_ts': lambda r: r['tt_ts'], 'land': lambda r: r['land'],
            'ocean': lambda r: r['ocean'], 'tt_s': lambda r: r['tt_s']}[SORT_BY]
    rows.sort(key=_KEY)
    return rows


def draw(ax, rows, xlim, tag, show_labels=True):
    n = len(rows)
    ax.set_xlim(*xlim)
    ax.set_ylim(-0.9, n - 1 + 0.9)
    ax.set_yticks(np.arange(n))
    ax.set_yticklabels(['%s (%g km)' % (r['sta'], r['land']) for r in rows]
                       if show_labels else [])
    ax.set_xlabel('Time from event origin (s)')
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    if tag:
        ax.text(0.004, 1.004, '(%s)' % tag, transform=ax.transAxes, ha='left', va='bottom', fontsize=9)

    for k, r in enumerate(rows):
        z = np.load(os.path.join(ENV_DIR, '%s_envelope_2-8Hz.npz' % r['sta']))
        t = z['time'].astype(float)
        w = z['waveform_t'].astype(float)
        m = (t >= xlim[0]) & (t <= xlim[1])
        if np.any(m):
            peak = float(np.abs(w[m]).max())
            if peak > 0:
                ax.plot(t[m], k + w[m] / peak * AMP_FRAC, color=COL, lw=LW, zorder=2)
        # 理论到时
        for tv, cv, ls in ((r['tt_tp'], C_TP, (0, (4, 3))),
                           (r['tt_ts'], C_TS, (0, (1.5, 1.5)))):
            if xlim[0] <= tv <= xlim[1]:
                ax.plot([tv, tv], [k - TICK_HALF, k + TICK_HALF], color=cv, lw=0.9, ls=ls, zorder=4)
        # 信号窗（橙）
        ax.plot([r['tt_ts'] - BAR_PAD, r['tt_ts'] + SIG_LEN + BAR_PAD], [k, k], color=C_SIG,
                lw=BAR_LW, solid_capstyle='butt', zorder=5)
        # 尾波窗（红）：按 CODA_MODE 取起点
        if CODA_MODE == '2S':
            c0 = 2.0 * r['tt_s']
        elif CODA_MODE == 'T-S+400':
            c0 = r['tt_ts'] + CODA_DT
        else:
            c0 = CODA_T0
        ax.plot([c0 - BAR_PAD, c0 + CODA_LEN + BAR_PAD], [k, k], color=C_CODA,
                lw=BAR_LW, solid_capstyle='butt', zorder=5)


rows = load_rows()
assert len(rows) == EXPECTED_N, '台站数 %d != %d' % (len(rows), EXPECTED_N)

if ONE_PANEL:
    fig, ax1 = plt.subplots(figsize=(FIG_W, FIG_H))
    draw(ax1, rows, XLIM, '')
else:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(FIG_W, FIG_H), sharey=True,
                                   gridspec_kw=dict(height_ratios=[1.0, 1.35], hspace=0.16))
    draw(ax1, rows, XLIM, 'a')
    draw(ax2, rows, XLIM, 'b')
ax1.set_ylabel('Station (land distance)' + (', sorted by T-S' if SORT_BY == 'tt_ts' else ''))

handles = [Line2D([], [], color=COL, lw=0.9, label='T component (2-8 Hz)'),
           Line2D([], [], color=C_SIG, lw=4.2, label='Signal window [T-S, T-S+20 s]'),
           Line2D([], [], color=C_CODA, lw=4.2, label={'2S': '2S coda window [2*tt_s, +20 s]',
                                                       'T-S+400': 'Coda window [T-S+%g, +20 s]' % CODA_DT,
                                                       'abs1200': 'Coda window [1200, 1220] s'}[CODA_MODE]),
           Line2D([], [], color=C_TP, lw=0.9, ls=(0, (4, 3)), label='T-P theory'),
           Line2D([], [], color=C_TS, lw=0.9, ls=(0, (1.5, 1.5)), label='T-S theory')]
_lg = ax1.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.0, 1.02), ncol=2,
           frameon=False, fontsize=6.8, handlelength=1.8, columnspacing=1.0,
           labelspacing=0.35, borderaxespad=0)

fig.subplots_adjust(left=0.235 if ONE_PANEL else 0.105, right=0.985,
                    top=0.845 if ONE_PANEL else 0.895, bottom=0.055)
fig.canvas.draw()
_r = fig.canvas.get_renderer()
_fw, _fh = fig.get_size_inches()
_bb = _lg.get_window_extent(renderer=_r).transformed(fig.transFigure.inverted())
_ab = ax1.get_window_extent(renderer=_r).transformed(fig.transFigure.inverted())
print('图幅 %.2f x %.2f in | 轴框 x %.3f-%.3f | 图例 x %.3f-%.3f, y %.3f-%.3f'
      % (_fw, _fh, _ab.x0, _ab.x1, _bb.x0, _bb.x1, _bb.y0, _bb.y1))
if _bb.x1 > 1.0 or _bb.x0 < 0.0:
    print('WARN: 图例横向出框')
if _bb.y1 > 1.0:
    print('WARN: 图例纵向出框')
fig.savefig(OUT, dpi=300)
print('saved: %s' % os.path.abspath(OUT))

_LAB = '2S' if CODA_MODE == '2S' else ('T-S+%g' % CODA_DT if CODA_MODE == 'T-S+400' else '[%g,%g]' % (CODA_T0, CODA_T0 + CODA_LEN))
print('\n[%s] %-5s %6s %9s %9s %10s' % (_LAB, 'sta', 'land', '尾波窗', 'T-S到时', '相对T-S(s)'))
for r in rows:
    c0 = 2.0 * r['tt_s'] if CODA_MODE == '2S' else (r['tt_ts'] + CODA_DT if CODA_MODE == 'T-S+400' else CODA_T0)
    print('%-5s %6.0f %6.0f-%-5.0f %9.1f %+10.0f' % (
        r['sta'], r['land'], c0, c0 + CODA_LEN, r['tt_ts'], c0 - r['tt_ts']))
