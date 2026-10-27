#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
step2_station_set_select.py   (2026-10-16, 第 2 版: 全子集枚举)

挑「最合适」的 Q 反演组合（台站池由 Q_POOL 传入）。

做法：先用 Q 反演引擎**自己的口径**把全池跑一遍（逐台逐频点 A_d/A_c 落在它的 DETAIL 里，
      滤波/旋转/窗位/几何扩散全部沿用引擎），再对每个 >=MIN_N 台的子集做纯算术的斜率拟合与指标统计。
      逐频点筛选口径 = 引擎那一条: A_d/A_c >= SNR_MIN(=1.0)（比值, 不是 dB）。

指标（全部来自数据本身，不拿论文的 Q 当靶子）：
  覆盖 : 台数 / 陆距范围 / 非空距离档数 / 最大相邻空档 / 方位跨度
  偏振 : 组合内各台在**定案信号窗** [该台理论 T-S, +20 s] 的 2.5-4.0 Hz frac_T 中位数, 不达标(<0.5)台数
  拟合 : 可用频点数（9 个频点里斜率为负的个数）/ 逐频点残差散布中位 / 各台平均残差最大绝对值
  拟合值: 幂律 Q0、alpha（误差 t*se/1.96, 导师 plot_alpha.m 口径）——只作报告, 不参与选台

输出：终端表 + <DATA>/台集选择_对比_全子集_<日期>.csv

【不改引擎调用链】只 import 它取 DETAIL（会顺带产出 q_pool*_*.csv / 图, 用完清掉）。
用法: Q_BASE=<process_flow> Q_POOL=A,B,C,... python select_station_set.py
"""
import csv, importlib.util, math, os, sys
from itertools import combinations
import numpy as np

# ---------------- ★ 参数区 ----------------
BASE = os.environ.get('Q_BASE', '').strip()
if not BASE:
    raise SystemExit('请用环境变量 Q_BASE 指定 process_flow 目录')
DATA = os.environ.get('Q_DATA') or (BASE + '/4_attenuation/数据')
QS5 = os.environ.get('Q_ENGINE') or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                'twave_q_inversion.py')
STAC = BASE + '/1_station_distribution/stations.csv'
SCAN = os.environ.get('Q_SCAN') or os.path.join(DATA, 'pol_window_scan.csv')
POOL = [x.strip() for x in os.environ.get('Q_POOL', '').split(',') if x.strip()]
if not POOL:
    raise SystemExit('请用环境变量 Q_POOL=A,B,C,... 指定候选台站池')
FCS = [2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0]
POL_F = [2.5, 3.0, 3.5, 4.0]       # 偏振判据用频点
POL_MIN = 0.5                      # frac_T 下限（>0.5 即切向大于径向）
SEL_L, SEL_S = 20, 0               # 定案信号窗
CM_MAIN, CM_ALT = 'T-S+400', 'abs1200'
MIN_N = 5                          # 子集最少台数
MAX_ANOM = 0.45                    # 「单台异常」上限 (dex): 各台平均残差绝对值 <= 此值 = 在允许范围内
TAGD = os.environ.get('Q_TAGDATE', '00000000')   # 输出文件名里的日期标记
BINS = [(0, 100), (100, 200), (200, 300), (300, 400), (400, 500), (500, 700)]   # 陆距档 (km)，按事件改
SPOT = ([('all %d' % len(POOL), POOL)] +
        # 逐个去掉一台、两两去掉（自动生成，不用写死台名）
        [('drop %s' % a, [x for x in POOL if x != a]) for a in POOL] +
        [('drop %s+%s' % (a, b), [x for x in POOL if x not in (a, b)])
         for i, a in enumerate(POOL) for b in POOL[i+1:]])
# -----------------------------------------

os.environ.update(Q_STAS=','.join(POOL), Q_TAG='pool%d' % len(POOL), Q_CODA_MODES='T-S+400,2S,abs1200',
                  Q_CODA_MAIN=CM_MAIN, Q_FIGNAME='tmp_pool%d' % len(POOL))
_spec = importlib.util.spec_from_file_location('q_s5mod', QS5)
_q = importlib.util.module_from_spec(_spec)
_stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')
_spec.loader.exec_module(_q)                      # 跑一次全池（产物 q_pool*_*）
sys.stdout = _stdout
print('反演引擎口径: 台站池 %d 台 | 信号窗 %.0f s (偏移 %+g) | 尾波口径 %s | v=%.1f km/s'
      % (len(_q.STAS), _q.SIG_WIN, _q.SIG_OFF, _q.CODA_MODES, _q.V_WAVE))

# ---- 台站表 + 定案窗偏振 ----
lines = [l for l in open(STAC) if l.strip() and not l.startswith('#')]
head = [c.strip() for c in lines[0].strip().split(',')]
meta = {d['sta']: d for d in (dict(zip(head, [c.strip() for c in l.strip().split(',')])) for l in lines[1:])}
pol = {}
for r in csv.DictReader(open(SCAN)):
    if int(r['L_s']) == SEL_L and int(r['s_s']) == SEL_S and float(r['fc_Hz']) in POL_F:
        pol.setdefault(r['sta'], []).append(float(r['frac_T']))
POLMED = {s: (float(np.median(v)) if v else float('nan')) for s, v in pol.items()}
DET = {(d['sta'], d['f']): d for d in _q.DETAIL}


def points(sub, cm, fc):
    """该子集该频点参与拟合的台 (x=r_km, A_d, A_c, sta), 条件同 q_s5.py"""
    P = []
    for s in sub:
        d = DET.get((s, fc))
        if d is None:
            continue
        Ad, Ac = d['A_theory'], d['A_coda_%s' % cm]
        if np.isfinite(Ad) and np.isfinite(Ac) and Ac > 0 and Ad > 0 and (Ad / Ac) >= _q.SNR_MIN:
            P.append((d['r_km'], Ad, Ac, s))
    return P


def line_of(sub, cm, fc):
    """拟合直线, 返回 (k, 逐台残差 dict dex, 点数)"""
    P = points(sub, cm, fc)
    if len(P) < 3:
        return np.nan, {}, 0
    x = np.array([p[0] for p in P])
    y = np.array([math.log(p[1] * _q.geom_factor(p[0]) / p[2]) for p in P])
    k, a = np.polyfit(x, y, 1)
    return k, {p[3]: v / math.log(10) for p, v in zip(P, y - (k * x + a))}, len(P)


def metrics(sub, cm):
    """可用频点数 / 逐频点残差散布中位 / 各台平均残差的绝对值最大者 / 可用频点 Q 中位"""
    n_us, allres, st_res, medQ = 0, [], {}, []
    for fc in FCS:
        k, rs, np_ = line_of(sub, cm, fc)
        if np_ < 3:
            continue
        allres += [abs(v) for v in rs.values()]
        for s, v in rs.items():
            st_res.setdefault(s, []).append(v)
        if k < 0:
            n_us += 1
            medQ.append(-math.pi * fc / (k * _q.V_WAVE))
    return dict(n_usable=n_us, rms_med=float(np.median(allres)) if allres else np.nan,
                st_mean={s: float(np.mean(v)) for s, v in st_res.items()}, medQ=float(np.median(medQ)) if medQ else np.nan)


def fit_q(sub, cm):
    """幂律拟合 Q0/alpha（误差 t*se/1.96）"""
    from scipy import stats
    f, Q = [], []
    for fc in FCS:
        k, rs, np_ = line_of(sub, cm, fc)
        if np_ < 3 or k >= 0:
            continue
        f.append(fc); Q.append(-math.pi * fc / (k * _q.V_WAVE))
    n = len(f)
    if n < 3:
        return dict(n=n, Q0=np.nan, sQ0=np.nan, al=np.nan, sal=np.nan)
    x, y = np.log(np.array(f)), np.log(np.array(Q))
    p = np.polyfit(x, y, 1)
    res = y - np.polyval(p, x)
    t = stats.t.ppf(0.975, n - 2) / 1.96
    # 斜率/截距标准误, 按定义算 (★ 2026-10-16 修正见下)
    s2 = float(res @ res) / (n - 2)                       # 残差方差
    Sxx = float(((x - x.mean()) ** 2).sum())
    se = np.array([math.sqrt(s2 / Sxx), math.sqrt(s2 * (1.0 / n + x.mean() ** 2 / Sxx))])
    return dict(n=n, Q0=math.exp(p[1]), sQ0=math.exp(p[1]) * t * se[1], al=p[0], sal=t * se[0])


def coverage(sub):
    r = sorted(float(meta[s]['land_km']) for s in sub)
    g = [r[i + 1] - r[i] for i in range(len(r) - 1)]
    return dict(rmin=r[0], rmax=r[-1], ratio=(r[-1] / r[0] if r[0] else np.nan),
                bins=sum(1 for a, b in BINS if any(a <= x < b for x in r)), gap=max(g) if g else 0.0,
                az=(max(float(meta[s]['az']) for s in sub) - min(float(meta[s]['az']) for s in sub)),
                polmed=float(np.median([POLMED[s] for s in sub])),
                polbad=sum(1 for s in sub if not (POLMED[s] > POL_MIN)))


# ---------- (一) 逐台异常量（全 10 台基准） ----------
base = metrics(POOL, CM_MAIN)
print('\n=== (一) 逐台异常量：全 10 台基准下, 各台在 9 个频点的残差（dex, ln(A_d*Z/A_c) 与拟合直线的差）===')
print('   正值 = 该台的 直达/尾波 比高于趋势; 负值 = 低于趋势（差 0.3 dex = 2 倍）')
print('%-6s %6s %7s | %s | %7s' % ('sta', 'land', 'fracT中位', ' '.join('%6.2f' % f for f in FCS), '平均'))
for s in sorted(POOL, key=lambda x: -abs(base['st_mean'].get(x, 0))):
    per = [line_of(POOL, CM_MAIN, fc)[1].get(s, np.nan) for fc in FCS]
    print('%-6s %6.0f %7.2f | %s | %+7.2f %s'
          % (s, float(meta[s]['land_km']), POLMED[s], ' '.join('%6.2f' % v for v in per),
             base['st_mean'].get(s, np.nan),
             '<== 超 %.2f dex' % MAX_ANOM if abs(base['st_mean'].get(s, 0)) > MAX_ANOM else ''))

# ---------- (二) 点名组合 ----------
print('\n=== (二) 点名组合（逐个去台 / 两两去台的对照）===')
print('%-20s %2s %5s %5s %3s %5s %4s %6s %6s | %5s %5s %7s %7s | %-22s %-22s'
      % ('组合', 'N', 'rmin', 'rmax', '档', '空档', '偏', '偏振中位', '方位', 'T可用', 'T残差', 'T台最大', 'T中位Q',
         'T-S+400 拟合', 'abs1200 拟合'))


def fitstr(q):
    return ('Q0=%.0f+-%.0f a=%+.2f+-%.2f(%d)' % (q['Q0'], q['sQ0'], q['al'], q['sal'], q['n'])
            if np.isfinite(q['al']) else '拟合不出(n=%d)' % q['n'])


for name, sub in SPOT:
    cv = coverage(sub)
    m, q, m2, q2 = metrics(sub, CM_MAIN), fit_q(sub, CM_MAIN), metrics(sub, CM_ALT), fit_q(sub, CM_ALT)
    stmax = max((abs(v) for v in m['st_mean'].values()), default=np.nan)
    print('%-20s %2d %5.0f %5.0f %3d %5.0f %4d %6.2f %6.1f | %5d %5.3f %7.3f %7.0f | %-22s %-22s'
          % (name, len(sub), cv['rmin'], cv['rmax'], cv['bins'], cv['gap'], cv['polbad'], cv['polmed'], cv['az'],
             m['n_usable'], m['rms_med'], stmax, m['medQ'], fitstr(q), fitstr(q2)))

# ---------- (三) 全子集枚举 ----------
rows = []
for n in range(MIN_N, len(POOL) + 1):
    for sub in combinations(POOL, n):
        sub = list(sub)
        cv = coverage(sub)
        m, q = metrics(sub, CM_MAIN), fit_q(sub, CM_MAIN)
        m2, q2 = metrics(sub, CM_ALT), fit_q(sub, CM_ALT)
        rows.append(dict(tag='+'.join(sorted(sub, key=lambda x: float(meta[x]['land_km']))), sub=sub, n=n, **cv,
                         n_usable=m['n_usable'], rms=m['rms_med'],
                         stmax=max((abs(v) for v in m['st_mean'].values()), default=np.nan), medQ=m['medQ'],
                         Q0=q['Q0'], sQ0=q['sQ0'], al=q['al'], sal=q['sal'], nf=q['n'],
                         n2=m2['n_usable'], rms2=m2['rms_med'],
                         stmax2=max((abs(v) for v in m2['st_mean'].values()), default=np.nan), medQ2=m2['medQ'],
                         Q02=q2['Q0'], sQ02=q2['sQ0'], al2=q2['al'], sal2=q2['sal'], nf2=q2['n']))

print('\n=== (三) 638 个子集里, 距离档 6/6 且各台平均残差最大 <= %.2f dex 的 ===' % MAX_ANOM)
cand = [r for r in rows if r['bins'] >= 6 and r['stmax'] <= MAX_ANOM]
print('共 %d 个; 按「可用频点多 -> 单台异常小 -> 残差小」排前 14:' % len(cand))
print('%-52s %2s %4s %5s %6s %6s %7s %7s | %-22s'
      % ('组合(按陆距)', 'N', '偏', 'T可用', 'T残差', 'T台最大', 'T中位Q', 'r范围', 'T-S+400 拟合'))
for r in sorted(cand, key=lambda r: (-r['n_usable'], r['stmax'], r['rms']))[:14]:
    print('%-52s %2d %4d %5d %6.3f %6.3f %7.0f %3.0f-%3.0f | %-22s'
          % (r['tag'], r['n'], r['polbad'], r['n_usable'], r['rms'], r['stmax'], r['medQ'], r['rmin'], r['rmax'],
             'Q0=%.0f a=%+.2f+-%.2f(%d)' % (r['Q0'], r['al'], r['sal'], r['nf']) if np.isfinite(r['al']) else '拟合不出'))

print('\n=== (四) 排序规则: 6/6 档 -> 偏振不达标台最少 -> 可用频点最多 -> 单台异常最小 -> 残差最小（前 10）===')
print('%-52s %2s %4s %5s %6s %6s %7s | %-22s %-22s'
      % ('组合(按陆距)', 'N', '偏', 'T可用', 'T台最大', 'T残差', 'T中位Q', 'T-S+400 拟合', 'abs1200 拟合'))
for r in sorted(rows, key=lambda r: (-r['bins'], r['polbad'], -r['n_usable'], r['stmax'], r['rms']))[:10]:
    print('%-52s %2d %4d %5d %6.3f %6.3f %7.0f | %-22s %-22s'
          % (r['tag'], r['n'], r['polbad'], r['n_usable'], r['stmax'], r['rms'], r['medQ'],
             'Q0=%.0f a=%+.2f+-%.2f(%d)' % (r['Q0'], r['al'], r['sal'], r['nf']) if np.isfinite(r['al']) else '拟合不出',
             'Q0=%.0f a=%+.2f+-%.2f(%d)' % (r['Q02'], r['al2'], r['sal2'], r['nf2']) if np.isfinite(r['al2']) else '拟合不出'))

with open(os.path.join(DATA, '台集选择_对比_全子集_%s.csv' % TAGD), 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['# 台站池 10 台的全子集枚举 (含 >=%d 台, 共 %d 个)' % (MIN_N, len(rows))])
    w.writerow(['# 信号窗 = 理论 T-S 起 20 s (定案); 逐频点筛选 = A_d/A_c >= 1 (q_s5.py SNR_MIN); 指标不含论文 Q'])
    w.writerow(['# 误差 = t(0.975,n-2)*se/1.96 (导师 plot_alpha.m 口径); CM_MAIN=%s, CM_ALT=%s' % (CM_MAIN, CM_ALT)])
    w.writerow(['stations_by_land_km', 'n_sta', 'r_min_km', 'r_max_km', 'bins_used', 'max_gap_km', 'az_span_deg',
                'fracT_median', 'pol_bad', 'T_n_usable', 'T_resid_med_dex', 'T_station_max_dex', 'T_median_Q',
                'T_Q0', 'T_sigma_Q0', 'T_alpha', 'T_sigma_alpha', 'T_n_freq',
                'A_n_usable', 'A_resid_med_dex', 'A_station_max_dex', 'A_median_Q',
                'A_Q0', 'A_sigma_Q0', 'A_alpha', 'A_sigma_alpha', 'A_n_freq'])
    g = lambda v: '%.4g' % v if np.isfinite(v) else ''
    for r in rows:
        w.writerow([r['tag'], r['n'], '%.0f' % r['rmin'], '%.0f' % r['rmax'], r['bins'], '%.0f' % r['gap'],
                    '%.1f' % r['az'], '%.3f' % r['polmed'], r['polbad'],
                    r['n_usable'], '%.4f' % r['rms'], '%.4f' % r['stmax'], g(r['medQ']),
                    g(r['Q0']), g(r['sQ0']), g(r['al']), g(r['sal']), r['nf'],
                    r['n2'], '%.4f' % r['rms2'], '%.4f' % r['stmax2'], g(r['medQ2']),
                    g(r['Q02']), g(r['sQ02']), g(r['al2']), g(r['sal2']), r['nf2']])
print('\n表 -> %s/台集选择_对比_全子集_%s.csv' % (DATA, TAGD))
