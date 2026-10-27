#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_timing_convention_params.py —— 换计时口径（绝对计时 / 陆地计时）时
F-K 参数要不要改的复算工具（只用标准库）。

【为什么需要它】用户在问"都改用陆地计时，是不是所有参数都要更新"时，期望的回答是"三类"：
  ① 随口径走 —— 窗的秒数、τ、读数含义；
  ② 不随口径走 —— nsp / Δf / 带内频点 / 合法 SMAX / SPTS / Δs（只要窗长落在同一个 2ⁿ 档）；
  ③ 只有跨了 2ⁿ 阶梯（T/dt 越过 128/256/512 s）才连带重算 ②。
数字一律现算，不手算、不凭印象：本脚本把"该不该沿用同一套参数"变成一行可核对输出。

【打印四类量】
  A 子阵内 ocean_km / land_km 跨度与"上岸时刻"跨度（= 两版窗长之差的上限）
  B 由 stations.csv 现算的陆地窗 [min(L_land/v)−30, max(L_land/v)+50]（v = 6.0 / 3.6 km/s）
  C 每个给定窗的 nsp / Δf / 带内频点 / f_hi / 合法 SMAX（可选用 --smax/--spts 核对既有取值）
  D 两口径 τ 判据：max|τ|/T 是**假 FAIL**（陆地计时 τ 全正、只差公共偏移），展布/T 才是判据

【示例】
python scripts/check_timing_convention_params.py \\
  --stations-csv <case>/process_flow/1_station_distribution/stations.csv \\
  --stas SA16,SA17,SA18,SA22,SA23,SA25,SA29 \\
  --band 0.008 0.0322 --abs-window 1173 1292 --land-window 74.2 191.2 --dt 0.5 \\
  --smax 17.1 --spts 99
"""

import argparse
import csv
import math
import sys

# 无混叠上限系数：SMAX <= A_LIMIT / f_hi（A_LIMIT = 111.195/(2·d_min)；G1 的 d_min≈101.05 km ⇒ 0.550231）
A_LIMIT_DEFAULT = 0.550231
OCEAN_V = 1.48          # 海段声速 km/s（论文式(1)）
V_P = 6.0               # 陆段 P 波速度（T-P 用）
V_S = 3.6               # 陆段 S 波速度（T-S 用）


def load_stations(path):
    with open(path) as f:
        return {r['sta']: r for r in csv.DictReader(l for l in f if not l.startswith('#'))}


def nsp_of(T, dt):
    """fk4d 的零填充点数：2^ceil(log2(T/dt))"""
    return 2 ** math.ceil(math.log2(T / dt))


def grid_report(name, T, band, dt, smax=None, spts=None):
    """打印某个窗的频率栅格：nsp/Δf/带内频点/f_hi/合法 SMAX；返回 (nsp, Δf, freqs, f_hi)"""
    n = nsp_of(T, dt)
    df = 1.0 / (n * dt)
    klo = max(1, math.ceil(band[0] / df))
    khi = math.floor(band[1] / df)
    kk = [k for k in range(klo, khi + 1) if k * df >= 1.0 / T]   # 低于 1/T 的档不算（本窗无此频点）
    freqs = [k * df for k in kk]
    head = f'  {name:<14} T={T:7.1f}s  nsp={n:5d}  df={df:.8f} Hz  带内频点={len(freqs)}'
    if not freqs:
        print(head + '  —— 本窗带内无合法频点')
        return n, df, [], None
    f_hi = freqs[-1]
    lim = A_LIMIT_DEFAULT / f_hi
    print(head + f'  [{ ", ".join(f"{x:.6f}" for x in freqs)}]  f_hi={f_hi:.6f}')
    print(f'{"":16} 合法 SMAX <= {lim:.2f}（A/f_hi）  窗长下限 T >= 1/f_hi = {1/f_hi:.1f} s'
          + (f'   脚本用 SMAX={smax} / SPTS={spts}' if smax else ''))
    if smax:
        ok_ali = smax * f_hi <= A_LIMIT_DEFAULT
        _msg = "PASS" if ok_ali else "FAIL（混叠：把 SMAX 收到 %.2f 或切带顶）" % (A_LIMIT_DEFAULT / f_hi)
        print(f'{"":16} 无混叠复核 f_hi·SMAX = {smax*f_hi:.4f} <= {A_LIMIT_DEFAULT} ? {_msg}')
        if spts:
            ds = 2.0 * smax / (spts - 1)
            print(f'{"":16} 步长 Δs = 2·SMAX/(SPTS−1) = {ds:.4f} s/deg；'
                  f'网格可容 |u| <= (1−k)·SMAX，目标慢度须 fit 在 {smax:.2f} 之内')
    return n, df, freqs, f_hi


def main():
    ap = argparse.ArgumentParser(description='换计时口径时 F-K 参数是否需重算（只用标准库）')
    ap.add_argument('--stations-csv', required=True, help='stations.csv（15 列，含 ocean_km/land_km）')
    ap.add_argument('--stas', required=True, help='子阵台站名，逗号分隔')
    ap.add_argument('--band', nargs=2, type=float, metavar=('FLO', 'FHI'), help='FRNG 频带 Hz')
    ap.add_argument('--abs-window', nargs=2, type=float, metavar=('T0', 'T1'), help='绝对计时窗（相对发震 s）')
    ap.add_argument('--land-window', nargs=2, type=float, metavar=('T0', 'T1'), help='陆地计时窗（0 = 上岸 s）')
    ap.add_argument('--dt', type=float, default=0.5, help='输入采样间隔 s（本项目包络 2 Hz ⇒ 0.5）')
    ap.add_argument('--smax', type=float, help='可选：脚本现用的 SMAX，做无混叠核对')
    ap.add_argument('--spts', type=int, help='可选：脚本现用的 SPTS，一并给 Δs')
    a = ap.parse_args()

    rows = load_stations(a.stations_csv)
    stas = [s.strip() for s in a.stas.split(',') if s.strip()]
    miss = [s for s in stas if s not in rows]
    if miss:
        sys.exit(f'stations.csv 里找不到这些台站：{miss}')

    ocean = [float(rows[s]['ocean_km']) for s in stas]
    land = [float(rows[s]['land_km']) for s in stas]
    tocean = [o / OCEAN_V for o in ocean]
    tp = [l / V_P for l in land]
    ts = [l / V_S for l in land]

    print('=' * 96)
    print(f'A. 几何：{len(stas)} 台  ocean_km {min(ocean):.1f}–{max(ocean):.1f}（跨度 {max(ocean)-min(ocean):.1f} km）'
          f'  land_km {min(land):.1f}–{max(land):.1f}（跨度 {max(land)-min(land):.1f} km）')
    print(f'   上岸时刻跨度 = ocean 跨度/1.48 = {max(tocean)-min(tocean):.1f} s'
          f'  ⇒ 两版同名窗的窗长之差最多这么大（决定会不会跳 2ⁿ 阶梯）')

    win = {}
    for tag, tau in (('T-P', tp), ('T-S', ts)):
        w = (min(tau) - 30.0, max(tau) + 50.0)
        win[tag] = w
        print(f'B. 陆地窗（{tag}, v = {V_P if tag == "T-P" else V_S} km/s）= '
              f'[min(tau)−30, max+50] = [{w[0]:.1f}, {w[1]:.1f}]  长 {w[1]-w[0]:.1f} s')

    print('D. τ 判据（同一窗长 T 下两口径各算一遍）')
    if a.land_window and a.abs_window:
        pairs = [('T-P', tp, (a.land_window[1] - a.land_window[0]), (a.abs_window[1] - a.abs_window[0])),
                 ('T-S', ts, (a.land_window[1] - a.land_window[0]), (a.abs_window[1] - a.abs_window[0]))]
    else:
        pairs = [('T-P', tp, win['T-P'][1] - win['T-P'][0], None),
                 ('T-S', ts, win['T-S'][1] - win['T-S'][0], None)]
    for tag, tau, Tland, Tabs in pairs:
        spread = max(tau) - min(tau)
        mx = max(abs(t) for t in tau)
        if Tabs:
            print(f'   {tag}: max|tau|/T_abs = {mx/Tabs:.3f}（假 FAIL，别用）'
                  f'  展布/T_abs = {spread/Tabs:.3f} | 陆地窗 T={Tland:.1f}s → max|tau|/T = {mx/Tland:.3f}（假 FAIL）'
                  f'  展布/T = {spread/Tland:.3f}')
        else:
            print(f'   {tag}: 陆地窗 T={Tland:.1f}s → max|tau|/T = {mx/Tland:.3f}（假 FAIL，τ 全正）'
                  f'  真判据 展布/T = {spread/Tland:.3f}')

    if not a.band:
        print('（未给 --band，跳过 C 频率栅格；给上即可核对两版窗的 nsp/Δf/频点是否相同）')
        return

    print('=' * 96)
    print('C. 频率栅格（同一 --band 跑两版窗）：')
    res = {}
    if a.abs_window:
        res['abs'] = grid_report('绝对计时', a.abs_window[1] - a.abs_window[0], a.band, a.dt, a.smax, a.spts)
    if a.land_window:
        res['land'] = grid_report('陆地计时', a.land_window[1] - a.land_window[0], a.band, a.dt, a.smax, a.spts)
    if a.abs_window and a.land_window:
        same = (res['abs'][0] == res['land'][0] and abs(res['abs'][1] - res['land'][1]) < 1e-12
                and [f'{x:.6f}' for x in res['abs'][2]] == [f'{x:.6f}' for x in res['land'][2]])
        print('=' * 96)
        if same:
            print('⇒ 两版 nsp/Δf/带内频点完全相同 ⇒ SMAX/SPTS/Δs/method/分量 一字不用改，只有窗在时间轴上的位置不同。')
        else:
            print('⚠ 两版 nsp 不同（窗长跨了 2ⁿ 阶梯）⇒ Δf 变、频点变、SMAX/SPTS/Δs 必须重算，'
                  '不能沿用另一口径的那套参数。')


if __name__ == '__main__':
    main()
