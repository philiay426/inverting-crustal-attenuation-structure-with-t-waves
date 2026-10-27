#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_grid_radius_choice.py —— 慢度网格半径「三档口径」复算（只用标准库）

用途：被问"网格取多少 / 为什么取整 / 不做慢度上限时怎么定"时先跑它，别手算。

三档口径互不相通，别互相套用：
  ① 受限网格（默认）：SMAX = round(抗混叠上限 / 5)·5           —— 用户口语"就近取整"
  ② 用户点名"不做慢度上限"：SMAX = ceil(抗混叠上限 / 5)·5      —— "就近向右取整"
     此时并复核"复制峰间距 111.195/(f_hi·d_min) 是否落在显示半径之外"
     ⇒ 在范围外则图面内不会出现混叠复制峰（脚本仍会打 ⚠️ 抗混叠违规，属预期）
  ③ 坐标轴半宽（纯观感，与 SMAX 无关，不在本脚本范围）

抗混叠上限 = CONST / f_hi，CONST = 111.195/(2·d_min_km)（每子阵自算，别跨子阵抄）

用法：
  python3 check_grid_radius_choice.py --window 60 170 --band 0.020 0.032 --dmin 95.594 --target-u 18.53
  python3 check_grid_radius_choice.py --window 142.5 293.1 --band 0.005 0.013 --const 0.581599 --step 0.35
  python3 check_grid_radius_choice.py --window 74.2 191.2 --band 0.008 0.0322 --const 0.550231   # G1
"""
import argparse
import math

KM_PER_DEG = 111.195


def odd_near(x):
    """取最接近的奇数（网格点数为奇数时中心格恰好落在 0）。"""
    n = int(round(x))
    if n % 2 == 0:
        n += 1
    return max(n, 3)


def main():
    ap = argparse.ArgumentParser(description="慢度网格半径三档口径复算（只用标准库）")
    ap.add_argument("--window", nargs=2, type=float, required=True, help="计时窗 起 止 (s)")
    ap.add_argument("--band", nargs=2, type=float, required=True, help="频带 下 上 (Hz)")
    ap.add_argument("--dt", type=float, default=0.5, help="采样间隔 (s)，包络链默认 0.5")
    ap.add_argument("--dmin", type=float, default=None, help="子阵最小台间距 (km) ⇒ 常数 = 111.195/(2·dmin)")
    ap.add_argument("--const", type=float, default=None, help="直接给抗混叠常数（G1 0.550231 / G2 0.581599）")
    ap.add_argument("--step", type=float, default=0.35, help="目标慢度步长 Δs（项目 0.35）")
    ap.add_argument("--round-to", type=float, default=5.0, help="取整到该倍数（项目 5）")
    ap.add_argument("--target-u", nargs="*", type=float, default=[],
                    help="目标慢度 s/°（T-P 18.53 / T-S 30.89），检查装不装得下")
    a = ap.parse_args()

    if a.const is None and a.dmin is None:
        raise SystemExit("给 --dmin 或 --const 至少一个")
    const = a.const if a.const else KM_PER_DEG / (2.0 * a.dmin)
    dmin = a.dmin if a.dmin else KM_PER_DEG / (2.0 * const)

    T = a.window[1] - a.window[0]
    nsp = 2 ** math.ceil(math.log2(T / a.dt))
    df = 1.0 / (nsp * a.dt)
    grid = [k * df for k in range(math.ceil(a.band[0] / df), math.floor(a.band[1] / df) + 1)]
    legal = [f for f in grid if f >= 1.0 / T]

    print(f"窗 {a.window[0]:g}–{a.window[1]:g} s（T={T:g} s，dt={a.dt:g}）⇒ nsp {nsp}、Δf {df:.8f} Hz")
    print(f"频带栅格点 : {['%.7f' % f for f in grid] or '无'}")
    print(f"合法频点   : {['%.7f' % f for f in legal] or '无（f < 1/T，该窗+带没有可用频点）'}")
    if not legal:
        return
    fhi = max(legal)
    lim = const / fhi
    print(f"f_hi {fhi:.8f} Hz ⇒ 抗混叠上限 = 常数/f_hi = {lim:.2f} s/°"
          f"（常数 {const:.6f} = 111.195/(2·{dmin:.3f} km)）")

    choices = {"① 受限网格·就近取整": round(lim / a.round_to) * a.round_to,
               "② 不做上限·就近向右取整": math.ceil(lim / a.round_to) * a.round_to}
    for lab, R in choices.items():
        R = max(R, a.round_to)
        choices[lab] = R
        n = odd_near(2.0 * R / a.step + 1.0)
        print(f"{lab}: SMAX {R:g} s/° · SPTS {n:g}（Δs {2.0 * R / (n - 1):.4f}，无混叠复核 f_hi·SMAX = {fhi * R:.3f}）")

    rep = KM_PER_DEG / (fhi * dmin)
    print(f"复制峰间距 111.195/(f_hi·d_min) = {rep:.1f} s/°")
    for lab, R in choices.items():
        ok = rep > R
        print(f"   {lab} 网格 {R:g}: 复制峰{'不在' if ok else '落在'}显示半径内 ⇒ "
              + ("图面内不会出现混叠复制峰，读数可用" if ok else "可能出现混叠复制峰，图上不能直接读数"))

    for u in a.target_u:
        R = choices["② 不做上限·就近向右取整"]
        print(f"   目标慢度 {u:g} s/°：{'装得下' if u <= R else '装不下'}（向右取整网格 {R:g}）；1.7×u = {1.7 * u:.1f} s/°")


if __name__ == "__main__":
    main()
