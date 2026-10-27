#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
F-K 功率图"有没有 ±u 镜像/中心对称"的复核工具（只用标准库）。

用途：被问"另外一个峰是不是镜像/对称影"时，别用互谱 Hermitian 口推，直接量。
  python check_power_map_symmetry.py <beam.txt> [<beam.txt> ...]

输入：fn 导出的 beam txt（列序 `x y beam_dB`，y 按行递减——按文件原序建网格即可，勿升序 sort）。
输出：
  · 全局峰 (u, BAZ)、反极点 (u, BAZ+180°) 的 dB 落差
  · 全网格 max|B(s) − B(−s)|（点对称量）与"差 > 1 dB 的格点占比"
判读：
  · max|Δ| ≈ 0 且占比 ≈ 0  ⇒ 这张图真的中心对称（几乎不会发生）
  · max|Δ| 接近峰值、反极点落差 ≫ 10 dB ⇒ 无镜像峰，反方向的峰只能是：
    ① 原点针尖（s≈0，方位由噪声定） ② 混叠复制像 u ± 111.195/(f·d_min) ③ 真实次级/泄漏震相
  （合成注入实测：max|Δ| = 0.96；注入 32.5 s/°@319° 的反极点只有 0.041，相对峰 −13.9 dB。）
"""
import math
import sys


def load(path):
    xs, ys, rows = set(), set(), []
    with open(path) as f:
        for ln in f:
            p = ln.split()
            if len(p) < 3:
                continue
            try:
                x, y, b = float(p[0]), float(p[1]), float(p[2])
            except ValueError:
                continue
            xs.add(x); ys.add(y); rows.append((x, y, b))
    return sorted(xs), sorted(ys), rows


def report(path):
    xs, ys, rows = load(path)
    if not rows:
        print(f"[skip] {path}: 没解析到数据")
        return
    X = {v: i for i, v in enumerate(xs)}
    Y = {v: i for i, v in enumerate(ys)}
    nx, ny = len(xs), len(ys)
    grid = [[None] * ny for _ in range(nx)]
    for x, y, b in rows:
        grid[X[x]][Y[y]] = b
    flat = [(b, x, y) for x, y, b in rows]
    bmax, x0, y0 = max(flat)
    u = math.hypot(x0, y0)
    baz = math.degrees(math.atan2(x0, y0)) % 360.0
    # 反极点（网格对称 ±SMAX 且步长相同 ⇒ 索引镜像成立）
    anti = grid[nx - 1 - X[x0]][ny - 1 - Y[y0]]
    # 点对称量：B(s) vs B(−s)
    diffs, pairs, big = [], 0, 0
    for i in range(nx):
        for j in range(ny):
            a, b = grid[i][j], grid[nx - 1 - i][ny - 1 - j]
            if a is None or b is None:
                continue
            pairs += 1
            d = abs(a - b)
            diffs.append(d)
            if d > 1.0:
                big += 1
    mx = max(diffs) if diffs else float("nan")
    frac = 100.0 * big / max(pairs, 1)
    print(f"== {path}")
    print(f"   网格 {nx} x {ny}，格点 {len(rows)}")
    print(f"   全局峰 u = {u:.3f} (~ 单位同网格；s/deg 或 s/km 看导出脚本)，"
          f"BAZ = {baz:.1f} deg，B = {bmax:.2f} dB")
    print(f"   反极点 (BAZ+180 = {(baz + 180) % 360:.1f} deg) B = {anti:.2f} dB"
          if anti is not None else "   反极点不在网格内")
    if anti is not None:
        print(f"   反极点相对峰落差 = {bmax - anti:.2f} dB")
    print(f"   点对称 max|B(s)-B(-s)| = {mx:.3f} dB；差 >1 dB 的格点 {frac:.1f}%")
    if mx < 0.01:
        print("   -> 该图为点对称（罕见）：反方向的峰可能真是等高对称影，需另找依据判源。")
    elif anti is not None and (bmax - anti) > 10:
        print("   -> 无镜像峰：反方向的峰按 ①原点针尖 ②混叠复制像 u±111.195/(f*d_min) ③真实次级震相 三条排查。")
    else:
        print("   -> 非对称、但反极点不算低：先查是否落在混叠复制像位置（换频点看峰是否随 u±111.195/(f*d_min) 移动）。")


if __name__ == "__main__":
    args = sys.argv[1:]
    if not args:
        print(__doc__)
        sys.exit(1)
    for a in args:
        report(a)
