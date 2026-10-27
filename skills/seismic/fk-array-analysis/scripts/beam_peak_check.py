#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
beam_peak_check.py —— 波束文件峰位体检（报数前跑一遍）：内部极大？贴边？源向锥内有没有独立的峰？

用途：F-K/波束慢度图定完 SMAX、或"换个网格复核"时，逐文件核"这个读数能不能用"。
判据出处：fk-parameter-selection.md §4b-1（换网格复核）与 §9 判读第 1 条（峰必须在网格内部）。

用法：
  python beam_peak_check.py <beam 文件或目录>... [--baz 137.3] [--cone 15] [--umin 20] [--umax 45]
  目录 = 该目录下全部 *.txt 逐文件一行；同一频点在两个 SMAX 下各跑一份，比 Δu / ΔBAZ / 内部性 / 离边距离。

输出列：u_peak(s/°) @ BAZ(°) | prom(dB) | 内部极大 | dx,dy（离 x/y 边距离）| 锥内最强及其相对全图峰(dB) | 判定
判定规则：非内部极大（含贴边）或 |u| > 内切圆 SMAX ⇒ 弃/换大网格复核；离边 < 15%×SMAX ⇒ 近边，也要复核。

⚠️ 网格约定：方格 ±SMAX（角点可达 √2·SMAX，角点区与贴边读数一律不可用）；
   x = E–W、y = N–S，单位 s/°；y 轴按 unique(...)[::-1] 重建（与 plot_landtime_map.py 同一读法）。
"""
import argparse
import glob
import math
import os
import re
import sys

import numpy as np

WARN_FRAC = 0.15          # 离边距离 < 15%×SMAX ⇒ 近边提示


def load(path):
    d = np.loadtxt(path, comments="#")
    x = np.unique(d[:, 0])
    y = np.unique(d[:, 1])[::-1]
    return x, y, d[:, 2].reshape(len(y), len(x))


def check(path, baz0, cone_half, u_lo, u_hi):
    x, y, b = load(path)
    j, i = np.unravel_index(int(np.argmax(b)), b.shape)
    xs, ys = float(x[i]), float(y[j])
    u = math.hypot(xs, ys)
    az = math.degrees(math.atan2(xs, ys)) % 360.0
    smax = float(max(abs(x.min()), abs(x.max())))            # 方格半宽（x、y 同范围）
    nb = b[max(j - 1, 0):j + 2, max(i - 1, 0):i + 2]
    interior = bool(nb.max() == b[j, i] and 0 < i < len(x) - 1 and 0 < j < len(y) - 1)
    dx, dy = smax - abs(xs), smax - abs(ys)
    XX, YY = np.meshgrid(x, y)
    U = np.hypot(XX, YY)
    AZ = np.degrees(np.arctan2(XX, YY)) % 360.0
    cone = (np.abs((AZ - baz0 + 180.0) % 360.0 - 180.0) <= cone_half) & (U >= u_lo) & (U <= u_hi)
    if cone.any():
        j2, i2 = np.unravel_index(int(np.argmax(np.where(cone, b, -1e9))), b.shape)
        cone_u, cone_baz, cone_db = float(U[j2, i2]), float(AZ[j2, i2]), float(b[j2, i2] - b.max())
    else:
        cone_u = cone_baz = cone_db = float("nan")
    m = re.search(r"SMAX=([0-9.]+)", open(path).readline())
    if not interior:
        verdict = "弃/复核: 非内部极大(贴边或边界)"
    elif u > smax + 1e-9:
        verdict = "弃/复核: 超出内切圆"
    elif min(dx, dy) < WARN_FRAC * smax:
        verdict = "近边: 换大网格复核"
    else:
        verdict = "OK"
    return dict(file=os.path.basename(path), smax_hdr=m.group(1) if m else "?", spts=len(x),
                smax=smax, u=u, az=az, prom=float(b.max() - np.median(b)),
                interior=interior, dx=dx, dy=dy,
                cone_u=cone_u, cone_baz=cone_baz, cone_db=cone_db, verdict=verdict)


def main(argv=None):
    ap = argparse.ArgumentParser(description="波束文件峰位体检（内部性/贴边/源向锥）")
    ap.add_argument("paths", nargs="+", help="beam .txt 文件或目录（目录 = 全部 *.txt）")
    ap.add_argument("--baz", type=float, default=137.3, help="理论 BAZ（子阵各自的值；G1 137.3 / G2 144.2）")
    ap.add_argument("--cone", type=float, default=15.0, help="源向锥半角(°)")
    ap.add_argument("--umin", type=float, default=20.0, help="锥内慢度下限(s/°)")
    ap.add_argument("--umax", type=float, default=45.0, help="锥内慢度上限(s/°)")
    a = ap.parse_args(argv)

    files = []
    for p in a.paths:
        if os.path.isdir(p):
            files += sorted(glob.glob(os.path.join(p, "*.txt")))
        else:
            files.append(p)
    if not files:
        print("没有找到 beam 文件", file=sys.stderr)
        return 2

    print(f"# 锥 = BAZ {a.baz:g}°±{a.cone:g}°、慢度 {a.umin:g}–{a.umax:g} s/°；判定见 fk-parameter-selection.md §4b-1")
    print(f"{'file':34s} {'SMAX(hdr/edge)':>16s} {'npts':>5s} {'u_peak':>8s} {'BAZ':>6s} {'prom':>6s} "
          f"{'内部':>4s} {'dx':>6s} {'dy':>6s} {'锥内u':>7s} {'锥内BAZ':>7s} {'锥内ΔdB':>8s}  判定")
    for f in files:
        try:
            r = check(f, a.baz, a.cone, a.umin, a.umax)
        except Exception as e:                              # noqa: BLE001 —— 逐文件容错，便于批量扫目录
            print(f"{os.path.basename(f):34s}  读取失败: {e}")
            continue
        print(f"{r['file']:34s} {r['smax_hdr']:>7s}/{r['smax']:<8.1f} {r['spts']:5d} {r['u']:8.2f} {r['az']:6.1f} "
              f"{r['prom']:6.2f} {('Y' if r['interior'] else 'N'):>4s} {r['dx']:6.2f} {r['dy']:6.2f} "
              f"{r['cone_u']:7.2f} {r['cone_baz']:7.1f} {r['cone_db']:+8.2f}  {r['verdict']}")
    print("\n读法：同一频点在两个 SMAX 下各跑一份 ⇒ 比 u/BAZ 差与内部性；"
          "锥内 ΔdB ≈ 0 说明源方向的峰就是全图峰；锥内最强钉在锥的内边界（≈umin）说明该方向没有独立于原点的峰。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
