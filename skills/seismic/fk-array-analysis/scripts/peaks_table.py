#!/usr/bin/env python3
"""波束多峰读数：列前 N 个局部极大（不只看 argmax），并给出与理论 BAZ 的差。

用法:
    python peaks_table.py <理论BAZ> <beam文件或目录> [前N名，默认5]

输入文件格式（本项目导出脚本与 SEIZMO 源码导出一致，列序固定）:
    # x(s/deg) y(s/deg) beam_dB
    30.0 -30.0 -7.3
    ...

读法约定（判据见 references/fk-parameter-selection.md §9「判读前四条」）:
  1) |s_peak| > SMAX（= 文件里 x 的最大值）的峰标 OUT —— 那是波束朝外上升被网格截断的假值；
  2) dB rel peak 只是相对最大峰；**判有没有信号要拿它跟峰下多少 dB + B(0) 一起看**，
     单频低频点上主瓣半宽可能大到罩住理论慢度（§8）；
  3) |dBAZ| = 实测 BAZ 与理论 BAZ 的角差（0–180°，取模 180；方向反 180° 时也算 0）。
"""
import sys
import math
import glob
import os

import numpy as np
from scipy.ndimage import maximum_filter, label

MAXFILT = 5          # 5x5 邻域取局部极大；网格很密/很疏时按 Δs 调
FLOOR_DB = 15.0      # 只看峰下 15 dB 以内的局部极大（再低的多半是噪声结构）


def load(path):
    d = np.loadtxt(path, comments="#")
    x, y, b = d[:, 0], d[:, 1], d[:, 2]
    xs = np.unique(x)
    ys = np.unique(y)[::-1]                 # ⚠️ y 必须降序还原，否则方位镜像
    if xs.size * ys.size != b.size:
        raise SystemExit(f"{path}: 网格不是规则矩形（或列序变了），无法 reshape")
    return xs, ys, b.reshape(ys.size, xs.size)


def read_peaks(path, baz_true, n=5):
    xs, ys, b = load(path)
    mx = float(b.max())
    med = float(np.median(b))
    smax_grid = float(np.max(np.abs(xs)))

    loc = (b == maximum_filter(b, size=MAXFILT)) & (b > mx - FLOOR_DB)
    lab, nlab = label(loc)
    rows = []
    for k in range(1, nlab + 1):
        m = lab == k
        iy, ix = np.unravel_index(int(np.argmax(np.where(m, b, -np.inf))), b.shape)
        s_deg = float(math.hypot(xs[ix], ys[iy]))
        baz = float(math.degrees(math.atan2(xs[ix], ys[iy])) % 360.0)
        dbaz = abs((baz - baz_true + 180.0) % 360.0 - 180.0)
        rows.append((float(b[iy, ix]), s_deg, baz, dbaz, s_deg > smax_grid * 1.0001))
    rows.sort(key=lambda r: -r[0])
    return mx, med, smax_grid, rows[:n]


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    baz_true = float(sys.argv[1])
    target = sys.argv[2]
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    files = sorted(glob.glob(os.path.join(target, "*.txt"))) if os.path.isdir(target) else [target]
    if not files:
        raise SystemExit(f"没有找到 beam 文件: {target}")

    for f in files:
        try:
            mx, med, smax_grid, rows = read_peaks(f, baz_true, n)
        except SystemExit as e:
            print(e)
            continue
        print(f"\n== {os.path.basename(f)}  (SMAX_grid={smax_grid:.2f} s/deg, "
              f"prom={mx - med:.2f} dB)")
        print(f"{'#':>2} {'dBrel':>7} {'s/deg':>8} {'BAZ':>7} {'|dBAZ|':>7}  flag")
        for i, (v, s_deg, baz, dbaz, out) in enumerate(rows, 1):
            print(f"{i:>2} {v - mx:>7.2f} {s_deg:>8.2f} {baz:>7.1f} {dbaz:>7.1f}  "
                  f"{'OUT(网格截断)' if out else ''}")


if __name__ == "__main__":
    main()
