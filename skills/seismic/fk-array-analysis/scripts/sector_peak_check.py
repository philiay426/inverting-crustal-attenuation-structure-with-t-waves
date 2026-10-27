#!/usr/bin/env python3
"""读一张慢度图 txt（三列：x(s/deg) y(s/deg) beam_dB），报四件套：
   全图峰 / 理论 BAZ 扇区内最强局部极大 / 全图中位 / 指定理论慢度点取值。

用法：
    python sector_peak_check.py <beam.txt> <理论BAZ> [扇区±deg=20] [最小|u| s/deg=5] [理论慢度,逗号分隔]

口径（与 references/peak-reading-and-method-comparison.md §5 一致）：
  · 局部极大 = 8 邻域最大值（等价 maximum_filter(b,5) 的去重前一步），排除 |u| < 最小值的原点针尖
  · "相对全图峰 dB" = 该值 − max(beam)；"相对中位 dB" = 该值 − median(beam)
  · 理论点坐标 = (u·sin(BAZ), u·cos(BAZ))，取最近网格点；s=0 参考 = 离原点最近的格点
  · 判据：慢度在理论量级 **且** |dBAZ| <= 15 deg 才算一个相；数值接近不等于看对相位
"""
import sys

import numpy as np


def load(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            p = line.split()
            if len(p) < 3:
                continue
            rows.append((float(p[0]), float(p[1]), float(p[2])))
    arr = np.asarray(rows, dtype=float)
    if arr.size == 0:
        raise SystemExit('no data rows in %s' % path)
    xs = np.unique(arr[:, 0])
    if arr.shape[0] % xs.size:
        raise SystemExit('not a full rectangular grid: %d cols, %d rows'
                         % (xs.size, arr.shape[0]))
    # ⚠ y 必须按文件里出现的顺序取: Seizmo fkmap 导出的 y 是降序(+smax -> -smax)。
    #   用 np.unique 的升序会把慢度图上下翻转、方位被镜像成 180-az
    #   (实测马达加斯加 ALL25 真值 33.14 s/deg @141.64 deg 被误读成 @38.36 deg)。
    ys = arr[::xs.size, 1]
    if ys.size * xs.size != arr.shape[0]:
        raise SystemExit('not a full rectangular grid: %d x %d != %d rows'
                         % (xs.size, ys.size, arr.shape[0]))
    return xs, ys, arr[:, 2].reshape(ys.size, xs.size)


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    path = sys.argv[1]
    baz0 = float(sys.argv[2])
    tol = float(sys.argv[3]) if len(sys.argv) > 3 else 20.0
    umin = float(sys.argv[4]) if len(sys.argv) > 4 else 5.0
    theos = [float(v) for v in sys.argv[5].replace(',', ' ').split()] if len(sys.argv) > 5 else []

    xs, ys, beam = load(path)
    vmax = beam.max()
    med = np.median(beam)
    i, j = np.unravel_index(np.argmax(beam), beam.shape)
    print('file        : %s' % path)
    print('grid        : %d x %d, |x|max %.2f |y|max %.2f s/deg'
          % (xs.size, ys.size, abs(xs).max(), abs(ys).max()))
    print('global peak : %.2f s/deg @ %.1f deg   (%.2f dB, median %.2f dB)'
          % (np.hypot(xs[j], ys[i]), np.degrees(np.arctan2(xs[j], ys[i])) % 360.0, vmax, med))

    best = None
    for a in range(1, beam.shape[0] - 1):
        for b in range(1, beam.shape[1] - 1):
            v = beam[a, b]
            if v < beam[a - 1:a + 2, b - 1:b + 2].max():
                continue
            u = np.hypot(xs[b], ys[a])
            if u < umin:
                continue
            az = np.degrees(np.arctan2(xs[b], ys[a])) % 360.0
            if abs((az - baz0 + 180.0) % 360.0 - 180.0) <= tol and (best is None or v > best[0]):
                best = (v, u, az)
    if best is None:
        print('sector %6.1f+/-%-4.0f : no local max with |u| >= %.1f s/deg' % (baz0, tol, umin))
    else:
        print('sector %6.1f+/-%-4.0f : %.2f s/deg @ %.1f deg   (%.2f dB rel. peak, %.2f dB rel. median)'
              % (baz0, tol, best[1], best[2], best[0] - vmax, best[0] - med))

    i0 = int(np.argmin(abs(ys)))
    j0 = int(np.argmin(abs(xs)))
    for u in theos:
        for az, tag in ((baz0, 'theo'), ((baz0 + 180.0) % 360.0, 'anti')):
            jj = int(np.argmin(abs(xs - u * np.sin(np.radians(az)))))
            ii = int(np.argmin(abs(ys - u * np.cos(np.radians(az)))))
            print('theo %6.2f s/deg %s @%6.1f : %7.2f dB rel. peak, %7.2f dB rel. median, %7.2f dB rel. s=0'
                  % (u, tag, az, beam[ii, jj] - vmax, beam[ii, jj] - med, beam[ii, jj] - beam[i0, j0]))


if __name__ == '__main__':
    main()
