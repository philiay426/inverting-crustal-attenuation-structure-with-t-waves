#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""图件出框自查：量 PNG 的墨迹范围与四边留白（免视觉模型）。

用法：
    python fig_ink_margins.py fig_a.png fig_b.png ...
    python fig_ink_margins.py <目录>          # 该目录下所有 png

输出：canvas 尺寸 / 墨迹 bbox / 四边留白(px)；任一边 <3 px 标 ⚠️。
判读：左右都贴边 ⇒ 横向溢出（页脚/标题行太长，拆行或加宽画布）；
      下边贴边 ⇒ 页脚 y 太低，抬高 y 并加大 subplots_adjust(bottom=...)。
⚠️ 只对 matplotlib 固定画布有效：pygmt/GMT 的 savefig 按墨迹自动裁边（psconvert -A），
   四边必然报 0 留白 —— 属假阳性。pygmt 图改用「像素尺寸 / DPI / 区域宽高比」核对，
   文字互压也换不出 matplotlib 的 text bbox，需靠布局时按估算框逐个避让。
阈值 250：白底图上 <250 的灰度算墨迹（浅灰填充/淡色包络也会被算进去）。
"""
import glob
import os
import sys

import numpy as np
from PIL import Image


def check(path, thresh=250):
    a = np.asarray(Image.open(path).convert("L"))
    ink = a < thresh
    if not ink.any():
        print(f"{os.path.basename(path)}: 全白？检查阈值/路径")
        return False
    ys, xs = np.where(ink)
    h, w = a.shape
    left, right = int(xs.min()), int(w - 1 - xs.max())
    top, bottom = int(ys.min()), int(h - 1 - ys.max())
    tight = [n for n, v in (("左", left), ("右", right), ("上", top), ("下", bottom)) if v < 3]
    ok = not tight
    print(f"{os.path.basename(path)}: canvas {w}x{h} | ink x[{xs.min()},{xs.max()}] y[{ys.min()},{ys.max()}] "
          f"| 留白 左{left} 右{right} 上{top} 下{bottom} px"
          + ("  ✓" if ok else f"  ⚠️ {'/'.join(tight)}贴边"))
    return ok


if __name__ == "__main__":
    args = sys.argv[1:]
    if not args:
        print(__doc__)
        sys.exit(1)
    paths = []
    for a in args:
        paths += sorted(glob.glob(os.path.join(a, "*.png"))) if os.path.isdir(a) else [a]
    bad = [p for p in paths if not check(p)]
    sys.exit(0 if not bad else 2)
