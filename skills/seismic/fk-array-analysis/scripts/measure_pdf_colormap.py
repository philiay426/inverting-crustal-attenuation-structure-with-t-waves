#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
measure_pdf_colormap.py —— 从论文插图的色度尺上**像素实测**配色，输出可直接用的控制点。

用法:
    python measure_pdf_colormap.py <paper.pdf> <page> [n_stops] [out_txt]
    python measure_pdf_colormap.py Zhou2021.pdf 7 13          # 默认 13 个控制点

为什么不能凭眼定色板：用户要求慢度/能量图"匹配论文色度尺"，而目测/截图描色会偏（且容易把
连续渐变做成离散分阶，被直接否）。实测做法 = 渲染该页 → 连通域找出"细长且颜色种类极多"的域
（色度尺的特征：h/w 很大、内部颜色数远高于任何数据团）→ 沿中线从下到上取样 → 取低→高控制点。

坑：
  * 阈值别取太松，否则会把彩色数据团也当候选；打印**每个候选**的 h/w、颜色数、bbox
    让人工确认（一篇论文常有多条色度尺，取"细长"且互相对齐的那几条平均更稳）。
  * 取样时**跳过上下各 2%**：色度尺两端有刻度线/边框。
  * 若是论文自己用离散分阶画的色度尺，取样结果会是台阶——此时交付也应该是离散版；
    "不要渐变"指的是**色序过渡别跳变**，不等于要把连续色板改成台阶。
"""

import os
import shutil
import subprocess
import sys

import numpy as np
from PIL import Image
from scipy import ndimage

DPI = 200


def run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    pdf, page = sys.argv[1], int(sys.argv[2])
    n_stops = int(sys.argv[3]) if len(sys.argv) > 3 else 13
    out_txt = sys.argv[4] if len(sys.argv) > 4 else None

    tmp = os.path.join('/tmp', f'cbar_{os.getpid()}')
    os.makedirs(tmp, exist_ok=True)
    try:
        run(['pdftoppm', '-r', str(DPI), '-f', str(page), '-l', str(page), '-png', pdf,
             os.path.join(tmp, 'pg')])
        png = [f for f in os.listdir(tmp) if f.endswith('.png')][0]
        im = Image.open(os.path.join(tmp, png)).convert('RGB')
        a = np.asarray(im).astype(int)

        sat = (a.max(axis=2) - a.min(axis=2) > 40) & (a.max(axis=2) > 60)
        lab, n = ndimage.label(sat)
        cands = []
        for k in range(1, n + 1):
            ys, xs = np.where(lab == k)
            if ys.size < 2000:
                continue
            h = ys.max() - ys.min() + 1
            w = xs.max() - xs.min() + 1
            hw = h / max(w, 1)
            ncol = len(set(map(tuple, a[ys[:3000], xs[:3000]])))
            cands.append((hw, ncol, xs.min(), ys.min(), xs.max(), ys.max(), h, w))
        cands.sort(key=lambda t: (-t[0], -t[1]))
        print(f'候选连通域（h/w 降序）：')
        for hw, ncol, x0, y0, x1, y1, h, w in cands[:10]:
            print(f'  h/w={hw:5.1f}  颜色数≈{ncol:5d}  bbox x[{x0},{x1}] y[{y0},{y1}]  h={h} w={w}')
        bars = [t for t in cands if t[0] >= 5 and t[1] >= 15][:4]
        if not bars:
            print('⚠️ 没找到"细长且颜色多"的域：放宽阈值，或该色度尺是离散分阶/它在另一页。')
            return
        print(f'\n取前 {len(bars)} 条当色度尺（应互相对齐；不对齐就手动指定）：')
        stops = []
        for hw, ncol, x0, y0, x1, y1, h, w in bars:
            xm = (x0 + x1) // 2
            ylo = y0 + int(0.02 * h)
            yhi = y1 - int(0.02 * h)
            ys = np.linspace(yhi, ylo, n_stops).astype(int)       # 下 → 上 = 低 → 高
            col = a[ys, xm]
            stops.append(col)
            print(f'  x={xm} 列：' + ' '.join('#%02x%02x%02x' % tuple(c) for c in col[::max(1, n_stops // 6)]))
        med = np.median(np.stack(stops), axis=0).astype(int)
        print('\n控制点（低 → 高）：')
        print('[' + ', '.join("'#%02x%02x%02x'" % tuple(c) for c in med) + ']')
        print('\n用法：')
        print("  from matplotlib.colors import LinearSegmentedColormap")
        print("  CMAP = LinearSegmentedColormap.from_list('paper', " +
              '[' + ','.join("'#%02x%02x%02x'" % tuple(c) for c in med) + '])')
        if out_txt:
            with open(out_txt, 'w', encoding='utf-8') as fh:
                for c in med:
                    fh.write('#%02x%02x%02x\n' % tuple(c))
            print(f'\n已写入 {out_txt}')
        print('\n交付前核验：出图后统计代码里那几个关键色（红/品红/紫/深端）的像素数都 >0，'
              '并说明"色板取自论文 Fig.X 色度尺实测"。')
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == '__main__':
    main()
