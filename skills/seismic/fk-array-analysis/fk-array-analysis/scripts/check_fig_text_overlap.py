#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""出图脚本的"文字互压"自查：monkey-patch Figure.savefig，每次保存前检查
   ① 任意两个文字/注释框是否相交（bbox.overlaps）  ② 是否有文字出画布。

用法（不修改任何出图脚本，脚本按原样执行）：
    python check_fig_text_overlap.py plot_a.py [plot_b.py ...]
    PDF=1 python check_fig_text_overlap.py plot_a.py        # 脚本自己认这个环境变量时生效

判据：相交面积 > MIN_AREA（默认 4 px²，100 dpi 坐标）才报 —— 相邻刻度标签 / 色标刻度的
      包围盒常相切 1 px²，属噪声；漏报比误报更贵，所以阈值只滤到 4。

为什么必须有这一条：fig_ink_margins.py 只量"墨迹是否贴画布边"，**查不出文字互相压**。
实测（5.8×5.4 in 单面板）：x 轴标注「E-W slowness (s/°)」y[22.0,35.7] 与页脚行 y[12.1,24.3]
相交 2.3 px（≈7 px @300 dpi）—— 纯像素留白检查完全看不出来。修法：页脚压成一行 +
subplots_adjust(bottom) 从 0.115 提到 0.145（0.115 时仍相压）；改完复跑本脚本应为 0 相交。

注意：本脚本只管「文本 vs 文本」，元素出画布/墨迹留白用 scripts/fig_ink_margins.py。
"""
import os
import runpy
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.figure as mf  # noqa: E402
import matplotlib.text as mtext  # noqa: E402

MIN_AREA = 4.0  # px²（100 dpi 坐标）
_orig_savefig = mf.Figure.savefig


def _patched(self, fname, *a, **k):
    try:
        self.canvas.draw()
        r = self.canvas.get_renderer()
        boxes = []
        for t in self.findobj(mtext.Text):
            s = t.get_text().strip()
            if not s:
                continue
            try:
                bb = t.get_window_extent(r)
            except Exception:
                continue
            boxes.append((s[:44].replace("\n", " / "), bb))
        hits = []
        for i in range(len(boxes)):
            for j in range(i + 1, len(boxes)):
                bi, bj = boxes[i][1], boxes[j][1]
                if not bi.overlaps(bj):
                    continue
                ov = (min(bi.x1, bj.x1) - max(bi.x0, bj.x0)) * (
                    min(bi.y1, bj.y1) - max(bi.y0, bj.y0)
                )
                if ov > MIN_AREA:
                    hits.append((ov, boxes[i][0], boxes[j][0]))
        W, H = self.canvas.get_width_height()
        out = [s for s, bb in boxes if bb.x0 < 0 or bb.y0 < 0 or bb.x1 > W or bb.y1 > H]
        hits.sort(reverse=True)
        print(
            f"[TEXT-CHECK] {os.path.basename(str(fname))}  文字 {len(boxes)} 条 | "
            f"相交 {len(hits)} | 出框 {len(out)}"
        )
        for ov, x, y in hits[:8]:
            print(f"    x 相交 {ov:7.0f} px² : 「{x}」x「{y}」")
        for s in out[:5]:
            print(f"    x 出框:「{s}」")
    except Exception as e:  # 检查本身不该阻断出图
        print(f"[TEXT-CHECK] {fname}: 检查失败 {e!r}")
    return _orig_savefig(self, fname, *a, **k)


mf.Figure.savefig = _patched


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    for scr in sys.argv[1:]:
        print(f"=== {os.path.basename(scr)} ===")
        runpy.run_path(scr, run_name="__main__")
    return 0


if __name__ == "__main__":
    sys.exit(main())
