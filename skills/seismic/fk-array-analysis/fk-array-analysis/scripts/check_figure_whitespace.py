#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""出图后的像素级布局自查（无需视觉模型）。

用法：
    python check_figure_whitespace.py 图1.png [图2.png ...] [--max-gap 2.0] [--min-ink 0.002]

做什么：PNG 转灰阶 → 按行统计"有墨"像素数（与背景中位灰阶相差 > 阈值）→ 报告
        顶部/底部留白，以及**内部空白带**（两块有墨区之间的连续空行）占图高的百分比。
        第一条内部空白带若位于图上部，标为"标题↔内容"（这是要压的）；其余多为
        面板之间的正常间距，仅作参考、不参与判定。

为什么：用户反复要求"表头简化、空白别太大"。最常中招的是 `set_aspect("equal")` 的
        方形子图被放进偏高的子图槽里整体居中 ⇒ 标题与面板之间留下 10%+ 空白。
        改 figsize / `ax.set_anchor("N")` / `subplots_adjust(top=...)` 后复测，
        "标题↔内容"空白带应 ≤ ~2%。

环境：/home/sy/miniforge3/envs/obspy/bin/python（numpy + PIL 齐备）。
退出码：全部通过 0，有图超标 1（便于串在出图脚本后面）。
"""
import argparse
import sys

import numpy as np
from PIL import Image


def ink_profile(path, ink_thresh=18):
    """返回 (每行有墨像素数, 高, 宽, 背景灰阶)。RGBA/P 先合成到白底。"""
    im = Image.open(path)
    if im.mode in ("RGBA", "LA", "P"):
        im = im.convert("RGBA")
        white = Image.new("RGBA", im.size, (255, 255, 255, 255))
        im = Image.alpha_composite(white, im)
    g = np.asarray(im.convert("L")).astype(np.int16)
    bg_level = float(np.median(g))          # 背景 = 中位灰阶（白底/浅灰底都适用）
    ink = (np.abs(g - bg_level) > ink_thresh).sum(axis=1)
    return ink, g.shape[0], g.shape[1], bg_level


def report(path, max_gap=2.0, min_ink_frac=0.002):
    try:
        ink, H, W, bg = ink_profile(path)
    except Exception as exc:                # 文件不存在 / 不是图片 / 被占用
        print(f"\n\u274c {path}: {exc}")
        return False
    row_thresh = max(1, int(min_ink_frac * W))   # 一行至少这么多像素算"有墨"
    rows = np.flatnonzero(ink > row_thresh)
    print(f"\n{path}  ({W}\u00d7{H}px, 背景灰阶 {bg:.0f}, 行阈 {row_thresh}px)")
    if rows.size == 0:
        print("  \u26a0\ufe0f 未检测到墨迹（阈值过松或图为纯色）")
        return False
    print(f"  顶部留白 {100 * rows[0] / H:.1f}%   底部留白 {100 * (H - 1 - rows[-1]) / H:.1f}%")

    # 内部空白带 = 两个有墨行区间之间的连续空行
    gaps, prev = [], rows[0]
    for r in rows[1:]:
        if r - prev > 1:
            gaps.append((prev + 1, r - 1))
        prev = r
    if not gaps:
        print("  内部空白带：无 \u2705")
        return True

    ok = True
    for i, (a, b) in enumerate(gaps):
        pct = 100 * (b - a + 1) / H
        head_gap = (i == 0 and a < 0.4 * H)          # 上部第一条 = 标题↔内容
        if head_gap:
            flag = ("\u26a0\ufe0f 过大 → 压 figsize 高度 / ax.set_anchor('N') / "
                    "subplots_adjust(top=…)") if pct > max_gap else "ok"
            ok = ok and pct <= max_gap
            role = "标题↔内容（判据）"
        else:
            flag, role = "info", "面板/块间（通常正常）"
        print(f"  内部空白带 行 {a}-{b}（{b - a + 1}px）: {pct:.1f}%  [{role}] {flag}")
    return ok


def main():
    ap = argparse.ArgumentParser(description="PNG 布局空白自查（标题↔内容 ≤ --max-gap）")
    ap.add_argument("images", nargs="+")
    ap.add_argument("--max-gap", type=float, default=2.0,
                    help="允许的最大标题↔内容空白带，占图高百分比（默认 2.0）")
    ap.add_argument("--min-ink", type=float, default=0.002,
                    help="判定有墨行的最小像素占比（默认 0.002 = 0.2%%）")
    args = ap.parse_args()
    results = [report(p, args.max_gap, args.min_ink) for p in args.images]
    ok = all(results)
    print("\n结论：" + (f"\u2705 全部通过（标题↔内容空白 ≤ {args.max_gap}%）"
                     if ok else "\u26a0\ufe0f 有图未通过，见上面标记"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
