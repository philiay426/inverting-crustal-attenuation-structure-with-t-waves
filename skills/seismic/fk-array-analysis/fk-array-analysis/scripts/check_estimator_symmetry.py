#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_estimator_symmetry.py —— 判定"这个 F-K 引擎在 ±u 上对不对称"的确定性探针

回答的问题：图上对面那一团是"镜像孪生"，还是另一团独立能量？
结论不能背，必须按估计量定义算：只用互谱实部（协谱）才严格 B(u)=B(−u)；保留虚部（正交谱=台间时延方向）
则方向唯一。本脚本一次给出三样证据：
  ① 合成注入（已知真值 |u|、BAZ）→ 五种估计量的峰与镜像点电平；含"只留互谱实部"对照组，
     它必须给 −0.00 dB / 0.00 dB，否则先怀疑注入或相位符号（见 §5 的坑）。
  ② 真实波束文件（可给多份）→ 全图 B(u)/B(−u) 差（中位/最大）与各局部极大与其镜像点的电平差。
  ③ 注入自检：台对互谱相位 + 2πf·(u0·Δr) ≈ 0（子波必须落在各台窗交集内，否则残差不收敛）。

用法（默认 = XA Zhou2021 G1 陆地计时案例）：
  python check_estimator_symmetry.py
  python check_estimator_symmetry.py --land-dir <dir> --stations SA31,SA32,... --comp T \
      --baz 144.2 --u 30.0 --win 142.5 293.1 --freq 0.0078125 --smax 42.2 --spts 243 \
      --beams <dir1> [<dir2> ...]

只依赖 numpy/obspy/scipy。读数口径：峰 = argmax；镜像点 = 同 |u|、方位 +180°（邻域 ±1 格取最大）。
"""
import argparse
import math

import numpy as np
from obspy import read
from scipy import ndimage

DEFAULT_LAND = "/mnt/e/sharedate/case/zhou_2021_test/1997-04-28_mw66-prince-edward/process_flow/4_fk_analysis/G1/04_landtime_fk/input_landtime"
DEFAULT_BEAMS = "/mnt/e/sharedate/case/zhou_2021_test/1997-04-28_mw66-prince-edward/process_flow/4_fk_analysis/G1/04_landtime_fk/beams"

ap = argparse.ArgumentParser()
ap.add_argument("--land-dir", default=DEFAULT_LAND)
ap.add_argument("--stations", default="SA16,SA17,SA18,SA22,SA23,SA25,SA29")
ap.add_argument("--comp", default="T")
ap.add_argument("--baz", type=float, default=137.3)          # 子阵理论 BAZ
ap.add_argument("--u", type=float, default=30.0)               # 注入真值 |u| (s/°)
ap.add_argument("--win", type=float, nargs=2, default=(120.0, 270.0))   # 陆地计时秒
ap.add_argument("--freq", type=float, default=0.0078125)
ap.add_argument("--smax", type=float, default=42.2)
ap.add_argument("--spts", type=int, default=243)
ap.add_argument("--beams", nargs="*", default=[DEFAULT_BEAMS])  # 真实波束目录或文件
ap.add_argument("--files", default="ts120_270_T,bin78g_T,s30_T")  # 目录模式下要查的文件（前缀）
a = ap.parse_args()
STAS = [s.strip() for s in a.stations.split(",") if s.strip()]

# ---------------------------------------------------------------- 几何
# 相对质心的局部 EN 度坐标（脚本自检用；与波束图的 E-W / N-S 慢度轴对应）
trs = [read(f"{a.land_dir}/{s}_land_{a.comp}.sac")[0] for s in STAS]
dt = trs[0].stats.delta
lat = np.array([t.stats.sac.stla for t in trs]); lon = np.array([t.stats.sac.stlo for t in trs])
la0, lo0 = lat.mean(), lon.mean()
dE = (lon - lo0) * math.cos(math.radians(la0)); dN = lat - la0
t_ax = (np.arange(trs[0].stats.npts) * dt) + float(trs[0].stats.sac.b)
u0x, u0y = a.u * math.sin(math.radians(a.baz)), a.u * math.cos(math.radians(a.baz))
dly = u0x * dE + u0y * dN

# ---------------------------------------------------------------- 注入
m = (t_ax >= a.win[0]) & (t_ax <= a.win[1])
tt = t_ax[m] - t_ax[m][0]
E = np.exp(-2j * np.pi * a.freq * tt)
# 子波要落在各台窗的交集内；半宽按需自缩（出现大残差时先调这里）
lo = a.win[0] - dly.max(); hi = a.win[1] - dly.min()
center = 0.5 * (lo + hi); half = 0.45 * (hi - lo) / 3.0
src = np.exp(-((t_ax - center) / half) ** 2) * np.cos(2 * np.pi * a.freq * (t_ax - center))
spec = np.array([np.interp(t_ax[m] - dly[k], t_ax, src) @ E for k in range(len(STAS))])
spec = spec / np.abs(spec)                                   # 白化（同引擎）

row, col = np.triu_indices(len(STAS), 1)
dEij, dNij = dE[row] - dE[col], dN[row] - dN[col]
dRij = np.stack([dEij, dNij], 1)
Cij = spec[row] * np.conj(spec[col])
r_ij = 2 * np.pi * a.freq * (u0x * dEij + u0y * dNij)
res = np.angle(np.exp(1j * (np.angle(Cij) + r_ij)))          # C_ij ∝ exp(-i2πf·dly) ⇒ 相位 + r 应 ≈ 0
print(f"注入自检：台对相位 + 预测 的最大残差 = {np.abs(res).max():.4f} rad（应 ≈ 0；不收敛先查子波是否落在窗交集 {lo:.1f}–{hi:.1f} s 内）")

G = np.linspace(-a.smax, a.smax, a.spts)
GX, GY = np.meshgrid(G, G); U = np.stack([GX.ravel(), GY.ravel()], 1)
W = np.exp(2j * np.pi * a.freq * (U @ dRij.T))


def report(beam, label):
    B = beam.reshape(a.spts, a.spts)
    j, i = np.unravel_index(np.argmax(B), B.shape)
    xp, yp = G[i], G[j]
    jm = int(np.argmin(abs(G + yp))); im = int(np.argmin(abs(G + xp)))
    nb = B[max(jm - 2, 0):jm + 3, max(im - 2, 0):im + 3].max()
    with np.errstate(divide="ignore"):
        d = 10 * np.log10(np.maximum(B[::-1, ::-1], 1e-300) / np.maximum(B, 1e-300))
    print(f"  {label:20s} 峰 {math.hypot(xp, yp):5.2f} s/° @ {math.degrees(math.atan2(xp, yp)) % 360:6.2f}° ｜ "
          f"镜像点 {10 * math.log10(nb / B.max()):+7.2f} dB ｜ 全图 B(u)/B(−u) 差 中位 {np.median(np.abs(d)):5.2f}、最大 {np.abs(d).max():6.2f} dB")


print(f"\n【① 合成注入】真值 |u| = {a.u} s/° @ BAZ {a.baz}°（{len(STAS)} 台，窗 {a.win[0]:g}–{a.win[1]:g} s，{a.freq} Hz，网格 ±{a.smax:g}/{a.spts}）")
report(np.abs(np.real(W @ Cij)), "coarray+real(引擎)")
report(np.abs(W @ Cij) ** 2, "coarray+abs2")
c2 = np.concatenate([Cij, np.conj(Cij)]); d2 = np.concatenate([dRij, -dRij])
report(np.abs(np.exp(2j * np.pi * a.freq * (U @ d2.T)) @ c2) ** 2, "全有序台对+abs2")
report(np.abs(np.exp(2j * np.pi * a.freq * (U @ np.stack([dE, dN], 1).T)) @ spec) ** 2, "延时叠加")
report(np.real(W @ np.real(Cij)), "只留互谱实部(对照)")
print("  ⇒ 对照组必须是 −0.00 dB / 0.00 / 0.00 dB（严格对称）；其余各行给出的是该估计量的真实方向性。")

# ---------------------------------------------------------------- 真实波束
print("\n【② 真实数据】B(u)/B(−u) 与各局部极大对其镜像点")
for bd in a.beams:
    import glob, os
    paths = sorted(glob.glob(os.path.join(bd, "*_coarray.txt"))) if os.path.isdir(bd) else [bd]
    for p in paths:
        name = os.path.basename(p).replace("_coarray.txt", "")
        if os.path.isdir(bd) and a.files and name not in a.files:
            continue
        d = np.loadtxt(p, comments="#")
        x = np.unique(d[:, 0]); y = np.unique(d[:, 1])[::-1]; b = d[:, 2].reshape(len(y), len(x))
        with np.errstate(divide="ignore"):
            dd = 10 * np.log10(np.maximum(b[::-1, ::-1], 1e-300) / np.maximum(b, 1e-300))
        mxf = ndimage.maximum_filter(b, size=5)
        lab, n = ndimage.label((b == mxf) & (b > np.median(b) + 2.0))
        pk = []
        for t in range(1, n + 1):
            jy, jx = np.where(lab == t); j, i = jy[0], jx[0]
            u = math.hypot(x[i], y[j]); az = math.degrees(math.atan2(x[i], y[j])) % 360
            jm = int(np.argmin(abs(y + y[j]))); im = int(np.argmin(abs(x + x[i])))
            nb = b[max(jm - 1, 0):jm + 2, max(im - 1, 0):im + 2].max()
            pk.append((b[j, i], u, az, nb))
        pk.sort(reverse=True)
        print(f"  {name}: 全图差 中位 {np.median(np.abs(dd)):5.2f}、最大 {np.abs(dd).max():6.2f} dB")
        for v, u, az, nb in pk[:3]:
            print(f"      极大 {v:7.2f} dB @ {u:5.1f} s/° {az:5.1f}°   其镜像点 {nb - v:+6.2f} dB")
print("\n判读：镜像点电平 ≈ 0 dB 且 |u| 相等 ⇒ 成对（近零慢度族常见）；−4 dB 以下 ⇒ 不对称，不是孪生。\n"
      "     详规见 references/estimator-direction-and-ambiguity.md。")
