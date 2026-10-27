#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""步 3 F-K 出图与判读（一个子阵一套；FK_TAG 区分多个子阵）

对应 skill：seismic/twave-crustal-attenuation  §5 步 3

流程：先跑 step3_fk_driver.m（Octave + SEIZMO fkmap 算波束，写 beam_<TAG>.txt）→ 本脚本读网格出图。

做的事：
  1) 读 beam_<TAG>.txt（x=E/W 慢度, y=N/S 慢度, beam dB）
  2) 换极坐标：|u| = hypot(x,y)，BAZ = atan2(x,y) % 360
  3) ★ 列出**全部**局部极大（> 峰 - 6 dB），不只 argmax —— 这是判读要求
  4) 画慢度图（等高线 + 峰标）+ 送进 fkmap 的波形（时窗落在哪些台上）
  5) 给出分辨率地板 Δs/u >= 1.28·d_min/D 与 BAZ 误差 Δu/|u|（算完才知道能不能报小数）

产物：<OUT>/fk_map_<TAG>.png、<OUT>/fk_peaks_<TAG>.csv

用法：
  T_BASE=<事件目录> FK_TAG=<子阵标签> FK_STAS=A,B,C /path/to/python step3_fk_demo.py
"""
import os
import csv
import numpy as np
from scipy.ndimage import maximum_filter
from obspy.geodetics import gps2dist_azimuth
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ==================== ★ 参数区 ====================
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.environ.get("T_OUTDIR") or HERE                    # 本步目录本身
SRC = os.environ.get("T_BASE", "").strip()                 # 事件根（必填）
if not SRC:
    raise SystemExit("请用环境变量 T_BASE 指定事件根目录")
STACSV = (os.environ.get("T_CSV")
          or os.path.join(SRC, "process_flow", "1_station_distribution", "stations.csv"))
EV = (float(os.environ.get("EVENT_LAT", "0.0")),
      float(os.environ.get("EVENT_LON", "0.0")))           # 震中（lat, lon）
# 子阵台站名单由 FK_STAS 传入（逗号分隔）；FK_TAG 只决定输入/输出文件名前缀
SUB = os.environ.get("FK_TAG", "G1")
STAS = [x.strip() for x in os.environ.get("FK_STAS", "").split(",") if x.strip()]
if not STAS:
    raise SystemExit("请用环境变量 FK_STAS=A,B,C 指定该子阵的台站")
PK_DROP_DB = 6.0                  # 只报"峰 - 6 dB"以内的局部极大
FK_WIN = os.environ.get("FK_WIN", "1200-2200")   # 仅用于图标题，须与驱动脚本 TWIN 一致
DPI = 300
# ===================================================


def read_beam(path):
    """读 SEIZMO 导出的网格。★ y 轴方向：文件里 N/S 慢度是**降序**写的，
    直接 reshape 会把图上下翻过来（峰 |u| 对、BAZ 错 180-差），所以按坐标建索引再填。"""
    rows = [tuple(float(v) for v in l.split()) for l in open(path) if not l.startswith("#")]
    a = np.array(rows)
    xs = np.unique(a[:, 0])
    ys = np.unique(a[:, 1])                     # 升序
    xi = {v: i for i, v in enumerate(xs)}
    yi = {v: i for i, v in enumerate(ys)}
    B = np.full((len(ys), len(xs)), np.nan)
    for xv, yv, bv in a:
        B[yi[yv], xi[xv]] = bv
    return xs, ys, B


def read_station_meta():
    meta = {}
    with open(STACSV) as f:
        for r in csv.DictReader(l for l in f if not l.startswith("#")):
            meta[r["sta"].strip()] = (float(r["lat"]), float(r["lon"]))
    return meta


def read_fkinput(path):
    """返回 (names, data[nsta, npts], delta, t_start)"""
    lines = [l for l in open(path) if not l.startswith("#")]
    names, data = [], []
    for l in lines:
        p = l.split()
        names.append(p[0])
        data.append(np.array([float(v) for v in p[1:]]))
    # 头里有两行 "# nsta npts delta t_start"（说明）与 "# 4 20001 0.05 1200"（数值）
    vals = None
    for l in open(path):
        if not l.startswith("#"):
            continue
        p = l.lstrip("#").split()
        if len(p) == 4 and all(v.replace(".", "", 1).replace("-", "", 1).isdigit() for v in p):
            vals = p
    if vals is None:
        raise SystemExit("fkinput 头里没找到 nsta npts delta t_start")
    _, npts, delta, t0 = vals
    return names, np.array(data), float(delta), float(t0)


def main():
    xs, ys, B = read_beam(os.path.join(OUT, f"beam_{SUB}.txt"))
    Bn = B - B.max()                                   # dB，峰归一到 0
    U = np.hypot(*np.meshgrid(xs, ys))                 # |u| 网格
    BAZ = np.degrees(np.arctan2(*np.meshgrid(xs, ys))) % 360

    # 局部极大（8 邻域），只留 峰-6 dB 以内
    mx = maximum_filter(Bn, size=9)
    cand = np.argwhere((Bn == mx) & (Bn >= -PK_DROP_DB))
    peaks = sorted([(Bn[i, j], U[i, j], BAZ[i, j], i, j) for i, j in cand], reverse=True)
    # 合并离得极近的（<2 s/deg）
    merged = []
    for b, u, baz, i, j in peaks:
        if all(abs(u - m[1]) > 2.0 or abs((baz - m[2] + 180) % 360 - 180) > 8 for m in merged):
            merged.append((b, u, baz, i, j))

    print(f"=== 步3 F-K | {SUB} {len(STAS)} 台子阵 ===")
    print(f"网格: {len(xs)}x{len(ys)}, |u| 最大 {U.max():.1f} s/deg；峰归一 0 dB")
    print(f"局部极大（峰 -{PK_DROP_DB:g} dB 内）:")
    for k, (b, u, baz, i, j) in enumerate(merged[:6], 1):
        print(f"  [{k}] {b:+6.2f} dB   |u| = {u:5.2f} s/deg   BAZ = {baz:6.1f} deg")

    # 分辨率地板与 BAZ 误差（几何决定，与频带无关）
    meta = read_station_meta()
    lat = np.array([meta[s][0] for s in STAS]); lon = np.array([meta[s][1] for s in STAS])
    dmat = np.zeros((len(STAS), len(STAS)))
    for i in range(len(STAS)):
        for j in range(len(STAS)):
            dmat[i, j] = gps2dist_azimuth(lat[i], lon[i], lat[j], lon[j])[0] / 1000.0
    dmin = dmat[dmat > 0].min()
    # 孔径：台站两两最大距离
    D = dmat.max()
    floor = 1.28 * dmin / D
    baz_mean = np.mean([gps2dist_azimuth(EV[0], EV[1], a, b)[2] for a, b in zip(lat, lon)])
    print(f"\n几何: d_min = {dmin:.0f} km, 孔径 D = {D:.0f} km")
    C = 111.195 / (2 * dmin)                       # 抗混叠常数 s/deg·Hz
    f_noalias = C / SMAX_                          # 带内上界必须 <= 该值
    print(f"抗混叠 C = 111.195/(2*d_min) = {C:.4f} s/deg·Hz -> 无混叠要求 f <= C/SMAX = {f_noalias:.4f} Hz")
    print(f"本次频带上界 1/25 = 0.0400 Hz {'OUT OF RANGE (> {:.4f})'.format(f_noalias) if 0.04 > f_noalias else '（在界内）'}"
          f"  -> 越界部分靠台数压制旁瓣，台数越少越危险（台太少会整体跑偏）")
    print(f"预测慢度 u_pred = 111.195/V_T = {SMAX_/1.4:.2f} s/deg（smax = 1.4*u_pred）")
    print(f"分辨率地板 Δs/u >= 1.28*d_min/D = {floor*100:.0f}%   （4 台比 7-9 台更差，峰更宽）")
    print(f"BAZ 误差 ~ Δu/|u| = {floor*100:.0f}% -> ±{floor*360/ (2*np.pi)*0 if False else floor*57.3:.0f} deg 量级（不能报小数）")
    print(f"理论 BAZ（{len(STAS)} 台均值，gps2dist_azimuth）= {baz_mean:.1f} deg")

    with open(os.path.join(OUT, f"fk_peaks_{SUB}.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([f"# {SUB} {len(STAS)} 台子阵 F-K 局部极大；dB 相对峰；BAZ = atan2(u_E, u_N) % 360"])
        w.writerow(["rank", "dB_rel_peak", "u_s_per_deg", "baz_deg"])
        for k, (b, u, baz, i, j) in enumerate(merged[:6], 1):
            w.writerow([k, f"{b:.2f}", f"{u:.2f}", f"{baz:.1f}"])

    # ---------------- 图 ----------------
    fig = plt.figure(figsize=(12, 5.2))
    ax = fig.add_subplot(1, 2, 1)
    pc = ax.pcolormesh(xs, ys, Bn, shading="auto", cmap="Spectral_r", vmin=-12, vmax=0)
    ax.contour(xs, ys, Bn, levels=[-6, -3, -2, -1], colors="k", linewidths=0.5)
    b0, u0, baz0, i0, j0 = merged[0]
    ax.plot(xs[j0], ys[i0], "k*", ms=13, mfc="none", mew=1.2)
    ax.annotate(f"peak {u0:.1f} s/deg, {baz0:.0f} deg", (xs[j0], ys[i0]),
                textcoords="offset points", xytext=(8, 6), fontsize=8)
    for k, (b, u, baz, i, j) in enumerate(merged[1:4], 2):
        ax.plot(xs[j], ys[i], "o", ms=4, mfc="none", mec="k", mew=0.8)
        ax.annotate(f"{u:.1f} s/deg, {baz:.0f} deg", (xs[j], ys[i]), textcoords="offset points",
                    xytext=(5, -10), fontsize=7)
    # 理论 BAZ 方向射线
    for uu in (10, 20, 30):
        ax.plot([0, uu * np.sin(np.radians(baz_mean))], [0, uu * np.cos(np.radians(baz_mean))],
                color="0.35", lw=0.7, ls="--")
    ax.text(0.2 * np.sin(np.radians(baz_mean)), 28, f"theory BAZ {baz_mean:.0f} deg",
            fontsize=7, color="0.3", rotation=-15)
    ax.set_xlim(-SMAX_, SMAX_)
    ax.set_ylim(-SMAX_, SMAX_)
    ax.set_aspect("equal")
    ax.set_xlabel("E-W slowness (s/deg)")
    ax.set_ylabel("N-S slowness (s/deg)")
    ax.set_title(f"F-K beam map, {SUB} {len(STAS)} stations, 2-8 Hz envelope, {FK_WIN} s", fontsize=9)
    fig.colorbar(pc, ax=ax, label="beam (dB re peak)", pad=0.02)

    # 右：送进 fkmap 的波形（时窗）
    names, data, delta, t0 = read_fkinput(os.path.join(OUT, f"fkinput_{SUB}.txt"))
    ax2 = fig.add_subplot(1, 2, 2)
    t = t0 + np.arange(data.shape[1]) * delta
    for i, nm in enumerate(names):
        y = data[i] / (np.abs(data[i]).max() or 1) * 0.4 + i
        ax2.plot(t, y, "k", lw=0.6)
        ax2.text(t[0] + 5, i + 0.45, nm, fontsize=7)
    ax2.axvline(t0, color="#d62728", lw=0.9, ls="--")
    ax2.axvline(t[-1], color="#d62728", lw=0.9, ls="--")
    ax2.set_yticks([])
    ax2.set_xlabel("Time relative to origin (s)")
    ax2.set_title("fkmap input: envelope after getspline / removemean / taper", fontsize=9)
    ax2.grid(axis="x", alpha=0.2, lw=0.4)

    fig.tight_layout()
    p = os.path.join(OUT, f"fk_map_{SUB}.png")
    fig.savefig(p, dpi=DPI)
    plt.close(fig)
    print(f"\n写: {p}")
    print(f"写: {os.path.join(OUT, f'fk_peaks_{SUB}.csv')}")
    print("步3 完成。")


# 网格半径：FK_SMAX 直接给，否则由预测慢度算（u_pred = 111.195/FK_VT）
_smax_env = os.environ.get("FK_SMAX", "").strip()
if _smax_env:
    SMAX_ = float(_smax_env)
else:
    _vt = os.environ.get("FK_VT", "").strip()
    if not _vt:
        raise SystemExit("缺 FK_VT（T 波视速度 km/s）或 FK_SMAX（网格半径 s/deg）；先跑 fk_geometry.py")
    SMAX_ = 1.4 * 111.195 / float(_vt)

if __name__ == "__main__":
    main()
