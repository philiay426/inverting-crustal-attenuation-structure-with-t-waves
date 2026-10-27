#!/usr/bin/env python3
"""步 3 前置：由台站几何现算 F-K 参数（预测慢度/网格/频带/分辨率地板）。

目的：换地区、换台组时**不需要人脑猜 smax**——用台站坐标与 T 波视速度算出来。
所有输入走环境变量，缺必填直接 SystemExit（不给案例默认值）。

环境变量
  FK_STAS   必填，逗号分隔台站名单（如 "AA,BB,CC"）
  FK_VT     必填，T 波视速度 km/s（由该区频散/文献/包络峰值组速度定）
  坐标来源（二选一）：
    T_CSV   台站表 csv（需 sta/lat/lon 列，可带 tt_tp/tt_ts）
    FK_SAC  SAC 目录（读 stla/stlo；文件名需含台名）
  FK_SMAX   可选，直接给网格半径 s/deg（给了就不按 1.4*u_pred）
  FK_SPTS   可选，网格点数，默认 301（奇数，0 在格点）
  FK_FMAX   可选，频带上界上限 Hz，默认 0.04
  EVENT_LAT / EVENT_LON  可选，给了就算理论 BAZ 与震中距
  T_OUTDIR  可选，写 fk_geometry.csv 的目录；不给只打印

输出：控制台表格 + 可选 fk_geometry.csv，含
  n_sta, d_min_km, aperture_km, C_sdeg_hz, u_pred_sdeg, smax_sdeg, spts,
  floor_pct, baz_deg, f_hi_hz, win_min_s
"""
import csv
import math
import os
import sys

import numpy as np
from obspy.geodetics import gps2dist_azimuth


def need(k):
    v = os.environ.get(k, "").strip()
    if not v:
        raise SystemExit(f"缺环境变量 {k}（本脚本无默认值，见文件头说明）")
    return v


def read_coords(stas):
    """按 FK_STAS 顺序取 (lat, lon)：优先台站表，其次 SAC 头。"""
    csv_path = os.environ.get("T_CSV", "").strip()
    if csv_path and os.path.exists(csv_path):
        rows = {}
        with open(csv_path, encoding="utf-8-sig") as f:
            for r in csv.DictReader(x for x in f if not x.startswith("#")):
                if r.get("sta"):
                    rows[r["sta"].strip()] = (float(r["lat"]), float(r["lon"]))
        miss = [s for s in stas if s not in rows]
        if miss:
            raise SystemExit(f"台站表里没有: {','.join(miss)}（T_CSV={csv_path}）")
        return [rows[s] for s in stas], csv_path

    sac_dir = os.environ.get("FK_SAC", "").strip()
    if sac_dir and os.path.isdir(sac_dir):
        from obspy import read
        coords = []
        for s in stas:
            hits = [os.path.join(sac_dir, f) for f in os.listdir(sac_dir) if s in f]
            if not hits:
                raise SystemExit(f"{sac_dir} 里找不到台站 {s} 的文件")
            h = read(sorted(hits)[0])[0].stats.sac
            coords.append((float(h.stla), float(h.stlo)))
        return coords, sac_dir

    raise SystemExit("缺坐标来源：设 T_CSV（台站表）或 FK_SAC（SAC 目录）")


def main():
    stas = [s.strip() for s in need("FK_STAS").split(",") if s.strip()]
    vt = float(need("FK_VT"))
    spts = int(os.environ.get("FK_SPTS", "301") or 301)
    fmax = float(os.environ.get("FK_FMAX", "0.04") or 0.04)
    if spts % 2 == 0:
        spts += 1  # 0 必须落在格点

    coords, src = read_coords(stas)
    lat = np.array([c[0] for c in coords])
    lon = np.array([c[1] for c in coords])

    n = len(stas)
    dmat = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            dmat[i, j] = gps2dist_azimuth(lat[i], lon[i], lat[j], lon[j])[0] / 1000.0
    d_min = float(dmat[dmat > 0].min())
    aperture = float(dmat.max())

    # ---- 核心三个量：预测慢度、网格半径、抗混叠上界 ----
    u_pred = 111.195 / vt                      # s/deg，由视速度直接换算
    smax_env = os.environ.get("FK_SMAX", "").strip()
    smax = float(smax_env) if smax_env else 1.4 * u_pred
    C = 111.195 / (2.0 * d_min)                # s/deg/Hz，台阵响应常数
    f_hi = min(fmax, C / smax)                 # 超过就混叠，峰会塌到零慢度
    floor_pct = 100.0 * 1.28 * d_min / aperture
    baz_deg = float("nan")
    if os.environ.get("EVENT_LAT", "").strip() and os.environ.get("EVENT_LON", "").strip():
        elat = float(os.environ["EVENT_LAT"])
        elon = float(os.environ["EVENT_LON"])
        bazs = [gps2dist_azimuth(elat, elon, a, b)[2] for a, b in zip(lat, lon)]
        baz_deg = float(np.mean(bazs))

    win_min = 1.0 / f_hi if f_hi > 0 else float("nan")   # 窗长下限（频点能分辨所需）
    print(f"台站来源: {src}    台数 n = {n}")
    print(f"{'sta':>6} {'lat':>8} {'lon':>8}")
    for s, (a, b) in zip(stas, coords):
        print(f"{s:>6} {a:8.4f} {b:8.4f}")
    print("-" * 46)
    print(f"d_min      = {d_min:8.2f} km      孔径 D = {aperture:8.2f} km")
    print(f"C          = {C:8.4f} s/deg/Hz   (= 111.195/(2*d_min))")
    print(f"u_pred     = {u_pred:8.3f} s/deg   (= 111.195/V_T, V_T = {vt:g} km/s)")
    print(f"smax       = {smax:8.2f} s/deg    (1.4*u_pred 或 FK_SMAX)")
    print(f"spts       = {spts:8d}")
    print(f"f_hi       = {f_hi:8.5f} Hz      (= min(FK_FMAX, C/smax))")
    print(f"窗长下限   = {win_min:8.1f} s       (1/f_hi)")
    print(f"分辨率地板 = {floor_pct:8.1f} %       (1.28*d_min/D) -> 慢度只能报到 10^0~10^1 s/deg 量级")
    if not math.isnan(baz_deg):
        print(f"理论 BAZ   = {baz_deg:8.2f} deg     (台组中心 -> 震中; 逐台散布可再算 std)")
    print("-" * 46)
    print("把 smax / f_hi 抄进 FK_SMAX / FK_FRNG，或在 Octave 里给 FK_C 让它自己算上界。")

    out = os.environ.get("T_OUTDIR", "").strip()
    if out:
        os.makedirs(out, exist_ok=True)
        p = os.path.join(out, "fk_geometry.csv")
        with open(p, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["n_sta", "d_min_km", "aperture_km", "C_sdeg_hz", "u_pred_sdeg",
                        "smax_sdeg", "spts", "floor_pct", "baz_deg", "f_hi_hz", "win_min_s"])
            w.writerow([n, f"{d_min:.2f}", f"{aperture:.2f}", f"{C:.4f}", f"{u_pred:.3f}",
                        f"{smax:.2f}", spts, f"{floor_pct:.1f}",
                        "" if math.isnan(baz_deg) else f"{baz_deg:.2f}",
                        f"{f_hi:.5f}", f"{win_min:.1f}"])
        print(f"写: {p}")


if __name__ == "__main__":
    main()
