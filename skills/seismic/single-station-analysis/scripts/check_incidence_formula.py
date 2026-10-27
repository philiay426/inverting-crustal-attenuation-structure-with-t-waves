#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_incidence_formula.py —— 入射角公式链现跑自检：合成真值 + 真实窗各算一遍

用途：用户贴一张公式图问"这个对吗 / 这个是 arctan 吗 / 这样吗"时，现跑一遍给数字背书，
不要只讲道理（公式链与读式陷阱见 references/incidence-angle-from-polarization.md §2/§4）。

自检项
  A 合成：已知 i、已知方位 az 的直线偏振波包 → 两个估计器都应原样复原真值（±0.05°）
  B 合成反例：Z–R 相差 90° 的椭圆偏振 → PCA 退化（报 90°），i_amp 仍给对角
  C 真实窗：去响应三分量 → 绕真实大圆 BAZ 转 R/T → 2–8 Hz 带通 → 窗内两估计器，
            结果应等于已存档值（默认 SA18 ±5 s = 57.459°）

用法
  python check_incidence_formula.py                              # 只跑 A/B（不需要数据）
  python check_incidence_formula.py --base <事件目录> --sta SA18 --pick 1195.15

口径（C 部分）
  时间轴 t = (b − o) + i/fs（SAC o 字段优先，见 SKILL.md 时间轴铁律）；baz 用 gps2dist_azimuth 的
  台站端真实大圆值，不用 stations.csv 的 baz 列（震源端口径，偏西约 11°）。
"""
import argparse
import csv
import glob
import math
import numpy as np
from scipy.signal import butter, filtfilt

EV = {'lat': -42.4134, 'lon': 42.6807}   # 默认事件（Zhou2021 1997-04-28 Mw6.6），--ev-lat/--ev-lon 可覆盖
BP = (2.0, 8.0)


def estimators(z, r, t):
    """返回 (i_PCA, i_amp, 直线度)。i 自垂直方向量，单位度。"""
    w, v = np.linalg.eigh(np.cov(np.c_[z, r, t].T))
    p = v[:, -1]                                        # 最大特征值对应的主轴（偏振方向）
    i_pca = math.degrees(math.acos(min(1.0, abs(p[0]))))  # arccos(|p_Z|/|p|)，|p|=1
    i_amp = math.degrees(math.atan2(math.hypot(np.std(r), np.std(t)), np.std(z)))
    return i_pca, i_amp, 1 - w[-2] / w[-1]


def synth(i_true, az, n=3000, noise=2e-9):
    """已知入射角 i_true 与方位 az 的直线偏振波包（窄高斯包络 × 正弦）。"""
    s = np.exp(-((np.arange(n) - n // 2) / 60.0) ** 2) * np.sin(2 * np.pi * np.arange(n) / 40.0)
    u = np.array([math.cos(math.radians(i_true)),
                  math.sin(math.radians(i_true)) * math.cos(math.radians(az)),
                  math.sin(math.radians(i_true)) * math.sin(math.radians(az))])
    X = np.outer(s, u) * 1e-6 + np.random.default_rng(0).normal(0, noise, (n, 3))
    return X[:, 0], X[:, 1], X[:, 2]


def part_ab():
    print('=== A 合成（直线偏振：两个估计器都应复原真值）===')
    for i_true, az in ((30, 50), (55.4, 20), (75, 0), (5, 0)):
        ip, ia, lin = estimators(*synth(i_true, az))
        print(f'  i_true={i_true:5.1f}  az={az:3d}  i_PCA={ip:6.2f}  i_amp={ia:6.2f}  直线度={lin:.3f}')
    print('=== B 合成反例（椭圆偏振，Z–R 相差 90°：PCA 退化）===')
    n = 3000
    tt = np.arange(n) / 100.0
    env = np.exp(-((np.arange(n) - 1500) / 200.0) ** 2)
    for i_true in (55.4,):
        z = env * math.cos(math.radians(i_true)) * np.sin(2 * np.pi * tt * 2)
        r = env * math.sin(math.radians(i_true)) * np.cos(2 * np.pi * tt * 2)
        ip, ia, lin = estimators(z, r, np.zeros(n))
        print(f'  i_true={i_true:5.1f}  i_amp={ia:6.2f}  i_PCA={ip:6.2f} (错)  直线度={lin:.3f}')
    print('  ⇒ 直线度低(<0.7)时只信 i_amp（见 reference §3）')


def part_c(base, sta, pick, half, expect):
    from obspy import read
    from obspy.signal.rotate import rotate_ne_rt
    from obspy.geodetics import gps2dist_azimuth

    rows = {r['sta']: r for r in csv.DictReader(
        l for l in open(f'{base}/process_flow/1_station_distribution/stations.csv') if not l.startswith('#'))}
    row = rows[sta]
    baz = gps2dist_azimuth(float(row['lat']), float(row['lon']), EV['lat'], EV['lon'])[1]

    trs = {c: read(glob.glob(f'{base}/sacprocess/remove_response/*{sta}..BH{c}*_resp.sac')[0])[0] for c in 'ZNE'}
    sac = trs['Z'].stats.sac
    o_val = getattr(sac, 'o', None)
    if o_val is not None and 0 < abs(float(o_val)) < 1e6:
        t_off = float(sac.b) - float(o_val)
    else:
        raise SystemExit('该事件 SAC 缺 o 字段：时间轴改用 starttime − EV_ORIGIN 口径（见 SKILL.md 时间轴铁律）')
    sr = trs['Z'].stats.sampling_rate
    t = t_off + np.arange(trs['Z'].stats.npts) / sr

    b, a = butter(4, list(BP), btype='band', fs=sr)
    Z = filtfilt(b, a, trs['Z'].data.astype(float))
    R, T = (filtfilt(b, a, x) for x in rotate_ne_rt(trs['N'].data.astype(float),
                                                    trs['E'].data.astype(float), baz))
    m = (t >= pick - half) & (t <= pick + half)
    ip, ia, lin = estimators(Z[m], R[m], T[m])

    print(f'=== C 真实窗（{sta}  ±{half:g} s  {BP[0]:g}-{BP[1]:g} Hz  baz={baz:.2f}°  样本 {int(m.sum())}）===')
    print(f'  i_PCA={ip:.3f}  i_amp={ia:.3f}  直线度={lin:.3f}')
    print('  口径：去响应速度、绕真实大圆 BAZ 旋转、带通 2–8 Hz、RMS 与 std 在本窗等价')
    if expect is not None:
        ok = abs(ip - expect) < 0.05 and abs(ia - expect) < 0.05
        print(f'  已存档值 {expect:.3f} → ' + ('一致' if ok else '不一致：先查读式/代码口径，不要改 README'))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--base', help='事件目录；不给则只跑合成 A/B')
    ap.add_argument('--sta', default='SA18')
    ap.add_argument('--pick', type=float, default=1195.15, help='T-P 窗中心（相对发震，s）')
    ap.add_argument('--half', type=float, default=5.0, help='窗半宽（s）')
    ap.add_argument('--expect', type=float, default=57.459, help='结果表/README 里已存档的 i，用于核对')
    ap.add_argument('--ev-lat', type=float, default=EV['lat'])
    ap.add_argument('--ev-lon', type=float, default=EV['lon'])
    args = ap.parse_args()
    EV.update(lat=args.ev_lat, lon=args.ev_lon)
    part_ab()
    if args.base:
        part_c(args.base, args.sta, args.pick, args.half, args.expect)
    else:
        print('（未给 --base：跳过真实窗核对）')
