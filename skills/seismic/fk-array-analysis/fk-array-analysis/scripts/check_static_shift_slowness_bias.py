#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""check_static_shift_slowness_bias.py —— 逐台时间平移的「等效慢度偏置」复算（只用标准库）

【回答什么问题】
  "按陆地计时对齐（0 = 上岸）和按震中距对齐再送 F-K，是有区别的吧？"
  "这个平移等效于多少 s/° 的慢度？F-K 的读数会被它改掉多少？"

【原理】F-K 只吃台间时延的平面波梯度。任何逐台平移 dt_i 都把读数改成
       u_读数 = u_真 - u(dt_i)
       u(dt_i) = 把 dt_i 对台站局部坐标 (xE, xN，单位度) 做最小二乘平面拟合得到的斜率矢量。
  拟合出的线性部分会被波束的平面波拟合整段吸收（= 直接改慢度读数）；拟合残差是吸收不掉的
  波前形变（rms 秒）。全局平移落进常数项、不进斜率 ⇒ 不改读数（与 SKILL.md 同一条）。

【三种平移模式】
  ocean : dt_i = L_ocean,i / v_ocean   （陆地计时对齐；v 默认 1.48 km/s）
  land  : dt_i = L_land,i  / v_land    （只扣陆段理论走时）
  dist  : dt_i = d_i       / v_dist    （震中距/走时对齐 —— 会扣掉要测的 u 本身！）
  实测锚点（XA 1997-04-28、G1 七台，见 references/land-time-fk.md §9）：
  ocean/1.48 = 1.12 s/°、dist/3.6 = 30.78、land/3.6 = 31.19 s/° ⇒ 跑完照这个核一遍。
  判据提醒：land（陆地计时）的偏置是小量且多横向 ⇒ u 基本保留；dist 的偏置 = 111.195/v_ref
  且严格沿传播向 ⇒ 扣掉的正是 u，v_ref 取真值时峰塌到原点，按它报慢度是循环论证。

【用法】
  python check_static_shift_slowness_bias.py \
    --stations-csv <case>/process_flow/1_station_distribution/stations.csv \
    --stas SA16,SA17,SA18,SA22,SA23,SA25,SA29 \
    --event-lat -42.4134 --event-lon 42.6807
  --modes ocean,dist,land 只跑其一；默认三种全跑。

【口径】BAZ = 台站→震中的大圆起始方位角（球面公式自算，与 obspy gps2dist_azimuth[1] 同定义）；
        传播方向 e = BAZ - 180（震中→台站）；dt_i 为正 = 该台时间轴零点被提前 dt_i。
        只用标准库（不依赖 obspy），自带两项自检：球面几何、平面拟合还原线性场。
"""

import argparse
import csv
import math

R_EARTH = 6371.0


# ----------------------------------------------------------------------------- 球面几何
def gc_dist_km(la1, lo1, la2, lo2):
    """两点大圆距离 (km)，haversine（球面，本项目 2000-3000 km 尺度误差 <0.3%）"""
    p1, p2 = math.radians(la1), math.radians(la2)
    dl = math.radians(lo2 - lo1)
    a = math.sin((p2 - p1) / 2.0) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2.0) ** 2
    return 2.0 * R_EARTH * math.asin(min(1.0, math.sqrt(a)))


def gc_baz(la1, lo1, la2, lo2):
    """点1 → 点2 的大圆起始方位角（度，北 = 0、顺时针）。台→震中即真实大圆 BAZ"""
    p1, p2 = math.radians(la1), math.radians(la2)
    dl = math.radians(lo2 - lo1)
    y = math.sin(dl) * math.cos(p2)
    x = math.cos(p1) * math.sin(p2) - math.sin(p1) * math.cos(p2) * math.cos(dl)
    return math.degrees(math.atan2(y, x)) % 360.0


# ----------------------------------------------------------------------------- 平面拟合
def fit_plane(xE, xN, t):
    """最小二乘 t = a + uE*xE + uN*xN（t 秒，x 度）→ [a, uE, uN]；正规方程 + 高斯消元"""
    n = len(t)
    A = [[1.0, xE[i], xN[i]] for i in range(n)]
    M = [[sum(A[k][i] * A[k][j] for k in range(n)) for j in range(3)]
         + [sum(A[k][i] * t[k] for k in range(n))] for i in range(3)]
    for c in range(3):
        piv = max(range(c, 3), key=lambda r: abs(M[r][c]))
        M[c], M[piv] = M[piv], M[c]
        for r in range(3):
            if r != c:
                f = M[r][c] / M[c][c]
                for k in range(c, 4):
                    M[r][k] -= f * M[c][k]
    return [M[i][3] / M[i][i] for i in range(3)]


# ----------------------------------------------------------------------------- 自检
def self_test():
    ok = True
    d1 = gc_baz(-28.9321, 23.2257, -42.4134, 42.6807)      # 正
    d2 = gc_baz(-42.4134, 42.6807, -28.9321, 23.2257)      # 反
    # ⚠ 收敛角 = |(d1-d2) - 180| 的偏离。旧代码写成 abs((d1-d2+180)%360-180):
    #   正反方位恰好差 180° 时该式恒等于 180 -> **永远报 FAIL**(实测 G1/MAHA5 都报 168.5,
    #   真收敛角 = 180-168.5 = 11.5°)。修法如下。
    conv = abs((d1 - d2) % 360.0 - 180.0)
    print(f'[自检1] 球面 BAZ 正反相差 180 的偏离 = {conv:.2f} deg（= 大圆收敛角，应 < 40）')
    ok &= conv < 40.0
    xE = [-0.5, 0.2, 1.1, 0.0]
    xN = [1.0, -0.3, 0.4, 0.7]
    t = [1.0 + 2.0 * xE[i] - 3.0 * xN[i] for i in range(4)]
    a = fit_plane(xE, xN, t)
    err = max(abs(a[0] - 1.0), abs(a[1] - 2.0), abs(a[2] + 3.0))
    print(f'[自检2] 平面拟合还原无噪声线性场，最大误差 = {err:.2e}（应 ~1e-15）')
    ok &= err < 1e-9
    print('[自检] ' + ('PASS' if ok else 'FAIL —— 下面的数字不可用'))
    print('=' * 96)
    return ok


# ----------------------------------------------------------------------------- 报告
def report(label, shift, baz, xE, xN, stas):
    a, uE, uN = fit_plane(xE, xN, shift)
    mag = math.hypot(uE, uN)
    az = math.degrees(math.atan2(uE, uN)) % 360.0
    e = [math.radians(b - 180.0) for b in baz]          # 传播方向（震中→台）
    n = len(stas)
    along = sum(uE * math.sin(e[i]) + uN * math.cos(e[i]) for i in range(n)) / n
    tran = sum(uE * math.cos(e[i]) - uN * math.sin(e[i]) for i in range(n)) / n
    res = [shift[i] - (a + uE * xE[i] + uN * xN[i]) for i in range(n)]
    rms = math.sqrt(sum(r * r for r in res) / n)
    print(f'  {label:<22} 平移量 {min(shift):8.2f}…{max(shift):8.2f} s（跨度 {max(shift) - min(shift):6.2f} s）')
    print(f'  {"":<22} 等效偏置 |u| = {mag:6.2f} s/°   方位 {az:6.1f}°   '
          f'沿传播向 {along:+6.2f}   横向 {tran:+6.2f}   平面拟合残差 rms {rms:.2f} s')
    return mag, az, along, tran, rms


def main():
    ap = argparse.ArgumentParser(description='逐台时间平移的等效慢度偏置复算（只用标准库）')
    ap.add_argument('--stations-csv', required=True, help='含 sta/ocean_km/land_km/lat/lon 列的 stations.csv')
    ap.add_argument('--stas', required=True, help='逗号分隔的子阵成员（用逗号，勿用空格）')
    ap.add_argument('--event-lat', type=float, required=True)
    ap.add_argument('--event-lon', type=float, required=True)
    ap.add_argument('--modes', default='ocean,dist,land', help='ocean / dist / land 的任意组合')
    ap.add_argument('--ocean-speed', type=float, default=1.48, help='海洋段速度 km/s（默认 1.48）')
    ap.add_argument('--dist-speeds', default='1.48,3.6,6.0')
    ap.add_argument('--land-speeds', default='3.6,6.0')
    args = ap.parse_args()

    self_test()

    rows = {r['sta']: r for r in csv.DictReader(l for l in open(args.stations_csv) if not l.startswith('#'))}
    stas = [s.strip() for s in args.stas.split(',') if s.strip()]
    missing = [s for s in stas if s not in rows]
    if missing:
        raise SystemExit(f'stations.csv 里缺这些台站：{missing}')

    lat = [float(rows[s]['lat']) for s in stas]
    lon = [float(rows[s]['lon']) for s in stas]
    la0, lo0 = sum(lat) / len(lat), sum(lon) / len(lon)
    xE = [(lo - lo0) * math.cos(math.radians(la0)) for lo in lon]      # 局部切平面（度）
    xN = [la - la0 for la in lat]
    baz = [gc_baz(la, lo, args.event_lat, args.event_lon) for la, lo in zip(lat, lon)]
    dist = [gc_dist_km(la, lo, args.event_lat, args.event_lon) for la, lo in zip(lat, lon)]
    ocean = [float(rows[s]['ocean_km']) for s in stas]
    land = [float(rows[s]['land_km']) for s in stas]

    print(f'子阵 {len(stas)} 台：真实大圆 BAZ {min(baz):.1f}–{max(baz):.1f}°（台→震中）')
    print(f'  震中距 {min(dist):.0f}–{max(dist):.0f} km（跨度 {max(dist) - min(dist):.0f} km）')
    print(f'  ocean_km 跨度 {max(ocean) - min(ocean):.1f} km ⇒ 上岸时刻跨度 '
          f'{(max(ocean) - min(ocean)) / args.ocean_speed:.1f} s（= 两版窗长之差的上界）')
    print(f'  land_km  跨度 {max(land) - min(land):.1f} km')
    print('-' * 96)

    modes = [m.strip() for m in args.modes.split(',') if m.strip()]
    out = []
    if 'ocean' in modes:
        out.append(report(f'ocean/{args.ocean_speed}（陆地计时）',
                          [o / args.ocean_speed for o in ocean], baz, xE, xN, stas))
    if 'dist' in modes:
        for v in [float(x) for x in args.dist_speeds.split(',')]:
            out.append(report(f'dist/{v}（震中距对齐）', [d / v for d in dist], baz, xE, xN, stas))
    if 'land' in modes:
        for v in [float(x) for x in args.land_speeds.split(',')]:
            out.append(report(f'land/{v}（只扣陆段）', [l / v for l in land], baz, xE, xN, stas))

    print('-' * 96)
    print('解读：读数会被改成 u_真 − u(Δt_i)。① ocean（陆地计时）的 |u| 是小量、方向多不沿传播向'
          '⇒ 相当于静校正，u 基本保留；')
    print('      ② dist/v 的 |u| ≈ 111.195/v 且横向≈0 ⇒ 它扣掉的就是要测的 u：v 取真值时峰塌到网格原点，'
          '只能读出参考速度误差，')
    print('         按它报慢度是循环论证（此时 SMAX 只需覆盖容许的速度误差，而非 111.195/v）。')
    print('      ③ 只改画图轴（不动数据）不进这张表 —— 那种做法读数一字不变。')
    print(f'  参考量级核对：111.195/3.6 = {111.195 / 3.6:.2f}、111.195/6.0 = {111.195 / 6.0:.2f}、'
          f'111.195/1.48 = {111.195 / 1.48:.2f} s/°')
    return out


if __name__ == '__main__':
    main()
