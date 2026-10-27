#!/usr/bin/env python3
# =============================================================================
# 台阵几何体检：指定的台站群能不能当 F-K 子阵用
# -----------------------------------------------------------------------------
# 回答"这 N 个台能不能凑一个子阵"，逐项给判据与 PASS/FAIL，不靠感觉。
# 加新子阵、或在窗口/网格定不下来之前，先跑这个（见 references/window-length-aperture-matching.md
# 的「子阵可行性体检」一节）。
#
# 判据与公式：
#   · 无混叠      f   <= C/SMAX        C = 111.195/(2*d_min) [s/deg/Hz]
#   · 分辨率      Du  = 0.61*111.195/(f*D)             [s/deg]
#   · 方位误差    sig_BAZ ~= Du/|u|                     [rad]
#   · 几何地板    Du/|u| >= 1.28*d_min/D               (与 f 无关)
#   · 频点数      N = (f2-f1)*T 落在 FFT 栅格上的点数 >= 10 (理想) / 5 (下限)
#   · 窗长        T >= 1.5*D/v_app,  v_app = 111.195/|u_pred|
#   · 二维性      PCA 长短轴比 —— 近共线阵只能约束一个方向
#   · 相关对数    coarray 用 N(N-1)/2 对
# 阈值 3/10/2 是本项目经验值，不是文献定理；对外报时要标明出处。
#
# 输入台站表：15 列口径 sta,az,dist,land_km,ocean_km,tt_p,tt_s,tt_tp,tt_ts,snr[,snr_r,snr_t],lon,lat[,baz]
#   一律 csv.DictReader 按列名读，禁止固定索引。只用到 sta/lon/lat/tt_tp/tt_ts/snr/land_km/dist。
#   BAZ 现场用 gps2dist_azimuth(sta, ev)[1] 重算，不用表里的 baz 列（可能是漏 cos 纬度的历史 bug）。
#
# 用法（解释器用 obspy 环境）：
#   /home/sy/miniforge3/envs/obspy/bin/python check_array_geometry.py MAHA,AMPY,FOMA,CPSM,BKTA
#   可选： --csv <stations.csv>  --ev <lat,lon>  --osync <秒>  --u <|u| s/deg>  --T <只算这个窗>
# =============================================================================
import csv, math, sys
import numpy as np
from obspy.geodetics import gps2dist_azimuth

# ============================== ★ 参数区（CLI 可覆盖）
CSV_DEF = ('/mnt/e/sharedate/case/Madagascar/2013-08-17-mw62-southwest-indian-ridge/'
           'process_flow/1_station_distribution/stations.csv')
EVLA_DEF, EVLO_DEF = -34.887, 54.128      # 事件震中
OSYNC_DEF = 599.98                          # 发震在 synchronize 系里的位置(秒)
SMAX_COPY, SPTS_COPY = 30.0, 300            # 照抄周2021 复现的 fkmap 口径
FRNG_COPY = (1/1000.0, 1/25.0)
DT = 0.05                                   # 20 sps


def argval(flag, cast, default):
    if flag in sys.argv:
        i = sys.argv.index(flag)
        if '=' in flag:
            return cast(flag.split('=', 1)[1])
        if i + 1 < len(sys.argv):
            return cast(sys.argv[i + 1])
    return default


def read_stations(path):
    raw = open(path).read().splitlines()
    i = [k for k, l in enumerate(raw) if l.startswith('sta,')][0]
    out = {}
    for d in csv.DictReader(raw[i:]):
        o = {}
        for k, v in d.items():
            try:
                o[k] = float(v)
            except (TypeError, ValueError):
                o[k] = v
        out[o['sta'].strip()] = o
    return out


def main():
    pos = [a for a in sys.argv[1:] if not a.startswith('--')]
    if not pos:
        raise SystemExit(__doc__)
    want = [s.strip().upper() for s in pos[0].split(',')]
    csv_path = argval('--csv', str, CSV_DEF)
    evla, evlo = EVLA_DEF, EVLO_DEF
    if '--ev' in sys.argv:
        a = sys.argv[sys.argv.index('--ev') + 1].split(',')
        evla, evlo = float(a[0]), float(a[1])
    o_sync = argval('--osync', float, OSYNC_DEF)
    u_ref = argval('--u', float, 40.0)
    T_org = argval('--T', float, None)

    S = read_stations(csv_path)
    miss = [s for s in want if s not in S]
    if miss:
        raise SystemExit(f'台站表里没有: {miss}')
    g = [S[s] for s in want]
    n = len(g)
    num = lambda v: (float(v) if v not in ('', None) else float('nan'))
    print(f'台站群 ({n} 台): ' + ' '.join(want))
    print(f'事件: {evla}, {evlo}   | 台站表: {csv_path}\n')

    print('--- 逐台台账 ---')
    print(f'{"sta":6s}{"lat":>9s}{"lon":>10s}{"snr":>7s}{"land_km":>9s}{"dist_km":>9s}'
          f'{"tt_tp":>8s}{"tt_ts":>8s}{"BAZ_zhen":>10s}')
    for s in g:
        baz = gps2dist_azimuth(s['lat'], s['lon'], evla, evlo)[1]
        print(f"{s['sta']:6s}{num(s['lat']):9.4f}{num(s['lon']):10.4f}{num(s.get('snr')):7.1f}"
              f"{num(s.get('land_km')):9.0f}{num(s.get('dist')):9.0f}"
              f"{num(s.get('tt_tp')):8.1f}{num(s.get('tt_ts')):8.1f}{baz:10.2f}")

    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = gps2dist_azimuth(g[i]['lat'], g[i]['lon'],
                                                 g[j]['lat'], g[j]['lon'])[0] / 1000.0
    off = D[D > 0]
    dmin, dmax, dmed = float(off.min()), float(off.max()), float(np.median(off))
    print('\n--- 两两大圆距离 (km) ---')
    print('      ' + ''.join(f'{s["sta"]:>8s}' for s in g))
    for i, s in enumerate(g):
        print(f'{s["sta"]:6s}' + ''.join(f'{D[i, j]:8.1f}' for j in range(n)))

    lat0 = float(np.mean([s['lat'] for s in g]))
    lon0 = float(np.mean([s['lon'] for s in g]))
    bazs = np.array([gps2dist_azimuth(s['lat'], s['lon'], evla, evlo)[1] for s in g])
    C = 111.195 / (2 * dmin)
    x = (np.array([s['lon'] for s in g]) - lon0) * 111.195 * math.cos(math.radians(lat0))
    y = (np.array([s['lat'] for s in g]) - lat0) * 111.195
    P = np.column_stack([x - x.mean(), y - y.mean()])
    ev = np.linalg.eigvalsh(P.T @ P / n)
    aspect = math.sqrt(ev[1] / max(ev[0], 1e-9))
    print('\n--- 几何 ---')
    print(f'  质心 {lat0:.4f}, {lon0:.4f} | d_min {dmin:.1f} km | D {dmax:.1f} km | 中位 {dmed:.1f} km')
    print(f'  孔径比 D/d_min = {dmax/dmin:.2f}    真实大圆 BAZ {bazs.min():.2f}..{bazs.max():.2f} (均值 {bazs.mean():.2f})')
    print(f'  PCA 特征值 {ev[1]:.1f} / {ev[0]:.1f} km^2 -> 长短轴比 {aspect:.2f}'
          f'  (>=3 视为二维覆盖良好 —— 本项目经验阈值, 非文献定理; <3 说明单向伸长)')

    print('\n--- 理论目标(走时场最小二乘平面拟合) ---')
    for key, tag in (('tt_tp', 'T-P'), ('tt_ts', 'T-S')):
        t = np.array([s[key] for s in g], float)
        A = np.column_stack([np.ones(n), x, y])
        coef, *_ = np.linalg.lstsq(A, t, rcond=None)
        b, c = coef[1], coef[2]
        u = math.hypot(b, c) * 111.195
        baz = math.degrees(math.atan2(-b, -c)) % 360
        rms = float(np.sqrt(np.mean((t - A @ coef) ** 2)))
        print(f'  [{tag}] u_pred {u:6.2f} s/deg  BAZ_pred {baz:6.1f} deg  rms {rms:5.2f} s'
              f'  | |grad t| 与真实大圆 BAZ 可差几十度(孔径大时更甚), 判峰用大圆 BAZ, 报 BAZ_pred 必带 rms')

    print(f'\n--- 分辨率 (D={dmax:.1f} km, |u| 参考 {u_ref:g} s/deg) ---')
    print(f'{"f (Hz)":>9s}{"Du (s/deg)":>12s}{"sig_BAZ (deg)":>14s}{"Du/|u|":>9s}')
    for f in (0.005, 0.01, 0.02, 0.04):
        du = 0.61 * 111.195 / (f * dmax)
        print(f'{f:9.4f}{du:12.1f}{math.degrees(du/u_ref):14.1f}{du/u_ref:9.2f}')
    floor = 1.28 * dmin / dmax
    print(f'  几何地板 Du/|u| >= 1.28*d_min/D = {floor:.2f} ({100*floor:.0f}%) —— 与 f 无关')

    tt_tp = [s['tt_tp'] for s in g]
    tt_ts = [s['tt_ts'] for s in g]
    lo = min(tt_tp) - 30 + o_sync
    hi = max(tt_ts) + 50 + o_sync
    print('\n--- 窗与频点 ---')
    rows = [('span(gongshi)', hi - lo), ('1000 s', 1000.0)]
    if T_org:
        rows = [(f'custom {T_org:g} s', T_org)]
    smax_geo = math.ceil(1.4 * u_ref / 5) * 5
    for name, T in rows:
        nsp = 2 ** int(math.ceil(math.log2(T / DT)))
        df = 1 / (nsp * DT)
        # copy: 照抄目标口径(smax/frng 字面不动); geo: smax=1.4*u_pred, f2=min(C/smax,0.25)
        for ptag, smax, f2v in (('copy', SMAX_COPY, FRNG_COPY[1]),
                                ('geo', smax_geo, min(C / smax_geo, 0.25))):
            k = len([kk for kk in range(1, int(f2v / df) + 1) if kk * df >= FRNG_COPY[0] - 1e-12])
            flag = 'OK' if k >= 10 else ('marginal' if k >= 5 else '**TOO FEW**')
            print(f'  {name:15s} T={T:6.1f} s nsp={nsp:6d} df={df:.6f} | {ptag:4s} smax={smax:5.1f} '
                  f'f2={f2v:.4f} -> in-band {k:3d} pts {flag}')
    v_app = 111.195 / u_ref
    print(f'  [REQUIRED] T >= 1.5*(D/v_app) = 1.5*({dmax:.0f}/{v_app:.2f}) = {1.5*dmax/v_app:.0f} s'
          f'  (v_app = 111.195/|u| = {v_app:.2f} km/s)  <-- 公式 span 窗可能比这短, 短了就换长窗')

    print('\n--- 判据汇总 ---')
    print(f'  无混叠(copy: f2={FRNG_COPY[1]:.4f} <= C/smax={C/SMAX_COPY:.4f}) : '
          f'{"PASS" if FRNG_COPY[1] <= C/SMAX_COPY else "FAIL(in-band high-freq aliased, known boundary)"}')
    pairs = n * (n - 1) // 2
    print(f'  台数 {n} -> coarray {pairs} pairs (>=10) : {"PASS" if pairs >= 10 else "lOW"}')
    print(f'  D/d_min = {dmax/dmin:.2f} (>=2) : {"PASS" if dmax/dmin >= 2 else "FAIL"}')
    print(f'  PCA aspect = {aspect:.2f} (>=3) : {"PASS" if aspect >= 3 else "weak (elongated, one direction only)"}')
    print('\nDONE_CHECK')


if __name__ == '__main__':
    main()
