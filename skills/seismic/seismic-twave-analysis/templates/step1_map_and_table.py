#!/usr/bin/env python3
"""第一步-台站表格与分布图：数据表格（xlsx/csv）+ 两张台站分布图。

输出:
- stations.xlsx（格式化台站表）
- stations.csv（后续分析输入）
- station_map_all_names.png（标注全部台站名）
- station_map_top3.png（标注 SNR 前3台站名）

用法:
    python step1_map_and_table.py
"""
import argparse
import csv
import glob
import math
import os
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pygmt
from geographiclib.geodesic import Geodesic
from obspy import read
from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side

CODE_ROOT = "/home/sy/codex"
if CODE_ROOT not in sys.path:
    sys.path.insert(0, CODE_ROOT)
from lib.arrivals import compute_theoretical_arrivals, dist_az

# ==================== 参数区（可修改） ====================
# ---- 路径：**全部由环境变量传入**（模板不含任何事件的路径）----
import os as _os
_DEMO = Path(__file__).resolve().parent                       # 本步目录本身
_T_BASE = _os.environ.get("T_BASE", "").strip()
if not _T_BASE:
    sys.exit("请用环境变量 T_BASE 指定事件根目录（其下应有 raw/ 或 SAC 文件，坐标从 SAC 头读）")
BASE = Path(_T_BASE).expanduser()
RESP_DIR = BASE / "sacprocess/remove_response"
OUT_DIR = Path(_os.environ.get("T_OUTDIR", _DEMO))               # stations.csv 写本步目录
MAP_DIR = Path(_os.environ.get("T_MAPDIR", _DEMO))                # 分布图也写本步目录

# 事件信息（★ 每个事件填一次；也可用 EVENT_* 环境变量覆盖）
EVENT_DATE = _os.environ.get("EVENT_DATE", "YYYY-MM-DD")
EVENT_MW = float(_os.environ.get("EVENT_MW", "0.0"))
_elat, _elon = _os.environ.get("EVENT_LAT"), _os.environ.get("EVENT_LON")
if not (_elat and _elon):
    raise SystemExit("必须给 EVENT_LAT / EVENT_LON：震中坐标决定震中距、到时与 SNR 窗，"
                     "用 (0,0) 兜底会算出 SNR=0 的废表（实测 39 台全部 SNR<6 dB）")
EVLA, EVLO = float(_elat), float(_elon)
EVDP = float(_os.environ.get("EVENT_DEP") or _os.environ.get("EVENT_DEPTH") or "10.0")

# 走时与SNR
V_OCEAN = 1.48   # SOFAR km/s
V_LOVE = 3.6     # Love km/s
SNR_THRESHOLD = 6.0   # 表格高亮阈值（与剖面图显示阈值一致）

# 台站筛选（★ 按事件填：EXCLUDE_STA=A,B 环境变量传入，或直接写在这里）
DEFAULT_EXCLUDE = set(x for x in _os.environ.get("EXCLUDE_STA", "").split(",") if x.strip())
#   剔除理由写进事件目录的 station_qc_notes.md：灵敏度错、全程嘈杂、无 T 波、偏出阵列……
# 排除框: [lon_min, lon_max, lat_min, lat_max]，落在框内的台站被剔除
EXCLUDE_BOXES = []
# 经度范围（只保留区域内台站）：None = 不裁剪；例 LON_RANGE=43.0,51.0
_LON = _os.environ.get("LON_RANGE", "").strip()
LON_MIN, LON_MAX = (tuple(float(v) for v in _LON.split(",")) if _LON else (None, None))

# 地图参数
STATION_SIZE = 0.5
STATION_PEN = "0.3p,black"
LABEL_FONT = "8p,Helvetica-Bold,black"
STROKE_FONT = LABEL_FONT.replace("black", "white")   # 描边字体 = 正文字体（只换色），字号必须一致
# 台站标签位置：auto(默认，右下→左下→正下候选避让) | right | left | bottom
# ★ 全站台名一起标时必须用 auto：固定 bottom 在横排密集处必然互相压住
LABEL_POS = "auto"
WHITE_STROKE = 0.03
# 只标注这些台站的名字（[] = 全部标注）—— 默认全部标，不单独突出某一台
LABEL_ONLY = []
# 红三角突出显示的台站（[] = 不突出任何台站）
HIGHLIGHT_STA = []
HL_SIZE = 0.60          # 红三角尺寸 (c)，略大于普通三角
HL_PEN = "0.4p,black"
# 图幅宽度：地图本体宽度（期刊双栏正文宽 19 cm；含右侧色标后总宽须 <=19 cm）
MAP_WIDTH = "15.5c"
# ★ 配色开关（默认 = 现行定稿；`MAP_STA_CMAP` / `MAP_TERRAIN_CMAP` 可临时换，不动代码）
#   台站三角：viridis（默认）| batlow | magma | plasma | turbo | cool | 或任何 GMT cpt 名
#   地形底图：geo（默认）| gray | earth | relief | etopo1
STATION_CMAP = _os.environ.get("MAP_STA_CMAP", "viridis")
TERRAIN_CMAP = _os.environ.get("MAP_TERRAIN_CMAP", "geo")
# 地形色标位置：right = 两条色标都竖排在右侧（省竖向空间，期刊版）
#               bottom = 原横向置底（总图更高）
ELEV_CBAR = "right"
DPI = int(_os.environ.get("MAP_DPI", 600))   # 组合图（矢量线条 + 栅格地形）：PNG 600 dpi（MAP_DPI 可降为预览用）
PDF_OUT = False         # 该图 82% 面积是栅格地形，PDF 会内嵌高分辨率底图（约 100 MB），
                        # 不适合投稿；600 dpi PNG 已满足期刊要求。需要 PDF 时改 True。
MARGIN = 6
GRID_RES = "@earth_relief_15s"

# 震源机制（★ 用事件自己的 GCMT 节面；环境变量 CMT_STRIKE / CMT_DIP / CMT_RAKE）
CMT_STRIKE = int(_os.environ.get("CMT_STRIKE", "0"))
CMT_DIP = int(_os.environ.get("CMT_DIP", "90"))
CMT_RAKE = int(_os.environ.get("CMT_RAKE", "0"))
# ========================================================


def _is_excluded(lon, lat):
    for lo, hi, la, ha in EXCLUDE_BOXES:
        if lo <= lon <= hi and la <= lat <= ha:
            return True
    return False


def compute_snr_welch(data, delta, tt_p, tt_tp, tt_ts, t_off=0.0):
    """T波 SNR (dB) — Welch PSD 2-8 Hz 能量比。

    信号窗: [tt_tp - 20, tt_ts + 100]；噪声窗: P 前与信号窗等长
    [tt_p - L, tt_p]（L = 信号窗长度；2026-09-02 用户定稿，替换原
    0.8*tt_p 动态窗，保证两窗 Welch 段数/频率分辨率一致）。
    t_off: 数据起点相对发震时刻的偏移（秒）。SAC 常从发震前 600s 开始，
           若不修正会把信号窗整体错位（2026-08-27 修复，原版 t 从 0 起算）。
    """
    from scipy.signal import welch
    sr = 1.0 / delta
    t = np.arange(len(data)) * delta + t_off
    s_start = int(max(0, tt_tp - 20 - t_off) / delta)
    s_end = int(min(t[-1], tt_ts + 100 - t_off) / delta)
    sig_len = (tt_ts + 100) - (tt_tp - 20)
    actual_noise = sig_len
    n_start = int(max(0, tt_p - actual_noise - t_off) / delta)
    n_end = int(max(0, tt_p - t_off) / delta)
    if s_end <= s_start or n_end <= n_start:
        return 0.0
    f_sig, psd_sig = welch(data[s_start:s_end].astype(np.float64), sr, nperseg=min(256, s_end - s_start))
    f_noise, psd_noise = welch(data[n_start:n_end].astype(np.float64), sr, nperseg=min(256, n_end - n_start))
    e_sig = np.sum(psd_sig[(f_sig >= 2) & (f_sig <= 8)])
    e_noise = np.sum(psd_noise[(f_noise >= 2) & (f_noise <= 8)])
    if e_noise > 0 and e_sig > 0:
        ratio = e_sig / e_noise
        if ratio > 0 and np.isfinite(ratio):
            return 10 * np.log10(ratio)
    return 0.0


def _get_t_off(tr):
    """数据起点相对发震时刻的偏移（秒）。SAC o 有效时用 b-o；否则从文件名解析。"""
    sac = tr.stats.sac
    o_val = getattr(sac, 'o', None)
    if o_val is not None and 0 < abs(o_val) < 1e6:
        return sac.b - o_val
    # fallback: 从文件名解析发震时刻
    import re
    m = re.search(r'\.(\d{4})\.(\d{3})\.(\d{6})', os.path.basename(tr.filename).upper())
    if m:
        from obspy import UTCDateTime
        y = int(m.group(1)); jd = int(m.group(2)); hm = m.group(3)
        origin = UTCDateTime(year=y, julday=jd, hour=int(hm[0:2]),
                             minute=int(hm[2:4]), second=int(hm[4:6]))
        return tr.stats.starttime - origin
    return 0.0


def _rotate_rt(n, e, baz):
    """N/E → R/T 旋转（obspy 约定，R 背离震源；2026-09-02 统一）。"""
    from obspy.signal.rotate import rotate_ne_rt
    r, t = rotate_ne_rt(np.asarray(n, dtype=float), np.asarray(e, dtype=float), float(baz))
    return r, t


def compute_station_rows():
    """读全部三分量 SAC → 理论到时 → Z/R/T SNR，返回台站行列表。"""
    bhz_files = sorted(glob.glob(f"{RESP_DIR}/*BHZ*_resp.sac"))
    if not bhz_files:
        bhz_files = sorted(glob.glob(f"{RESP_DIR}/*.BHZ*"))
    print(f"找到 {len(bhz_files)} 个 BHZ 文件")

    rows, seen_sta = [], set()
    for fn in bhz_files:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                tr = read(fn)[0]
            sac = tr.stats.sac
            stla = getattr(sac, 'stla', 0)
            stlo = getattr(sac, 'stlo', 0)
            if stla == 0 or stlo == 0:
                continue
            sta = os.path.basename(fn).split('.')[1]
            if sta in seen_sta or sta in DEFAULT_EXCLUDE:
                continue
            if _is_excluded(stlo, stla):
                continue
            if LON_MIN is not None and not (LON_MIN <= stlo <= LON_MAX):
                continue
            seen_sta.add(sta)

            dist_km, az, baz = dist_az(EVLO, EVLA, stlo, stla)
            arrs = compute_theoretical_arrivals(EVLO, EVLA, EVDP, stlo, stla)
            tt_p = arrs.get('P', 0)
            tt_s = arrs.get('S', 0)
            tt_tp = arrs.get('T-P', 0)
            tt_ts = arrs.get('T-S', 0)
            ocean_km = arrs.get('ocean_km', 0)
            land_km = arrs.get('land_km', 0)

            t_off = _get_t_off(tr)
            # Z 分量 SNR
            if tt_tp > 0:
                snr = compute_snr_welch(tr.data, tr.stats.delta,
                                        tt_p=tt_p, tt_tp=tt_tp, tt_ts=tt_ts,
                                        t_off=t_off)
            else:
                snr = -999

            # R/T 分量 SNR（需 BHN/BHE；缺水平分量台站置空）
            # 2026-09-02: 未带通 + 按分量起点样本差对齐（同一天同网络 reference 偏差一致，
            #             差值可信），修复 trim 后仍用 BHZ t_off 导致的 R/T 窗错位（VOI 曾低估~35dB）。
            snr_r, snr_t = "", ""
            net = os.path.basename(fn).split('.')[0]
            tr_n = tr_e = None
            m_n = sorted(glob.glob(f"{RESP_DIR}/{net}.{sta}..BHN*_resp.sac"))
            m_e = sorted(glob.glob(f"{RESP_DIR}/{net}.{sta}..BHE*_resp.sac"))
            if m_n and m_e:
                tr_n = read(m_n[0])[0]
                tr_e = read(m_e[0])[0]
                dt = tr.stats.delta
                dn = int(round((tr_n.stats.starttime - tr.stats.starttime) / dt))
                de = int(round((tr_e.stats.starttime - tr.stats.starttime) / dt))
                imax = max(0, dn, de)
                # 公共段长度（各分量可用样本数约束，BHZ 以其全长 npts 为界）
                cmin = min(tr.stats.npts - imax,
                           tr_n.stats.npts - imax + dn,
                           tr_e.stats.npts - imax + de)
                if cmin > 0 and tt_tp > 0:
                    i1 = imax + cmin
                    n = tr_n.data[imax - dn: i1 - dn].astype(np.float64)
                    e = tr_e.data[imax - de: i1 - de].astype(np.float64)
                    r, tc = _rotate_rt(n, e, baz)
                    # 公共段第0样本相对发震 = t_off + imax*dt
                    t_off_seg = t_off + imax * dt
                    snr_r = round(compute_snr_welch(r, dt,
                                                    tt_p=tt_p, tt_tp=tt_tp, tt_ts=tt_ts,
                                                    t_off=t_off_seg), 1)
                    snr_t = round(compute_snr_welch(tc, dt,
                                                    tt_p=tt_p, tt_tp=tt_tp, tt_ts=tt_ts,
                                                    t_off=t_off_seg), 1)

            rows.append({
                'sta': sta, 'az': round(az, 1), 'dist': round(dist_km, 0),
                'land_km': round(land_km, 0), 'ocean_km': round(ocean_km, 0),
                'tt_p': round(tt_p, 1), 'tt_s': round(tt_s, 1),
                'tt_tp': round(tt_tp, 1), 'tt_ts': round(tt_ts, 1),
                'snr': round(snr, 1), 'snr_r': snr_r, 'snr_t': snr_t,
                'lon': round(stlo, 4), 'lat': round(stla, 4),
                'baz': round(baz, 1),
            })
        except Exception:
            continue

    rows.sort(key=lambda x: x['snr'], reverse=True)
    snr_ok = sum(1 for r in rows if r['snr'] >= SNR_THRESHOLD)
    print(f"总台站: {len(rows)} | SNR≥{SNR_THRESHOLD:g}dB: {snr_ok}站")
    return rows


def write_station_table(rows):
    """生成 stations.xlsx（微软雅黑格式化，输出 maps_info/）和 stations.csv（根目录，下游依赖）。"""
    out_dir = MAP_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    snr_ok = sum(1 for r in rows if r['snr'] >= SNR_THRESHOLD)

    # ---- xlsx ----
    wb = Workbook()
    ws = wb.active
    ws.title = 'Sheet'
    ws.merge_cells('A1:O1')
    info = ws.cell(row=1, column=1,
                   value=f"事件: {EVENT_DATE} Mw {EVENT_MW}  |  {EVLA}°N, {EVLO}°E  |  "
                         f"深度 {EVDP} km  |  总台站: {len(rows)}  |  "
                         f"SNR≥{SNR_THRESHOLD:g}dB: {snr_ok}站")
    info.font = Font(name='微软雅黑', bold=True, size=16)
    info.alignment = Alignment(horizontal='center', vertical='center')

    hdrs_cn = ['台站', '方位角', '震中距', '陆地距', '海洋距',
               '信噪比', 'R信噪比', 'T信噪比', 'P波', 'S波', 'T-P', 'T-S',
               '经度', '纬度', '反方位角']
    for ci, h in enumerate(hdrs_cn, 1):
        c = ws.cell(row=2, column=ci, value=h)
        c.font = Font(name='微软雅黑', bold=True, size=14, color='000000')
        c.alignment = Alignment(horizontal='center', vertical='center')

    hdrs_en = ['STA', 'AZ(°)', 'DIST(km)', 'LAND(km)', 'OCEAN(km)',
               'SNR(dB)', 'SNR_R(dB)', 'SNR_T(dB)', 'P(s)', 'S(s)', 'T-P(s)', 'T-S(s)',
               'LON', 'LAT', 'BAZ(°)']
    for ci, h in enumerate(hdrs_en, 1):
        c = ws.cell(row=3, column=ci, value=h)
        c.font = Font(name='微软雅黑', bold=True, size=14, color='FFFFFF')
        c.fill = PatternFill('solid', fgColor='4472C4')
        c.alignment = Alignment(horizontal='center', vertical='center')

    for ri, r in enumerate(rows, 4):
        vals = [r['sta'], r['az'], r['dist'], r['land_km'], r['ocean_km'],
                r['snr'], r['snr_r'] if r['snr_r'] != '' else '', r['snr_t'] if r['snr_t'] != '' else '',
                r['tt_p'], r['tt_s'], r['tt_tp'], r['tt_ts'],
                r['lon'], r['lat'], r['baz']]
        for ci, v in enumerate(vals, 1):
            cell = ws.cell(row=ri, column=ci, value=v)
            cell.alignment = Alignment(horizontal='center', vertical='center')
            cell.font = Font(name='微软雅黑', size=14)
            if ci == 6 and isinstance(v, (int, float)) and v >= SNR_THRESHOLD:
                cell.fill = PatternFill('solid', fgColor='C6EFCE')

    widths = [12, 10, 12, 12, 12, 12, 12, 12, 10, 10, 12, 12, 14, 14, 12]
    for i, w in enumerate(widths, 1):
        col_letter = chr(64 + i) if i <= 26 else 'A' + chr(64 + i - 26)
        ws.column_dimensions[col_letter].width = w

    thin = Side(style='thin', color='000000')
    border = Border(left=thin, right=thin, top=thin, bottom=thin)
    for row in ws.iter_rows(min_row=1, max_row=len(rows) + 3, max_col=15):
        for cell in row:
            cell.border = border

    xlsx_path = out_dir / "stations.xlsx"
    try:
        wb.save(xlsx_path)
    except PermissionError:
        xlsx_path = out_dir / "stations_v2.xlsx"
        wb.save(xlsx_path)
    print(f"✓ {xlsx_path}")

    # ---- csv（必须留在根目录，下游脚本硬依赖）----
    csv_path = OUT_DIR / "stations.csv"
    lines = [f"# Event: {EVENT_DATE} Mw {EVENT_MW} | {EVLA}N, {EVLO}E | Depth {EVDP}km\n",
             "sta,az,dist,land_km,ocean_km,tt_p,tt_s,tt_tp,tt_ts,snr,snr_r,snr_t,lon,lat,baz\n"]
    for r in rows:
        lines.append(f"{r['sta']},{r['az']:.1f},{r['dist']:.0f},{r['land_km']:.0f},{r['ocean_km']:.0f},"
                     f"{r['tt_p']:.1f},{r['tt_s']:.1f},{r['tt_tp']:.1f},{r['tt_ts']:.1f},{r['snr']:.1f},"
                     f"{r['snr_r'] if r['snr_r'] != '' else ''},{r['snr_t'] if r['snr_t'] != '' else ''},"
                     f"{r['lon']:.4f},{r['lat']:.4f},{r['baz']:.1f}\n")
    csv_path.write_text("".join(lines))
    print(f"✓ {csv_path}")


def plot_station_map(rows, label_mode="all", out_name="station_map_all_names.png", label_only=None):
    """台站分布图（期刊级）：地形底图 + SNR 着色 + 震源机制球 + 标注。"""
    stations = [dict(sta=r["sta"], lon=r["lon"], lat=r["lat"],
                     snr=r["snr"], land_km=r["land_km"]) for r in rows]
    print(f"台站: {len(stations)}")

    snrs = [s["snr"] for s in stations]
    lons_all = [s["lon"] for s in stations] + [EVLO]
    lats_all = [s["lat"] for s in stations] + [EVLA]
    R = [min(lons_all) - MARGIN, max(lons_all) + MARGIN,
         min(lats_all) - MARGIN, max(lats_all) + MARGIN]

    pygmt.config(MAP_GRID_PEN_PRIMARY="0.3p,gray40,--", PS_MEDIA="A0")
    fig = pygmt.Figure()

    elev_cpt = os.path.join(tempfile.gettempdir(), "elev_relief.cpt")
    pygmt.makecpt(cmap=TERRAIN_CMAP, series=[-8000, 5000, 100], continuous=True, output=elev_cpt)
    fig.grdimage(grid=GRID_RES, region=R, projection=f"M{MAP_WIDTH}",
                 cmap=elev_cpt, shading="+a315+nt0.8+m0")
    fig.coast(shorelines="0.3p,gray30", borders=["1/0.3p,gray50", "2/0.2p,gray70"],
              resolution="h", water="lightblue@60", land="white@80")
    fig.basemap(frame=["xa10f5", "ya10f5", "WSne"])
    fig.grdcontour(grid=GRID_RES, region=R, levels=[-1000], pen="0.8p,black,--")

    snr_valid = [s for s in snrs if np.isfinite(s) and s > 0]
    snr_min, snr_max = max(0, min(snr_valid)), max(snr_valid) if snr_valid else 10
    pygmt.makecpt(cmap=STATION_CMAP, series=[math.floor(snr_min), math.ceil(snr_max), 1])

    tf = tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False)
    for s in stations:
        tf.write(f"{s['lon']} {s['lat']} {s['snr']}\n")
    tf.close()
    fig.plot(data=tf.name, style=f"t{STATION_SIZE}c", cmap=True, pen=STATION_PEN)
    os.unlink(tf.name)

    # 红三角突出关键台站（HIGHLIGHT_STA 非空时才画）
    if HIGHLIGHT_STA:
        hl = [s for s in stations if s["sta"] in set(HIGHLIGHT_STA)]
        if hl:
            th = tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False)
            for s in hl:
                th.write(f"{s['lon']} {s['lat']}\n")
            th.close()
            fig.plot(data=th.name, style=f"t{HL_SIZE}c", fill="red", pen=HL_PEN)
            os.unlink(th.name)

    # 台站名标注：位置可选（右下/左下/正下），auto 时按候选避让避免重叠（无引线，贴台站）
    label_stations = sorted(stations, key=lambda x: x["snr"], reverse=True)
    if label_only:
        # 只标注指定台站（LABEL_ONLY = [] 时走全标分支，不进入这里）
        want = set(label_only)
        label_stations = [s for s in label_stations if s["sta"] in want]
    elif label_mode != "all":
        label_stations = label_stations[:3]

    # 候选锚点偏移：一律以 cm 定义再换算成度（图幅宽度变化时自动适配）
    # 8p 粗体台名（4 字符）约 0.62 cm 宽、0.28 cm 高；符号半宽/半高 = 符号尺寸/2 cm
    _cm_w = float(MAP_WIDTH.rstrip("c"))
    _dpc_x = (R[1] - R[0]) / _cm_w                                  # 每 cm 对应经度(度)
    _dpc_y = _dpc_x / math.cos(math.radians(0.5 * (R[2] + R[3])))   # 纬度方向(Mercator 近似)
    TXT_W_CM, TXT_H_CM = 0.62, 0.28

    def _dgx(cm):
        """cm -> 经度增量(度)"""
        return cm * _dpc_x

    def _dgy(cm):
        """cm -> 纬度增量(度)"""
        return cm * _dpc_y

    LABEL_W, LABEL_H = 1.15 * _dgx(TXT_W_CM), 1.2 * _dgy(TXT_H_CM)  # 避让 box（略大于文字）
    _HALF = HL_SIZE / 2                                             # 三角半高/半宽 (cm)
    _RIGHT = (_dgx(_HALF + 0.06), -_dgy(_HALF))                     # 右侧（与三角同底）
    _LEFT = (-_dgx(_HALF + 0.06) - 2 * _dgx(TXT_W_CM), -_dgy(_HALF))
    _BOTTOM = (-_dgx(TXT_W_CM / 2), -_dgy(_HALF + 0.08 + TXT_H_CM))  # 正下方居中，贴三角底边
    _TOP = (-_dgx(TXT_W_CM / 2), _dgy(_HALF + 0.08))

    def _fan(r_cm):
        """以三角中心为圆心、半径 r_cm 的一圈 8 方向候选（返回标签框左下角偏移）"""
        out = []
        for ang in (0, 45, 90, 135, 180, 225, 270, 315):
            a = math.radians(ang)
            out.append((_dgx(r_cm * math.cos(a)) - LABEL_W / 2,
                        _dgy(r_cm * math.sin(a)) - LABEL_H / 2))
        return out

    CAND_MAP = {"right": [_RIGHT], "left": [_LEFT], "bottom": [_BOTTOM],
                # auto：先贴符号四向，再向外扩两圈（0.95 cm / 1.55 cm）逐个试
                # 全阵列 37 个台名一起标时必需，否则密集区会有 ~10 个标签互相压住
                "auto": [_RIGHT, _LEFT, _BOTTOM, _TOP] + _fan(0.95) + _fan(1.55)}
    cands = CAND_MAP.get(LABEL_POS, CAND_MAP["auto"])

    def _lbox(lon, lat, dx, dy):
        return (lon + dx, lon + dx + LABEL_W, lat + dy, lat + dy + LABEL_H)

    def _overlap(a, b, gap=0.10):
        return not (a[1] + gap <= b[0] or b[1] + gap <= a[0]
                    or a[3] + gap <= b[2] or b[3] + gap <= a[2])

    placed = []
    for s in label_stations:
        chosen = cands[0]
        for dx, dy in cands:
            bx = _lbox(s["lon"], s["lat"], dx, dy)
            if not any(_overlap(bx, pb) for pb in placed):
                chosen = (dx, dy)
                placed.append(bx)
                break
        lx = s["lon"] + chosen[0]
        ly = s["lat"] + chosen[1]
        # 白色细描边（4 方向白字偏移 WHITE_STROKE cm）+ 黑字
        # ★ 描边字号必须与正文字号严格一致（LABEL_FONT 换色而来）；用 6p 描边 + 8p 正文时，
        #   白字比黑字小一圈，白边缩在字里面 -> 远看是"空心字 / 边缘脏白"，不是描边
        for ddx, ddy in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
            fig.text(x=lx, y=ly, text=s["sta"], font=STROKE_FONT,
                     justify="BL", offset=f"{ddx * WHITE_STROKE}c/{ddy * WHITE_STROKE}c")
        fig.text(x=lx, y=ly, text=s["sta"], font=LABEL_FONT, justify="BL")

    geod = Geodesic.WGS84

    # 最近/最远陆地台站两条方位角线
    by_land = sorted(stations, key=lambda x: x["land_km"])
    for s in (by_land[0], by_land[-1]):
        g = geod.Inverse(EVLA, EVLO, s["lat"], s["lon"])
        az = (g["azi1"] + 360) % 360
        g2 = geod.Direct(EVLA, EVLO, az, g["s12"] + 3.0 * 111000)
        fig.plot(x=[EVLO, g2["lon2"]], y=[EVLA, g2["lat2"]], pen="0.5p,gray30,--")

    # 震源机制球
    fig.meca(spec={'strike': CMT_STRIKE, 'dip': CMT_DIP, 'rake': CMT_RAKE,
                   'magnitude': EVENT_MW},
             longitude=EVLO, latitude=EVLA, depth=EVDP,
             scale='0.60c', compression_fill='black', extension_fill='white',
             outline='0.55p,black', convention='aki')

    # 事件信息（白描边）：12p 全串约 4.2 cm 宽，向右会顶出图框（右边界距震中仅 ~3 cm）
    #   -> 锚在沙滩球左上方空海区，BR 对齐让文字向左延伸
    lx, ly, jt = EVLO + _dgx(0.30), EVLA + _dgy(1.30), "BR"
    etxt = f"{EVENT_DATE}  Mw {EVENT_MW}"
    for ddx, ddy in [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]:
        fig.text(x=lx, y=ly, text=etxt, font="12p,Helvetica-Bold,white",
                 justify=jt, offset=f"{ddx * WHITE_STROKE}c/{ddy * WHITE_STROKE}c")
    fig.text(x=lx, y=ly, text=etxt, font="12p,Helvetica-Bold,black", justify=jt)

    # 色标：两条都竖排右侧时总图更矮（期刊限高 23.7 cm）
    if ELEV_CBAR == "bottom":
        fig.colorbar(frame="xa10f1+lSNR (dB)", position="JMR+w7c/0.3c+o1c/0c+v")
        fig.colorbar(cmap=elev_cpt, frame="xa2000f500+lElevation (m)",
                     position="JBC+w10c/0.3c+h+o0c/1.6c")
    else:
        fig.colorbar(frame="xa10f1+lSNR (dB)", position="JMR+w6c/0.3c+o0.9c/-2.4c+v")
        fig.colorbar(cmap=elev_cpt, frame="xa2000f500+lElevation (m)",
                     position="JMR+w6c/0.3c+o0.9c/4.2c+v")

    # 图例符号必须与图上实际符号一致：震中在图上由沙滩球表示，故删去原“红星 Epicenter”条
    legend_spec = ("G 0.1c\n"
                   f"S 0.3c t {STATION_SIZE}c gray30 {STATION_PEN} 0.8c Station\n"
                   "G 0.15c\n")
    # 红三角条目只在真的有突出台站时才写（HIGHLIGHT_STA 为空时图例里不该出现台站名）
    if HIGHLIGHT_STA:
        legend_spec += (f"S 0.3c t {HL_SIZE}c red {HL_PEN} 0.8c "
                        f"Station {','.join(HIGHLIGHT_STA)}\n"
                        "G 0.15c\n")
    legend_spec += "S 0.3c - 0.6c - 0.8p,black,-- 0.8c -1000 m depth\n"
    lf = tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False)
    lf.write(legend_spec)
    lf.close()
    fig.legend(spec=lf.name, position="jBR+w4.3c+o0.55c/0.55c", box="+gwhite+p0.5p,black")
    os.unlink(lf.name)

    fig.text(x=EVLO + 1, y=EVLA + 1, text="-1000 m", font="8p,Helvetica,black", angle=-30)

    out_file = MAP_DIR / out_name
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=DPI)
    print(f"OK {out_file}  (PNG {DPI} dpi)")
    if PDF_OUT:
        pdf_file = out_file.with_suffix(".pdf")
        fig.savefig(pdf_file)
        print(f"OK {pdf_file}  (vector PDF)")


def main():
    parser = argparse.ArgumentParser(description="第一步-台站表格与分布图")
    parser.add_argument("--map-only", action="store_true", help="只画分布图（读已有 stations.csv）")
    parser.add_argument("--table-only", action="store_true", help="只更新 stations.csv/xlsx（不画分布图）")
    args = parser.parse_args()

    if args.map_only:
        rows = []
        with (OUT_DIR / "stations.csv").open() as f:
            for r in csv.DictReader(line for line in f if not line.startswith("#")):
                rows.append({k: (float(v) if v not in ('', None) and k not in ('sta',) else v)
                             for k, v in r.items()})
    else:
        rows = compute_station_rows()
        write_station_table(rows)

    if args.table_only:
        print("表格已更新（未画分布图）")
        return

    # 默认只出一版：全站标注（LABEL_ONLY = [] 表示全标，HIGHLIGHT_STA = [] 不突出某台）
    plot_station_map(rows, "all", "station_map_all_names.png", label_only=LABEL_ONLY)
    # 需要"只标 SNR 前 3 台"那版时：MAP_TOPN=1 python step1_map_and_table.py
    if _os.environ.get("MAP_TOPN", "0") == "1":
        plot_station_map(rows, "top3", "station_map_top3.png", label_only=None)
    print("第一步完成：表格 + 分布图（全站标注）")


if __name__ == "__main__":
    main()
