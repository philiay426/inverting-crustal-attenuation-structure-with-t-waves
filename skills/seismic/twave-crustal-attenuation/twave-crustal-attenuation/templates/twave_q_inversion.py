#!/usr/bin/env python3
# =============================================================================
# T 波地壳衰减 Q0/alpha 反演引擎（尾波归一化 CNM）
# -----------------------------------------------------------------------------
# 适用：事件目录已完成
#         <process_flow>/0_data_processing      (去响应速度 SAC, m/s)
#         <process_flow>/1_station_distribution  (stations.csv)
#       本脚本只读这两处 + 写 4_attenuation/ 的 CSV/PNG。
#
# 公式链（待求量在左，出处见 references/q-inversion-math.md）：
#        Q(f)  = -pi*f / (V * k)                                    ...(1)
#        k     = d/dr { ln[ A_d * Z(r)^-1 / A_c ] }   最小二乘斜率   ...(2)
#        Z(r)  = r            (r <= h_m)                            几何扩散倒数
#              = sqrt(r*h_m)  (r >  h_m)     [GEOM_MODE='paper']
#              = sqrt(r*70)                  [GEOM_MODE='his' 导师口径]
#        A_d   = 信号窗 [T_ref, T_ref+20] 内的 rms
#        A_c   = 尾波窗 [t_coda, t_coda+20] 内的 rms
#        Q(f)  = Q0 * f^alpha          log-log 不加权回归
#        误差  = 回归系数 95% 置信区间 / 1.96  (= 1 sigma)
#        r     = stations.csv 的 land_km（陆段距离，不是震中距）
#
# 环境变量（Q_CASE 必填，其余可选）：
#   Q_CASE=<事件目录>          默认见下方 CASE
#   Q_STAS="A,B,C,..."         台站名单（逗号分隔；优先于 Q_STATION_SET）
#   Q_EXCLUDE="BITY,MMBE"      剔除台站（仅振幅异常/坏台，见 references/qc-criteria.md）
#   Q_WAVE=S|P                 S=Qs(切向 T 分量) | P=Qp(垂直 Z 分量)
#   Q_VS=3.6  Q_VP=6.0         陆段速度 km/s（同时决定到时与式(1) 的 V）
#   Q_CODA_MODES="2S,T-S+400"  尾波口径列表（拟合前写定，不是挑结果）
#   Q_CODA_MAIN="2S"           主口径（决定 A_d/A_c 列与剖面窗）
#   Q_CODA_FACTOR=2.0          '2S' 模式倍率: t_coda = FACTOR * tt_s
#   Q_SIGN=abs|neg             abs(定案口径) = Q 用 |k|；neg(严格) = 斜率必须为负
#   Q_KMIN=1.0                 频点判据 |k|/sigma_k >= KMIN 才进幂律拟合（0=只报告不剔除）
#   Q_FREQ_EXCL="2.25,2.75"    人工剔除频点 Hz（须有判据支撑，见 qc-criteria.md）
#   Q_BANDS="2.0-3.0,3.0-4.0"  宽带替代逐点窄带（默认逐点 fc +-0.1 Hz）
#   Q_TAG=south13              数据标签 -> q_south13_detail.csv / _byfreq.csv
#   Q_FIGNAME=Q随频率_south13   图名（不含扩展名）
#   Q_OUT=<目录>               产物目录（默认 <process_flow>/4_attenuation）
#
# 用法：
#   Q_CASE=<事件目录> Q_STAS="A,B,C,..." Q_CODA_MODES="2S" Q_CODA_MAIN="2S" \\
#   Q_FREQ_EXCL="2.25,2.75" Q_TAG=<标签> python twave_q_inversion.py
#   产物：<Q_OUT>/数据/q_<标签>_detail.csv、q_<标签>_byfreq.csv + 终端打印 Q0/alpha
# =============================================================================
import csv, glob, os, math, re
import numpy as np
from obspy import read
from obspy.geodetics import gps2dist_azimuth
from obspy.signal.rotate import rotate_ne_rt
from scipy.signal import butter, filtfilt, hilbert
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================== 路径 ==============================
CASE = os.environ.get("Q_CASE", "").strip()
if not CASE:
    raise SystemExit("请用环境变量 Q_CASE 指定事件目录（其下有 process_flow/ 与 0_data_processing/）")
PROC = os.path.join(CASE, "process_flow")
STACV = os.path.join(PROC, "1_station_distribution", "stations.csv")
# 速度(m/s) 去响应数据目录（论文口径）；两个候选路径按存在性自动探测
_RESP_CANDS = [os.path.join(PROC, "0_data_processing", "waveform_image", "remove_response_vel"),
               os.path.join(PROC, "0_data_processing", "remove_response_vel"),
               os.path.join(PROC, "0_data_processing", "remove_response"),
               os.path.join(CASE, "sacprocess", "remove_response")]   # 旧布局（2026-09-28 前）
RESP = os.environ.get("Q_RESP") or next((p for p in _RESP_CANDS if os.path.isdir(p)), _RESP_CANDS[0])
OUT = os.environ.get("Q_OUT") or os.path.join(PROC, "4_attenuation")
OUT_DATA = os.path.join(OUT, "数据")          # 表放 数据/，图放组目录（用户约定）
print("数据目录: %s" % RESP)
print("产物目录: %s" % OUT)

# ---------------------------- 读台站表 ----------------------------
# stations.csv 前三行是注释（第 1 行含事件坐标），正文按列名取值，禁固定索引
if not os.path.isfile(STACV):
    raise SystemExit("缺 stations.csv: %s（先跑 1_station_distribution）" % STACV)
_lines = open(STACV).read().splitlines()
_hdr = [l for l in _lines if l.startswith("sta,")]
if not _hdr:
    raise SystemExit("stations.csv 缺 'sta,...' 表头行")
_i0 = _lines.index(_hdr[0])
rows_all = {r["sta"]: r for r in csv.DictReader(_lines[_i0:]) if r.get("sta")}
_head = "\n".join(_lines[:_i0])

# 事件坐标：从注释行 "# Event: <日期> Mw <震级> | <lat>N, <lon>E | Depth <d>km | ..." 解析；
# 也支持 12.3S / 45.6W 这类半球后缀写法。解析失败则用 Q_EVLAT / Q_EVLON。
def _ev_from_head(txt):
    m = re.search(r"(-?\d+(?:\.\d+)?)\s*([NS])\s*,\s*(-?\d+(?:\.\d+)?)\s*([EW])", txt)
    if not m:
        return None
    la, ls, lo, lw = m.groups()
    # 数字带显式负号时以数字为准（表头常把 -34.887 写成 "-34.887N"），
    # 只有数字本身无符号时才用 N/S/E/W 定半球
    lat = float(la) if la.startswith("-") else abs(float(la)) * (1.0 if ls == "N" else -1.0)
    lon = float(lo) if lo.startswith("-") else abs(float(lo)) * (1.0 if lw == "E" else -1.0)
    return lat, lon

_ev = _ev_from_head(_head)
if _ev is None:
    if "Q_EVLAT" in os.environ and "Q_EVLON" in os.environ:
        _ev = (float(os.environ["Q_EVLAT"]), float(os.environ["Q_EVLON"]))
    else:
        raise SystemExit("无法从 stations.csv 注释行解析事件坐标，请设 Q_EVLAT/Q_EVLON")
EV_LAT, EV_LON = _ev
print("事件坐标: %.4f, %.4f（来源: %s）"
      % (EV_LAT, EV_LON, "stations.csv 注释行" if _ev_from_head(_head) else "环境变量"))

# ============================== ★ 参数区 ==============================
GROUP = "grp"                             # 日志标签（仅打印用）
STATION_SET = os.environ.get("Q_STATION_SET", "all")
# 台站组解析优先级: Q_STAS 名单 > Q_STATION_SET ('all')；两者都没有则报错
EXCLUDE = [s.strip() for s in os.environ.get("Q_EXCLUDE", "").split(",") if s.strip()]
# 台站池怎么来：先按 SNR + 波形目视筛（见 qc-criteria.md），再用 select_station_set.py 枚举子集看稳定性。
# 劣质台的剔除理由写进事件目录的 station_qc_notes.md，别只写在代码里。
V_OCEAN, V_P, V_S = 1.48, 6.0, 3.6        # km/s（海段 SOFAR / 陆段 P / 陆段 S）
V_S = float(os.environ.get("Q_VS", V_S))  # 视速度检验用；同时改到时与式(1) 的 V
V_P = float(os.environ.get("Q_VP", V_P))
WAVE_TYPE = os.environ.get("Q_WAVE", "S").upper()
assert WAVE_TYPE in ("S", "P"), "Q_WAVE 只能是 S 或 P"
V_WAVE = V_P if WAVE_TYPE == "P" else V_S

F_C = np.arange(2.0, 4.001, 0.25)         # 9 个中心频率 2.00..4.00（Zhou2021 SI）
HALF_BW = 0.1                             # +-0.1 Hz 窄带
_BANDS = [b.strip() for b in os.environ.get("Q_BANDS", "").split(",") if b.strip()]
BANDS = [tuple(float(v) for v in b.split("-")) for b in _BANDS]
FC_LIST = [0.5 * (f1 + f2) for f1, f2 in BANDS] if BANDS else list(F_C)
BAND_OF = {0.5 * (f1 + f2): (f1, f2) for f1, f2 in BANDS}
FREQ_EXCL = [float(v) for v in os.environ.get("Q_FREQ_EXCL", "").split(",") if v.strip()]

WIN = 20.0                                # 窗长 s（信号窗 = 尾波窗 = 20 s）
SIG_WIN = float(os.environ.get("Q_SIG_WIN", WIN))    # 信号窗长（定案 20）
SIG_OFF = float(os.environ.get("Q_SIG_OFF", 0.0))    # 信号窗起点偏移（定案 0）
CODA_WIN = WIN
MODE_WIN = {}                             # 逐口径窗长覆盖（可选）
SNR_MIN = float(os.environ.get("Q_SNR_MIN", 1.0))    # A_d/A_c 下限
K_MIN = float(os.environ.get("Q_KMIN", 0.0))         # |k|/sigma_k 下限（0=只报告不剔除）
K_SIGN = os.environ.get("Q_SIGN", "abs")             # 'abs'=定案口径(用 |k|) | 'neg'=要求 k<0
assert K_SIGN in ("abs", "neg"), "Q_SIGN 只能是 abs 或 neg"
H_M = 80.0                                # 几何扩散特征长度 = 2 x 35 km（CRUST1.0）
ORIGIN_FALLBACK = 600.0                   # SAC o 字段异常时的兜底（数据从 -600 s 起）
CODA_FACTOR = float(os.environ.get("Q_CODA_FACTOR", 2.0))
CODA_ANCHOR = "T-S"                       # 锚点: 'T-S'|'S'|'T-Sm'|'Sm'（逐台/中位数）
CODA_DT = 400.0                           # '<锚点>+<延迟>' 口径的默认延迟
CODA_ABS_WINS = [1200.0]                  # 'abs<t0>' 口径候选
CODA_MODES = [m.strip() for m in os.environ.get("Q_CODA_MODES", "2S").split(",") if m.strip()]
CODA_MAIN = os.environ.get("Q_CODA_MAIN", CODA_MODES[0]).strip()
DIRECT_MODES = ["theory"]                 # 'theory' = [理论 T-S, +WIN]（定案）
SEARCH_A, SEARCH_B = -30.0, 60.0          # 'peak' 模式的峰位搜索区间
PEAK_HALF = 10.0
AMP_STAT = os.environ.get("Q_AMP", "rms") # 'rms'(导师口径) | 'max'
NB_ORDER = 4                              # 4 阶 Butterworth + filtfilt（SAC bp n 4 p 2）
GEOM_MODE = os.environ.get("Q_GEOM", "his")   # 'his' = sqrt(r*70) | 'paper' = 分段
T4_OFFSET, T4_V = 12.60, 4.518            # 导师 T-S_local = T5+12.60+land/4.518（'his' 直达窗）
DPI_FIG = 300
# 量纲: 目录名带 _vel / 文件名带 _resp_vel = 已是速度(m/s)；否则视为位移(m)，需频域微分
if "Q_TO_VEL" in os.environ:
    TO_VELOCITY = os.environ["Q_TO_VEL"] == "1"
else:
    TO_VELOCITY = "vel" not in os.path.basename(RESP.rstrip("/"))
CODA_T0 = 340.0                           # 'fixed' 口径起点
PRE_BAND = None                           # 窄带前预滤波（导师流程无，故 None）
FITS_MODES = []                           # 非空则额外出逐频点 3x3 斜率拟合图
QF_SHOW = []                              # Q(f) 图要画的组合（下方按 CODA_MODES 生成）
PLOT_CODA_QC = os.environ.get("Q_PLOT_QC", "0") == "1"  # 尾波窗落点 QC 图
QF_BY_DIRECT = True
# 导师口径 T-S 直达窗偏移（'his' 专用），与 coda_t0_of 的 'land2ts' 配套
_ = T4_OFFSET, T4_V

# ---------------------------- 台站组解析 ----------------------------
if os.environ.get("Q_STAS"):
    STAS = [s.strip() for s in os.environ["Q_STAS"].split(",") if s.strip()]
    _miss = [s for s in STAS if s not in rows_all]
    if _miss:
        raise SystemExit("Q_STAS 里有 stations.csv 中没有的台站: %s" % _miss)
    STATION_SET = "list"
else:                                     # 'all' = stations.csv 全部 - EXCLUDE
    STAS = [s for s in rows_all if s not in EXCLUDE]
STAS = [s for s in STAS if s not in EXCLUDE]
if not STAS:
    raise SystemExit("台站组为空（检查 stations.csv / Q_STAS / Q_EXCLUDE）")

RUN_TAG = os.environ.get("Q_TAG") or (STATION_SET if STATION_SET != "list" else "list")
if WAVE_TYPE == "P":
    RUN_TAG = RUN_TAG + "_p"
FIG_STEM = (os.environ.get("Q_FIGNAME") or f"Q随频率_{RUN_TAG}").strip()

os.makedirs(OUT, exist_ok=True)
os.makedirs(OUT_DATA, exist_ok=True)

rows = {s: rows_all[s] for s in STAS}
print(f"=== 尾波归一化 Q{WAVE_TYPE} 反演 | 台集 {STATION_SET} ({len(rows)} 台) ===")
print("台站: " + ", ".join(f"{s}(land={float(rows[s]['land_km']):g}km)" for s in STAS))
R_KM = {s: float(rows[s]["land_km"]) for s in STAS}
print("衰减距离 r = land_km: 跨度 %.0f-%.0f km" % (min(R_KM.values()), max(R_KM.values())))

tt_s_med = float(np.median([float(rows[s]["tt_s"]) for s in STAS]))
TT_TS = {s: float(rows[s]["ocean_km"]) / V_OCEAN + float(rows[s]["land_km"]) / V_WAVE for s in STAS}
tt_ts_med = float(np.median([TT_TS[s] for s in STAS]))
print("理论 T-S 到时 %.0f ~ %.0f s | S 到时中位 %.1f s"
      % (min(TT_TS.values()), max(TT_TS.values()), tt_s_med))
print("尾波口径 %s（主 = %s）| 信号窗 %g s (偏移 %+g) | 频点 %s"
      % (CODA_MODES, CODA_MAIN, SIG_WIN, SIG_OFF,
         ("宽带 " + str(BANDS)) if BANDS else "fc+-%g Hz x %d" % (HALF_BW, len(FC_LIST))))
print("振幅量 %s | 几何扩散 %s | H_M=%g | A_d/A_c>=%g | 斜率符号 %s | 频点判据 |k|/sigma_k>=%g%s"
      % (AMP_STAT, GEOM_MODE, H_M, SNR_MIN, K_SIGN, K_MIN,
         (" | 人工剔除 %s Hz" % FREQ_EXCL) if FREQ_EXCL else ""))

# 真实大圆 BAZ（stations.csv 的 baz 列漏 cos(纬)，偏西可达 11 度，不用）
BAZ_TRUE = {}
for s in STAS:
    # gps2dist_azimuth 返回 (dist, az=台站->事件, baz=事件->台站)；旋转要 az=[1]
    _, baz, _ = gps2dist_azimuth(float(rows[s]["lat"]), float(rows[s]["lon"]), EV_LAT, EV_LON)
    BAZ_TRUE[s] = baz
QF_SHOW = [(cm, m) for cm in CODA_MODES for m in DIRECT_MODES]

# 事件坐标自检：旋转方位与 stations.csv 的 baz 列应差不超过几十度；
# 若接近 180 度 = 事件坐标半球/符号写反了（曾因表头 "-34.887N" 解析错踩过）
_d = [abs((BAZ_TRUE[s] - float(rows[s]["baz"]) + 180.0) % 360.0 - 180.0) for s in STAS]
if np.median(_d) > 30.0:
    print("⚠ 旋转方位与 stations.csv baz 列中位差 %.1f 度 -> 检查事件坐标符号/Q_EVLAT/Q_EVLON" % np.median(_d))
else:
    print("旋转方位自检: 与 stations.csv baz 列中位差 %.1f 度 (ok)" % np.median(_d))

def geom_factor(r):
    """几何扩散的**倒数**（把扩散补回去）。
    'his'   : 导师 amp2Q_ecnm.m L28-29 —— hm = 2*thk_crust = 70 km, 一律 sqrt(r*hm)
    'paper' : 用户推导文档第 2/5/6 步 —— r<=h_m 用 r; r>h_m 用 sqrt(r*h_m) (h_m=80)"""
    if GEOM_MODE == 'his':
        return math.sqrt(r * 70.0)
    return r if r <= H_M else math.sqrt(r * H_M)


def coda_t0_of(mode, ocean_km, land_km, tt_s=None, tt_p=None):
    """尾波窗起点（s, 绝对时间 = 相对发震）以及在陆地计时轴上的位置。

    mode 写法:
      '<锚点>+<延迟>' 锚点 = 该台理论到时, 延迟 s: t_coda = ANCHOR + delay
                      (★ 参数化口径: 换事件/换台阵只需改锚点与延迟, 不硬编码绝对时刻)
      'land2ts'       t_coda = L_ocean/V_OCEAN + 2*(L_land/V_S)
                      (已判定不可用, 代码保留以便复现, 见 coda_window_land2ts_trial.md)
      'fixed'         t_coda = CODA_T0 (全台同一绝对窗, 硬编码, 仅作对照)
      'abs<t0>'       t_coda = t0 (全台同一绝对窗 [t0, t0+WIN]; 试口径, 见参数区 CODA_ABS_WINS)
      '2S'            t_coda = 2 x 该台 tt_s (逐台窗; 依导师源码 case '2S' ts2 = O + 2*(T2-O); 用户 2026-10-13 要求用作分母)
      'preP'          t_coda = 该台 tt_p - WIN, 即窗 [tt_p - WIN, tt_p] (直达 P 波之前; 用户 2026-10-13 试)
      'preP-snr'      t_coda = 该台 tt_p - 60, 窗 [tt_p-60, tt_p-10] (长 50 s; = 项目 SNR 噪声窗口径, 2026-10-14 加)

    锚点取值（★ 参数区 CODA_ANCHOR）:
      'T-S'  = L_ocean/V_OCEAN + L_land/V_S  该台理论 T-S 到时（恒速口径: 海段 1.48, 陆段 V_S = 3.6 km/s; 与 stations.csv 的 tt_ts 不同, 后者早 0~27.5 s）
      'S'    = 该台 tt_s（投影走时理论 S 到时, 调用时需传 tt_s）
      'T-Sm' = 台站集 T-S 到时的中位数（全局常量, 与传入台站无关）
      'Sm'   = 台站集 tt_s 的中位数（全局常量）

    返回 (t_coda, t_coda 在陆地计时轴 t - L_ocean/1.48 上的值; 'fixed' 口径返回 None)
    """
    if mode == 'fixed':
        return (CODA_T0 if CODA_T0 is not None else CODA_FACTOR * tt_s_med), None
    tl_ocean = ocean_km / V_OCEAN                       # 海段耗时 = 陆地计时轴的原点偏移
    if mode.startswith('abs'):                          # 全台同一绝对窗 [t0, t0+WIN] (试口径; 由参数区 CODA_ABS_WINS 生成)
        t0 = float(mode[3:])
        return t0, t0 - tl_ocean
    if mode == '2S':                                    # 该台 2 x tt_s 逐台窗 (导师源码 case '2S')
        t0 = CODA_FACTOR * float(tt_s)
        return t0, t0 - tl_ocean
    if mode == 'preP':                                  # 该台直达 P 波之前 [tt_p - WIN, tt_p]
        t0 = float(tt_p) - WIN
        return t0, t0 - tl_ocean
    if mode == 'preP-snr':                              # ★ 项目 SNR 噪声窗口径 [tt_p-60, tt_p-10] (窗长 50 s)
        t0 = float(tt_p) - 60.0                         #   出处: 1_station_distribution/step1_map_and_table.py
        return t0, t0 - tl_ocean                        #   (SNR_NOISE_LEN=50, SNR_NOISE_GAP=10)
    if mode == 'land2ts':
        tl = 2.0 * land_km / V_S                        # 陆地计时轴上的窗起点
        return tl_ocean + tl, tl
    anchor, dt = mode.split('+', 1)                     # 例: 'Sm+10' -> ('Sm', '10')
    if anchor == 'Sm':
        t0 = tt_s_med                                   # 全台共用（中位数锚点）
    elif anchor == 'T-Sm':
        t0 = tt_ts_med                                  # 全台共用（中位数锚点）
    elif anchor == 'S':
        t0 = float(tt_s)
    else:                                               # 'T-S'
        t0 = tl_ocean + land_km / V_S
    return t0 + float(dt), t0 + float(dt) - tl_ocean


_starts = [coda_t0_of(CODA_MAIN, float(rows[s]["ocean_km"]), float(rows[s]["land_km"]),
                      float(rows[s]["tt_s"]), float(rows[s]["tt_p"]))[0] for s in STAS]
print(f"归一化窗 (主口径) = {CODA_MAIN} -> 逐台起点 [{min(_starts):.1f}, {max(_starts):.1f}] s, 窗长 {WIN:g} s")
print(f"  锚点参考: {len(STAS)} 台 S 到时中位 = {tt_s_med:.1f} s; T-S 到时范围 "
      f"{min(TT_TS.values()):.1f} ~ {max(TT_TS.values()):.1f} s")
print(f"直达窗模式: {DIRECT_MODES}（peak = 包络峰 +/-{PEAK_HALF:g}s, 搜索 "
      f"{SEARCH_A:+.0f}..{SEARCH_B:+.0f} s; theory = 理论 T-S 后 {WIN:g}s）")
print(f"量纲 {('速度 m/s' if not TO_VELOCITY else '位移 m + 频域微分')} | 预滤波 {PRE_BAND} | H_M={H_M:g} km")
print(f"窗口径: 信号窗 {SIG_WIN:g} s (偏移 {SIG_OFF:+g} s) | 尾波窗 {CODA_WIN:g} s"
      + (f" | 频带合并 {BANDS}" if BANDS else " | 逐点窄带 fc+-{HALF_BW:g}"))
print(f"理论 T-S 到时约 {min(float(rows[s]['ocean_km'])/V_OCEAN + float(rows[s]['land_km'])/V_S for s in STAS):.0f}"
      f"~{max(float(rows[s]['ocean_km'])/V_OCEAN + float(rows[s]['land_km'])/V_S for s in STAS):.0f} s")

# 真实大圆 BAZ（stations.csv 的 baz 列偏西约 11 度，不用）
BAZ_TRUE = {}
for s in STAS:
    # gps2dist_azimuth 返回 (dist, az=台站->事件, baz=事件->台站);
    # 旋转要的是台站看向事件的方位 = az = [1]（用 [2] 会整整差 180 度）
    _, baz, _ = gps2dist_azimuth(float(rows[s]["lat"]), float(rows[s]["lon"]), EV_LAT, EV_LON)
    BAZ_TRUE[s] = baz
print("真实大圆 BAZ: " + ", ".join(f"{s}={BAZ_TRUE[s]:.1f}" for s in STAS)
      + " | stations.csv baz 列: " + ", ".join(f"{s}={float(rows[s]['baz']):.1f}" for s in STAS))


def load_three(sta):
    """读该台 Z/N/E 去响应波形 -> (dt, Z, N, E)，长度对齐"""
    got = {}
    for c in ("BHZ", "BHN", "BHE"):
        m = sorted(glob.glob(os.path.join(RESP, f"*.{sta}.*{c}*_resp*.sac")))
        if not m:
            return None
        got[c] = read(m[0])[0]
    n = min(len(t.data) for t in got.values())
    return (got["BHZ"].stats.delta, got["BHZ"].data[:n].astype(float),
            got["BHN"].data[:n].astype(float), got["BHE"].data[:n].astype(float),
            got["BHZ"])


def origin_of(tr):
    """时间基准：SAC o 字段（= 数据起点相对发震的秒数, 数据常从 -600 s 起）"""
    o = tr.stats.sac.get("o", None)
    return abs(float(o)) if (o is not None and 0 < abs(float(o)) < 1e6) else ORIGIN_FALLBACK


def spec_deriv(x, dt):
    """频域微分（位移 -> 速度）: x i*w。论文口径见文件头第 6 条。"""
    n = len(x)
    X = np.fft.rfft(x)
    w = 2.0 * np.pi * np.fft.rfftfreq(n, dt)
    return np.fft.irfft(1j * w * X, n)


def prefilter(data, sr):
    """窄带前预滤波（默认 2-8 Hz）: 窄带完全落在带内, 结果基本不变,
    但可压掉带外（长周期）内容, 减少窄带滤波的带外泄漏。"""
    b, a = butter(3, list(PRE_BAND), btype="band", fs=sr)
    return filtfilt(b, a, data)


def narrowband(data, sr, fc):
    lo, hi = BAND_OF.get(fc, (fc - HALF_BW, fc + HALF_BW))   # ★ 频带合并时用整段带宽
    b, a = butter(NB_ORDER, [lo, hi], btype="band", fs=sr)
    return filtfilt(b, a, data)


def win_amp(data, t0, t1, origin, sr):
    """窗内振幅统计量 (★ 参数区 AMP_STAT)。
    'rms' : 导师 measure_amplitude_ecnm.m L113 sqrt(sum(x^2))/sqrt(sum(y^2)) —— 等长窗下即 RMS 之比
    'max' : 绝对值最大值 (本校本默认口径)"""
    i0 = max(0, int(round((t0 + origin) * sr)))
    i1 = min(len(data), int(round((t1 + origin) * sr)))
    if i1 - i0 < 5:
        return np.nan
    seg = data[i0:i1]
    if AMP_STAT == 'rms':
        return float(np.sqrt(np.mean(seg ** 2)))
    return float(np.max(np.abs(seg)))


def peak_in(data, t0, t1, origin, sr):
    i0 = max(0, int(round((t0 + origin) * sr)))
    i1 = min(len(data), int(round((t1 + origin) * sr)))
    if i1 - i0 < 5:
        return np.nan
    return float(np.max(np.abs(data[i0:i1])))


def direct_amp(Tb, sr, tt_ts, origin, mode, t_sig=None):
    """直达波振幅 (振幅, 用的窗心时刻 t_used)。
    his    : 导师口径 —— T-S_local 固定窗 [t_sig, t_sig+WIN], t_sig = ocean/1.48 + 12.60 + land/4.518 (其 SAC 头段 T4 同构式)
    theory : [tt_ts, tt_ts+WIN]（本项目 T-S 到时口径）
    peak   : 窄带包络峰（搜 [tt-30, tt+60]）为窗心, +/-PEAK_HALF 取最大绝对值"""
    if mode == 'his':
        return win_amp(Tb, t_sig, t_sig + WIN, origin, sr), t_sig
    if mode == 'theory':
        t0 = tt_ts + SIG_OFF
        return win_amp(Tb, t0, t0 + SIG_WIN, origin, sr), t0
    ei = np.abs(hilbert(Tb))
    i0 = max(0, int(round((tt_ts + SEARCH_A + origin) * sr)))
    i1 = min(len(Tb), int(round((tt_ts + SEARCH_B + origin) * sr)))
    if i1 - i0 < 5:
        return np.nan, np.nan
    it = i0 + int(np.argmax(ei[i0:i1]))
    tpk = it / sr - origin
    return peak_in(Tb, tpk - PEAK_HALF, tpk + PEAK_HALF, origin, sr), tpk


# ---------------------------------------------------------------- 主循环
DETAIL = []            # 每台每频点明细
per_freq = {(cm, m): [] for cm in CODA_MODES for m in DIRECT_MODES}   # ★ 键 = (尾波模式, 直达模式)

for fc in FC_LIST:
    st_pts = {(cm, m): [] for cm in CODA_MODES for m in DIRECT_MODES}   # (land, A_d, A_c, t_used, snr, sta)
    for sta in STAS:
        r = rows[sta]
        land = float(r["land_km"]); ocean = float(r["ocean_km"])
        rd = R_KM[sta]                                  # 衰减距离 r = land_km
        three = load_three(sta)
        if three is None:
            continue
        dt, Z, N, E, hdr = three
        sr = 1.0 / dt
        origin = origin_of(hdr)
        if WAVE_TYPE == 'P':                            # ★ Qp: 用垂直分量 Z（P 波近垂直入射）
            T = Z.copy()
        else:                                           # Qs: 用切向分量 T（obspy 旋转约定）
            _, T = rotate_ne_rt(N, E, BAZ_TRUE[sta])
        if TO_VELOCITY:                                 # ★ 位移(m) -> 速度(m/s)（论文口径）
            T = spec_deriv(T, dt)
        if PRE_BAND is not None:                        # ★ 2-8 Hz 预滤波（T 波优势频带）
            T = prefilter(T, sr)
        tt_ts = ocean / V_OCEAN + land / V_WAVE         # 理论到时（Qp = T-P 用 V_P; Qs = T-S 用 V_S）
        t_sig_his = ocean / V_OCEAN + T4_OFFSET + land / T4_V   # ★ 导师口径 T-S_local: T5 + 12.60 + land/4.518
        Tb = narrowband(T, sr, fc)
        rec = dict(sta=sta, f=fc, wave=WAVE_TYPE, land_km=land, r_km=rd, tt_ts=tt_ts)
        for cm in CODA_MODES:                           # ★ 逐尾波口径各取一个 A_c
            tc, tl = coda_t0_of(cm, ocean, land, float(r["tt_s"]), float(r["tt_p"]))   # 窗起点: 锚点+延迟 / fixed 
            _w = MODE_WIN.get(cm, CODA_WIN)                 # ★ 逐口径窗长 (默认 20 s)
            A_c = win_amp(Tb, tc, tc + _w, origin, sr)
            rec[f"A_coda_{cm}"] = A_c; rec[f"t_coda_{cm}"] = tc
            if cm == CODA_MODES[0]:
                rec["A_coda"] = A_c                     # 兼容旧列名（= 主口径）; A_d/A_c 只记主口径
            for m in DIRECT_MODES:
                A_d, t_used = direct_amp(Tb, sr, tt_ts, origin, m, t_sig_his)
                snr = (A_d / A_c) if (np.isfinite(A_d) and np.isfinite(A_c) and A_c > 0) else np.nan
                if cm == CODA_MAIN:
                    rec[f"A_{m}"] = A_d; rec[f"t_{m}"] = t_used; rec[f"snr_{m}"] = snr
                if np.isfinite(A_d) and np.isfinite(A_c) and A_c > 0 and snr >= SNR_MIN:
                    st_pts[(cm, m)].append((rd, A_d, A_c, t_used, snr, sta))   # x = 衰减距离 rd
        DETAIL.append(rec)

    for cm in CODA_MODES:
        for m in DIRECT_MODES:
            pts = st_pts[(cm, m)]
            _used = ""
            if len(pts) >= 3:
                x = np.array([p[0] for p in pts]); y = np.array([math.log(p[1] * geom_factor(p[0]) / p[2]) for p in pts])
                _used = " ".join(p[5] for p in pts)
                k, a = np.polyfit(x, y, 1)
                # 斜率标准误 -> 频点判据 |k|/sigma_k（与 0 分不开的频点不带衰减信息）
                _res = y - np.polyval([k, a], x)
                _s2 = float(_res @ _res) / max(len(x) - 2, 1)
                _sxx = float(((x - x.mean()) ** 2).sum())
                se_k = math.sqrt(_s2 / _sxx) if _sxx > 0 else float("nan")
                # 斜率符号口径: 'abs' = 定案口径 Q = pi*f/(V*|k|)（不要求 k<0）
                #              'neg' = 严格口径, 只接受 k<0 的频点
                if K_SIGN == "abs":
                    Q_doc = math.pi * fc / (abs(k) * V_WAVE) if k != 0 else float("nan")
                else:
                    Q_doc = -math.pi * fc / (k * V_WAVE) if k < 0 else float("nan")
            else:
                x = y = None; k = Q_doc = se_k = float("nan")
            per_freq[(cm, m)].append(dict(fc=fc, Q_doc=Q_doc, k=k, se_k=se_k,
                                          n=len(pts), used=_used, x=x, y=y))
            if m == DIRECT_MODES[0] and cm == CODA_MAIN:
                print(f"\n--- f={fc:.2f} Hz   coda={cm} 直通台站 ---")
                if not pts:
                    print("    (无台站过 SNR 判据)")
                for p in pts:
                    print(f"    {p[5]:5s} r={p[0]:4.0f} A_dir={p[1]:.3e} A_coda={p[2]:.3e} "
                          f"SNR={p[4]:6.2f} t_win={p[3]:.1f}s")


def _freq_ok(c, apply_qc=True):
    """频点入选判据（拟合前写定，不是按结果挑）:
       1) Q 有定义且为正
       2) 不在人工剔除名单 FREQ_EXCL 里
       3) apply_qc 且 K_MIN>0 时要求 |k|/sigma_k >= K_MIN
          (K_MIN=0 表示不自动剔除, 只把判据打印出来由人决定)"""
    if not (np.isfinite(c["Q_doc"]) and c["Q_doc"] > 0):
        return False
    if any(abs(c["fc"] - fe) < 1e-9 for fe in FREQ_EXCL):
        return False
    if apply_qc and K_MIN > 0:
        sk = c.get("se_k", float("nan"))
        if not (np.isfinite(sk) and sk > 0 and abs(c["k"]) / sk >= K_MIN):
            return False
    return True


def fit_q(ms, apply_qc=True):
    """Q(f) = Q0 f^alpha 的不加权 log-log 最小二乘; 返回 (Q0, alpha, f 数组, Q 数组)"""
    sel = [c for c in ms if _freq_ok(c, apply_qc)]
    f = np.array([c["fc"] for c in sel])
    Q = np.array([c["Q_doc"] for c in sel])
    if len(f) >= 3:
        al, b = np.polyfit(np.log(f), np.log(Q), 1)
        return math.exp(b), al, f, Q
    return float("nan"), float("nan"), f, Q


def fit_q_report(ms, apply_qc=True):
    """拟合 + 1 sigma 误差: 回归系数 95% 置信区间 / 1.96 (导师 plot_alpha.m 口径)"""
    Q0, al, f, Q = fit_q(ms, apply_qc=apply_qc)
    if not (np.isfinite(Q0) and np.isfinite(al) and len(f) >= 3):
        return Q0, al, float("nan"), float("nan"), f, Q
    _p, _cov = np.polyfit(np.log(f), np.log(Q), 1, cov=True)
    _se = np.sqrt(np.diag(_cov))                       # [se_slope, se_intercept]
    _t = stats.t.ppf(0.975, len(f) - 2) / 1.96
    return Q0, al, Q0 * _t * _se[1], _t * _se[0], f, Q


print("\n================ 逐频点结果 ================")
for _cm in (FITS_MODES or CODA_MODES):
    for m in DIRECT_MODES:
        print(f"[coda={_cm}  direct={m}]")
        for c in per_freq[(_cm, m)]:
            qs = f"Q={c['Q_doc']:7.1f}" if np.isfinite(c["Q_doc"]) else "Q=    nan"
            _sk = c.get("se_k", float("nan"))
            _ks = (abs(c["k"]) / _sk) if (np.isfinite(_sk) and _sk > 0) else float("nan")
            _tag = "" if _freq_ok(c, True) else "  <-- 不入选"
            print(f"  f={c['fc']:.2f} Hz  N={c['n']:2d}  {qs}  k={c['k']:+.3e}"
                  f"  se_k={_sk:.2e}  |k|/se_k={_ks:4.2f}{_tag}")

print("\n================ Q(f) = Q0 * f^alpha (f0 = 1 Hz) ================")
FIT = {}
for _cm in (FITS_MODES or CODA_MODES):
    for m in DIRECT_MODES:
        for _qc in (True, False):
            Q0, al, sQ0, sal, ff, Q = fit_q_report(per_freq[(_cm, m)], apply_qc=_qc)
            if _qc:
                FIT[(_cm, m)] = (Q0, al)
            _lab = "过判据" if _qc else "全频点"
            tag = f"[coda={_cm:8s} direct={m:6s} {_lab}]"
            if np.isfinite(al):
                print(f"  {tag} Q_{WAVE_TYPE} = ({Q0:.0f} +- {sQ0:.0f}) f^({al:+.2f} +- {sal:.2f})"
                      f"   (拟合 {len(ff)} 个频点, f = {[float(v) for v in np.round(ff, 2)]})")
            else:
                print(f"  {tag} 可拟合频点不足 (<3 个合格频点)")
_bad = []
for _cm in (FITS_MODES or CODA_MODES):
    for m in DIRECT_MODES:
        for c in per_freq[(_cm, m)]:
            _sk = c.get("se_k", float("nan"))
            if np.isfinite(_sk) and _sk > 0 and abs(c["k"]) / _sk < 1.0 and np.isfinite(c["Q_doc"]):
                _bad.append("%.2f Hz(%.2f)" % (c["fc"], abs(c["k"]) / _sk))
if _bad:
    print("  * |k|/sigma_k < 1 (衰减率与 0 不可区分): " + ", ".join(sorted(set(_bad)))
          + "  -> 用 Q_FREQ_EXCL 显式剔除并在结果里照实报告, 或设 Q_KMIN=1 自动剔除")
print("  参考: Zhou et al. (2021) Q_P = 204 f^1.48, Q_S = 685 f^0.53")

# ---------------------------------------------------------------- 输出
csvp = os.path.join(OUT_DATA, f"q_{RUN_TAG}_detail.csv")
with open(csvp, "w", newline="") as f:
    w = csv.writer(f)
    head = ["sta", "f_hz", "land_km", "r_km", "tt_ts_s"]
    for cm in CODA_MODES:
        head += [f"A_coda_{cm}", f"t_coda_{cm}_s"]
    for m in DIRECT_MODES:
        head += [f"A_{m}", f"t_{m}_s", f"snr_{m}"]
    w.writerow(head)
    for d in DETAIL:
        row = [d["sta"], f"{d['f']:.2f}", f"{d['land_km']:g}", f"{d['r_km']:.1f}", f"{d['tt_ts']:.1f}"]
        for cm in CODA_MODES:
            row += [f"{d['A_coda_' + cm]:.6g}", f"{d['t_coda_' + cm]:.1f}"]
        for m in DIRECT_MODES:
            row += [f"{d['A_' + m]:.6g}", f"{d['t_' + m]:.1f}", f"{d['snr_' + m]:.3f}"]
        w.writerow(row)
csvq = os.path.join(OUT_DATA, f"q_{RUN_TAG}_byfreq.csv")
with open(csvq, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow([f"# A_<mode>/snr_<mode> 列均按主口径尾波窗 coda={CODA_MAIN} 计算; 模式名见列名"])
    w.writerow(["coda_mode", "direct_mode", "f_hz", "Q",
                "slope_k", "se_k", "abs_k_over_se_k", "used_in_fit", "n_used", "stations_used"])
    for cm in CODA_MODES:
        for m in DIRECT_MODES:
            for c in per_freq[(cm, m)]:
                _sk = c.get("se_k", float("nan"))
                w.writerow([cm, m, f"{c['fc']:.2f}", f"{c['Q_doc']:.2f}",
                            f"{c['k']:.6e}", f"{_sk:.6e}",
                            (f"{abs(c['k']) / _sk:.3f}" if (np.isfinite(_sk) and _sk > 0) else ""),
                            int(_freq_ok(c, True)), c["n"], c["used"]])
print(f"\n写: {csvp}\n写: {csvq}")

# ---------------------------------------------------------------- 图
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.linewidth": 0.6,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.labelsize": 6, "ytick.labelsize": 6,
                     "axes.labelsize": 6.5, "axes.titlesize": 6.5})
# 配色: 按 CODA_MODES 顺序循环取色; 实线 = 主直达口径, 虚线 = 其他直达口径
_CYC = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e"]
COLS = {(cm, m): _CYC[(i * len(DIRECT_MODES) + j) % len(_CYC)]
        for i, cm in enumerate(CODA_MODES) for j, m in enumerate(DIRECT_MODES)}


def _slug(s):
    return s.replace("+", "p").replace("-", "m").replace(" ", "")

# (1) FITS_MODES 里每个尾波口径 x 主直达口径 一张 3x3 斜率拟合图
for cm in FITS_MODES:
    ms = per_freq[(cm, DIRECT_MODES[0])]
    fig, axes = plt.subplots(3, 3, figsize=(6.6, 6.0))
    for i, c in enumerate(ms):
        ax = axes.ravel()[i]
        if c["x"] is not None and len(c["x"]):
            yy = c["y"]
            ax.plot(c["x"], yy, "o", ms=3.5, mfc="none", mec="#1f77b4", mew=0.9)
            if np.isfinite(c["k"]):
                xs = np.linspace(c["x"].min(), c["x"].max(), 10)
                ax.plot(xs, np.polyfit(c["x"], yy, 1)[1] + c["k"] * xs, "-", lw=0.8, color="k")
        ttl = f"{c['fc']:.2f} Hz   Q={c['Q_doc']:.0f}" if np.isfinite(c["Q_doc"]) else f"{c['fc']:.2f} Hz   Q=nan"
        ax.set_title(ttl)
        ax.tick_params(top=True, right=True, width=0.6, length=2.5)
        for s in ax.spines.values():
            s.set_linewidth(0.6)
        if i == 7:
            ax.set_xlabel("land distance r (km)")
        if i == 3:
            ax.set_ylabel(r"ln[ $A_{dir}\cdot Z(r)^{-1}/A_{coda}$ ]")
    fig.suptitle(f"{STATION_SET} $Q_{WAVE_TYPE}$ coda-normalization fits (coda = {cm}, direct = {DIRECT_MODES[0]})", fontsize=7)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fp1 = os.path.join(OUT, f"{FIG_STEM}_逐频点拟合.png")
    fig.savefig(fp1, dpi=DPI_FIG); print(f"写: {fp1}")

# (2) Q(f) 对照: QF_BY_DIRECT=True 时按直达口径分左右面板 (左 = DIRECT_MODES[0] = 理论到时; 右 = peak)
#     每个面板内 3 条 = CODA_MODES 三档尾波窗 (2S / S+10s / 1200s), 均配论文 Q_S 参考线
_fl = np.array([2.0, 4.0])
# ★ 参考线按波型: Qp 对 Q_P = 204 f^1.48, Qs 对 Q_S = 685 f^0.53 (Zhou et al. 2021)
REF_Q0, REF_AL = (204.0, 1.48) if WAVE_TYPE == 'P' else (685.0, 0.53)
REF_LABEL = (r"Zhou2021: $204f^{1.48}$" if WAVE_TYPE == 'P' else r"Zhou2021: $685f^{0.53}$")


def _q_curve(ax, cm, m):
    """一个 (尾波口径, 直达口径) 组合: 散点 + 幂律拟合线 + 图例"""
    ff = np.array([c["fc"] for c in per_freq[(cm, m)] if np.isfinite(c["Q_doc"]) and c["Q_doc"] > 0])
    QQ = np.array([c["Q_doc"] for c in per_freq[(cm, m)] if np.isfinite(c["Q_doc"]) and c["Q_doc"] > 0])
    if len(ff) == 0:
        print(f"  提示: coda={cm} direct={m} 无有效频点, 不进 Q(f) 图")
        return
    Q0, al, _, _ = fit_q(per_freq[(cm, m)])
    ax.loglog(ff, QQ, "o", ms=3.5, mfc="none", mec=COLS[(cm, m)], mew=0.9)
    ffl = np.linspace(ff.min(), ff.max(), 50)
    ax.loglog(ffl, Q0 * ffl ** al, "-", lw=0.9, color=COLS[(cm, m)],
              label=fr"{cm}: $Q_{WAVE_TYPE}={Q0:.0f}f^{{{al:.2f}}}$")


if QF_BY_DIRECT and len(DIRECT_MODES) > 1:         # ★ 只剩 1 个直达口径时走下面单面板
    fig2, axes2 = plt.subplots(1, len(DIRECT_MODES), figsize=(7.8, 3.1), sharey=True)
    for j, m in enumerate(DIRECT_MODES):
        ax2 = axes2[j]
        for cm in CODA_MODES:
            _q_curve(ax2, cm, m)
        ax2.loglog(_fl, REF_Q0 * _fl ** REF_AL, ":", lw=1.0, color="0.45", label=REF_LABEL)
        ax2.set_xlabel("frequency (Hz)")
        ax2.set_title(f"direct window: {m}", fontsize=7)
        ax2.legend(fontsize=5.0, loc="upper left", frameon=False)
        ax2.tick_params(top=True, right=True, width=0.6, length=2.5)
        for s in ax2.spines.values():
            s.set_linewidth(0.6)
    axes2[0].set_ylabel(f"$Q_{WAVE_TYPE}$")
    fig2.tight_layout()
else:                                              # (旧) 单面板
    fig2, ax2 = plt.subplots(figsize=(3.9, 3.0))
    for cm, m in QF_SHOW:
        _q_curve(ax2, cm, m)
    ax2.loglog(_fl, REF_Q0 * _fl ** REF_AL, ":", lw=1.0, color="0.45", label=REF_LABEL)
    ax2.set_xlabel("frequency (Hz)"); ax2.set_ylabel(f"$Q_{WAVE_TYPE}$")
    ax2.legend(fontsize=5.0, loc="upper left", frameon=False)
    ax2.tick_params(top=True, right=True, width=0.6, length=2.5)
    for s in ax2.spines.values():
        s.set_linewidth(0.6)
    fig2.tight_layout()
fp2 = os.path.join(OUT, f"{FIG_STEM}.png")
fig2.savefig(fp2, dpi=DPI_FIG); print(f"写: {fp2}")


if PLOT_CODA_QC:
    # (3) 时窗落点 QC 图: 归一化窗 + 直达峰窗画在波形上
    fig3, ax3 = plt.subplots(len(STAS), 1, figsize=(11.0, max(7.6, 0.55 * len(STAS))), sharex=True)
    plt.subplots_adjust(hspace=0.12, left=0.085, right=0.985, top=0.905, bottom=0.075)
    C_MAIN, C_PK, C_ANC, C_TS = "#ff7f0e", "#d62728", "#9467bd", "#2ca02c"
    sta_ord = sorted(STAS, key=lambda s: float(rows[s]["land_km"]))
    _anc_main = CODA_ANCHOR if CODA_ANCHOR in ("Sm", "T-Sm") else None   # 仅中位数锚点口径画锚点线
    _xmax = max(1300.0, max(TT_TS.values()) + 80.0, max(_starts) + WIN + 80.0)
    for i, sta in enumerate(sta_ord):
        r3 = rows[sta]
        land3, ocean3 = float(r3["land_km"]), float(r3["ocean_km"])
        dt3, Z3, N3, E3, hdr3 = load_three(sta)
        sr3 = 1.0 / dt3
        o3 = origin_of(hdr3)
        _, T3 = rotate_ne_rt(N3, E3, BAZ_TRUE[sta])
        if TO_VELOCITY:
            T3 = spec_deriv(T3, dt3)
        bb, aa = butter(3, [2.0, 4.0], btype="band", fs=sr3)     # 2-4 Hz（展示用, 计算用逐频点窄带）
        Td = filtfilt(bb, aa, T3)
        tt3 = np.arange(len(Td)) / sr3 - o3
        m3 = np.abs(Td).max()
        ax3[i].plot(tt3, Td / m3, "-", lw=0.3, color="k")
        tc3, tl3 = coda_t0_of(CODA_MAIN, ocean3, land3, float(r3["tt_s"]), float(r3["tt_p"]))   # 主口径窗起点
        t_ts3 = coda_t0_of("T-S+0", ocean3, land3, float(r3["tt_s"]))[0]              # 该台理论 T-S（参考）
        ax3[i].axvspan(tc3, tc3 + WIN, color=C_MAIN, alpha=0.40, zorder=0,
                       label=(f"normalization window {CODA_MAIN} = [{tc3:.1f}, {tc3 + WIN:.1f}] s" if i == 0 else None))
        _d0 = DIRECT_MODES[0]
        tp3 = next((d["t_" + _d0] for d in DETAIL if d["sta"] == sta), float("nan"))   # 主口径直达窗起点 (f=2 Hz)
        if np.isfinite(tp3):
            ax3[i].axvspan(tp3, tp3 + WIN, color=C_PK, alpha=0.22, zorder=0,
                           label=(f"direct window ({_d0}, 20 s)" if i == 0 else None))
        if _anc_main is not None:
            t_ref = coda_t0_of(f"{_anc_main}+0", ocean3, land3, float(r3["tt_s"]), float(r3["tt_p"]))[0]   # 锚点时刻
            ax3[i].axvline(t_ref, color=C_ANC, ls=(0, (1, 2)), lw=0.7, zorder=3,
                           label=(f"anchor {_anc_main} = {t_ref:.1f} s" if i == 0 else None))
        ax3[i].axvline(t_ts3, color=C_TS, ls=(0, (4, 2)), lw=0.7, zorder=3,
                       label=("this station's theoretical T-S (reference)" if i == 0 else None))
        ax3[i].set_xlim(-100, _xmax); ax3[i].set_ylim(-1.15, 1.15)
        ax3[i].set_ylabel(f"{sta}\nland={land3:.0f} km", fontsize=7)
        ax3[i].tick_params(direction="in", top=True, right=True, labelsize=7, length=2.5)
        for s in ax3[i].spines.values():
            s.set_linewidth(0.6)
        print(f"  {sta:5s} land={land3:4.0f}  主口径窗=[{tc3:.1f},{tc3 + WIN:.1f}]"
              f"  该台 T-S={t_ts3:6.1f}  直达窗={tp3:.0f}")
    ax3[0].legend(fontsize=5.6, loc="upper left", frameon=False, ncol=3,
                  handlelength=1.6, columnspacing=0.9)
    ax3[-1].set_xlabel("time from origin (s)")
    fig3.suptitle(f"{STATION_SET} normalization window on T component (2-4 Hz, per-trace normalized)"
                  f"  |  main = {CODA_MAIN}", fontsize=8)
    fp3 = os.path.join(OUT, f"{FIG_STEM}_尾波窗落点.png")
    fig3.savefig(fp3, dpi=DPI_FIG); print(f"写: {fp3}")

print("DONE_Q")
