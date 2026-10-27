#!/usr/bin/env python3
"""
T波全时段偏振扫描（多台站批量版；一视同仁，不预设 T-P/T-S 窗口）
==================================================================
对每站：在 [T-P理论-100s, T-S理论+500s] 整段T波波列上，以 0.25s 步长逐窗口
做 ±3s PCA（1s 子窗 / 75% 重叠 / analytic / ZNE→ZRT），不区分震相，
记录 偏振主方向 / 径向残差 / 横向残差 / 线性度 / 平面度 / 扇区占比 / 2-8Hz SNR。

图（三面板，共享 x；竖线**只画理论到时**，蓝/红虚线 + 数值标注）:
  (a) 偏振主方向 R-T 角（左轴, deg）+ 径向/横向扇区占比（右轴, %）
  (b) 线性度 rectilinearity（左轴）
  (c) 平面度 planarity（左轴）
⚠️ 线性度/平面度定义与 templates/polarization_attributes.py 完全一致（同一 pca() 输出）:
     linearity = (l1 - l2) / l1, planarity = (l2 - l3) / l1  (l1 >= l2 >= l3, clip 到 [0,1])

方法开关 --method: cov | amp | count（默认 count 采样点等权/投票）
  优先级: --method > 环境变量 POLAR_METHOD > 默认 count。方向+权重见 polar_methods.py
  · cov   = ±3s 窗内 pca() 子窗主轴 R-T 方位, 权重 λ1（与旧版数值一致）
  · amp   = 2-8Hz R/T 逐采样点 atan2(T,R), 权重 R²+T²
  · count = 2-8Hz R/T 逐采样点 atan2(T,R), 每点 1 票
  线性度/平面度三种方法都来自同一次 pca()（属性口径不随方法变）。

用法（改顶部 EV_ORIGIN / BASE 适配其他事件）:
    python full_scan_polarization.py                      # 默认 10 个高SNR台站
    python full_scan_polarization.py VOI BATG KIRI        # 指定台站
    python full_scan_polarization.py AMPY --method cov    # 换偏振方法
    python full_scan_polarization.py AMPY --out-dir /tmp/x --flat   # 输出重定向（烟测）
输出:
    {OUT_ROOT}/{STA}/{STA}_polarization_fullscan.{png,csv}      (默认, 每站一子目录)
    {OUT_ROOT}/{STA}_polarization_fullscan.{png,csv}            (--flat, 单台平铺)

CSV 缓存：若 CSV 已存在、行数/首中心一致、表头含新列且 method 相同则直接读回
→ 改图不重扫，秒出；换方法或改扫描范围/步长会自动重扫。
判读用**连续簇**（≥1s、residual<10°、lin>0.99）而非全局最小窗口，详见
references/window-velocity-sensitivity-scan.md。
"""
from pathlib import Path
import sys, os, glob, csv, argparse
import multiprocessing as mp
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from obspy import read, Stream, UTCDateTime
from obspy.signal.filter import bandpass
from obspy.signal.util import smooth
from obspy.signal.rotate import rotate_ne_rt
from scipy.signal import welch, hilbert

for _codex in ('/home/sy/codex/t_phase_analysis', '/home/sy/codex'):
    if os.path.isdir(_codex):
        sys.path.insert(0, _codex)
from t_phase_analysis.pca import pca

# 共享模块 polar_methods.py 与本脚本同目录(驱动脚本会一起拷到台站目录)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import polar_methods

# ── 每事件需改: EV_ORIGIN / BASE ──
EV_ORIGIN = UTCDateTime(2013, 8, 17, 16, 22, 31)
BASE = Path('/mnt/e/sharedate/case/Madagascar/2013-08-17-mw62-southwest-indian-ridge')
RESP = BASE / 'sacprocess/remove_response'
ROOT = BASE / 'process_flow/2_single_station_analysis'
STATIONS_CSV = BASE / 'process_flow/1_station_distribution/stations.csv'

HALF = 3.0            # 窗口半宽，与搜索方案一致
STEP = 0.25           # 扫描步长
PAD_BEFORE, PAD_AFTER = 100.0, 500.0   # T-P 前 100s ~ T-S 后 500s
DEFAULT_STATIONS = ['AMPY', 'BAND', 'BATG', 'BKTA', 'CPSM', 'FOMA', 'KIRI', 'SOLA', 'VINA', 'VOI']

# CSV 列名（纯文字表头, 不含 -> ｜ @ Δ ° 等符号）
CSV_COLS = ['center_s', 'rt_angle_med_deg', 'linearity_med', 'radial_res_deg', 'transv_res_deg',
            'snr_db', 'n_subwin', 'planarity_med', 'frac_radial', 'frac_transverse', 'method']
# 旧版(无 planarity/扇区/方法列)缓存判定用
LEGACY_COLS = ['center_s', 'rt_angle_med_deg', 'linearity_med', 'radial_res_deg',
               'transv_res_deg', 'snr_db', 'n_subwin']

# 运行时全局(由 main() 填充, 子进程继承)
METHOD = polar_methods.DEFAULT_METHOD
METHOD_LBL = polar_methods.method_label(METHOD)
OUT_ROOT = ROOT
FLAT = False


def station_params(sta):
    with open(STATIONS_CSV) as f:
        for row in csv.DictReader((l for l in f if not l.startswith('#'))):
            if row['sta'] == sta:
                return float(row['tt_tp']), float(row['tt_ts']), float(row['baz']), float(row['ocean_km'])
    raise KeyError(f'{sta} not in stations.csv')


def packet_snr(z, sr, t, center, half=HALF):
    m = (t >= center - half) & (t <= center + half)
    n = (t >= -300) & (t <= -60)
    if m.sum() < 20 or n.sum() < 20:
        return np.nan
    fs, ps = welch(z[m], sr, nperseg=min(256, m.sum()))
    fn, pn = welch(z[n], sr, nperseg=min(256, n.sum()))
    a = ps[(fs >= 2) & (fs <= 8)].sum(); b = pn[(fn >= 2) & (fn <= 8)].sum()
    return 10 * np.log10(a / b) if a > 0 and b > 0 else np.nan


def process_station(sta):
    OUT = OUT_ROOT if FLAT else (OUT_ROOT / sta)
    OUT.mkdir(parents=True, exist_ok=True)
    method = METHOD                      # 子进程里用传进来的全局
    ttp, tts, baz, ocean_km = station_params(sta)
    print(f'[{sta}] start: T-P={ttp:.1f}s T-S={tts:.1f}s baz={baz:.1f} method={method}')

    trs = {}
    for ch in ('BHZ', 'BHN', 'BHE'):
        m = glob.glob(str(RESP / f'*.{sta}.*.{ch}*_resp.sac'))
        if m:
            trs[ch] = read(m[0])[0]
    if len(trs) < 3:
        print(f'[{sta}] missing channel, skip')
        return sta, -1
    t0 = max(tr.stats.starttime for tr in trs.values())
    t1 = min(tr.stats.endtime for tr in trs.values())
    for ch in trs:
        trs[ch].trim(t0, t1)
    delta = trs['BHZ'].stats.delta
    sr = 1 / delta
    sac = trs['BHZ'].stats.sac
    o_val = getattr(sac, 'o', None)
    t_off = (sac.b - o_val) if (o_val is not None and 0 < abs(o_val) < 1e6) else (trs['BHZ'].stats.starttime - EV_ORIGIN)
    t = np.arange(trs['BHZ'].stats.npts) * delta + t_off

    # 三分量 2-8Hz 带通（Z/R/T, obspy 约定）—— 只用于逐采样点方向/权重
    zf = bandpass(trs['BHZ'].data.astype(float), 2, 8, df=sr, corners=4, zerophase=True)
    r_raw, t_raw = rotate_ne_rt(trs['BHN'].data.astype(float), trs['BHE'].data.astype(float), baz)
    rf = bandpass(r_raw, 2, 8, df=sr, corners=4, zerophase=True)
    tf = bandpass(t_raw, 2, 8, df=sr, corners=4, zerophase=True)
    st = Stream([trs['BHZ'], trs['BHN'], trs['BHE']])

    def win_metrics(c):
        """返回一行 CSV 字段; 失败返回 None。线性度/平面度一律来自 pca()（属性口径固定）。"""
        stw = st.slice(starttime=t0 + (c - HALF - t_off), endtime=t0 + (c + HALF - t_off))
        if stw[0].stats.npts < 4 * sr:
            return None
        try:
            r = pca(stw, window_length_s=1.0, window_step_s=0.25,
                    back_azimuth=baz, covariance_method='analytic')
        except Exception:
            return None
        lin = float(np.median(r['linearity']))
        pla = float(np.median(r['planarity']))
        if method == polar_methods.METHOD_COV:
            v = np.asarray(r['principal_vector'])
            if len(v) == 0:
                return None
            ang_deg = np.degrees(np.arctan2(v[:, 2], v[:, 1]))
            w_sub = np.asarray(r['eigenvalues'])[:, 0]
        else:
            m = (t >= c - HALF) & (t <= c + HALF)
            ang_rad, w_sub = polar_methods.sample_direction_weights(method, rf[m], tf[m])
            ang_deg = np.degrees(ang_rad)
        if len(ang_deg) == 0:
            return None
        ang_med = float(np.median(ang_deg))
        d_rad = abs(((ang_med + 90) % 180) - 90)
        sec = polar_methods.sector_fractions(np.radians(ang_deg), w_sub)
        return (c, ang_med, lin, d_rad, abs(d_rad - 90), packet_snr(zf, sr, t, c), len(ang_deg),
                pla, sec['frac_radial'], sec['frac_transverse'], method)

    c_lo, c_hi = ttp - PAD_BEFORE + HALF, tts + PAD_AFTER - HALF
    centers = np.arange(c_lo, c_hi + 1e-6, STEP)

    csv_path = OUT / f'{sta}_polarization_fullscan.csv'
    rows = None
    if csv_path.exists():
        cached = list(csv.DictReader(open(csv_path)))
        hdr = set(cached[0].keys()) if cached else set()
        if not cached:
            rows = None
        elif len(cached) != len(centers) or abs(float(cached[0]['center_s']) - centers[0]) > 1e-6:
            print(f'[{sta}] cache mismatch (rows/center) -> rescan')
        elif not set(CSV_COLS).issubset(hdr):
            print(f'[{sta}] cached CSV has legacy columns -> rescan')
        elif cached[0].get('method', method) != method:
            print(f"[{sta}] cached CSV method={cached[0].get('method')} != {method} -> rescan")
        else:
            print(f'[{sta}] reuse cached CSV ({len(cached)} rows, method={method})')
            def _f(x):
                try:
                    return float(x)
                except (TypeError, ValueError):
                    return np.nan
            rows = [(float(r['center_s']), _f(r['rt_angle_med_deg']), _f(r['linearity_med']),
                     _f(r['radial_res_deg']), _f(r['transv_res_deg']), _f(r['snr_db']),
                     int(float(r['n_subwin'])), _f(r['planarity_med']),
                     _f(r['frac_radial']), _f(r['frac_transverse']), r['method']) for r in cached]
    if rows is None:
        rows = []
        for i, c in enumerate(centers):
            row = win_metrics(c)
            if row is None:
                rows.append((c, np.nan, np.nan, np.nan, np.nan, np.nan, 0,
                             np.nan, np.nan, np.nan, method)); continue
            rows.append(row)
            if (i + 1) % 500 == 0:
                print(f'[{sta}] progress {i+1}/{len(centers)}')
        with open(csv_path, 'w', newline='') as f:
            w = csv.writer(f); w.writerow(CSV_COLS); w.writerows(rows)
        print(f'[{sta}] saved CSV ({len(rows)} rows, method={method})')

    # ── 图: 3 行共享 x —— (a) 主方向+扇区占比  (b) 线性度  (c) 平面度 ──
    plt.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'DejaVu Sans'],
        'font.size': 9, 'axes.labelsize': 10, 'axes.titlesize': 10,
        'xtick.labelsize': 8, 'ytick.labelsize': 8,
        'xtick.direction': 'in', 'ytick.direction': 'in',
    })
    fig, axs = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    X0, X1 = ttp - PAD_BEFORE, tts + PAD_AFTER

    centers_a = np.array([r[0] for r in rows])
    ang_a = np.array([r[1] for r in rows]); lin_a = np.array([r[2] for r in rows])
    pla_a = np.array([r[7] for r in rows])
    fr_a = np.array([r[8] for r in rows]); ft_a = np.array([r[9] for r in rows])
    ok = ~np.isnan(ang_a)

    # (a) 偏振主方向（左轴, R-T 角）+ 扇区占比（右轴, %）
    ax0b = axs[0].twinx()
    ax0b.plot(centers_a[ok], fr_a[ok] * 100, color='#1f77b4', lw=0.7, alpha=0.55, label='R sector')
    ax0b.plot(centers_a[ok], ft_a[ok] * 100, color='#D55E00', lw=0.7, alpha=0.55, label='T sector')
    ax0b.set_ylabel('Sector fraction (%)')
    ax0b.set_ylim(0, 100)
    axs[0].plot(centers_a[ok], ang_a[ok], color='#222222', lw=0.9, label='median R-T angle')
    axs[0].axhline(0, c='green', ls='--', lw=0.8)
    axs[0].axhline(90, c='red', ls='--', lw=0.8)
    axs[0].axhline(-90, c='red', ls='--', lw=0.8)
    axs[0].text(X0 + 3, 4, 'radial (0 deg)', color='green', fontsize=8, fontweight='bold', va='bottom')
    axs[0].text(X0 + 3, 88, 'transverse (+/-90 deg)', color='red', fontsize=8, fontweight='bold', va='bottom')
    axs[0].set_ylim(-180, 180)
    axs[0].set_ylabel('R-T angle (deg)')
    axs[0].set_title(f'{sta} T-wave polarization scan (2-8 Hz, \u00b1{HALF:.0f}s win, method: {METHOD_LBL})',
                     fontsize=11)
    h1, l1 = axs[0].get_legend_handles_labels()
    h2, l2 = ax0b.get_legend_handles_labels()
    axs[0].legend(h1 + h2, l1 + l2, loc='upper right', fontsize=7, ncol=3, framealpha=0.9)

    # (b) 线性度 (rectilinearity) —— 全程 0.99+：纵轴必须放大，否则是一条贴顶直线
    lmin = np.nanmin(lin_a[ok])
    axs[1].plot(centers_a[ok], lin_a[ok], color='#555555', lw=0.9)
    axs[1].set_ylabel('Linearity')
    axs[1].set_ylim(lmin - 0.003, 1.001)
    axs[1].text(X0 + 3, lmin - 0.0025, f'min = {lmin:.3f}', color='gray', fontsize=8, va='bottom')

    # (c) 平面度 (planarity)
    pmin = np.nanmin(pla_a[ok])
    axs[2].plot(centers_a[ok], pla_a[ok], color='#7b3294', lw=0.9)
    axs[2].set_ylabel('Planarity')
    axs[2].set_ylim(pmin - 0.003, 1.001)
    axs[2].text(X0 + 3, pmin - 0.0025, f'min = {pmin:.3f}', color='gray', fontsize=8, va='bottom')
    axs[2].set_xlabel('Time from origin (s)')

    # 竖线只画理论到时（用户定稿：不画未标注的辅助线）
    for ax in axs:
        ax.axvline(ttp, c='#4d9fff', ls='--', lw=0.9)
        ax.axvline(tts, c='#ff4d4d', ls='--', lw=0.9)
    axs[0].text(ttp - 5, 168, f'T-P = {ttp:.1f} s', color='#4d9fff',
                fontsize=9, ha='right', va='center', fontweight='bold')
    axs[0].text(tts + 5, 168, f'T-S = {tts:.1f} s', color='#ff4d4d',
                fontsize=9, ha='left', va='center', fontweight='bold')

    for ax in axs:
        ax.set_xlim(X0, X1)
        ax.grid(alpha=0.08)
    fig.align_ylabels(axs)
    fig.tight_layout()
    fn = OUT / f'{sta}_polarization_fullscan.png'
    try:
        fig.savefig(fn, dpi=200, bbox_inches='tight')
    except OSError as e:
        print(f'[{sta}] PNG save failed ({e}); retry after unlink')
        try:
            fn.unlink()
        except Exception:
            pass
        fig.savefig(fn, dpi=200, bbox_inches='tight')
    plt.close(fig)

    # 理论窗扇区占比（判读用；跨方法只比占比，不比权重和）
    for lbl, tc_ in (('T-P', ttp), ('T-S', tts)):
        mw = np.abs(centers_a - tc_) <= HALF
        if mw.any():
            print(f'[{sta}] {lbl} theory win [{method}]: '
                  f'radial={np.nanmedian(fr_a[mw])*100:.0f}% transverse={np.nanmedian(ft_a[mw])*100:.0f}%')
    print(f'[{sta}] FINISHED -> {fn}')
    return sta, len(rows)


def _worker(job):
    """子进程入口：显式接收 out_root/flat/method。

    🔴 不能用“main() 里设全局变量再让子进程读”的写法：Python 3.14 在 Linux 上默认
    multiprocessing start method 已从 fork 改为 **forkserver**（子进程重新 import 模块，
    拿不到 main() 里赋的全局量）→ 输出目录会静默退回模块级 ROOT（= 事件目录里），
    图被写进 case 目录。2026-09-15 烟测踩到：--out-dir /tmp/... 但文件落在
    <case>/process_flow/2_single_station_analysis/{STA}/。
    """
    global METHOD, METHOD_LBL, OUT_ROOT, FLAT
    sta, out_root, flat, method, method_lbl = job
    OUT_ROOT, FLAT = Path(out_root), bool(flat)
    METHOD, METHOD_LBL = method, method_lbl
    return process_station(sta)


def main():
    global METHOD, METHOD_LBL, OUT_ROOT, FLAT
    p = argparse.ArgumentParser(description='T-wave full-scan polarization (3-panel figure + CSV)')
    p.add_argument('stations', nargs='*', help='台站名（默认 10 个高SNR台站）')
    polar_methods.add_common_args(p)
    p.add_argument('--flat', action='store_true',
                   help='输出直接平铺在 --out-dir 下（单台用；默认每站一个 {STA}/ 子目录）')
    args = p.parse_args()

    METHOD = polar_methods.resolve_method(args.method)
    METHOD_LBL = polar_methods.method_label(METHOD)
    OUT_ROOT = Path(polar_methods.resolve_out_dir(args.out_dir, str(ROOT)))
    FLAT = args.flat

    stations = args.stations or DEFAULT_STATIONS
    nproc = min(len(stations), 8)
    print(f'batch scan: {len(stations)} stations [{", ".join(stations)}], {nproc} processes, '
          f'method={METHOD_LBL}, out_root={OUT_ROOT}{" (flat)" if FLAT else ""}')
    jobs = [(s, str(OUT_ROOT), FLAT, METHOD, METHOD_LBL) for s in stations]
    with mp.Pool(nproc) as pool:
        for sta, n in pool.imap_unordered(_worker, jobs):
            print(f'>>> [{sta}] done, windows={n}')


if __name__ == '__main__':
    main()
