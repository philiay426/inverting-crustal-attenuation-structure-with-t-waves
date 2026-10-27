#!/usr/bin/env python3
"""
单台分析一键驱动 (v5, 2026-09-15)
用法:
  python run_single_station.py --event sw_indian_ridge_2013 --stations B1665,AMPY,W15SS
  python run_single_station.py --event sw_indian_ridge_2013 --stations AMPY --method cov
  python run_single_station.py --event sw_indian_ridge_2013 --stations AMPY --all

默认 (v5) 只出 **3 张图** (输出到 process_flow/2_single_station_analysis/{STA}/):
  0. copy 预处理图 (spectrogram + bandpass)     — 拷贝行为保持现状
  1. polarization_analysis.py     → {STA}_PS_polar.png          (P/S 理论到时 ±2s)
  2. pca_theory_polarization.py   → {STA}_pca_theory_TP_TS.png  (T-P/T-S 理论到时 ±2s)
  3. full_scan_polarization.py    → {STA}_polarization_fullscan.png + .csv
                                    (三面板: 主方向+扇区占比 / 线性度 / 平面度)

--all 再跑旧图: plot_waveform_windows (Twave_3comp) / pca_polarization (pca_TP_TS,
偏振引导拾取) / polarization_attributes (linearity_planarity_3D)

--method cov|amp|count (默认 count 采样点等权) 传给三个生产脚本;
优先级: --method > 环境变量 POLAR_METHOD > 默认 count。方向+权重见 templates/polar_methods.py
(polar_methods.py 会随模板一起拷到台站目录, 三个脚本都 import 它)。

--out-dir / --base-dir: 输出根目录 / 事件数据根目录覆盖 (烟测重定向、同事件多 case 分支副本)。
时间轴铁律：`sac.o` 有效时用 `b-o`；缺 o 时从文件名解析发震时刻并用 `starttime-EV_ORIGIN` fallback。
"""
import os, sys, re, shutil, subprocess, argparse, glob

SKILL = os.path.expanduser("~/.hermes/skills/seismic/single-station-analysis")
sys.path.insert(0, os.path.join(SKILL, "scripts"))
from single_station_common import load_event, CONFIG_PATH

# 默认 3 张图的生产模板 (fullscan 单独调用, 见 FULLSCAN)
TEMPLATES_DEFAULT = ['polarization_analysis', 'pca_theory_polarization']
# --all 才跑的旧图模板
TEMPLATES_LEGACY = ['plot_waveform_windows', 'pca_polarization', 'polarization_attributes']
# 共享模块: 必须跟模板一起拷到台站目录(三个脚本都 import 它)
SHARED_MODULES = ['polar_methods']
FULLSCAN = 'full_scan_polarization'

PYTHON = os.path.expanduser("~/miniforge3/envs/obspy/bin/python")
if not os.path.exists(PYTHON):
    PYTHON = sys.executable


def _stage0_dir(base):
    """步0 数据目录（2026-09-28 规范）：
       新事件 = process_flow/0_data_processing/ ；旧事件（有 sacprocess/）= <event>/sacprocess/"""
    new = os.path.join(base, 'process_flow', '0_data_processing')
    old = os.path.join(base, 'sacprocess')
    return new if (os.path.isdir(new) or not os.path.isdir(old)) else old


def render_template(tpl_name, ev_cfg, sta, net, out_dir):
    """复制模板到台站目录并注入事件/台站参数"""
    src = os.path.join(SKILL, "templates", f"{tpl_name}.py")
    dst = os.path.join(out_dir, f"{tpl_name}.py")
    if not os.path.exists(src):
        print(f"  ⚠ 模板不存在: {tpl_name}.py"); return None
    s = open(src).read()

    # 通用替换: EV / BASE / STA / NET / OUT / CSV
    s = re.sub(r'EVLA, EVLO, EVDP = [^\n]+',
               f"EVLA, EVLO, EVDP = {ev_cfg['lat']}, {ev_cfg['lon']}, {ev_cfg['depth_km']}", s)
    s = re.sub(r'EV_ORIGIN = UTCDateTime\([^)]*\)',
               f"EV_ORIGIN = UTCDateTime({_origin_args(ev_cfg)})", s)
    s = re.sub(r'BASE = "[^"]+"', f'BASE = "{_base_dir(ev_cfg)}"', s)
    # full_scan_polarization.py 用 BASE = Path('...')
    s = re.sub(r"BASE = Path\('[^']+'\)", f"BASE = Path('{_base_dir(ev_cfg)}')", s)
    s = re.sub(r"STA, NET = '[^']+', '[^']+'", f"STA, NET = '{sta}', '{net}'", s)
    # OUT 可能是 f-string 形式
    s = re.sub(r"OUT = f'\{BASE\}/process_flow/2_single_station_analysis/\w+'",
               f"OUT = f'{{BASE}}/process_flow/2_single_station_analysis/{sta}'", s)
    s = re.sub(r'OUT = f"\{BASE\}/process_flow/2_single_station_analysis/\w+"',
               f'OUT = f"{{BASE}}/process_flow/2_single_station_analysis/{sta}"', s)
    s = re.sub(r"OUT = '[^']*2_single_station_analysis/\w+'",
               f"OUT = '{_base_dir(ev_cfg)}/process_flow/2_single_station_analysis/{sta}'", s)
    s = re.sub(r'OUT = "[^"]*2_single_station_analysis/\w+"',
               f'OUT = "{_base_dir(ev_cfg)}/process_flow/2_single_station_analysis/{sta}"', s)
    # CSV 路径
    s = re.sub(r"CSV = '[^']*stations\.csv'",
               f"CSV = '{_base_dir(ev_cfg)}/process_flow/1_station_distribution/stations.csv'", s)
    s = re.sub(r'CSV = "[^"]*stations\.csv"',
               f'CSV = "{_base_dir(ev_cfg)}/process_flow/1_station_distribution/stations.csv"', s)

    # ★ 步0 数据目录（2026-09-28 规范）：新事件 → process_flow/0_data_processing，
    #   旧事件（有 sacprocess/）→ 原样；模板里两种写法都改写
    if _stage0_dir(_base_dir(ev_cfg)) != os.path.join(_base_dir(ev_cfg), 'sacprocess'):
        s = (s.replace('/sacprocess/', '/process_flow/0_data_processing/')
              .replace('"sacprocess"', '"process_flow/0_data_processing"')
              .replace("'sacprocess'", "'process_flow/0_data_processing'"))

    open(dst, 'w').write(s)
    return dst


def copy_shared_modules(out_dir):
    """把 polar_methods.py 拷到台站目录(模板脚本 import 它, 漏拷会 ModuleNotFoundError)"""
    for mod in SHARED_MODULES:
        src = os.path.join(SKILL, "templates", f"{mod}.py")
        if not os.path.exists(src):
            print(f"  ⚠ 共享模块不存在: {mod}.py"); continue
        shutil.copy(src, os.path.join(out_dir, f"{mod}.py"))


def _origin_args(ev_cfg):
    """发震时刻: 优先 config origin_time, 否则从 raw SAC 文件名解析"""
    if ev_cfg.get('origin_time'):
        t = str(ev_cfg['origin_time']).split(':')
        y, mo, da = [int(x) for x in ev_cfg['date'].split('-')]
        return f"{y}, {mo}, {da}, {int(t[0])}, {int(t[1])}, {int(t[2])}"
    # 从 raw 目录 SAC 文件名解析 (year.julday.HHMMSS)
    from single_station_common import parse_origin_from_filename
    import glob as _g
    raw_dir = os.path.join(_base_dir(ev_cfg), "raw")
    for f in _g.glob(f"{raw_dir}/*.BHZ*.SAC")[:20]:
        o = parse_origin_from_filename(os.path.basename(f))
        if o:
            return f"{o.year}, {o.month}, {o.day}, {o.hour}, {o.minute}, {o.second}"
    return "2000, 1, 1, 0, 0, 0"


def _base_dir(ev_cfg):
    """事件实际目录；支持同事件在不同 case 分支的目录覆盖。"""
    if ev_cfg.get('_base_dir_override'):
        return ev_cfg['_base_dir_override']
    return f"/mnt/f/sharedate/case/{ev_cfg['case_dir']}"


def resolve_network(ev_cfg, sta):
    """从 raw 文件名推断台站网络 (XK.B1665... → XK)"""
    raw_dir = os.path.join(_base_dir(ev_cfg), "raw")
    if os.path.isdir(raw_dir):
        import glob
        m = glob.glob(f"{raw_dir}/*.{sta}.*.SAC")
        if m:
            return os.path.basename(m[0]).split('.')[0]
    return 'XX'


def _run(script, args=(), timeout=None):
    r = subprocess.run([PYTHON, script] + list(args), capture_output=True, text=True, timeout=timeout)
    out = (r.stdout + r.stderr).strip().split('\n')
    key = [l for l in out if any(k in l for k in ('=>', '✓', 'OK ', 'Error', 'error', 'Traceback'))]
    for l in key[-4:]:
        print(f"  {l.strip()}")
    return r.returncode


def main():
    p = argparse.ArgumentParser(description='单台分析 v5 一键驱动 (默认 3 张图)')
    p.add_argument('--event', required=True, help='events.yaml 中的事件 key')
    p.add_argument('--stations', required=True, help='逗号分隔台站名')
    p.add_argument('--snr-min', type=float, default=3.0)
    p.add_argument('--method', default=None, choices=['cov', 'amp', 'count'],
                   help='偏振方法: cov | amp | count (默认 count 采样点等权; env POLAR_METHOD 亦可)')
    p.add_argument('--all', action='store_true',
                   help='额外跑旧图: Twave_3comp / pca_TP_TS(引导拾取) / linearity_planarity_3D')
    p.add_argument('--base-dir', default=None, help='覆盖 events.yaml 的事件 BASE (同事件多 case 分支)')
    p.add_argument('--out-dir', default=None,
                   help='输出根目录 (默认 <BASE>/process_flow/2_single_station_analysis)')
    args = p.parse_args()

    ev_cfg = load_event(args.event)
    if args.base_dir:
        ev_cfg = dict(ev_cfg); ev_cfg['_base_dir_override'] = args.base_dir
    base = _base_dir(ev_cfg)
    out_root = args.out_dir or f"{base}/process_flow/2_single_station_analysis"
    stas = [s.strip() for s in args.stations.split(',') if s.strip()]

    method_args = ['--method', args.method] if args.method else []
    method_show = args.method or os.environ.get('POLAR_METHOD') or 'count'

    print(f"事件: {ev_cfg['date']} Mw {ev_cfg['mw']} | 台站: {', '.join(stas)} | "
          f"SNR≥{args.snr_min} | method={method_show}")
    print(f"BASE: {base}\nOUT ROOT: {out_root}" + ("   (--all: 额外旧图)" if args.all else ""))
    print(f"步0 数据目录: {_stage0_dir(base)}   (新规范 process_flow/0_data_processing；旧事件 sacprocess)")

    for sta in stas:
        net = resolve_network(ev_cfg, sta)
        out_dir = f"{out_root}/{sta}"
        os.makedirs(out_dir, exist_ok=True)
        print(f"\n========== {net}.{sta} ==========")

        # 0. copy 预处理图（默认不出：用户 2026-09-15 定「默认只出三张偏振图」；--all 时恢复拷贝）
        if args.all:
            wim = os.path.join(_stage0_dir(base), "waveform_image", sta)
            for src_n, dst_n in [('3_spectrogram', 'spectrogram'), ('bandpass_2-8Hz', 'bandpass_2-8Hz')]:
                s, d = f"{wim}/{sta}_{src_n}.png", f"{out_dir}/{sta}_{dst_n}.png"
                if os.path.exists(s):
                    shutil.copy(s, d); print(f"  ✓ {dst_n}")
                else:
                    print(f"  ⚠ 预处理图不存在: {s}")

        # 共享模块 (polar_methods.py) 必须先拷贝
        copy_shared_modules(out_dir)

        # 1-2. 模板脚本 (默认 2 个 + --all 时的旧图)
        tpls = TEMPLATES_DEFAULT + (TEMPLATES_LEGACY if args.all else [])
        for tpl in tpls:
            script = render_template(tpl, ev_cfg, sta, net, out_dir)
            if not script:
                continue
            # 只有支持 --method 的生产脚本才带方法参数
            use_method = tpl in ('polarization_analysis', 'pca_theory_polarization')
            print(f"  [{tpl}] ...")
            # 🔴 必须显式传 --out-dir：否则模板用脚本内 OUT（= case 目录），
            #    会把图写到事件目录里（2026-09-15 烟测踩到：AMPY 三张图被写到 case 目录）
            code = _run(script, ['--out-dir', out_dir] + (method_args if use_method else []))
            if code != 0:
                print(f"  ⚠ {tpl} 退出码 {code}")

        # 3. 全时段扫描 (三面板 + CSV), 输出直接落在台站目录 (--flat)
        script = render_template(FULLSCAN, ev_cfg, sta, net, out_dir)
        if script:
            print(f"  [{FULLSCAN}] ...")
            code = _run(script, [sta, '--out-dir', out_dir, '--flat'] + method_args, timeout=None)
            if code != 0:
                print(f"  ⚠ {FULLSCAN} 退出码 {code}")

        # 重命名 PS_polar (p_s_polar_P*_S* → PS_polar)
        for f in glob.glob(f"{out_dir}/*p_s_polar_P*_S*.png"):
            os.rename(f, f"{out_dir}/{sta}_PS_polar.png")
        # 清理中间图
        for pat in ['*p_s_fullwindow*', '*p_s_onset10s*', '*tp_ts_polar*']:
            for f in glob.glob(f"{out_dir}/{pat}"):
                os.remove(f)

        pngs = sorted(os.path.basename(f) for f in glob.glob(f"{out_dir}/*.png"))
        print(f"  ✓ {sta} 完成 ({len(pngs)} 张图): {', '.join(pngs)}")

    print("\n全部完成")


if __name__ == '__main__':
    main()
