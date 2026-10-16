#!/usr/bin/env python3
"""步骤1: 数据预处理 — 去仪器响应+2-8Hz滤波+频谱图(仅垂直分量)"""
import sys, os, glob, argparse, yaml

CODE_ROOT = "/home/philiay/shared/codex"
for p in [CODE_ROOT, CODE_ROOT + '/t_phase_analysis']:
    if p not in sys.path: sys.path.insert(0, p)

from lib.processing import process_sac_file
from obspy.io.sac import SACTrace

CONFIG = yaml.safe_load(open(f"{CODE_ROOT}/config/events.yaml"))

def main():
    p = argparse.ArgumentParser()
    p.add_argument('event', help='事件key')
    p.add_argument('--skip-existing', action='store_true')
    args = p.parse_args()

    cfg = CONFIG[args.event]
    case_dir = cfg.get('case_dir', cfg['label'])
    raw_dir = cfg['data_dir']
    sacpz_dir = os.path.join(raw_dir, 'SACPZ')
    out_dir = f"/home/philiay/shared/case/{case_dir}/sacprocess"
    os.makedirs(out_dir, exist_ok=True)

    sac_files = sorted(glob.glob(f"{raw_dir}/*.SAC"))
    print(f"事件: {cfg['label']} | {len(sac_files)} SAC文件")
    print(f"  去响应 → {out_dir}/\n")

    ok = 0
    for fn in sac_files:
        base = os.path.basename(fn).replace('.SAC','').replace('.sac','')
        if args.skip_existing and os.path.exists(f"{out_dir}/{base}_proc.sac"):
            continue
        try:
            process_sac_file(fn, sacpz_dir, out_dir, only_bhz=True)
            ok += 1
        except Exception as e:
            print(f"  ❌ {os.path.basename(fn)}: {e}")

    print(f"\n完成: {ok}/{len(sac_files)} 成功")

if __name__ == '__main__':
    main()
