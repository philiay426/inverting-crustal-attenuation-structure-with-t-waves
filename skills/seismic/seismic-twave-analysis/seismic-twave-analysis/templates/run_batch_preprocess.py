#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
步 0 批量预处理：对台站清单逐个调用 preprocess_one_station.py，输出六张图 + 一个 npz。

用法（务必用 obspy 环境）：
    /home/sy/miniforge3/envs/obspy/bin/python run_batch.py                 # 全部台站
    STA_LIST=<A>,<B>,<C> python run_batch_preprocess.py                   # 只跑指定台站
    OUT=/tmp/xx python run_batch.py                                        # 换输出根目录

产物落盘（照案例 sacprocess/ 布局）：
    <OUT>/waveform_image/<STA>/<STA>_1_raw.png
    <OUT>/waveform_image/<STA>/<STA>_2_remove_response.png
    <OUT>/waveform_image/<STA>/<STA>_3_spectrogram.png
    <OUT>/waveform_image/<STA>/<STA>_4_envelope.png
    <OUT>/waveform_image/<STA>/<STA>_bandpass_2-8Hz.png
    <OUT>/waveform_image/<STA>/<STA>_packet.png
    <OUT>/envelope/<STA>_envelope_2-8Hz.npz
日志：<OUT>/batch_log.txt（每台的 exit code、耗时、SNR 行、失败时的 stderr 尾部）

为什么用 subprocess 而不是循环调用 main()：
  单台脚本的参数在 import 时就固化了（STA/OUT 都是模块级常量），
  反复 import 同一模块不会重新执行顶层代码；起子进程最干净，且单台崩溃不影响整批。
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable                                   # 当前解释器 = obspy 环境

# 案例路径（与 preprocess_one_station.py 的 SRC 保持一致）
SRC = os.environ.get("SRC", "").strip()                       # 事件根（必填）
if not SRC:
    raise SystemExit("请用环境变量 SRC 指定事件根目录")
STACSV = os.path.join(SRC, "process_flow", "1_station_distribution", "stations.csv")
OUT = os.environ.get("OUT") or HERE                        # 本步目录本身
FIG_OUT = os.path.join(OUT, "waveform_image")
ENV_OUT = os.path.join(OUT, "envelope")


def station_list():
    """台站清单：默认取 stations.csv（跳过 # 注释行），可用 STA_LIST 覆盖"""
    env_list = os.environ.get("STA_LIST", "").replace(" ", "")
    if env_list:
        return [s for s in env_list.split(",") if s]
    stas = []
    with open(STACSV) as f:
        for line in f:
            if line.startswith("#") or line.startswith("sta,"):
                continue
            sta = line.split(",")[0].strip()
            if sta:
                stas.append(sta)
    return stas


def main():
    stas = station_list()
    n = len(stas)
    os.makedirs(FIG_OUT, exist_ok=True)
    os.makedirs(ENV_OUT, exist_ok=True)
    log = open(os.path.join(OUT, "batch_log.txt"), "w")
    log.write(f"# 批量预处理 {n} 台 | 解释器 {PY}\n")
    log.write(f"# 图 -> {FIG_OUT}\n# 数据 -> {ENV_OUT}\n\n")
    print(f"共 {n} 台，输出到 {FIG_OUT} 和 {ENV_OUT}")
    t_all = time.time()
    ok, bad = [], []
    for i, sta in enumerate(stas, 1):
        t0 = time.time()
        env = dict(os.environ, STA=sta, OUT=OUT, FIG_OUT=FIG_OUT, ENV_OUT=ENV_OUT)
        r = subprocess.run([PY, os.path.join(HERE, "preprocess_one_station.py")],
                           cwd=HERE, env=env, capture_output=True, text=True)
        dt = time.time() - t0
        snr_line = ""
        for ln in r.stdout.splitlines():
            if ln.strip().startswith("SNR("):
                snr_line = ln.strip()
        if r.returncode == 0:
            ok.append(sta)
            print(f"[{i}/{n}] {sta} ok  {dt:5.1f}s  {snr_line}")
        else:
            bad.append(sta)
            print(f"[{i}/{n}] {sta} FAIL {dt:5.1f}s")
        log.write(f"----- {sta}  exit={r.returncode}  {dt:.1f}s\n")
        log.write(r.stdout)
        if r.returncode != 0:
            log.write("STDERR (tail):\n" + r.stderr[-3000:] + "\n")
        log.flush()

    log.write(f"\n# 成功 {len(ok)}: {','.join(ok)}\n# 失败 {len(bad)}: {','.join(bad)}\n")
    log.close()
    print(f"\n完成 {len(ok)}/{n} 台，总耗时 {(time.time() - t_all) / 60:.1f} min")
    if bad:
        print("失败台站：", ", ".join(bad), "（详见 batch_log.txt）")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
