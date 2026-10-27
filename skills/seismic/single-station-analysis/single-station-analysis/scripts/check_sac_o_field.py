#!/usr/bin/env python3
"""
检查事件目录 SAC 头是否含 o 字段 — 预测 t_off 时间轴偏移风险

用法: python3 check_sac_o_field.py <SAC目录> [--limit N]
例:   python3 check_sac_o_field.py /mnt/f/sharedate/case/Southern_Africa/2019-01-22-mww67-prince-edward-islands-station_low/raw

背景 (2026-08 南非 Mw6.7 实测):
部分 D 型日数据 SAC 缺 o 字段(发震时刻偏移)。若脚本用
`t_off = sac.b - getattr(sac, 'o', 0)`, 缺 o 时 t_off=0, 时间轴从
记录开始时刻起算而非发震时刻, 整体偏移(南非 2019-01-22 偏移 1903s,
T 波窗口画成全噪声, 振幅差 100 倍)。这是静默失败 — 不报错, 图能画出,
但内容全错。

实测结果:
- 2019-01-22-mww67-prince-edward-islands-station_low: 缺 o ❌
- 2023-05-21-mww68-prince-edward:                           缺 o ❌
- 2013-08-17 / 2018-04-19 / 2020-08-06 / 2020-09-26:       有 o ✅
规律: 与 D/M 类型无必然关系, 取决于下载数据中心/批次, 无法预先判断。

输出: 每个台站 BHZ 的 o 字段值。全部缺失 → 该事件的所有单台分析脚本
必须用 `t_off = starttime - EV_ORIGIN`(发震时刻) 方案, 禁止 getattr 兜底。
"""
import sys, glob, os
from obspy import read

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    d = sys.argv[1]
    limit = None
    if '--limit' in sys.argv:
        limit = int(sys.argv[sys.argv.index('--limit') + 1])

    files = sorted(glob.glob(f"{d}/*.SAC") + glob.glob(f"{d}/*.sac"))
    if not files:
        print(f"目录无 SAC 文件: {d}")
        sys.exit(2)

    # 按台站+分量去重, 只检查 BHZ (每个台站一个代表即可)
    seen = set()
    missing = 0
    checked = 0
    for f in files:
        base = os.path.basename(f)
        # 提取台站名 (如 AF.SOE..BHZ.D... → SOE)
        parts = base.split('.')
        sta = parts[1] if len(parts) > 1 else base
        if not any(x in base.upper() for x in ['BHZ', 'HHZ', 'EHZ']):
            continue
        key = sta
        if key in seen:
            continue
        seen.add(key)
        try:
            tr = read(f)[0]
            sac = tr.stats.sac
            o = getattr(sac, 'o', None)
            o_str = f"{o:.1f}" if o is not None else "MISSING ❌"
            if o is None:
                missing += 1
            print(f"{sta:8s} o={o_str}  starttime={tr.stats.starttime}")
            checked += 1
        except Exception as e:
            print(f"{sta:8s} 读取失败: {e}")
        if limit and checked >= limit:
            break

    print(f"\n检查 {checked} 个台站, {missing} 个缺 o 字段")
    if missing == checked and checked > 0:
        print(">>> 该事件全部 SAC 缺 o 字段 → 单台分析必须用 t_off = starttime - EV_ORIGIN")
        sys.exit(1)
    elif missing > 0:
        print(">>> 部分台站缺 o → 脚本需逐台判断 (有 o 用 sac.b-o, 缺 o 用 starttime-EV_ORIGIN)")
        sys.exit(1)
    else:
        print(">>> 所有台站都有 o 字段, t_off = sac.b - sac.o 安全")
        sys.exit(0)

if __name__ == '__main__':
    main()
