#!/usr/bin/env python3
"""
单台分析公共模块 (v4, 2026-08-04)
- 事件参数: events.yaml
- 发震时刻: 从 SAC 文件名解析 (year.julday.HHMMSS) —— 不依赖 SAC o 字段
- 时间轴：`o` 有效时 t=b-o；缺 o 时 t=starttime-origin+i*delta，t=0=发震时刻

⚠️ 时间轴铁律：数据常从发震前 600s 开始。SAC `o` 有效（如 o=600,b=0）时，
`t_off=b-o=-600`；不能仅看 ObsPy starttime（reference time 可能被写成发震时刻）。
"""
import os, re, glob, csv, sys
import numpy as np
from obspy import read, UTCDateTime

try:
    import yaml
except ModuleNotFoundError:      # obspy 环境常缺 PyYAML -> 用内置兜底解析(见 _parse_events_minimal)
    yaml = None

CODE_ROOT = "/home/sy/codex"
CONFIG_PATH = os.path.join(CODE_ROOT, "config", "events.yaml")


def _parse_events_minimal(path):
    """PyYAML 缺失时的兜底解析: 只支持 events.yaml 的扁平两级结构(标量/数字列表)。
    仅在没有 PyYAML 的解释器里生效, 有 PyYAML 时走 yaml.safe_load。"""
    import ast as _ast
    events, cur = {}, None
    for raw in open(path):
        line = raw.split('#')[0].rstrip()
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip())
        m = re.match(r'([^:\s]+):\s*(.*)$', line.strip())
        if not m:
            continue
        key, val = m.group(1), m.group(2).strip()
        if indent == 0:
            cur = key
            events[cur] = {}
        elif cur is not None:
            if val == '':
                events[cur][key] = None
            else:
                try:
                    events[cur][key] = _ast.literal_eval(val)
                except (ValueError, SyntaxError):
                    events[cur][key] = val
    return events


def _load_all_events(config_path):
    if yaml is not None:
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    return _parse_events_minimal(config_path)


def load_event(event_name, config_path=None):
    """从 events.yaml 读事件参数"""
    if config_path is None:
        config_path = CONFIG_PATH
    all_events = _load_all_events(config_path)
    if event_name not in all_events:
        raise KeyError(f"事件 '{event_name}' 不在配置中: {', '.join(all_events)}")
    return all_events[event_name]


def parse_origin_from_filename(bname):
    """
    从 SAC 文件名解析发震时刻: ...2013.229.162231_resp.sac → UTCDateTime
    返回 UTCDateTime；解析失败返回 None
    """
    m = re.search(r'\.(\d{4})\.(\d{3})\.(\d{6})', bname.upper())
    if not m:
        return None
    y = int(m.group(1)); jd = int(m.group(2)); hm = m.group(3)
    try:
        return UTCDateTime(year=y, julday=jd,
                           hour=int(hm[0:2]), minute=int(hm[2:4]), second=int(hm[4:6]))
    except Exception:
        return None


def load_station_csv(csv_path, sta=None):
    """读 stations.csv，返回 dict[sta] = row 或单站 dict"""
    rows = {}
    with open(csv_path) as f:
        for r in csv.DictReader((l for l in f if not l.startswith('#'))):
            rows[r['sta']] = r
    if sta is not None:
        return rows.get(sta)
    return rows


def find_channel(resp_dir, net, sta, ch):
    """在 remove_response 目录找单分量文件"""
    for pat in [f"{resp_dir}/{net}.{sta}.*.{ch}*_resp.sac",
                f"{resp_dir}/{net}.{sta}..{ch}*_resp.sac",
                f"{resp_dir}/{net}.{sta}.{ch}*_resp.sac"]:
        m = glob.glob(pat)
        if m:
            return m[0]
    return None


def load_3comp(resp_dir, net, sta, origin=None):
    """
    加载三分量 ZNE，统一时间轴
    返回 dict: {trs, delta, sr, t, origin, t_off}
    t 轴: t=0 = 发震时刻 (文件名解析, 不依赖 SAC o)
    """
    trs = {}
    for ch in ['BHZ', 'BHN', 'BHE']:
        fn = find_channel(resp_dir, net, sta, ch)
        if fn:
            trs[ch] = read(fn)[0]
    if len(trs) < 3:
        return None
    t0 = max(tr.stats.starttime for tr in trs.values())
    t1 = min(tr.stats.endtime for tr in trs.values())
    for ch in trs:
        trs[ch].trim(t0, t1)
    delta = trs['BHZ'].stats.delta
    sr = 1 / delta

    # 发震时刻: 文件名解析优先, fallback trace 起点
    if origin is None:
        bhz_fn = find_channel(resp_dir, net, sta, 'BHZ')
        origin = parse_origin_from_filename(os.path.basename(bhz_fn)) if bhz_fn else None
        if origin is None:
            origin = trs['BHZ'].stats.starttime

    # 时间轴: SAC o 字段可靠时用 b-o (数据常从发震前开始记录, 2013南非 o=600)
    # ObsPy starttime 可能被 SAC reference time 误导 (reference 写成了发震时刻)
    sac = trs['BHZ'].stats.sac
    o_val = getattr(sac, 'o', None)
    if o_val is not None and 0 < abs(o_val) < 1e6:
        t_off = sac.b - o_val
    else:
        t_off = trs['BHZ'].stats.starttime - origin
    t = np.arange(trs['BHZ'].stats.npts) * delta + t_off

    return {'trs': trs, 'delta': delta, 'sr': sr, 't': t, 'origin': origin, 't_off': t_off}


def _first_file(resp_dir, net, sta):
    f = find_channel(resp_dir, net, sta, 'BHZ')
    return f or ''
