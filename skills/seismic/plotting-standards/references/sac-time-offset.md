# SAC `o` header: 有效则用 b-o，缺失才用 event_origin - starttime

⚠️ **规则（2026-08-17 修正，勿用旧版"一律不信任 sac.o"）**：

| 情况 | 做法 |
|------|------|
| `o` 有效（0<|o|<1e6） | **`t_off = sac.b - sac.o`**（如 2013-08-17 南非 `o=600` → `t_off=-600`，数据从发震前 600s 开始，实测正确） |
| `o` 缺失 | `t_off = tr.stats.starttime - EV_ORIGIN`（EV_ORIGIN 从文件名解析） |
| `o=0` 或超界 | 视为无效，走缺失分支 |

```python
def get_t0_offset(tr, event_origin_utc):
    """发震偏移（秒，负=数据从发震前开始）"""
    sac = tr.stats.sac
    o = getattr(sac, 'o', None)
    if o is not None and 0 < abs(o) < 1e6:
        return sac.b - o
    return event_origin_utc - tr.stats.starttime
```

Then slice data correctly:

```python
t0_offset = get_t0_offset(tr, UTCDateTime(2019, 1, 22, 19, 1, 43))
i1 = int((t_min + t0_offset) * sr)
i2 = int((t_max + t0_offset) * sr)
t = np.arange(i2 - i1) / sr + (i1 / sr - t0_offset)  # origin-relative time
```

## Common failure pattern

```
AttributeError: 'o'  # SAC header has no 'o' field at all
ValueError: could not convert string to float  # CSV header row not skipped
```

## How to diagnose

```python
print(f"o={getattr(sac, 'o', 'N/A')}  b={sac.b}")
print(f"starttime={tr.stats.starttime}")
print(f"origin - start = {origin - tr.stats.starttime:.0f}s")
```

⚠️ 不要仅凭 `starttime == EV_ORIGIN` 判断 trace 从发震开始：ObsPy starttime 读 SAC reference time，
该字段可能被写成发震时刻而数据实际从发震前开始（2013-08-17 南非 `o=600` 案例）。

自检：窗口内应有纺锤形能量包络；若只见平坦噪声+包络平直线 → 先查 t_off（曾致南非 2019 偏移 1903s）。
