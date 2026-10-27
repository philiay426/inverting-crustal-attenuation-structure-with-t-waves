# 波形剖面图（waveform_proflie.py）时间轴与归一化陷阱

## ⚠️ 时间轴 bug（2026-08 南非/马达加斯加事件实测）

**症状**：单台三分量图 T 波很明显，但波形剖面图里 T 波"看不见"（全是平直噪声）。

**根因**：`waveform_proflie.py` 模板旧版时间轴逻辑：
```python
sac_o = getattr(tr.stats.sac, 'o', 0)
if sac_o and 0 < abs(sac_o) < 1e6:
    t = sac_b - sac_o + np.arange(len(data)) * delta
else:
    t = np.arange(len(data)) * delta   # ← 错误！SAC 缺 o 字段时时间轴从数据开始起算
```
当 SAC 头缺 `o` 字段（D 型日数据常见，如南非 2019/2023、马岛 2013），else 分支把 t=0
设在**数据开始时刻**（常为发震前 600s），而 T-P/T-S 理论到时是相对**发震时刻** →
整体偏移 → T 波窗口画的全是噪声。

**修复**（已合入 skill 模板 2026-08-12）：缺失 `o` 时用文件名解析发震时刻：
```python
_origin = _parse_origin_from_filename(bname)   # NET.STA..BHZ.D.2013.203.070100 → 2013-07-22T07:01:00
if _origin is not None:
    t = (tr.stats.starttime - _origin) + np.arange(len(data)) * delta
else:
    t = np.arange(len(data)) * delta
```
模板已含 `_parse_origin_from_filename()` 函数。

**自检**：剖面图 T 波应在 T-P/T-S 虚线附近有纺锤形能量团。若见平直噪声 →
先查时间轴（`sac_o` 是否缺失、`_parse_origin_from_filename` 是否解析成功）。

## 归一化陷阱：强振幅台站压扁全局 P98

马达加斯加 2013-08-17 事件：CRAIL 振幅比其他台大 79 倍 → 全局 P98 归一化被它主导，
其余 23 条线全被压平，看起来像"噪声"。**不是信号弱，是归一化被绑架**。

处置：检查每站 T 波窗口包络峰 vs 全图 P98 的相对量级（>20x 需警惕），
把异常台站加入 `PROFILE_EXCLUDE`（skill 模板已有该集合）。

## 台站振幅诊断脚本模式

对比各台站去响应后全局 std / max：
- 正常台站 std 同量级（如 ~1.3-2.3e-05）
- 异常台站 std 大 10-50 倍（如 CRAIL 5.99e-04）→ 响应文件或增益问题

## 频段选择

- 2-8 Hz：T 波主频，但陆地远距衰减快（马岛 ~1000km 外 12 站仅 4 站可见）
- 0.5-2 Hz：T 波低频成分衰减慢，远处仍可见，但受微震噪声污染更重
- 模板默认 0.5-2，T 波标准分析用 `--band 2-8`
