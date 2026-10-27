# 批量生成单台分析脚本：参数替换遗漏检查（2026-08-05）

## 症状

用 python 字符串替换批量生成 `{STA}_*.py`（从 skill 模板复制）后，
单台分析结果异常：P=0s、T-P=6488s（应为 1278s）、输出跑到 SOE 目录、
`pca() requires exactly 3 traces` 等。

## 根因：批量替换时的 4 类遗漏

| # | 遗漏 | 后果 | 检查 |
|---|------|------|------|
| 1 | **EV 坐标**仍为模板旧事件（`19.0428,-66.8096` 加勒比 2014） | P/S 到时全错（0s / 6488s） | `grep "EVLA"` 确认 = 事件坐标 |
| 2 | **EV_ORIGIN** 仍为模板旧事件（`UTCDateTime(2013,8,17,...)`） | 时间轴偏移 → slice 切空 → pca 报错 | `grep "EV_ORIGIN"` 确认 = 本事件发震 |
| 3 | **OUT 硬编码**（`2_single_station_analysis/SOE`） | 输出写到别的台站目录 | `grep "OUT"` 确认 = `{STA}` 展开 |
| 4 | **`rotate_stream_to_zrt` 后传给 `pca()`** | Non-ZNE 报错（pca 收到 ZRT 不认） | 传原始 ZNE 流 + `back_azimuth`，让 pca 内部旋转 |

## 另一坑：Python `glob` 不支持花括号

```python
glob.glob('./{WTBG,MARK,PUMBA}_*.py')   # ❌ 返回空列表！循环静默不执行
```
要用 `glob.glob('./*_*.py')` + 前缀过滤，或显式列表 `['WTBG','MARK','PUMBA']`。

## 推荐做法

模板参数替换统一用**正则**（而非精确字符串，模板里参数行带注释时精确匹配失败）：
```python
import re
c = re.sub(r'BASE = "[^"]*"', NEW_BASE, c)
c = re.sub(r"STA, NET = '[^']*', '[^']*'", f"STA, NET = '{sta}', 'ZT'", c)
c = re.sub(r"EVLA, EVLO, EVDP = [\d.,\s-]+", "EVLA, EVLO, EVDP = -43.1219, 42.3568, 13.0", c)
```

## 自检

跑完单台分析后检查：
1. print 行 `{STA}: P=xxx S=xxx T-P=xxx` 与 stations.csv 的 tt_p/tt_tp 一致
2. 输出文件名在 `2_single_station_analysis/{STA}/` 下（不是 SOE）
3. 时间窗合理（T-P 相对发震，不是 6000+s）
