# SAC 时间轴调试（reference time vs o vs ObsPy starttime）

2026-08-04 南非 2013-08-17 Mw6.2 事件，多次误判的完整复盘。

## 症状

偏振分析窗口"对上"理论到时（T-P=879s），但窗口内波形不像 T 波，反而像 P/S 波偏振。
用户："波形对不上，到时是对的"；"这像将 P/S 波做了个偏振"。

## 数据事实（用户确认，最终定稿）

- 波形数据**从发震前 10 分钟（600s）开始记录**。
- SAC 头：`nzyear..nzsec = 16:22:31`（= 发震时刻！），`b = 0`，`o = 600`。
- 数据实际起点 = 16:12:31（发震前 600s）。
- **`o = 600` 是正确的**：发震在数据起点后 600s。
- ObsPy `tr.stats.starttime` 被 SAC reference time 误导 = 16:22:31（发震时刻），
  因为 reference time 字段被写成了发震时刻而不是数据起点。

## 正确时间轴

```python
o_val = getattr(sac, 'o', None)
if o_val is not None and 0 < abs(o_val) < 1e6:
    t_off = sac.b - o_val          # = 0 - 600 = -600  ← 正确
else:
    t_off = tr.stats.starttime - EV_ORIGIN   # 仅缺 o 时 fallback
t = t_off + np.arange(npts) * delta
```

- `t[0] = -600`：数据从发震前 600s 开始。
- T-P=879s → idx = (879 - (-600))/delta。

## 我犯过的错（不要重蹈）

1. 看到 `starttime == EV_ORIGIN` 就以为 trace 从发震开始 → 错。
   `starttime` 读的是 SAC reference time，该字段可能被写成发震时刻。
2. 据此把 `o=600` 判为"假字段"，改写成 `starttime - EV_ORIGIN = 0` →
   全部图偏移 600s，P 波（实际 188s）在 t=0 坐标系下看起来像 788s 的"T 波"。
3. 中间又反复（b-o → starttime → b-o），每次都是只改一处模板，副本不同步。

## 验证方法

- 读 SAC 头原始字段：`sac.nzyear..nzsec`、`sac.b`、`sac.o`（不是只看 ObsPy starttime）。
- `t[0]` 应为负值（发震前开始记录）。若 SAC 有 o 字段但 t[0]=0 → 坐标系错了。
- 对比全 trace 包络峰值：2013 南非 AMPY 全 trace 包络最大在 P 波（188s），
  T 波在 850-1000s。若"T 波窗口"里看到的是全局最大峰 → 多半错位到了 P 波。

## 模板同步

改时间轴逻辑后必须同步所有副本：
- `single-station-analysis/templates/*.py`（plot_waveform_windows / polarization_analysis /
  pca_polarization / pca_theory_polarization / polarization_attributes）
- `seismic-twave-analysis/templates/waveform_proflie.py` 与 `scripts/step1_station_table_and_map.py`
- `run_single_station.py` 用 `--event --stations` 重新渲染副本
