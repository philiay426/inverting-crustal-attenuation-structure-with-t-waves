# 波形图纵轴稳健定标（防体波尖峰顶爆）

日期：2026-09-02（马达加斯加 2013-08-17 Mw6.2，AMPY 案例）

## 症状（用户反馈）
单台 `{STA}_bandpass_2-8Hz.png`（`process_image.py` 预处理图，会复制进
`2_single_station_analysis/{STA}/`）"纵轴上下 10 了，太大了，T 波非常不明显"。

## 根因
`plot_bandpass` 旧定轴规则取**窗口内最大值**：
```python
all_wave = np.concatenate([np.abs(filt[l][i1:i2]) for l, _ in comps])  # 全窗 0..T-S+500
ep98 = np.percentile(all_wave, 98) or 1.0
y_max = 1.05 * max(max(np.abs(filt[l][i1:i2]/ep98)), max(env[l][i1:i2]/ep98), 1.0)
```
AMPY 数值：2–8Hz 带通后 **P 波（t≈187.5s，理论 P=185s）振幅达 11.7×ep98**
（Z 分量；R/T 亦 4.5/9.7×）→ y_max=**±12.26**；T 波窗口（T-P−100~T-S+200）
包络峰值只有 0.57–0.97×ep98 → 占轴高 ~5%，肉眼不可见。
尖峰是**真实 P 波**（非毛刺），只是短时高振幅；ep98 被 T 波/背景主体占住，故 P 远超 P98。

## 修复（process_image.py::plot_bandpass，2026-09-02）
改用稳健分位数定轴，排除 P/S 短时尖峰：
```python
all_norm = np.concatenate(
    [np.abs(filt[l][i1:i2]) / ep98 for l, _ in comps] +
    [env[l][i1:i2] / ep98 for l, _ in comps])
y_max = 1.05 * max(float(np.percentile(all_norm, 99)), 1.0)
```
AMPY 结果：y_max=±2.04（旧 ±12.26）；T 波包络占轴高 ~48%。
P/S 顶部被截断是预期行为（到时参考线仍在），这正是 T 波视图想要的效果。

## 验证方式（不依赖视觉模型）
复算 y_max 与 T 波窗口包络占轴高百分比，确认目标量可见度：
```python
tw = (t >= tt_tp - 100) & (t <= tt_ts + 200)
tw_env_max = max(env[l][i1:i2][tw].max() / ep98 for l in comps)
print(y_max, tw_env_max / y_max * 100)   # AMPY: 2.04, ~48%
```

## 遗留
- `plot_raw`（1_raw）与 `plot_decon`（2_remove_response）仍是同一 max 定轴逻辑，
  同样会被 P/S 顶爆（±10 以上）。这两张是数据质量检查图，用户未要求改；
  若用户要求 T 波在这些图也可读，套用同一 p99 规则即可。
- 重生成后若文件在 2_single_station_analysis 下有副本（md5 相同），
  需 `cp` 同步（process_image.py 只写 `sacprocess/waveform_image/{STA}/`）。
- 重跑 `process_image.py` 会全量重做 25 站去响应（约数分钟，有逐站进度打印），
  用后台 + notify_on_complete 即可。
