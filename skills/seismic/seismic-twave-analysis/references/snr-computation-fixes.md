# SNR 计算注意事项

> ⚠️ **2026-10-14 起已被取代**：现行 SNR 口径 =「统一噪声基准 SNR」——2–8 Hz **整段带通**后取窗内 RMS，
> 除以**阵列噪声中位数**；完整参数与两个坑见 `references/snr-unified-noise-baseline.md`。
> 本文（Welch 频段积分 + P 前动态噪声窗 + 旧脚本名 `step1_station_table_and_map.py`）仅作历史存档。

## 动态噪声窗

固定 300s 噪声窗对近台（P < 300s 到达）会截断。
`step1_station_table_and_map.py` 已改为动态窗：

```python
available_pre = noise_center  # P 波到达时间
actual_noise = max(20, min(300, available_pre * 0.8))
```

即取 P 前可用时长的 80%，最少 20s，最多 300s。
`noise_center` 直接用 `tt_p`（不再减 100s）。

## 数值溢出

`scipy.signal.welch` 默认输出 float32，去响应后的速度数据值可能非常大（>1e24），
导致 PSD 溢出 → SNR=0。修复：显式 `astype(np.float64)`。

```python
f_sig, psd_sig = welch(data[s_start:s_end].astype(np.float64), sr, ...)
f_noise, psd_noise = welch(data[n_start:n_end].astype(np.float64), sr, ...)
```

## 异常处理

```python
if e_noise > 0 and e_sig > 0:
    ratio = e_sig / e_noise
    if ratio > 0 and np.isfinite(ratio):
        return 10 * np.log10(ratio)
return 0.0  # 回退为 0（非 NaN），不阻塞表格生成
```
