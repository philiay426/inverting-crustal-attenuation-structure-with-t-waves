# 包络阈值选择与文献边界

## 结论

固定的 30% 或 50% 峰值比例都不是普适的“初至拾取标准”。在 T 波单台偏振分析中应区分三种量：

1. **包络峰值**：搜索窗内 `argmax(envelope)`，用于定位主要波列中心；
2. **偏振分析窗口**：可用 `0.30 * envelope_peak` 作为经验边界，以保留起始—峰值—尾部的完整波列；
3. **观测初至**：应使用噪声统计阈值，例如 `mu_noise + 3*sigma_noise`，并要求连续超过阈值 0.5–1.0 s 或至少 2–3 个主周期。

因此当前推荐：30% 作为偏振窗口的主方案，同时做 20%/30%/50% 敏感性分析；不能把 30% 写成有统一文献依据的固定标准。

## 为什么不直接把 50% 当标准

50% 窗口更集中于主能量，较少包含前导噪声，但可能截断 T 波起始部分；30% 窗口更适合 T 波持续波列的 PCA/偏振分析，但低 SNR 时更容易吸收噪声。应根据窗口中心、窗口长度、偏振角和线性度的稳定性决定是否可靠。

## 推荐的验证指标

对高 SNR 台站分别运行 20%、30%、50%：

- 窗口中心差异：目标小于约 2–3 s；
- R-T 中值角差异：目标小于约 10°；
- 线性度/平面度趋势不发生定性改变；
- T-P 与 T-S 窗口不重叠，且不明显覆盖 P/S 或其他波列。

若不满足，应保留“拾取不稳定”标记，不强行选一个比例。

## 文献边界

- Al-Mashhor, A. A., Al-Shuhail, A. A., Hanafy, S. M., & Mousa, W. A. (2019). *First Arrival Picking of Seismic Data Based on Trace Envelope*. IEEE Access, 7, 128806–128815. DOI: 10.1109/ACCESS.2019.2939320. 支持使用地震道包络进行初至识别，但不规定通用 30%/50% 峰值比例。
- Bargees, A., & Al-Shuhail, A. A. (2022). *First Arrival Picking of Zero-Phase Seismic Data by Hilbert Envelope Empirical Half Window (HEEH) Method*. Sensors, 22, 7580. DOI: 10.3390/s22197580. 使用 Hilbert 包络和噪声异常检测/三倍标准差思想，说明严格初至更适合统计噪声阈值，而不是固定峰值比例。
- Flinn (1965) 与 Jurkevics (1988) 可支持协方差偏振和子窗覆盖若干主频周期，但不能直接作为 30% 或 50% 包络阈值的依据。

## 实现边界

如果代码写成：

```python
i0 = np.argmax(envelope_segment)
emax = envelope_segment[i0]
above = local_envelope > 0.30 * emax
```

应在方法中称为“包络峰值定位 + 相对峰值阈值窗口提取”，不能称为“30% 阈值首次到时拾取”。
