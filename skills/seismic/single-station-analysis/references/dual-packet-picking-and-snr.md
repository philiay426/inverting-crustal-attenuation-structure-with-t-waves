# 双波包独立拾取与 SNR 质量控制（2026-08-18）

> ⚠️ 本文的「单波包 SNR」是**偏振窗口内部的局部质量指标**（包络功率比，判断波形包可用性），
> 与台站表/剖面图选台用的 **统一噪声基准 SNR**（20log10 信号窗 RMS / 阵列噪声中位数，
> 2–8 Hz 整段带通后切窗，见 seismic-twave-analysis `references/snr-unified-noise-baseline.md`）
> **不是同一个量**，不要互相套用或拿本文的数值当台站筛选依据。

## 问题

T-P 与 T-S 是两个需要分别分析的小波包。若在整个 T 波时间段内只做一次包络 `argmax`，算法会偏向最大波包，导致另一个波包漏拾、重复拾取，或把强波包误标成另一个波包。

## 推荐拾取流程

对 T-P、T-S 各自建立理论到时中心的搜索窗，例如：

```python
from scipy.signal import find_peaks

def candidate_peaks(t, env, t0, t1, sr):
    m = (t >= t0) & (t <= t1)
    ts, es = t[m], env[m]
    if len(es) == 0:
        return []
    peaks, props = find_peaks(
        es,
        distance=max(1, int(5 * sr)),
        prominence=0.15 * np.max(es),
    )
    return [
        {
            "time": float(ts[i]),
            "height": float(es[i]),
            "prominence": float(props["prominences"][j]),
            "dt_theory": float(ts[i] - (t0 + t1) / 2),
        }
        for j, i in enumerate(peaks)
    ]
```

候选峰按照峰高、prominence、与理论到时的一致性及局部 SNR 综合评价；不要仅按全局最大峰选择。T-P、T-S 禁止共用一个峰，偏振窗口不能重叠。没有合格候选峰时，波包标记为不可用。

## 相对阈值的正确定位

30% 或 50% 相对峰值阈值只用于定义已确认波包的边界（峰值两侧下降到相对幅度阈值处），不是严格的初至判据，也不能单独证明波包存在。建议对 20%/30%/50% 做敏感性分析，比较峰值时间、窗口长度、偏振角和线性度。

## 单波包 SNR

对 T-P 和 T-S 分开计算频带能量比，而不是只用整个 T 波段总 SNR：

```python
snr_db = 10 * np.log10(power(packet_window) / power(noise_window))
```

噪声窗应位于 P 波前、避开前驱 T 波；信号窗应为该波包独立窗口或固定的局部波包窗，2–8 Hz 频带与当前流程一致。

## 建议分级

| 单个波包 SNR | 用途 |
|---:|---|
| <3 dB | 不可用 |
| 3–6 dB | 不作可靠偏振结论 |
| 6–10 dB | 探索性结果，需额外稳定性检查 |
| ≥10 dB | 正式双波包偏振分析最低门槛 |
| ≥15 dB | 高可信度核心台站 |

正式结论使用 `min(SNR_TP, SNR_TS) >= 10 dB`，不能用总 SNR 替代。

## 偏振稳定性附加门槛

- 每个波包中位线性度 `>= 0.6`；
- 偏振角圆形标准差建议 `<= 20°`，或至少 60% 子窗落在中位角 ±20°；
- T-P、T-S 窗口不重叠；无法分离时标记为不可用，不强行后移窗口。

## 文献边界

Flinn (1965) 和 Jurkevics (1988) 支持协方差/特征值偏振和滑动窗口，但不规定 30% 或 50% 固定包络比例。包络初至研究（如 Al-Mashhor et al., 2019, DOI: 10.1109/ACCESS.2019.2939320；Bargees & Al-Shuhail, 2022, DOI: 10.3390/s22197580）更支持包络结合局部峰值、prominence 和噪声统计的自适应判据，而非宣称固定百分比为普适标准。
