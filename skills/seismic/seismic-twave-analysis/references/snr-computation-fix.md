# SNR 计算修复 (2026-08-03)

> ⚠️ **2026-10-14 起已被取代**：现行 SNR 口径 =「统一噪声基准 SNR」——2–8 Hz **整段带通**后取窗内 RMS，
> 除以**阵列噪声中位数**；完整参数与两个坑见 `references/snr-unified-noise-baseline.md`。
> 本文（Welch 频段积分 + P 前动态噪声窗）仅作历史存档，勿再据此实现。

## 问题

1. **NaN SNR**：近台（如 FOMA, land=54km, P=166s）噪声窗不足导致 PSD 能量为零 → log(0) = NaN
2. **float32 溢出**：大振幅数据 (10^24) 平方后溢出 → PSD inf/NaN
3. **噪声窗太小**：固定 300s 窗口对近台不合理

## 修复

### 1. 动态噪声窗
```python
available_pre = noise_center  # P 波前的全部可用时间
actual_noise = max(20, min(noise_width, available_pre * 0.8))
```
- 近台（P=166s）→ 噪声窗 ~133s ✓
- 远台（P=400s+）→ 噪声窗 300s（上限）✓

### 2. 噪声中心用 `tt_p`
改用 `noise_center=tt_p`（而非 `tt_p - 100`），充分利用全部可用噪声。

### 3. float64 强制转换
```python
data[s_start:s_end].astype(np.float64)  # 修复溢出
```

### 4. 回退值
```python
if s_end <= s_start or n_end <= n_start:
    return 0.0  # 非 NaN 或 -999
```

## 验证

FOMA: 修复前 NaN → 修复后 30.6 dB
