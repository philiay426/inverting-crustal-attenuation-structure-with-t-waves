# SNR 计算方法（2026-07-24 定稿）

## 频谱 dB 能量比 (2-8Hz 频段)

```python
from scipy.signal import welch

# 信号窗: T-P 到 T-S+80s（T波能量区）
m_sig = (t >= tt_TP) & (t <= tt_TS + 80)

# 噪声窗: P波前 P-300 到 P-60（用 P 波前而非 T-P 前）
m_noise = (t >= tt_P - 300) & (t <= tt_P - 60)

# Welch PSD → 2-8 Hz 积分
f_sig, psd_sig = welch(data[m_sig], fs=sr, nperseg=256)
f_noise, psd_noise = welch(data[m_noise], fs=sr, nperseg=256)

energy_sig = np.sum(psd_sig[(f_sig>=2)&(f_sig<=8)])
energy_noise = np.sum(psd_noise[(f_noise>=2)&(f_noise<=8)])

SNR = 10 * log10(energy_sig / energy_noise)  # dB
```

## 包络计算（tpef 方法）

```python
from obspy.signal.filter import envelope
from obspy.signal.util import smooth

env = envelope(data)                       # sqrt(data² + hilbert(data)²)
env_smooth = smooth(env, max(1, int(3*sr)))  # 3秒滑动窗口
```

**对比旧方法：** `np.abs(hilbert(data))` + `uniform_filter1d(size=2s)` → 废弃

## 关键说明

- 噪声窗用 P 波前而非 T-P 前：P 波前更干净，避免 T-P 前可能含有的前驱 T 波信号
- 对于高 SNR 台站（>10 dB），两种噪声窗口差异很小（1-2 dB）
- 对于中等 SNR 台站，差异可能较大（最高 7 dB），排序基本不变
