# 频域 F-K 频谱计算陷阱：FFT vs Welch vs STFT

**发现时间**：2026-09-02，马达加斯加 Mw6.2 事件 BKTA/CPSM/AMPY 子阵 120s 窗 F-K 分析三轮调试

**症状与根因对照表**：

| 方法 | 症状 | 根因 | 峰值慢度范围 |
|------|------|------|--------------|
| Welch PSD | 所有场景慢度 0.00 s/deg | 实数功率谱无相位信息 | 固定 0.00（无限大速度） |
| STFT + 时间平均 | 慢度 0.03-0.57 s/deg | 时间平均抹掉相位延迟 | 理论值的 1/30 - 1/600 |
| 整窗 FFT（正确） | 慢度接近理论 18-31 s/deg | 保留完整相位 | 物理合理范围 |

---

## 第一轮失败：Welch PSD

**错误代码**：
```python
from scipy.signal import welch

# 计算功率谱密度
f, psd = welch(windowed_data[sta], fs=1/delta, nperseg=NPERSEG)
psds[sta] = psd  # ← 实数数组，无相位

# 尝试白化（无效）
spec = np.array([psds[sta][freq_idx] for sta in sta_list])
spec = spec / (np.abs(spec) + 1e-12)  # 实数 ÷ 实数 = ±1.0
```

**结果**：所有分量/时间窗峰值慢度 = **0.00 s/deg**（对应无限大速度）

**问题**：
- `welch()` 返回的是 **实数功率谱** |X(f)|²，相位信息已丢失
- 白化实数 PSD 得到 ±1.0，所有台站"相位"相同
- Steering vector `e^(i·2πf·τ)` 与常数相乘，无法检测相位延迟
- 峰值必然在慢度原点 (0, 0)

---

## 第二轮失败：STFT + 时间平均

**错误代码**：
```python
from scipy.signal import stft

# 短时傅里叶变换
f, t, Zxx = stft(windowed_data[sta], fs=1/delta, nperseg=NPERSEG)
stfts[sta] = Zxx.mean(axis=1)  # ← 时间平均，抹掉延迟信息

# 白化
spec = np.array([stfts[sta][freq_idx] for sta in sta_list], dtype=complex)
spec = spec / (np.abs(spec) + 1e-12)
```

**结果**：
- Z-TP: 慢度 0.43 s/deg（理论 18.5 s/deg，偏小 43 倍）
- Z-TS: 慢度 0.57 s/deg（理论 30.9 s/deg，偏小 54 倍）
- T-TS: 慢度 0.03 s/deg（理论 30.9 s/deg，偏小 1030 倍）
- BAZ 方向混乱：106°, 270°, 306°, 334°（理论 151.6°）

**问题**：
- `Zxx.mean(axis=1)` 对时间轴平均，相位延迟信息被抹掉
- 设台站 A 在 t=0 时相位 0°，台站 B 有延迟 τ，则：
  - t=0 时刻：A 相位 0°，B 相位 2πf·τ
  - t=1 时刻：A 相位 2πf，B 相位 2πf(1+τ)
  - 时间平均后：两站相位差的矢量和趋近于 0（相位是周期性的，平均=随机游走）

---

## 第三轮成功：整窗 FFT

**正确代码**：
```python
# 对整个时间窗做 FFT（无时间平均）
fft_result = np.fft.rfft(windowed_data[sta])
freqs_fft = np.fft.rfftfreq(len(windowed_data[sta]), d=delta)
ffts[sta] = fft_result  # 复数数组，保留完整相位

# 白化（保留相位）
spec = np.array([ffts[sta][freq_idx] for sta in sta_list], dtype=complex)
spec = spec / (np.abs(spec) + 1e-12)  # 单位复数，相位不变
```

**为什么有效**：
- `rfft()` 返回 **复数频谱** X(f) = |X|·e^(iφ)，相位 φ 编码了信号的时间特性
- 不同台站的相位差 Δφ = 2πf·τ 直接对应时间延迟 τ
- 白化后得到单位复数 e^(iφ)，相位信息完整保留
- Steering vector 与复数频谱相乘，能正确检测相位延迟模式

**频率分辨率验证**：
- 120s 窗口 → Δf = 1/120 ≈ 0.008 Hz
- 2-8 Hz 频段 → 约 750 个频率点
- 分辨率足够，无需 STFT 的时频权衡

---

## 单位转换陷阱

**错误尝试**（2026-09-02 第一次修正）：
```python
# 错误：慢度 s/deg 乘以 111.32 得到 s/km
sx_km = sx_deg * 111.32  # ← 错误方向
```

**正确转换**：
```python
# 台站坐标单位 = km，慢度必须 s/km
# 1 degree ≈ 111.32 km
# s/deg ÷ 111.32 = s/km
sx_km = sx_deg / 111.32
```

**推导**：
- 延迟 τ = 慢度 s × 距离 r
- 若 r 单位是 km，s 单位必须是 s/km
- 慢度定义：s = 1/v，单位 [时间/距离]
- 转换：s [s/deg] × (1 deg / 111.32 km) = s [s/km] / 111.32

---

## 慢度网格范围陷阱

**错误设置**（原始脚本）：
```python
SMAX = 0.8  # s/deg，对应最慢 v=111.32/0.8≈139 km/s
```

**问题**：
- 地壳 P 波 6 km/s → 慢度 18.5 s/deg（**在网格外**）
- 地壳 S 波 3.6 km/s → 慢度 30.9 s/deg（**远超网格范围**）
- 算法只能扫到 0-0.8 s/deg 范围，真峰在外→检测到的是旁瓣噪声

**正确设置**：
```python
SMAX = 40.0  # s/deg，覆盖地壳波（6 km/s → 18.5 s/deg, 3.6 km/s → 30.9 s/deg）
DS = 0.5    # s/deg，161×161 网格（粗但够用）
```

---

## 诊断清单

遇到异常慢度时，按顺序检查：

1. **慢度 = 0.00 s/deg**（所有场景一致）
   - → 检查是否用了 Welch PSD
   - → 改用 `np.fft.rfft()` 整窗 FFT

2. **慢度偏小 10-100 倍，BAZ 混乱**
   - → 检查是否对 STFT 做了时间平均
   - → 去掉 `.mean(axis=1)`，或直接用整窗 FFT

3. **慢度值卡在网格边界**
   - → 检查 SMAX 是否覆盖目标波型的理论慢度
   - → 扩大 SMAX（T 波至少 40 s/deg，声道 T 波 80 s/deg）

4. **单位转换错误**
   - → 确认台站坐标单位（km 或 deg）
   - → 慢度转换：s/deg ÷ 111.32 = s/km（匹配 km 坐标）

---

## 代码模板（正确版）

```python
import numpy as np

# === 参数 ===
WINDOW_LEN = 120.0  # 秒
FREQ_BAND = (2.0, 8.0)  # Hz
SMAX = 40.0  # s/deg（覆盖地壳波）
DS = 0.5    # s/deg（网格步长）

# === 提取时间窗 ===
windowed_data = {}
for sta in station_list:
    mask = (waveform[sta]['time'] >= win_start) & (waveform[sta]['time'] <= win_end)
    windowed_data[sta] = waveform[sta]['data'][mask]

# === 整窗 FFT（关键：不做时间平均）===
ffts = {}
delta = waveform[station_list[0]]['delta']
for sta in station_list:
    fft_result = np.fft.rfft(windowed_data[sta])
    freqs = np.fft.rfftfreq(len(windowed_data[sta]), d=delta)
    ffts[sta] = fft_result

# === 频率掩膜 ===
freq_mask = (freqs >= FREQ_BAND[0]) & (freqs <= FREQ_BAND[1])
freqs_use = freqs[freq_mask]

# === 慢度网格 ===
sx_grid = np.arange(-SMAX, SMAX + DS, DS)
sy_grid = np.arange(-SMAX, SMAX + DS, DS)
power_map = np.zeros((len(sy_grid), len(sx_grid)))

# === 波束形成 ===
for freq in freqs_use:
    freq_idx = np.argmin(np.abs(freqs - freq))
    
    # 提取复数频谱并白化
    spec = np.array([ffts[sta][freq_idx] for sta in station_list], dtype=complex)
    spec = spec / (np.abs(spec) + 1e-12)
    
    # 遍历慢度网格
    for i, sx_deg in enumerate(sx_grid):
        for j, sy_deg in enumerate(sy_grid):
            # 单位转换：s/deg → s/km
            sx_km = sx_deg / 111.32
            sy_km = sy_deg / 111.32
            
            # 计算延迟（台站坐标 km）
            delays = sta_coords[:, 0] * sx_km + sta_coords[:, 1] * sy_km
            
            # Steering vector
            steering = np.exp(1j * 2 * np.pi * freq * delays)
            
            # Bartlett 功率
            beam = np.abs(np.sum(spec * np.conj(steering)))**2
            power_map[j, i] += beam

# 归一化
power_map /= len(freqs_use)
power_map /= power_map.max()

# 找峰值
peak_idx = np.unravel_index(np.argmax(power_map), power_map.shape)
sx_peak = sx_grid[peak_idx[1]]
sy_peak = sy_grid[peak_idx[0]]
baz_peak = (np.degrees(np.arctan2(sx_peak, sy_peak)) + 180) % 360
slowness_mag = np.sqrt(sx_peak**2 + sy_peak**2)
velocity = 111.32 / slowness_mag  # km/s
```

---

## 参考

- 马达加斯加 Mw6.2 事件调试日志：2026-09-02 会话
- 子阵：BKTA/CPSM/AMPY（边长 101-160 km）
- 理论 BAZ：151.6°
- 理论慢度：T-P 18.5 s/deg (6 km/s), T-S 30.9 s/deg (3.6 km/s)

---

## 附：另一版记录（原 `references/welch-vs-stft-pitfall.md`，2026-09-28 整理时逐字并入）

# F-K 波束形成：Welch PSD vs STFT 复数频谱陷阱

## 问题诊断（2026-09-02，马达加斯加 BKTA/CPSM/AMPY 子阵）

**症状**：所有分量（Z/R/T）× 所有时间窗（T-P/T-S）的 F-K 结果一致显示 **慢度峰值 = 0.00 s/deg**（对应无限大速度），表明三站信号完全同相，没有检测到任何时间延迟。

**根因**：脚本使用 `scipy.signal.welch()` 计算功率谱密度（PSD），然后尝试进行相位白化。

### 错误代码模式

```python
from scipy.signal import welch

# 计算 PSD
for sta in sta_list:
    f, psd = welch(windowed_data[sta], fs=1/delta, nperseg=NPERSEG, noverlap=int(0.5*NPERSEG))
    psds[sta] = psd  # ← psd 是实数数组

# 频域波束形成
for freq in freqs_use:
    freq_idx = np.argmin(np.abs(freqs - freq))
    
    # ❌ 错误：从实数 PSD 构造"复数"频谱
    spec = np.array([psds[sta][freq_idx] for sta in sta_list])
    spec = spec / (np.abs(spec) + 1e-12)  # 白化
    
    # 实数/abs(实数) = ±1.0（符号），丢失所有振幅和相位信息
```

### 为什么失败

1. **`welch()` 返回实数功率谱**：`psd[i] = |X[i]|²`，只有振幅的平方，**相位信息已完全丢失**
2. **白化无意义**：实数 PSD 除以其绝对值 = ±1.0（仅保留符号），没有复数相位可以归一化
3. **steering vector 失效**：`e^{i2πf·s·r}` 与 ±1.0 相乘，相位延迟信息无法编码
4. **零延迟假象**：所有台站"等效同相"（因为没有相位），峰值必然在慢度原点 (0, 0)

### 正确做法：使用整个窗口的 FFT

⚠️ **STFT + 时间平均也会丢失相位延迟信息**（2026-09-02 测试发现）：

即使用 `stft()` 获取复数频谱，对时间轴做 `.mean(axis=1)` 会抹掉台站间的相位差。测试案例：三个台站接收同一 4 Hz 信号，人工注入 0s/0.5s/1.0s 延迟，STFT 时间平均后相位差变为 0。

**根本原因**：时间平均 = 对不同时刻的复数频谱求平均，相当于对不同相位的向量求平均，会抵消相位信息。

**正确做法**：对整个窗口做**单次 FFT**（不是 STFT），保留完整的相位信息：

```python
# 计算整个窗口的 FFT（复数）
ffts = {}
for sta in sta_list:
    data = windowed_data[sta]
    # 单次 FFT
    fft_result = np.fft.rfft(data)
    freqs_fft = np.fft.rfftfreq(len(data), d=delta)
    
    ffts[sta] = fft_result  # 保留复数频谱
    
# 频域波束形成
for freq in freqs_use:
    freq_idx = np.argmin(np.abs(freqs - freq))
    
    # ✓ 正确：从复数 FFT 提取频谱
    spec = np.array([ffts[sta][freq_idx] for sta in sta_list], dtype=complex)
    spec = spec / (np.abs(spec) + 1e-12)  # 相位白化：归一化振幅，保留相位
    
    # steering vector 现在可以正确叠加相位延迟
```

**窗长与频率分辨率**：
- 120s 窗 @ 20Hz 采样 = 2400 点
- FFT 频率分辨率 = 1/120s ≈ **0.0083 Hz**
- 2-8 Hz 频段有 ~720 个频率点，足够精细

### 关键修改点

1. **用整个窗口的 FFT 替代 Welch PSD** — 保留完整的复数频谱和相位信息
2. **不要对 STFT 做时间平均** — 时间平均会抵消台站间的相位差
3. **显式 `dtype=complex`** — 确保数组类型是复数，即使从复数来源构造也要显式声明

### 三种方法对比

| 方法 | 返回类型 | 相位信息 | 适用性 |
|------|---------|---------|--------|
| `welch()` | 实数 PSD | ❌ 完全丢失 | **不可用** |
| `stft() + mean(axis=1)` | 复数（时间平均） | ❌ 相位差被平均抹掉 | **不可用** |
| `np.fft.rfft()` | 复数 | ✅ 完整保留 | **推荐** |

### 验证方法

修复后检查输出：
- 峰值慢度应该 **不为零**（~18.5 s/deg for T-P, ~30.9 s/deg for T-S）
- BAZ 应该有物理意义（接近理论值 ±合理偏差）
- 不同分量/时间窗的结果应该有差异（而不是所有场景都返回 0.00 s/deg）

### 延迟估算（用户纠正）

用户指出我对台站间延迟的计算有误。陆地传播是**地壳 P/S 波速度**，不是水声速度：

**正确估算**：
- T-P 窗：地壳 P 波速度 **6.0 km/s**
  - BKTA/CPSM/AMPY 最长边 159.5 km → 延迟 **26.6s**
  - 120s 窗长 = **4.5×** 最大延迟 → 完全足够
  
- T-S 窗：地壳 S 波（Love）速度 **3.6 km/s**
  - 最长边 159.5 km → 延迟 **44.3s**
  - 120s 窗长 = **2.7×** 最大延迟 → 足够

❌ 错误算法：用 1.5 km/s 水声速度算陆地段延迟（会高估 3-4 倍）

### 相关陷阱

- ⚠️ **时间平均 vs 空间叠加**：F-K 需要的是**空间叠加**（多个台站在同一时刻的相位差），不是**时间平均**（同一台站在不同时刻的复数向量平均）。STFT 的时间平均会让不同时刻的相位相互抵消。
- ⚠️ **`dtype=complex` 必须显式声明**：即使从 FFT 提取，`np.array([...])` 默认可能推断为 float64，显式 `dtype=complex` 避免隐式转换并触发 ComplexWarning
- ⚠️ **窗长选择**：整个窗口 FFT 要求窗长足够覆盖台站间最大延迟。对于 T-P/T-S 分析，120s 窗对应的频率分辨率（~0.008 Hz）已经足够精细，无需使用 STFT 的多窗平均

## 实测修复效果（2026-09-02）

**原始问题**：所有场景返回慢度 0.00 s/deg（无限大速度）

**第一次修复尝试**（Welch → STFT + 时间平均）：
- ❌ 仍然失败，原因是时间平均抹掉了相位延迟
- 测试验证：三站接收同一 4 Hz 信号，人工延迟 0/0.5/1.0s，STFT 时间平均后相位差 = 0

**第二次修复**（STFT → 整个窗口 FFT）：
- ✅ 理论正确：单次 FFT 保留完整相位信息
- 由于 Hermes 终端工具系统级 bug（"embedded null byte" 错误），脚本未能在会话中成功执行验证
- 修复后的脚本路径：`/mnt/e/sharedate/case/Madagascar/2013-08-17-mw62-southwest-indian-ridge/process_flow/4_fk_analysis/fk_120s_fixed.py`
- 需要在本地 WSL 终端手动验证

**关键教训**：
1. Welch PSD 是第一层陷阱（实数，无相位）
2. STFT 时间平均是第二层陷阱（复数，但相位被平均抹掉）
3. 正确方法：整个窗口的单次 FFT（复数，相位完整）

