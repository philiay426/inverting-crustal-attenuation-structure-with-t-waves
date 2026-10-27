# SNR 数据源必须用去响应未带通波形（2026-09-02 定，勿双重带通）

> ⚠️ **2026-10-14 起已被取代**：现行口径 =「统一噪声基准 SNR」——2–8 Hz **整段带通**后取窗内 RMS，
> 除以**阵列噪声中位数**（见 seismic-twave-analysis `references/snr-unified-noise-baseline.md`）。
> 本文的「未带通 + Welch 频段积分」不再适用；其中**唯一仍有效的教训**是：
> npz 里 `waveform_*` 已是 2–8 Hz 缓存，**不得再滤一遍**（双重带通扭曲 SNR）。

## 结论

SNR 应基于**去响应后、未带通**的速度波形计算（`sacprocess/remove_response/*_resp.sac`），
用 Welch PSD 后只积分 2-8 Hz 频段。**频段积分本身就是频域带通**——
`SNR(dB)=10log10(∫₂⁸P_sig df / ∫₂⁸P_noise df)`。

先做时域 2-8 Hz 带通、再 Welch 积 2-8 Hz = 双重带通，会引入滤波器滚降/旁瓣/振铃，
改变频段内功率比例 → SNR 数值被扭曲且方向不定。

## 实测差异（马达加斯加 2013-08-17，同站同窗）

| 台站 | 未带通 resp | 带通 npz | 差值 |
|---|---|---|---|
| VOI | 39.6 | 39.2 | −0.4 |
| AMPY | 29.5 | 32.9 | +3.4 |
| KIRI | 18.7 | 13.0 | −5.7 |
| SOLA | 31.9 | 32.5 | +0.6 |
| BATG | 26.5 | 25.5 | −1.0 |
| CPSM | 12.3 | 10.8 | −1.5 |

−5.7 ~ +3.4 dB、无系统性方向 → 在 3/6/8 dB 筛选阈值边缘的台站会改变筛选结果。

## 项目现状（不一致，需统一）

- ✅ **stations.csv（step1_map_and_table.py）**：读 `remove_response/*_resp.sac` 未带通 →
  标准口径，正确。复现验证：VOI 39.6 vs csv 39.8、CPSM 12.3 vs csv 12.3 吻合。
- ❌ **剖面/频谱 SNR 筛选（plot_unbinned_waveform_profiles.py）**：Z 分量从
  `sacprocess/envelope/*_envelope_2-8Hz.npz` 取 `waveform_z`（**2-8 Hz 带通后缓存**）→
  双重带通。⚠️ npz 缓存（waveform_*/envelope_*）只可用于画波形/包络，**不得当 SNR 数据源**。
- ❌ 南非台站表（write_sa_station_table）内 **Z 走 npz（带通）、R/T 读 resp 未带通** →
  同一张表内部口径不一致。

统一方案：所有 SNR（Z/R/T）改读 remove_response resp.sac；npz 保留画图用途。
修改点：plot_unbinned_waveform_profiles.py 的 compute_snr_welch_time 调用处
（Z 分量改为与 R/T 相同的 resp.sac 读入 + rotate 流程）。

## 验证脚本思路

两种口径都按同一窗口定义（信号 [tt_tp-20, tt_ts+100]、噪声 P 前
max(20,min(300,0.8·tt_p))）算，时间轴各自用 SAC o 的 t_off（resp）或 npz time 数组；
nperseg=min(256, nwin)。差异 ±6 dB 量级说明口径必须统一，不能混用。

## 注意

FOMA/ABPO 的 snr_r/snr_t 列为空 ≠ 缺水平通道——resp 文件存在且数据正常，
R/T 可算（2026-09-02 核查过）。空值只是当时没算，不是物理缺失。
