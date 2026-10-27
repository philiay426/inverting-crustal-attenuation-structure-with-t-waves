# 偏振玫瑰图：统计方法与绘制规范（SA18 单台实战，2026-09-08）

Zhou et al. (2021) Fig 2e/f 的 polarization rose 复现/对照规范。用户会问
"这是按协方差算的吗？论文好像是按振幅算的" —— 先如实交代方法，再给对比图。

## 论文原文事实（查证过）

- 正文 3.1 + SI 只说 "segmented ... at 25 s before and after the theoretical arrival
  time and measured the polarization of the waveforms (Figures 2e and 2f)"，**未给出
  协方差 / 振幅的具体公式**；SI 只有 back-projection (Text S1) 与 coda 衰减 (Text S2)。
- 因此"论文用哪种统计"无法从文献确认 —— 不要替论文断言；正确做法是做**多方法对比图**
  让用户自己判断哪种更贴论文图。
- 论文窗口径：**理论 T-P / T-S 到时 ±25 s**（理论值用 Equation 1，不是 stations.csv 的
  TauP 到时；TauP 偏早 ~20 s，用它当窗中心会罩不住真实 T-S 能量峰）。

## 三种统计方法（同一张图三行对比，本项目定稿产物）

脚本：`process_flow/2_single_station_analysis/SA18/plot_SA18_rose_methods.py`
（输出 `SA18_rose_methods_compare.png`，3 行 × 2 列，左 T-P 右 T-S）

| 行 | 统计 | 实现 |
|---|---|---|
| Covariance (PCA) | 1 s / 0.25 s 滑窗 N/E 协方差 `np.linalg.eigh`，取最大特征向量方向，λ1 加权 | 对每子窗算主轴，一窗一票 |
| Amplitude (R²+T²) | 逐样本瞬时方向 `atan2(E,N)`，权重 N²+E² | 能量加权 |
| Sample count | 逐样本瞬时方向 `atan2(E,N)`，每样本 1 票 | 等权计数（用户 2026-09-08 最终要的形态）|

要点：
- 三分量的协方差版要用 t_phase_analysis 的 analytic（Hilbert 复协方差）；二维版直接用 N/E。
- 协方差法对 **±25 s 长窗** 会被窗内混合能量带偏（SA18 ±25 s T-S 中值 1.3°，而 ±2 s
  拾取窗是 95° 横向）→ 长窗请用振幅/计数法，扇区占比判别更鲁棒。
- 方向判别用**扇区能量/计数占比**（相对 BAZ 轴 ±22.5° 的 radial vs BAZ±90 的 transverse），
  比主轴角中值稳：SA18 振幅法 T-P radial 38% vs trans 14%；T-S trans 43% vs radial 8%
  （与论文"T-P 径向、T-S 横向"一致）。
- 占比计算时注意：`(w*mask).sum()` 是加权和（~1e-16），除以 `w.sum()` 才是百分比
  （曾直接 f"{er:.1%}" 打印出 0.0% 的假 bug）。

## 绘制规范（用户 2026-09-08 逐轮定稿）

1. **恢复旋转前方向**：玫瑰画在原始 N/E 地理方位（0°=北），不要在 R/T 旋转坐标系里画；
   R/T 版会把能量"钉"在径向轴上，N/E 版能看出实际到达方向相对理论 BAZ 的偏离。
2. **全 360° 全窗**：`phi = (degrees(atan2(E,N)) % 360)`，bins `linspace(0,360,37)`；
   不要只画有向半圈。
3. **理论方向线**：T-P 子图只标**理论径向**（BAZ）、T-S 子图只标**理论横向**（BAZ+90）；
   用**黑色实线**，且必须是**过圆心的完整直径**（`plot([rad(a), rad(a+180)], [mx, mx])`），
   不要只画半轴（用户原话"黑色实线不要只画一半"）。
4. **文字精简**：总标题一行（`SA18 polarization (N/E, theory ±25 s, 2-8 Hz)`）；子图标题只两行
   `T-P` / `radial 19%`（每张只报与该相位相关的那个占比）；极坐标刻度只标 N/E/S/W，
   不要 45° 全标。用户明确嫌"字太多了"。
5. 台站 BAZ 从 `stations.csv` 的 `baz` 列取（SA18=126.8°），不要硬编码。

## 处理过程对比图（env chain）规范

`plot_env_chain_subs.py` → `env_chain_G1_southern.png` / `env_chain_G2_northern.png`：
- 3 行 = 2-8 Hz bandpass（灰）→ |Hilbert| envelope（蓝）→ 0.25 Hz low-pass（红，送 F-K）；
- **台站排列 = F-K 分析输入顺序**（论文子阵名单顺序：G1=SA16/17/18/22/23/25/29，
  G2=SA31/32/37/38/39/40/45/46），**不要按 land_km 排序**（用户 2026-09-08 要求）；
- 每台标理论 T-P（蓝 #1f77b4）/ T-S（紫 #9467bd）**短竖虚线**（行内 ±0.45），三行都标；
- 横轴陆地计时（`t − ocean_km/1.48`），到时用论文公式换算成同一轴。

## 相关

- 慢度图样式：`references/slowness-figure-paper-style.md`
- 单台偏振（T-P/T-S 窗 PCA、引导拾取）：user-owned `single-station-analysis`
- SEIZMO / Octave 环境：`references/octave-seizmo-fk4d.md`
