# 偏振引导的 T-P/T-S 自动拾取（高 SNR 台站）

## 目的
理论走时仅作宽先验，不把固定理论中心或固定后移量当作实际到时。对高 SNR 三分量记录，在候选搜索窗内滑动扫描，以偏振分类和信号质量共同选择 T-P/T-S 窗口；输出 CSV 供后续统计、走时拟合和人工复核。

## 适用门槛
- 仅处理完整 Z/N/E 三分量台站；
- 总体 2–8 Hz SNR >= 10 dB 进入扫描；
- 最终 T-P、T-S 各自单波包 SNR >= 10 dB；
- 不达标写入 CSV 并标记 rejected，不强行给出拾取。

## 预处理与滑窗
- 2–8 Hz 零相位带通；
- 窗口宽度：默认 ±3 s（总长 6 s）；滑动中心步长 0.25 s；
- PCA 子窗：1 s，75% 重叠（步长 0.25 s）；
- 波形展示：Z/R/T 共用全局 P98；灰色波形 + 分量彩色 Hilbert 包络，包络 4 s 平滑。

## 搜索逻辑
For this user's routine batch processing, theory defines a fixed, symmetric prior first:
- T-P: theoretical T-P ±15 s;
- T-S: theoretical T-S ±15 s.

Within each prior, scan fixed ±3 s analysis windows at 0.25 s center steps. Select the qualifying candidate with the smallest radial residual for T-P or transverse residual for T-S. Do not let a broad search wander into a neighboring packet merely because its polarization is cleaner. A candidate at the ±15 s boundary is a review flag, not an automatic physical onset; refine the travel-time model only after cross-station residual fitting.

每个滑动中心计算：
- R–T 中位角、径向残差、横向残差；
- 线性度及角度离散度；
- R/T 能量比；
- 单波包 2–8 Hz SNR；
- 包络起跳/局部峰显著性；
- 与理论到时残差。

**窗口选择原则**：不能把“径向/横向残差的全局最小窗口”直接视为初动。应先找满足质量阈值的连续偏振簇，再选**第一个连续簇的起始窗口**作为实际初动候选；簇内角度最优窗口仅记录为主偏振峰，不替代初动。这样避免 AMPY 中较晚、更纯径向的峰值窗口取代较早且已合格的 T-P 起始能量。

### T-P 选择
在 T-P 搜索窗内选择连续稳定的径向候选，优先：
1. 径向轴向残差最小（目标 0/180°）；
2. 线性度高；
3. R 能量相对 T 更强；
4. 单波包 SNR 高；
5. 相邻滑窗连续、方向稳定；
6. 走时残差较小。

### T-S 选择
在 T-S 搜索窗内选择连续稳定的横向候选，优先：
1. 横向残差最小（目标 ±90°）；
2. 线性度高；
3. T 能量相对 R 更强；
4. 单波包 SNR 高；
5. 相邻滑窗连续；
6. 严格晚于选定 T-P 至少 3 s。

## 质量判据
推荐 accepted 条件：
- T-P/T-S 单波包 SNR >= 10 dB；
- 线性度 >= 0.8；
- 至少 3 个连续的 0.25 s 滑动中心满足分类阈值；
- T-P 径向残差 <= 20°；
- T-S 横向残差 <= 20°；
- 角度离散度或 IQR <= 20°；
- T-S 在 T-P 后 >= 3 s。

## CSV 输出字段
每台站一行；必要字段：
`sta, snr_total, tp_theory, ts_theory, tp_pick_center, tp_window_start, tp_window_end, tp_snr, tp_rt_angle, tp_radial_residual, tp_linearity, tp_rt_energy_ratio, ts_pick_center, ts_window_start, ts_window_end, ts_snr, ts_rt_angle, ts_transverse_residual, ts_linearity, ts_tr_energy_ratio, tp_score, ts_score, accepted_tp, accepted_ts, reject_reason`

## 重要限制
- 偏振引导拾取可用于 T-P/T-S 波型和实际波包窗口识别，不能把“在理论 BAZ 旋转后的 R–T 径向性”当作独立反方位角反演；反方位角需另用未旋转 N–E 偏振轴并与转换点几何比较。
- 不允许为每站手工固定后移秒数作为通用模型；后移量应由扫描结果、CSV 质量标志与跨站稳健走时拟合约束。
