# 方位/慢度的“互相关滞后拟合”是什么、出处、怎么定位（用户问过“有文献出处吗”）

## 1. 它不是论文的步骤，是本项目的附加估计器
Zhou2021 只报 F-K 慢度图（Fig.4，三分量 Z/R/T）。带符号台对互相关滞后拟合（`xcorr_direction_check.py` 一族）
是本项目加的，定位 = **方位/慢度的独立估计 + 误差条**，不替代 F-K，也不当“对得上/对不上”的唯一判据。
用户质疑过“这是你自己加的检查”→ 必须这样直说，并给出下面的出处；他不想要就不当判据。

## 2. 文献出处（都是标准方法，可查）
| 用途 | 出处 |
|---|---|
| 互相关定台间时差 + 最小二乘 → 平面波拟合（= 滞后拟合本体） | **Schweitzer et al. (2012), “Seismic Arrays”, NMSOP-2 Ch.9**, doi:10.2312/GFZ.NMSOP-2_ch9 —— §9.8 “Slowness estimation using seismic arrays” 内含 *plane wave fitting – slowness estimation by time picks* 与 *the correlation method used at the UKAEA arrays* 两小节；并写明 “slowness resolution … improves with increasing the array aperture” |
| 多道互相关 + 最小二乘求精确相对到时（MCC） | **VanDecar & Crosson (1990), BSSA 80(1), 150–159** |
| 互相关定背方位与视速度的成型算法（PMCC） | **Cansi (1995), GRL 22(9), 1021–1024**, doi:10.1029/95GL00468 —— “progressive study of the interstation correlation functions … consistent set of time-delays” |
| 为什么 F-K argmax 的方位不准（主瓣宽度量级） | **Rost & Thomas (2002), Rev. Geophys. 40(3)**, doi:10.1029/2000RG000100 —— 式(30)/(34) `E(k−k0)=|A(k−k0)|²`；§[47]–[52]、Fig.13；§[51]：“the maximum of the power spectral density is shifted to the corresponding values, **without changing the form of the ARF**”；§[82]：“sidelobes due to spatial aliasing should not be within the wavenumber window of interest” |

## 3. 主瓣半宽系数 C 必须写明口径（用户抓过这点）
`Δu = C·111.195/(f·D)`，无量纲量 `X = f·D·u`（`B(x)=|mean_i exp(−i2π f x·rr_i)|`）：

| 口径 | X（理论，连续矩形孔径） | G1 七站实测 C | 说明 |
|---|---|---|---|
| 第一零点（主瓣边界） | 1.0000 | — | 主瓣全宽的 1/2 |
| 半幅 `|B|=0.5` | 0.6034 | **0.6096**（5 个频点一致，与 f 无关） | 项目通用值取 **0.64**（保守 +5%） |
| 半功率 `|B|²=0.5`（−3 dB） | 0.4429 | 0.4502 | 报误差条常用 |

各向异性：每 5° 扫一遍 C ∈ [0.556, 0.775]（中位 0.619）。
G1（D=312.4 km、d_min=101.0 km）例：f=0.0117 Hz → 半幅 19.5 / 半功率 13.5 / 零点 30.4 s/°；
f=0.018 Hz → 12.7 / 8.8 / 19.8 s/°。**报数时必须写“半幅”还是“−3 dB”**，否则等于没给不确定度。
⚠️ 论文里没有这个系数（Zhou2021 全文无 resolution/mainlobe/sinc 类说法）⇒ 报数、写稿一律标“本项目自检项，非论文参数”，
概念出处写 Rost & Thomas (2002) 的 ARF。

## 4. 用户已定的事项（别再回头质疑）
- 三分量都要做：论文 Fig.4 就是 Z/R/T 三分量（**T-S 显现在 T 分量、T-P 显现在 R 分量**，Z 作对照）。
- 子阵名单用户已定（= 论文所用，靠台站位置即可对比）⇒ **不要再反复强调“名单不确定”“要注明名单”**。
- “预测慢度由射线理论独立算出（用户自己会推）”⇒ **不要把它叫作循环论证**；只有“拿论文实测值回头调参数”才算循环。
- **滞后拟合不进判据**（用户 2026-09-23 定）：主线只需“理论预测 |u|（前向，射线理论）↔ F-K 峰的 |u| + 峰显著性（prom）”，
  物理逻辑闭环即可，**不做多方法交叉验证**。滞后拟合仅在“要把方位（BAZ）当数字报/画方位线”时才补（F-K argmax 方位
  误差 = 主瓣半宽/|u| ≈ ±30° 量级），且是事后计算、不进主流程。因“平图/无峰”而需归因时，改用合成注入标定，
  不用滞后拟合。
