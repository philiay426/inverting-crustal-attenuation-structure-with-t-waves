# Flinn (1965) 协方差矩阵极化分析法

## 原理

在分析窗口内滑动子窗口（默认1s，75%重叠），每个子窗口算3×3协方差矩阵并特征分解，用主方向做玫瑰图，线性度做权重。

## 步骤

1. **切子窗**：`[ws, we]` 内以 `sub_win=1s`、`step=0.25s` 滑动
2. **去均值**：`X = [Z-μz, R-μr, T-μt]` → 3×N
3. **协方差**：`C = X·Xᵀ / N` → 3×3
4. **特征分解**：`λ₁≥λ₂≥λ₃`, 对应 `v₁, v₂, v₃`
5. **线性度**：`Lin = 1 - (λ₂+λ₃)/(2λ₁)`  (Flinn 1965)
   - Lin≈1：质点走直线(P波/强偏振)
   - Lin≈0：质点乱晃(噪声)
6. **方向**：仅 R-T 平面 `θ = atan2(v₁_T, v₁_R)`（v3，移除 Z-R/Z-T）
   - θ ≈ 0°：径向（P/SV 波）
   - θ ≈ ±90°：横向（SH/Love 波）
7. **玫瑰图**：所有子窗的角度加权叠加（权重=Lin），max归一化

## 与旧方法对比

| 旧 (逐采样点arctan2) | 新 (Flinn协方差) |
|:--------------------|:-----------------|
| 每采样点一个方向 | 每1s子窗一个主方向 |
| 振幅²为权重 | 线性度为权重 |
| 噪声主导分布 | 主导偏振方向 |

## 实现

代码见 `templates/polarization_analysis.py`:
- `covariance_polarization(t, z, r, tc, ws, we)` — 子窗滑动，仅返回 R-T 角度和权重
- `rose_from_covar(angles, weights)` — 角度→玫瑰图直方图

## 已知坑

- SAC 文件可能缺 `sac.o`，且 **o 存在时优先用 b-o**（2013 事件 `o=600` 有效，trace 从发震前 600s 开始，`t_off = b - o = -600`）；缺 o 时用 `t_off = tr.stats.starttime - EV_ORIGIN`（发震时刻从文件名解析）。⚠️ 勿用 `getattr(sac,'o',0)`（缺 o 返回 0 导致时间轴偏移）
- SAC 文件可能缺 `stla/stlo`，需从 inventory 获取

## 参考

Flinn, E. A. (1965). Signal analysis using rectilinearity and direction of particle motion. Proceedings of the IEEE, 53(12), 1874-1876.
