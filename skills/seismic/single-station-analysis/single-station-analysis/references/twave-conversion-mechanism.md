# T 波接收端波型识别（P/SV vs Love）

2026-09-15 由独立 skill `twave-conversion-mechanism` **收缩并入**本条参考（原 skill 已归档至 `.unused/seismic/`）。

## 关键结论

- **T-P 起震为径向（P/SV）**；Love 波（若存在）在尾波中逐步发育，无尖锐初至；三要素闭合验证：**到时 + 偏振 + 波列位置**（缺一不可）。
- 判据集中在 **R-T 平面**：R-T 0° = 径向 = P/SV，R-T ±90° = 横向 = SH/Love；Z-R 平面已移除（R-T 足以区分波型，且 Z 分量故障时仍可分析）。

## 区域差异

| 区域 | T-P | T-S / T-Love | 原因 |
|---|---|---|---|
| PR2014（波多黎各） | ✓ P/SV | ✓ Love 波 | 大陆坡缓 → 发生 Love 转换 |
| 南非（Prince Edward） | ✓ P/SV | ✗ **无 Love** | 大陆坡陡 → 仅 P/SV |

- 与 De Groot-Hedlin (2020) 一致。文献：AELUMA 台阵、Love 波识别、Hawaii 2006 对照、De Groot-Hedlin & Orcutt (2001)。
- 南非 Mw6.7 (2019-01-22) 实测：33 站 AZ 302.6–329.6°（连续单一群），仅 T-P 无 T-S → `south-africa-mw67-2019-polarization.md`。
- ⇒ 南非这类事件的 T 波**只需确认径向（R-T ≈ 0°）**，不必区分两种波型。

## 方法（口径演进）

| 代次 | 方法 | 说明 |
|---|---|---|
| 历史 v3 | Flinn 1965 协方差，窗 ±5 s、子窗 5 s/75% 重叠 | 玫瑰图权重 = 线性度 |
| **现行默认** | **采样点等权（count，每点 1 票）** | 逐采样点 `atan2(E,N)`，无幅度权重、无窗参数；三张默认输出图用它 |
| 可选 | 振幅加权（amp，权重 N²+E²）/ 协方差（cov，λ₁） | 三方法都保留，脚本开关 `--method count|amp|cov` |
| 拾取 | 偏振引导（理论 ±15 s 内扫子窗，取最径向/最横向，半宽 ±2 s） | → `polarization-guided-picking.md` |

## 输出（2026-09-15 用户定：默认只出三张）

1. `{STA}_PS_polar.png` — P、S 理论到时偏振
2. `{STA}_pca_theory_TP_TS.png` — T-P、T-S 理论到时偏振
3. `{STA}_polarization_fullscan.png` + `.csv` — 全时段扫描：偏振方向、线性度、平面度

## F-K 与波型识别的关系（非洲南部实测）

**33 站（~100 km 间距）任何 F-K 方法都不可行**（传统 PMCC / 包络 PMCC / 三站延迟拟合 / 延迟扫描全部失败）；
唯一可靠的反方位角验证是**偏振 R-T 平面**（单台三分量即可，不需要台阵）。
Zhou (2021) 能做 F-K 是因为 82 站 + ~400 km 孔径的统计堆叠，不是单组 F-K。
→ `fk-array-analysis/references/south-africa-33stn-fk-feasibility.md` 与 `zhou2021-envelope-fk-method.md`。
