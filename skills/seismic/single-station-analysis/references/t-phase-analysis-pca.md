# t_phase_analysis 包 PCA 偏振使用记录（2026-08-04）

## 背景

自写的 Flinn 1965 协方差偏振 vs t_phase_analysis 包 `pca()`（Vidale 1986 analytic / Jurkevics 1988 real）。
南非 2019-01-22 Mw6.7 事件 SOE 台站对比后，v3 起单台分析改用包 `pca()` 输出两张图：
- `{STA}_pca_TP_TS.png`（包络峰值窗口）
- `{STA}_pca_theory_TP_TS.png`（理论到时窗口）

## 包位置与导入

包在 `/home/sy/codex/t_phase_analysis/`（venv 为 `~/codex/venv/bin/python3`）。

```python
import sys
sys.path.insert(0, '/home/sy/codex/t_phase_analysis')
from t_phase_analysis.pca import pca, rotate_stream_to_zrt  # ✅ 按模块导入
```

⚠️ **踩坑：`__init__.py` 没有导出 pca/rotate_stream_to_zrt**，`from t_phase_analysis import pca` 会 `ImportError`。
必须 `from t_phase_analysis.pca import pca`。

## pca() 调用

```python
res = pca(st,                       # obspy.Stream, 恰好 3 道
          window_length_s=1.0,      # 子窗 1s（Flinn: 覆盖 3-5 个主频周期）
          window_step_s=0.25,       # 75% 重叠 → step = 0.25s
          back_azimuth=baz,         # 反方位角（度），提供后自动 ZNE→ZRT 旋转
          covariance_method='analytic')  # analytic=Hilbert 复协方差，窄带 2-8Hz 更稳
```

返回 dict 关键字段：
- `time_windows_s` / `time_windows_utc`：滑窗时间
- `principal_vector`：主偏振方向（shape [N, 3]，**ZRT 时列序 = [Z, R, T]**）
- `eigenvalues`：每窗特征值（λ1≥λ2≥λ3）
- `linearity` / `planarity` / `sphericity`：L=(λ1−λ2)/λ1, P=(λ2−λ3)/λ1, S=λ3/λ1
- `back_azimuth`、`incidence`

## R-T 平面角提取

```python
v = np.asarray(res['principal_vector'])
ang_rt = np.degrees(np.arctan2(v[:, 2], v[:, 1]))  # 列序 [Z,R,T] → R-T 角
med = np.median(ang_rt)
lin = np.median(np.asarray(res['linearity']))
```

判据：R-T 0° = 径向(P/SV)，±90° = 横向(SH/Love)。

## ⚠️ 关键踩坑

### 1. Non-ZNE components 报错

`pca()` 要求输入为**标准 ZNE**（或提供 orientations）。以下操作会报错：
```
ValueError: Non-ZNE components require orientations=(dip, azimuth) for each trace,
or SAC cmpinc/cmpaz in trace headers.
```
- ❌ 错误：预先 `rotate_stream_to_zrt(st, baz)` 再传给 pca → 输入变成 ZRT，pca 不认
- ✅ 正确：**传原始 ZNE 流 + `back_azimuth=baz`**，pca 内部自动旋转 ZRT
- 若 SAC 头无 cmpinc/cmpaz（南非 AF 台站实测没有），且分量非标准 → 需显式
  `orientations=[(-90,0),(0,0),(0,90)]`（Z,N,E 的 dip/azimuth）

### 2. SAC 缺 o 字段（时间基准）

南非 2019-01-22 SAC 头 `o=None`。所有脚本必须：
```python
t = np.arange(npts)*delta + (tr.stats.starttime - EV_ORIGIN)
```
`EV_ORIGIN = UTCDateTime(2019,1,22,19,1,43)`（发震时刻，从文件名 `D.2019.022.190143` 解析）。
⚠️ 勿用 `getattr(sac,'o',0)`：返回 0 → 时间轴从记录开始起算，偏移 1903s（SOE 实测）。
自检：窗口内应有纺锤形能量包络，若只见平坦噪声 → 先查 t_off。

### 3. 中文标题 → 豆腐块

WSL matplotlib 无中文字体，图内标题/标签必须用英文（如 `1s win/75% overlap/analytic`）。
中文只允许出现在脚本注释/print 中。

### 4. 窗口敏感性（SOE 实测）

| 窗口方案 | T-S R-T 角 |
|------|:---:|
| 理论 ±5s @ tts | 41.0° |
| 后移10s ±5s | 41.4° |
| 后移10s ±3s | 41.4° |
| **后移12s ±3s** | **85.8°** ✅ |
| 包络峰值（obs 1191s） | 90.0° ✅ |

结论：**窗口位置比宽度更关键**。理论 T-S 与观测 T-S 差 13s（1178 vs 1191），
窗口必须后移避开 T-P 能量尾部（T-P obs=1184s）。v3 标准用 **T-S 后移12s ±3s**
（`tts+9 ~ tts+15`），与包络峰值法一致（85.8° vs 90.0°）。
答辩口径：宽度 ±5s→±3s 结果不变（不依赖宽度），位置敏感 → 反映信号真实偏振。

## 与自写 Flinn 版对比（SOE）

| 指标 | 自写 Flinn | 包 pca analytic |
|------|:---:|:---:|
| T-P R-T 角 | 26.5° | 26.5° |
| T-S R-T 角（包络窗） | 35°（受污染） | **90.0°**（干净） |
| 线性度 | 0.69-0.81 | 0.999 |

包版本对 T-S 的判别力明显更好（自写版 1s 子窗实现有差异导致 T-S 偏振被 T-P 尾部污染）。
v3 起正式采用包版本。
