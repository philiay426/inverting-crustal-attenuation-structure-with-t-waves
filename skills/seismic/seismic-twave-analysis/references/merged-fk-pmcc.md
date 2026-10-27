# FK/PMCC 分析笔记

## 原理简述（面向非技术读者）

**核心问题**：怎么知道波从哪个方向来？

**类比**：人有两只耳朵，声音先到左耳→声源在左边，同时到→在正前方。地震台阵也是同样道理——用多个"耳朵"（台站）接收波的到达时间差，反推波的来向。

**PMCC 三步法**：
1. **测时间差**：3 个台站两两做互相关，得到 A-B、A-C、B-C 三组延迟
2. **网格搜索**：假设平面波，扫描所有可能的方向 (0°–360°) 和速度 (1200–1800 m/s)，对每个组合计算理论延迟：`dt = −(台站间距 / 速度) × cos(台站连线方位 − 波来向)`
3. **找最佳匹配**：哪个（方向, 速度）能让理论延迟与实测延迟最接近，就是答案。所有有效窗口取众数得最终结果。

**名称澄清**：代码中叫 "FK"，实际是 **PMCC**（时域波束形成），不是传统频率-波数域的 F-K 分析。代码位置：`/home/sy/codex/t_phase_analysis/t_phase_analysis/internals/PMCC/`。

**⚠️ 解释风格**：对该用户解释 PMCC/FK 原理时，用物理类比（耳朵判方向），不要展示代码。

---

# TA 台阵 FK/PMCC 局限性 — PR2014 实测验证

## 测试矩阵

对 PR2014 350° 群，使用 K62A 子台阵（K62A+K63A+L61B+J62A+L63A, 间距 44-150 km）验证多种台阵方法：

| 方法 | 配置 | 结果 | 反方位角偏差 |
|:---|:---|:---|:---|
| FK (ObsPy array_processing) | 3站, 2-8 Hz | 无有效窗口 | — |
| PMCC (triad) | 3站, Z分量, 10s窗 | 0% 有效窗 | — |
| 包络互相关 | 3站, T分量, 100s全窗 | dt=0.00s (CC>0.95) | 168° |
| 包络互相关 | 3站, 10s滑动窗 | dt=0.00s (18窗) | 168° |
| 包络互相关 | 5站, 10s滑动窗 | dt=0.00s (10对) | 171° |

**所有方法均失败**。反方位角偏差 168-171°，视速度 999 km/s（无解）。

## 失败原因

1. **包络过度相干**：2-8 Hz + 2s 平滑后，T 波包络在 50-70 km 站间距下几乎完全相同（CC 0.94-0.99），互相关峰值始终在滞后 0 处
2. **延迟小于分辨率**：预期台站间延迟 ~5-20s (v=3.6 km/s)，但全窗互相关被包络整体形状主导，细微时移被淹没
3. **不同转换点**：每个台站的大圆路径不同，T 波到达时经历了略微不同的海洋路径和转换点，波前不满足平面波假设

## DGH (2020) 为何成功

DGH 的 AELUMA 方法与我们的尝试有本质区别：

| | AELUMA (DGH 2020) | 我们的尝试 |
|:---|:---|:---|
| 三元组数 | 数百个重叠 | 1-10 个 |
| 测量量 | 沿阵传播速度 | 台站对延迟 → 慢度矢量 |
| 统计方法 | 数百三元组堆叠 | 单次最小二乘 |
| 滤波器 | 1-3 Hz | 2-8 Hz |

AELUMA 测的是 T 波波前从南到北的整体传播速度（大尺度测量），不依赖单对互相关精度。我们试图从台站对延迟反推反方位角，属于空间混叠严重的小尺度测量。

## 台站间距阈值

基于 PR2014 实测，TA 台网 FK/包络互相关的可行间距上限：

| 测量目标 | 最小间距 | TA 平均间距 | 可行性 |
|:---|:---|:---|:---|
| 沿阵传播速度 (AELUMA) | <125 km | 70 km | ✅ |
| 反方位角 (单个 triad) | <10 km | 70 km | ❌ |

## 论文表述建议

_"While traditional array methods (FK, PMCC) proved unsuitable for back-azimuth estimation at USArray station spacing (~70 km), the geometric conversion-point method provides a physics-based alternative that is fully determined by the great-circle path and bathymetry, requiring no free parameters."_

## 参考文献
- De Groot-Hedlin (2020), SRL 91(6), 3444–3453
- Rost & Thomas (2002), Rev. Geophys. 40(3), 1008
- De Groot-Hedlin & Hedlin (2015), BSSA 105(2A), 784-794
- Zhou et al. (2021), GRL 48, e2021GL094410 — XA阵列 100 km 间距 F-K（见下节）

## Zhou et al. 2021 GRL — 100 km 间距为何能做 F-K

用户质疑：Zhou et al. 2021 用 XA 阵列（间距 ~100 km）做 F-K 验证 T-P/T-S，2-8 Hz 为何可行？

**四个原因**：
1. **阵列总孔径大**（1800 km）→ 角分辨率极高，补偿混叠
2. **82 站统计堆叠**→ 假峰（栅瓣）在不同台站对位置不同，叠不出；真峰一致
3. **震中先验约束**→ 已知反方位角，可排除错误峰
4. **内陆信号偏低频**→ 高频衰减后实际信号 <1 Hz，波长更长

**论文用 SEIZMO 做频率域 F-K，不是 PMCC**。PMCC（你的代码）是 3 站时域互相关；SEIZMO F-K 是 82 站频率域波束形成。原理不同。

# FK/PMCC for TA Array — Limitations

## Conclusion
TA station spacing (~50-70 km) is too large for triad-based beamforming (FK, PMCC, envelope CC).

## Evidence (PR2014)
- PMCC K62A+K63A+L61B: 0% valid windows
- Envelope CC 10s windows: all delays 0.00s (CC 0.94-0.99)
- 5-station LS slowness: BAZ=0.0, v=999 km/s

## Why DGH (2020) succeeded
- AELUMA: hundreds of overlapping triads statistically stacked
- Measured along-array velocity, not back-azimuth
- N-S line array: good N-S resolution, poor E-W

## References
- De Groot-Hedlin (2020), SRL 91(6)
- Rost & Thomas (2002), Rev. Geophys. 40(3)

# FK/PMCC 分析集成

## 环境
- PMCC 脚本位于 `t_phase_analysis/internals/PMCC/`
- `__init__.py` 自动在 `sys.path` 中注入 PMCC 路径
- 依赖 `future` 包：`pip install future`
- 统一重采样到 40Hz（`run_fk.py` 自动处理）

## 关键修复（azimuth_speed_array.py）
1. `plt.subplots(5, figsize=(14, 16))` → 加 figsize
2. `plt.subplots_adjust(hspace=0.05)` → 减少子图间距
3. 删 `plt.show()` → 非交互模式会卡死
4. 删 `../../` 路径 → `savefig` 存当前目录，加 `bbox_inches='tight'`

## 通用 FK 工具
```bash
python /home/sy/codex/t_phase_analysis/t_phase_analysis/internals/PMCC/run_fk.py \
  --proc-dir <sacprocess> --out-dir <output> \
  --stations S1,S2,S3 --t0-times t1,t2,t3
```

## PMCC 要求
- 必须恰好 3 个台站
- 台站间距不要太大（<50km，否则 T 波不相干）
- T0 用 T-P 理论到时（从 CSV 读 tt_tp）
- 参数：`--start 0 --slice 20 --slide 10 --total 100`

## 结果解读
- Consistency：时延三边闭合度，越低越不可靠
- 有效解：azimuth ≠ -12345, speed ≠ -12345
- 速度 1.5-1.8 km/s = T 波 SOFAR 通道

# PMCC FK Analysis — Integration & Usage

## Package Integration

PMCC internals live in `t_phase_analysis/internals/PMCC/`:

```
internals/PMCC/
├── __init__.py                   # NEW: enables import
├── azimuth_speed_array.py        # Core FK engine (fixed: plt.show(), ../../ path)
├── cross_correlation_WENBO.py    # CC engine (requires `future` package)
├── global_parameters_PMCC.py     # t0 = reference_time['Pwave']
├── PMCC_auxiliary.py
├── run_fk.py                     # NEW: generic CLI utility
```

## Fixes Applied (July 2026)

1. **plt.show() removed** — hangs in non-interactive mode
2. **Save path fixed** — was `../../` now saves to current directory
3. **Auto-resampling added** — all traces resampled to 40Hz to fix "Sampling rates do not match"
4. **PMCC dir added to sys.path** — via `t_phase_analysis/__init__.py`
5. **`future` package installed** — required by cross_correlation_WENBO

## Usage

```bash
python /path/to/run_fk.py \
  --proc-dir <sacprocess_dir> \
  --out-dir <output_dir> \
  --stations S1,S2,S3 \
  --t0-times 959,960,960
```

Required: exactly 3 stations. Station spacing should be ~1-20km for T-wave frequencies.

## Output Files

- `n_t0_*_*.png` — PMCC summary plot (Consistency, CC, Azimuth, Speed, Waveforms)
- `n_*_*.npz` — Raw numpy arrays
- `file_list.txt` — SAC file list
- `*.BHZ.sac` — Resampled SAC files with t0 set

## Interpretation

| Panel | Meaning |
|-------|---------|
| Consistency | Delay closure residual. **Lower = better** (< 1 is good) |
| Avg CC | Average cross-correlation across 3 station pairs |
| Azimuth | Back-azimuth (°). -12345 = invalid |
| Speed | Apparent velocity (m/s). -12345 = invalid |

## Pitfalls

- **CTAO+BW1H+BW2S failed**: stations too far apart (~200km), T-wave incoherent
- **RK2S+AUNRC+RK1H partially worked**: only 1/10 windows valid, array nearly linear
- PMCC requires **40Hz** input; auto-resampling handles conversion
- Data must be **deconvolved velocity** (proc.sac files from process_obs.py)

## 实测案例 (2026-07-23)

### 澳大利亚 Mw7.0 事件
| 台站组 | 结果 |
|--------|------|
| RK2S + AUNRC + RK1H | ✅ 1/10窗有效：BAZ=284.7°, v=1799 m/s |
| CTAO + BW1H + BW2S | ❌ 全部无效（台站间距~200km，不相干） |

### 南非 Mw6.3 事件
RK2S/AUNRC/RK1H 三站间距仅 ~8km，线状排列，FK分辨力有限。

## 剖面图技巧（本次会话新增）
- **转换点对齐**：横轴 `t - tt_ocean`（t=0=波上岸）
- **理论到时虚线**：T-P=land_km/6.0, T-S=land_km/3.6（不用 tt_t 算的对齐后到时）
- **包络剖面模式**：去掉波形只留包络线 + fill_between

## 环境配置
- 在 tphase_env 中安装 `future`：`pip install future`
- PMCC 目录：`t_phase_analysis/internals/PMCC/`
- 已在 `t_phase_analysis/__init__.py` 中注入 sys.path
- 已在 `internals/PMCC/` 添加 `__init__.py`

## 修复记录
### azimuth_speed_array.py
- 移除 `plt.show()`（非交互环境下卡死）
- 修复保存路径：`../../` → 当前目录
- 增大画布：`plt.subplots(5, figsize=(14,16))` + `subplots_adjust(hspace=0.05)`
- bbox_inches='tight'

### run_fk.py（CLI 入口）
```bash
python /home/sy/codex/t_phase_analysis/t_phase_analysis/internals/PMCC/run_fk.py \
  --proc-dir <sacprocess> \
  --out-dir <output> \
  --stations S1,S2,S3 \
  --t0-times t1,t2,t3
```
- 自动重采样到 40Hz
- 直接输出 PNG（`n_t0_*.png`）
- 同时保存 NPZ 数据

## 结果解读
- **Consistency**：时延闭合程度，越小越好（<1 可靠）
- **Azimuth**：有效值 = 波到达方位角；-12345 = 无效
- **Speed**：有效值 = 视速度 (m/s)；-12345 = 无效
- 三站需紧凑排列（间距 1-20km），线状排列 FK 分辨力差
- T 波有效速度约 1.5-1.8 km/s（SOFAR 通道）

## requirements.txt
```
future
numpy
scipy
matplotlib
obspy
```
