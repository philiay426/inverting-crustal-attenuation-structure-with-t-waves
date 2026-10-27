---
name: twave-crustal-attenuation
category: seismic
description: >
  Use when 跑 T 波分析全流程并测地壳衰减（Q0/alpha）：预处理 -> 台站表与分布图 -> 波形剖面 ->
  偏振 -> F-K -> 尾波归一化反演。通用参数化模板 + 每步算法/判据/出图规范；换事件只改环境变量。
---

# T 波分析全流程 → 地壳衰减 Q0 / α

一条链、五个环节，终点是 Q(f) 与 Q = Q0·f^α。**方法固定、判据前置、每步一组标准图**。
本文只写"怎么做、怎么判、怎么画"；**不含任何事件的路径、台站名单、测值或结论** ——
那些属于具体事件目录。

**测的是什么**：陆段路径平均的有效 Q（含 T 波在陆壳内的转换与散射），不是体波 t* 的深部 Q；
只与同样用 T 波的研究对比。**不适用**：单台、陆段跨度 < 300 km、尾波发育不起来、台站近垂直于传播方向。

---

## 模块化调用（只做其中一步时看这一节）

本 skill 有两种用法：**全套衰减反演**（按 §1 顺序跑五步），或**只取其中一段**。
每一段都是一个独立可用的模块 skill，本 skill 只负责编排与步间契约，不重复模块内部的算法细节。

| 要做什么 | 加载 | 入口 | 输入 | 输出（交给下游的契约） |
|---|---|---|---|---|
| 只要数据能用来做包络 | `seismic-twave-analysis` | `<skill>/templates/preprocess_one_station.py` | 原始 SAC + SACPZ | `envelope/<STA>_envelope_2-8Hz.npz`（`time, waveform_z/r/t, envelope_z/r/t`）+ 六图 |
| 只要能算走时/陆海段/剖面 | `seismic-twave-analysis` | `<skill>/templates/step1_map_and_table.py` | SAC 头 + 等深线网格 | `stations.csv`（`sta,lat,lon,dist_km,land_km,ocean_km,tt_p,tt_s,tt_tp,tt_ts,snr_db,baz`） |
| 只要判断单台可不可用 | `single-station-analysis` | `<skill>/templates/polarization_rose_geo.py` | 三分量速度 + `tt_tp/tt_ts` | 可用台清单 + 逐台剔除理由 |
| 只要来向/慢度 | `fk-array-analysis` | `<skill>/templates/fk_geometry.py` → `step3_fk_driver.m` | `envelope` 或 2–8 Hz SAC + 台站几何 | `beam_<TAG>.txt`、`fk_peaks_<TAG>.csv`、慢度图 |
| 只要 Q(f) | 本 skill | `templates/twave_q_inversion.py` | 去响应**速度** SAC + `stations.csv`（含 `tt_tp/tt_ts`） | 逐台明细 + 逐频点 `f_hz,Q,sigma_Q,k,se_k,n_used` |
| 全套衰减结构反演 | 本 skill | §1 五步 | 事件目录 | Q0/α + 四张图 + 定案文档 |

**跨步铁律（模块间必须一致，否则结果对不上）**

1. **SNR 口径**：步 0 定的噪声窗/信号窗，步 1 的表、步 2 的判据、步 3 的筛台全部沿用；换口径要三处同时改。
2. **时间轴**：`t = b - o`，任何模块都不得自己另起一套时间基准；跨模块传时间一律传「相对发震」的秒数。
3. **单位**：波形速度 m/s（或注明 count）；`A_d/A_c` 无量纲；`land_km` 为投影陆段，不是震中距。
4. **台站名单**：以 `stations.csv` 的 `sta` 列为准（大小写不敏感），模块间传名单用逗号串环境变量，不写死在代码里。
5. **成品图不跨模块改**：上游出的图下游只引用不覆盖；同一图不同口径加后缀另存。

---


## 0. 通用约定（每一步都适用，先读这一节）

### 0.1 参数与代码约定

- **模板不含事件内容**：所有路径、台站名单、事件参数、频点、窗长一律走环境变量；
  必填项缺失时脚本 `raise SystemExit("缺 Q_CASE / FK_STAS ...")` **报错退出**，不许给"看起来合理"的默认值。
- 单位统一：距离 km、时间 s（相对发震，除非注明）、速度 km/s、慢度 s/deg、频率 Hz。
- **时间轴铁律**：`t = b - o`（SAC 头）。`stats.starttime` 可能被写成发震时刻，`o` 才是事件在文件内的位置；
  一次只看起点不看 `o`，整条时间轴平移几百秒，时窗取到噪声段。
- 目录：`NN_step/` 下**脚本 + 产物平铺**，不要多层 `code/`、`out/`；数据表可进 `数据/`。
  历史版本移 `_superseded/<日期>_清理/` 并写 `README_去向.md`（归档不硬删、删前先问）。
- **改图绝不原地覆盖**：先 `mv` 旧图到 /tmp 或归档目录，再写新文件；共享盘被 Windows 占用会 `Errno 22`。

### 0.2 制图规范（指针，规范本体不在这里）

**全部制图规定的唯一来源 = `seismic:plotting-standards` 的「§0 定版基线」**（栏宽、字号、线宽、图例位置、
轴标签格式、命名、rcParams 代码块都在那里，逐条是数字，不靠举例）。本 skill 只额外规定下面 4 条：

1. **Q 图**：双对数，纵轴一律取倒数 `Q^{-1}`；x/y 刻度显式 `set_xticks/set_yticklabels(["$10^{-4}$",...])`，禁 offset 文本；y 范围要覆盖一个整十倍程（如 `1e-4 ~ 1.2e-2`），数据带落图中部；**不画误差棒**。
2. **衰减散点图**：9 个频点 3×3 面板，事件用空心圆区分、黑色实线拟合；子图注记只写 `k`、`Q`、`n`。
3. **记录剖面**：逐道按自身峰值归一（保留相对振幅的图另出，不混）；到时线只用 `ARRIVAL_COLORS`；陆地计时与绝对计时**各出一版**，文件名带 `_landtime` / `_abs`。
4. **改图不覆盖**：新图另起文件名；确实要替换时先 `mv` 旧图到归档目录，写完 `ls -la --time-style=full-iso` 核对时间戳。

### 0.3 自检约定（脚本级，不看图）

模板末尾自带断言，跑完必须打印并通过：

- 输入完整性：台站数 = 名单长度、频率数 = 设定值、无 NaN/Inf、无空数组；
- 几何/数值：`sum(|v_i|)` 或窗边界与表逐项一致（`max|diff| = 0`）；
- 图面：轴范围覆盖全部数据、刻度数合理、文字越界由布局约束保证（`tight_layout` + 固定 `left/right`）。

**代码里不塞检测逻辑**（不写"自检模式"开关），项目原有自检保持不动；数值结论要能指到脚本或日志。

### 0.4 数值复现顺序（有新链时）

1. 先用**成品/参照数据**跑通公式，复现基准值；
2. 再用自研预处理跑同一事件，差多少记下来；
3. 差在步 0/1（窗、量纲、旋转、距离），**不要改公式**；
4. 报结果时写明是哪种数据、哪个口径。

---

## 1. 全链一览

| 步 | 干什么 | 关键判据 | 本步特有图 | 模板 |
|---|---|---|---|---|
| 0 | 去响应、带通、包络、逐站六图 + npz | 无缺口 / 无巨值 / 频谱无窄带 / 包络单峰 | `_1_raw`、`_2_remove_response`、`_3_spectrogram`、`_4_envelope`、`_bandpass`、`_packet`（全 1:1） | `preprocess_one_station.py`、`run_batch_preprocess.py` |
| 1 | 台站表（陆/海段、走时、SNR）+ 分布图 + 剖面 | `land_km` 必须投影；SNR 口径全表统一 | `station_map_all_names.png`、`fig3_record_sections_3x2{,_proj}.png` | `step1_map_and_table.py`、`step1_record_sections_3x2.py` |
| 2 | 单台分析（波形/SNR/偏振/质点运动四类证据） | 连续簇 + 目视 | `<STA>_rose_geo_TP_TS.png`、`<STA>_polarization_fullscan.png` | 外调 `single-station-analysis`；本 skill 两个偏振模板 |
| 3 | 子阵 F-K（慢度/BAZ） | 主次峰差 > 1 dB；BAZ 差 < 5° | 时窗选择图、分布图、处理对比图、慢度图 | `step3_fk_*`、`plot_fk_*` |
| 4 | 尾波归一化反演 Q(f) → Q0、α | `A_d/A_c ≥ 1`；`\|k\|/σ_k ≥ 1`；`\|α\|/σ_α ≥ 2` | Fig5 型 Q⁻¹ 图、9 频点衰减图、时窗选择图 | `twave_q_inversion.py`、`plot_q_inverse.py`、`plot_*`（论文复现类脚本按事件自建，不进模板） |

---

## 2. 步 0 · 预处理

> 细节与模板归 `seismic-twave-analysis`；本节只列编排需要的必做项、判据与图件，冲突时以模块 skill 为准。

**输入**：原始 SAC（+SACPZ）或 mseed + StationXML；事件参数（时刻/经纬度/深度）。

1. 读三分量，`t_off = b - o`（`0 < |o| < 1e6` 时），否则 `starttime - origin`。
2. 去响应 → **速度**（`output='VEL'`）；SACPZ 用 `read_inventory(f, format="SACPZ")`，自己解析时**只取第一个 ZEROS/POLES 块**。
3. **整段**带通（默认 2–8 Hz，corners=4、zerophase）再切窗；先切后滤会因窗边瞬变虚高噪声。
4. 包络 = `uniform_filter1d(|hilbert(x)|, size=int(smooth*fs))`，`smooth` 默认 6 s。

**本步图规范**（六张全 1:1，7.2×7.2 in @300 dpi）：

| 图 | 内容 | 要点 |
|---|---|---|
| `_1_raw` | 原始三分量 | 只标发震时刻，不标到时；看空道/缺口/毛刺 |
| `_2_remove_response` | 去响应后三分量 | 四条到时线（`ARRIVAL_COLORS`） |
| `_3_spectrogram` | Z/R/T 三面板 STFT，0–10 Hz | 共用 `Power (dB)` 色标、`viridis`、`vmin/vmax=P2/P95`、`nperseg≈13 s`、`noverlap 93.75%` |
| `_4_envelope` | 三分量包络 | 只画包络不画波形、不镜像 |
| `_bandpass` | 带通后按 BAZ 旋成 Z/R/T | 面板内右上角只写 Z / R / T，无标题 |
| `_packet` | SNR 计算图（信号窗/噪声窗色块） | 标题写 `SNR = ... dB (Welch 2-8 Hz)` |

横轴 `0 → tt_ts + tail`（`tail` 默认 500 s）；SNR 口径：`20log10(信号窗 RMS / 噪声窗 RMS)`，
信号窗 `[T-S±25]`（T）/ `[T-P±25]`（Z、R），噪声窗在 P 之前留 gap。**全表只能选一种噪声口径**并在表头写明。

**判据**：无缺口、无 1e6 级巨值、频谱无窄带干扰、T 波到时附近包络单峰。不合格写 `station_qc_notes.md`（哪张图什么现象）。

---

## 3. 步 1 · 台站表、分布图、记录剖面

> 细节与模板归 `seismic-twave-analysis`；本节只列步间契约与判据。

**算法（全部现算，不许抄表）**

```python
# 1) 走时：TauP 给远震 P/S（tt_p, tt_s）
# 2) 分段：沿震中->台站大圆，用区域地形/等深线网格找最后一个深于 REF_DEPTH 的点
#    REF_DEPTH 默认 -1000 m（SOFAR 声道轴深度，也是 T 波转换区）
# 3) 转换相到时：ocean/1.48 + land/V_P 与 ocean/1.48 + land/V_S
# 4) SNR：与步 0 同一口径
```

`land_km` 必须是**投影陆段**；用震中距会毁掉后面的斜率。`baz` 列不可信（常见漏 `cos(lat)`），
用 `gps2dist_azimuth` 现算。**负纬度**要检查半球字母是否把负号吃掉（纬度反号 → BAZ 反 180° → 振幅差几倍）。

**本步图规范**

- 分布图：岸线 + 候选台（三角按 SNR 着色，`viridis`）+ 震中星号 + 台站表栏；全站标名（8p 粗体 + 同字号白描边）；
  `REF_DEPTH` 等深线画黑虚线；`DPI=600`；底图是栅格时不出 PDF。
- 记录剖面（3×2）：上排纵轴震中距、下排纵轴陆段距离；横轴 0=发震 / 0=上岸；三列 Z/R/T；
  每道按自身峰值归一；到时用短竖虚线。`TT_MODE=proj`（投影走时）与 `const`（恒速公式）各出一版，
  两版横向错位量就是"该用哪套口径"的证据。

---

## 4. 步 2 · 单台分析

> 细节归 `single-station-analysis`（含三方法偏振对比、双波包拾取、图件自查）；本节只列"进反演前必须有哪几条证据"。

回答"这个台能不能用"，偏振只是四类证据之一：① 步 0 六图 ② 到时与 SNR ③ 偏振 ④ 三分量形态与质点运动。
最终产出 **可用台清单 + 每个落选台的理由**。

**本步图规范（定版两张）**

| 图 | 内容 | 固定参数 |
|---|---|---|
| `<STA>_rose_geo_TP_TS.png` | 两个地理方位玫瑰（T-P / T-S）+ Z/R/T 包络时窗 | PCA 滑窗（协方差），子窗 1 s / 步 0.25 s，窗 = 理论到时 ±25 s；扇区半宽 30°、权重 λ1；10° 分箱、环形 3 点平滑 |
| `<STA>_polarization_fullscan.png` | 竖排三栏：R 占比 / T 占比 / 线性度+平面度 | 步长 0.25 s，范围 `[T-P-100, T-S+500]`，不预设震相窗 |

版式要点：玫瑰画整圆、0°=N、顺时针、N/E/S/W；**理论线只画半段**（T-P 沿径向 `baz`、T-S 沿切向 `baz+90`），
落在样本多的一侧；柱宽 = 分箱宽、`edgecolor` 同色 `lw=0`（白描边会漏缝）；图内文字只留台名与距离；
`rt_angle_med_deg` 在 45° 附近会乱跳，判读看**扇区占比的对比**，别把单点中位角当来向；**PCA 有 180° 模糊**，只能报轴向角 `%180`。

---

## 5. 步 3 · 子阵 F-K

> 引擎、假象诊断、子阵选择判据归 `fk-array-analysis`；本节只列本流程用到的参数口径与判据。

**铁律：不要拿到事件就直接跑。** ① 选台 → ② 分布图 + 台站信息 + 处理对比 → ③ 时窗选择（参数在此定稿）→ ④ 结果图。
每组一个目录，产物不混。

### 5.1 几何量与预测值（先跑 `fk-array-analysis/templates/fk_geometry.py`，不许抄别处）

```python
C     = 111.195 / (2 * d_min)        # s/deg/Hz，台阵响应常数（d_min 为最近台间距, km）
u_pred = 111.195 / V_T               # s/deg；V_T 由该区频散/文献/波形先验定（km/s）
smax  = 1.4 * u_pred                 # 网格半径
spts  = 301                          # 奇数，0 在格点
f_hi  = min(1.0/25.0, C / smax)      # 频带上界：超过就会混叠塌到零慢度
f_lo  = 1.0 / win
floor = 1.28 * d_min / D             # 分辨率地板（D = 孔径），决定能报几位
```

`fk-array-analysis/templates/fk_geometry.py` 吃 `FK_STAS` + `FK_VT`（+ 台站表或 SAC 头），直接打印上面五个量并写 `fk_geometry.csv`；
Octave 驱动用同一套规则（`FK_VT` 或 `FK_SMAX`，`FK_C` 给频带上界），换地区只改 `FK_VT`。

台数下限按台间相干定（少于 4 台压不住旁瓣、会翻向）；是否允许两组共用过渡台由使用人决定，
报告写清本组为什么选它（台数、`d_min`、孔径、与上一组重叠几台）。

### 5.2 时窗（参数在这一步定稿）

```python
W0 = min(tt_tp) - x;  W1 = max(tt_ts) + x;  win = W1 - W0
tau_span = max(tt_ts) - min(tt_tp)
# 约束：tau_span / win <= 0.6，否则 T 波峰与零慢度假峰等强、主次分不出
```

时窗选择图：横轴从事件 0 画到 `W1 + 200 s`（全时段，不截近窗）；窗用阴影 + 两条竖虚线；
脚注写清波形链与 synchronize 起点。**判据**：窗内段与送检输入逐点比，去掉两端 taper 后 `max|diff| = 0`。

### 5.3 结果图与判据

- 慢度图：`flipud(fire)` + 固定 `CLim`（如 `[-0.7, 0] dB`）；半径按 `smax/3`、`2*smax/3`、`smax` 画灰色虚线圆（最下层）+ 每 30° 辐条；
  理论 BAZ 黑实线（最上层 zorder 10）；峰值靠色斑表达，另用绿色十字标理论慢度位；色斑只画 CLim 内。
- 标题只写组名 + 窗长；下一行 `measured <|u|> s/deg, BAZ <deg> deg`。
- 判据：只报峰不够，列 −6 dB 内**全部**局部极大并说明各对应哪个相位；主次峰差 > 1 dB 才算分得开；
  BAZ 与理论差 < 5° 视为一致、≥ 30° 判不可用；`|u|` 的小数位由 `floor` 决定，argmax 不报小数。

---

## 6. 步 4 · 尾波归一化反演 Q(f)

> 本步是**本 skill 的核心内容**（其它 skill 没有）；只想算 Q 而不做 F-K 时，从 §6.2 直接开始。

### 6.1 先定事件

多事件合并是**方法要求**，不是可选：单事件台数少、距离跨度窄，斜率被个别近台支配，Q0 与 α 都不可信。
把选中的事件写进 `events.csv`（时刻/经纬度/深度/震级/台数/数据目录），并说明为什么选它（同一源区、同一台阵、T-P/T-S 可辨）。

### 6.2 公式

```
Q(f) = pi*f / (V*|k|)          k = d/dr[ ln(A_d*sqrt(r*h_m)/A_c) ]
A_d = 信号窗 [T_ref, T_ref+WIN] rms     A_c = 尾波窗 [2*tt_s, +WIN] rms（lapse 从发震算）
r = land_km    h_m = 几何扩散特征长度（默认 70 km）    WIN 默认 20 s
ln Q = alpha*ln f + ln Q0（不加权；1 sigma = 95% CI / 1.96）
```

分量：**Qp 用垂直 Z（SAC 头 T3 为参考）、Qs 用切向 T（T4 为参考）**；尾波用同一分量。
`V = 6.0`(P) / `3.6`(S) km/s；理论到时用**步 1 算的投影到时**（`tt_tp/tt_ts`），不要在现场用 `ocean/1.48+land/V` 重算。

### 6.3 取窗的两个硬约束

1. **索引必须扣 `b`**：`i0 = int((t - b) / dt)`。不扣时整窗偏移，P 相位最敏感，Q0 能差 30%、α 差 0.3。
2. **多事件合并不扣事件常数**：直接对全部"台-路径"做一次 OLS（尾波归一化已扣掉源与场地）。
   按事件减均值会让斜率只由事件内距离跨度决定，近台变成高杠杆点，Q 可被一个台拉偏几倍。

### 6.4 判据（跑拟合前写定）

- 台站：`A_d/A_c ≥ 1`（不够 = 信噪不足），**不剔点、不做 3σ 离群剔除**；
- 频点：`|k|/σ_k ≥ 1`，且该频点台数 ≥ 8；
- 整体：频点数 ≥ 6 且 `|α|/σ_α ≥ 2` 才叫显著；否则只报"与频率无关、量级 X"。
- 被剔的频点/台站照实列出，理由写口径而不是"剔了好看"。

### 6.5 本步图规范（四张）

| 图 | 内容 | 要点 |
|---|---|---|
| `q_paper_fig5.png` | (a) Q⁻¹ 对频率 (b) Q⁻¹ 对波数 κ=f/V | **双对数、纵轴取倒数、无误差棒**；点 + 同色拟合线；图内文字写 `Qx = Q0 f^alpha`（纯文字）；y 轴范围拉到能看见一个整十倍程（如 1e-4 ~ 1e-2），数据带落在中部 |
| `q_decay_QP.png` | 9 频点：`ln[A_d·sqrt(r·h_m)/A_c]` 对陆段距离 | 3×3 面板，三事件三色（空心圆）、黑拟合线，子图注 `k / Q / n` |
| `q_decay_QS.png` | 同上，S 相位 | 同上 |
| `q_window_section.png` | 时窗选择：上排 Z（Qp）、下排 T（Qs） | 左列全时段、右列窗附近放大；信号窗橙条、尾波窗红条；37 台按陆距升序、逐道归一 |

`q_paper_fig5` 与 decay 图的数据全部来自同一份逐台明细表，改口径重跑即可，不改绘图代码。

---

## 7. 复现的落地顺序（换新项目照抄）

模板按模块分家：步 0/1 在 `seismic-twave-analysis/templates/`，步 2 在 `single-station-analysis/templates/`，
步 3 在 `fk-array-analysis/templates/`，步 4 在本 skill 的 `templates/`。

```bash
PY=/home/sy/miniforge3/envs/obspy/bin/python
SK0=~/.hermes/skills/seismic/seismic-twave-analysis/templates
SK2=~/.hermes/skills/seismic/single-station-analysis/templates
SK3=~/.hermes/skills/seismic/fk-array-analysis/templates
SK4=~/.hermes/skills/seismic/twave-crustal-attenuation/templates
EV=<事件目录>; OUT=<本步输出目录>

# 步 0 预处理
SRC=$EV OUT=$OUT STA=<台站>        $PY $SK0/preprocess_one_station.py
SRC=$EV OUT=$OUT STA_LIST=<A>,<B>  $PY $SK0/run_batch_preprocess.py

# 步 1 台站表 / 剖面
T_BASE=$EV T_OUTDIR=$OUT                          $PY $SK0/step1_map_and_table.py
T_BASE=$EV T_OUT=$OUT/fig3_record_sections_3x2.png TT_MODE=proj $PY $SK0/step1_record_sections_3x2.py

# 步 2 单台（只做偏振时在这一步停）
SRC=$EV/remove_response_vel OUT=$OUT STA=<台站>   $PY $SK2/polarization_rose_geo.py

# 步 3 F-K：先算几何，再跑引擎，再出图
FK_STAS=<A>,<B>,... FK_VT=<视速度 km/s> T_CSV=$OUT/stations.csv \
  EVENT_LAT=<..> EVENT_LON=<..> T_OUTDIR=$OUT     $PY $SK3/fk_geometry.py
FK_SEIZMO=<seizmo> FK_DATA=<子阵 SAC 目录> FK_TWIN="<W0>,<W1>" FK_VT=<视速度> \
  T_OUTDIR=$OUT conda run -n octave octave-cli $SK3/step3_fk_driver.m
T_BASE=$EV FK_STAS=<...> FK_TAG=<tag> T_OUTDIR=$OUT $PY $SK3/step3_fk_demo.py

# 步 4 衰减
Q_CASE=$EV Q_STAS="<A>,<B>,..." Q_CODA_MODES=2S Q_AMP=rms Q_HM=70 Q_PCOMP=z \
  $PY $SK4/twave_q_inversion.py
```

报告里必须写明本次用的 `V_T`（或 `FK_SMAX`）与 `FK_TWIN` —— 这两个是唯一需要人判断的量，其余都由算法给。
所有脚本的必填变量在各自 `templates/README.md`，缺一个就报错退出。

---

## 8. 全链陷阱（每条都踩过）

1. **`o` 字段**：用 `b-o` 定时间轴；别用 `starttime == origin` 判断起点。
2. **量纲**：`remove_response` 默认可能给位移；速度要么 `output='VEL'` 要么频域微分，并在图上注明。
3. **先滤波后切窗**（反了噪声虚高 20 dB 量级）。
4. **`land_km` 必须投影**，用震中距会毁掉 Q 的斜率。
5. **`baz` 现算**：漏 `cos(lat)` 会偏十几度；纬度负号被半球字母吃掉会反 180°。
6. **F-K 必须用包络**；`dB` 平均会把峰拽到另一相位；频带上界超 `C/smax` 或采样降到 2 sps，峰会塌到零慢度。
7. **窗太短**（`tau_span/win > 0.7`）时主次峰只差零点几 dB，判不可用而不是"也还行"。
8. **尾波 lapse 从发震算**（`2S`），不能用陆地计时取 2 倍。
9. **取窗扣 `b`**；不扣 = Q0/α 全偏。
10. **多事件不扣事件常数**；扣了近台成杠杆点。
11. **不剔点、不按结果挑台**（挑台 = 循环论证）。
12. **新链先复现再改**：成品数据对上基准值之前，不动公式。

---

## 9. 定案文档结构（每个事件照抄）

```
一、定案口径（表：事件/台组/尾波窗/信号窗/频点/振幅/速度/剔除）
二、定案结果（Q0±σ、α±σ、Q(3 Hz)、1σ 区间）+ 逐频点表（f, n, |k|, σ_k, |k|/σ_k, Q）
三、为什么能这样定（判据 + 独立性证据：换口径收敛 / 剔台救不回来 / 坏频点不是某一台带的）
四、台站剔除原则（逐台理由）
五、对比对象（其它事件、文献，给 Q(3 Hz) 差百分比）
六、产物清单        七、作废口径与理由        八、复算命令
```

---

## 10. 参考

- `templates/` —— **只剩步 4** 的 6 个脚本 + `README.md`（环境变量表）。步 0/1/2/3 的模板已归各自模块 skill。
- `references/q-inversion-math.md` —— 公式链与逐项出处（文献 / 源码行号 / 推导 / 自创）。
- `references/qc-criteria.md` —— 三层判据、剔除红线、敏感性检验。
- `references/figure-specs.md` —— 步 0–4 标准图版式、命名、返工点。
- 旁支 skill：`seismic-twave-analysis`（预处理）、`single-station-analysis`（偏振）、`fk-array-analysis`（F-K）、
  `record-section-plotting` / `plotting-standards`（制图）、`seismic-attenuation-inversion`（理论）。
