---
name: single-station-analysis
category: seismic
description: >
  单台偏振分析：T-P/T-S窗口滑动PCA，仅R-T平面；偏振引导拾取(理论±15s扫最径向/横向)与理论窗
  均默认 ±2s（2026-09-02 前十站敏感性定稿），1s子窗75%重叠。预处理见 seismic-twave-analysis skill。
---

# 单台分析

> **本流程定版模板（`templates/`）**：`polarization_rose_geo.py`（地理方位玫瑰 T-P/T-S + 包络时窗）、
> `polarization_fullscan_3panel.py`（R 占比 / T 占比 / 线性度+平面度 三栏扫描）；PCA 滑窗口径 = `cov`。
> 变量表见 `templates/README.md`，必填缺失即报错。

单台分析默认只做**台站分布图上标注的那 3 个台站**（最近陆 / 最远陆 / SNR 最高；用户 2026-09-28 定）；用户明确指定台站时以用户名单为准。预处理由 `seismic-twave-analysis` 的 `process_image.py` 完成，单台分析在此基础上做偏振。

**目录编号（2026-09-28 定稿）**：单台分析统一用 `process_flow/2_single_station_analysis/{STA}/`（**步2**）；不要再生成旧的 `3_single_station_analysis/`。
**0 号目录（数据处理）改名**：新规范 = `process_flow/0_data_processing/{remove_response,envelope,waveform_image}/`；本文下面出现的 `sacprocess/...` 是**旧事件的实际路径**，处理新事件时按 0 号目录写。

## 输出（v5，2026-09-15 定稿 3 张）⚠️ 默认方法 = 采样点等权（count，每点 1 票）

| # | 文件 | 脚本 | 内容 |
|---|------|------|------|
| 1 | `{STA}_PS_polar.png` | `polarization_analysis.py` | P、S 理论到时窗偏振方向 |
| 2 | `{STA}_pca_theory_TP_TS.png` | `pca_theory_polarization.py` | T-P、T-S 理论到时窗偏振方向（文件名沿用历史名；图内标题写实际方法） |
| 3 | `{STA}_polarization_fullscan.png` + `.csv` | `full_scan_polarization.py` | 全时段三面板：偏振方向 / 线性度 / 平面度 |

**方法开关**：三张图统一支持 `--method count|amp|cov`（默认 `count`；环境变量 `POLAR_METHOD` 同名生效，优先级 命令行 > 环境变量 > 默认）；输出目录可用 `--out-dir`（或环境变量 `OUT_DIR`）重定向，用于烟测不污染事件目录。
三方法全部保留：

| 方法 | 方向 | 权重 |
|---|---|---|
| `count`（**默认**） | 逐采样点 `atan2(E,N)` | 每点 1 票（无偏、不看幅度） |
| `amp` | 逐采样点 `atan2(E,N)` | `N²+E²`（瞬时能量） |
| `cov` | 子窗 3×3 协方差最大特征向量方位 | 代码实际用 `lin·√λ₁`（历史口径，`--method cov` 时数值不变） |

- **其余图默认不出**（需要时加 `--all` 或单独跑脚本）：`{STA}_Twave_3comp_*.png`、`{STA}_pca_TP_TS.png`（引导拾取版）、`{STA}_linearity_planarity_3D.png`、三方法对比图、窗口敏感性图等。
- ⚠️ **更正（2026-09-28）**：此处原先写着"既有旧图一律保留不删"，**用户指出这不是他的要求**（本文件其余用户口径均带"用户/原话"标注，唯此句没有 → 系 agent 自行归纳）。事实只有两条：
  ① **脚本只写不删**（重跑不会清掉已有文件）；
  ② **要不要清理目录由用户决定**，清时按定稿口径（每台 3 图 + fullscan csv + 脚本）归档到 `<目录>/_superseded/<日期>_清理/` 并写 README_去向.md，**别硬删**（删前先问）。
- 波型识别（P/SV vs Love、区域差异、R-T 判据）→ `references/twave-conversion-mechanism.md`。

数据源: 预处理输出 `sacprocess/remove_response/*_resp.sac`。

> 注意: 波形剖面图已并入步1 (`step1_station_table_and_map.py` → `profile.png`)，不再有独立步2。
> v3 起不再输出：`{STA}_envelope_polar.png`、`{STA}_theory_polar.png`、`{STA}_TP_TS_polar.png`、`{STA}_overview.png`。

**三合一总览图（用户常要的排版）**：`{STA}_3panel.png` = 原始波形 / 2–8 Hz+8 s 包络 / 频谱 0–10 Hz 三块**横排一行、每块 1:1**（`figsize=(3W, W)`，如 (18, 6)；块内 Z/R/T 三行共用横轴 0–2000 s）。到时竖虚线**只标用户要的口径**（常见要求："只留恒速计算的到时"，且表头不再重复写到时数值）；图例放面板内左上、三列各一个（轴外 `bbox_to_anchor` 会与列标题重叠）；三块共用横轴时列标题里写台站名，避免单独截图时看不出是哪台。

## 单台分析流程（v4，2026-08-04 一键驱动）

> ⚠️ 2026-09-15 起驱动默认只跑三张图（见上方“输出 v5”）；旧七张图需 `--all`。

**推荐方式**：用通用驱动脚本，不需要手改模板参数：
```bash
python ~/.hermes/skills/seismic/single-station-analysis/scripts/run_single_station.py \
    --event <EVENT_KEY> --stations STA1,STA2,STA3 [--snr-min 3.0]
```
- 自动从 events.yaml 读事件参数、从 raw 文件名解析发震时刻、从 raw 推断台站网络
- 自动 copy 预处理图 + 跑 5 个模板脚本 + 重命名/清理
- 输出到 `process_flow/2_single_station_analysis/{STA}/`

## 偏振反演观测 BAZ（只用于验证 F-K 来向）

`obs_baz = (baz + tp_rt_angle) % 360`（`pca()` 的 ZRT `principal_vector` 里 `a = atan2(T,R)`，径向偏振轴沿传播方向）。
只取震源侧分支，不画 +180 分支；轴向角差用 `%180` 不用 `%360`。
它只作**事后诊断**（与 F-K 的 BAZ 比残差、看是否同号同量级），不用它反过来挑时窗或挑台。
完全不以理论 BAZ 参与挑窗的独立流程见 `references/independent-baz-energy-window-workflow.md`。

**SNR 口径（2026-10-14 起，取代 2026-09-02 的“未带通 + Welch 频段积分”）**：统一噪声基准 SNR = 20log10( 信号窗 RMS / **阵列噪声中位数** )；波形**先 2–8 Hz 整段带通、再切窗算 RMS**。`stations.csv` 主列 `snr_u` 即此口径（`snr/snr_z/snr_r` 是自噪对照列）。`sacprocess/envelope/*.npz` 的 `waveform_*` 已是 2–8 Hz 缓存，**不得再滤一遍**（双重带通扭曲 SNR），只可画图或 `band=()` 直接取窗。完整参数、两个坑（整段滤波、噪声窗避开 P 波）、下游接法 → seismic-twave-analysis `references/snr-unified-noise-baseline.md`。

**偏振原理说明的表达**：用户只问原理时，连续、直接地回答；**默认给 3–6 行要点版**（用户嫌长会直说"太复杂了"），只在其要求时才展开成多段；除非明确要求，不贴代码，也不要拆成很多细小段落或清单。必须明确区分“无向偏振轴”和“有向 BAZ”。

## ⚠️ 时间轴铁律（2026-08-04，必须遵守）

**SAC `o` 字段优先用 `b-o` 计算时间偏移**，数据常从发震前开始记录：
- 2013 南非事件：SAC `o=600, b=0` → `t_off = b - o = -600`，数据从发震前 600s 开始 ✓
- ObsPy `starttime` 可能被 SAC reference time 误导：reference time 被写成了发震时刻
  （nzyear..nzsec = 16:22:31 = 发震），但数据实际从 16:12:31 开始。
  此时 `starttime - EV_ORIGIN = 0` 是**错误**的，必须用 `b - o = -600`。

**统一做法**：
1. `o_val = getattr(sac, 'o', None)`
2. `if o_val is not None and 0 < abs(o_val) < 1e6: t_off = sac.b - o_val`
3. `else: t_off = trs['BHZ'].stats.starttime - EV_ORIGIN`（无 o 字段时 fallback）
4. `t = np.arange(npts)*delta + t_off`，t=0 = 发震时刻
5. slice 窗口：`starttime + (w_abs - t_off)`，w_abs 是相对发震的绝对坐标

**判断技巧**：`t[0]` 应为负值（发震前开始记录），如 -600；若 t[0]=0 但 SAC 有 o 字段，说明用错了坐标系。
完整调试复盘（reference time vs o vs ObsPy starttime，误判教训）→ `references/sac-time-axis-b-o.md`。

公共模块：`scripts/single_station_common.py`（load_event / parse_origin_from_filename / load_3comp）

按顺序执行，输出到 `process_flow/2_single_station_analysis/{STA}/`：

1. **copy 预处理三分量图（最先做）**：从 `sacprocess/waveform_image/{STA}/` 复制
   - `{STA}_3_spectrogram.png` → `{STA}_spectrogram.png`（频谱图）
   - `{STA}_bandpass_2-8Hz.png` → `{STA}_bandpass_2-8Hz.png`（带通 2-8Hz）
   - ⚠️ 预处理**不生成** `{STA}_three_filter.png`（run_single_station.py 第1步复制的是不存在文件），直接复制上面两个
2. **T波三分量窗口**（不再画 P/S 窗口）：`plot_waveform_windows.py` → `{STA}_Twave_3comp_{t0}-{t1}s.png`
   - 时间范围：T-P-100s ~ T-S+200s
   - 三分量 Z/R/T 三行并排，波形灰色 (`color='0.55'`)，包络红色平滑 **8s**
   - **归一化（2026-08-24 定稿）**：Z/R/T 共用一个**全局波形 P98**（`np.percentile(np.abs(concatenate(z_bp,r_bp,t_bp[m])), 98)`），不分量独立归一化——独立 P98 会抹掉三分量相对振幅，无法比较分量强弱。纵轴标签 `Norm. Amp.`
   - 删除旧单分量 `{STA}_Twave_window_*.png` 和 `{STA}_PS_window_*.png`
3. P/S 偏振：`polarization_analysis.py` → `{STA}_PS_polar.png`（重命名自 `{STA}_p_s_polar_P*_S*.png`）
4. **偏振引导的实际 T-P/T-S 拾取与 PCA**：理论到时仅作 ±15 s 宽先验；对总体 SNR≥10 dB 的完整三分量台站扫描 2–8 Hz 滑窗（中心步长0.25 s、默认窗口±2 s，2026-09-02 前十站敏感性定稿）。弱短时窗数据**不强制连续簇**：在满足单波包 SNR、线性度、分类角度门槛的候选中，联合选择径向 T-P 与横向 T-S，且必须 `T-S ≥ T-P + 3 s`。CSV 同时记录实际窗口、偏振/SNR/质量标志，供走时拟合、转换点和人工复核；偏振最优短窗不是严格初至，须在结果中标注。详见 `references/polarization-guided-picking.md` 与 `references/short-window-polarization-picking.md`。

**搜索方案窗口回填单台主图**（2026-08-23 VOI 实测）：`run_polarization_guided_picks.py` 输出的 CSV 窗口是**相对发震的绝对时间**；`pca_theory_polarization.py` 模板内部用**陆地走时坐标**（`w_abs = w + TT_OCEAN`）。回填时按使用方区分坐标：
- 给 PCA slice 用：`w = PICK_ABS - TT_OCEAN`（保持模板 run_pca 内部 +TT_OCEAN 逻辑不变）；
- 给绝对时间轴图面板（ax3 的 axvspan/axvline）用：**直接用 `PICK_ABS`**，不能把陆地走时值画上去（VOI 原模板把 ~35/63s 的陆地走时画到 900–1000s 信号轴，色块完全错位）；
- 图标题/注释标注窗口来源（guided-pick）与拾取偏移（如 VOI T-S pick −14s vs 理论），便于与理论窗版本对比；旧图先备份再覆盖。
- 搜索方案拾取偏移大（VOI T-S −14s）说明理论窗（含后移12s）严重偏离实际波包，主图必须用搜索窗重画，不能只依赖理论窗。
详见 `references/guided-pick-window-application.md`。

**T-P 前峰包偏振检查**（2026-08-23 VOI 实测）：用户指出 T-P 拾取窗口前有峰包时，在 **T-P 窗口起点之前**的段（VOI: 940~956s）取 8s 平滑包络 argmax 定位峰值 tpk，对 tpk 分别做 ±1s 与 ±3s 的 R-T PCA。判读：
- ±1s 峰值窗角度干净（VOI: 91.0°，lin 0.995，p25/p75=80.5/95.0）→ 峰核单相偏振；±3s 角度跨大半圆（−103~80.5°）说明两侧混入不同偏振能量，峰核才是纯相；
- R-T≈90°=横向（Love/SH 型），与 T-P 径向（≈0°）相反 → 独立横向到达而非 T-P 前驱；距 T-S 拾取仍远（VOI: 24.9s）则与 T-S 无关；
- **实现位置（用户要求 2026-08-23）**：❌ 不新写脚本、不新增图排版（`plot_preTP_peak_polarization.py` 独立 2×2 图曾被用户明确否决并删除）。✅ 峰包偏振并入搜索方案脚本 `run_polarization_guided_picks.py`，逐站输出 CSV 新列 `pre_tp_peak_time / pre_tp_peak_dt / pre_tp_peak_angle / pre_tp_peak_lin / pre_tp_peak_snr`（8s 平滑包络在 `[tp_pick−20, tp_pick−HALF−0.25]` 内 argmax，±1s 窗 PCA，复用 packet_snr）。"看一下"= 报 CSV 行/打印数值；临时数值查询直接终端 heredoc，不落盘脚本。
调试坑：去响应后速度单位振幅 ~1e-8，print 用 `:10.1f` 显示 0.0 误判无峰——打印包络先归一化或用科学计数法。详见 `references/pre-tp-packet-polarization.md`。
4. **PCA 观测到时偏振（偏振引导版，2026-08-26 定稿）**：`pca_polarization.py` → `{STA}_pca_TP_TS.png`
   - **拾取逻辑（用户定稿，替代旧包络峰+30%阈值法）**：在**理论到时 ±15s** 内扫描 1s 子窗（75% 重叠），T-P 选最接近**径向 R**（0°/±180°）的时窗中心、T-S 选最接近**横向 T**（±90°）的时窗中心；最终窗口 = 中心 ±2s（与理论版口径一致，2026-09-02 定稿）。角度距离按折叠到 [0,90] 计算，线性度作为同距离下的次要判据。
   - ⚠️ **不要用 ±100s 宽搜索**：会选到远处偶然径向/横向的时段（VINA 实测 T-P 偏到理论+99s）；±15s 搜索范围是定稿值。
   - **run_pca 坐标坑（2026-08-26 VINA 复现）**：`run_pca` 内部已执行 `w_abs = w_l + TT_OCEAN`，**必须传陆地走时坐标**（`w_tp0_l`），不能传绝对坐标（`w_tp0 = w_tp0_l + TT_OCEAN` 仅供 axvspan/打印）——否则双重偏移 +936s，最终 R-T 角与拾取角矛盾（如拾取 98° 横向但输出 31.5°）。自检：最终 R-T 角中值应 ≈ 拾取角，T-P 径向、T-S 横向。
   - 内核：t_phase_analysis 包 `pca()`（1s 子窗 75% 重叠，analytic 协方差，ZNE→ZRT 自动旋转）
   - 布局 2×2：T-P R-T 玫瑰 | T-S R-T 玫瑰 | 三分量 Z/R/T 波形 + 拾取窗 | 线性度
     - **下方波形格三分量（2026-09-07 SA18 迭代定稿）**：Z/R/T 三条**浅灰 #999999、lw≈0.5、不分颜色**，各按自身窗口 max 归一，纵向错开**行距 2.0**（Z 基线 0 / R +2 / T +4，ylim 留 -1.2~5.6）；左侧 **Z/R/T 黑色字体**标注放**各行波形顶部上方的空隙**（基线 +1.25，如 y=1.25/3.25/5.25），不与波形重叠；Z 包络灰点线叠 Z 轨作参考；拾取窗色块与理论虚线保留（早期彩色版 Z 黑/R 蓝/T 红、行距 1.7、标签压轨中已被用户否决）
   - **图标注（用户偏好 2026-08-23 定稿）**：表头只标理论→实际时间（`T-P th 964.4s → pick 956.7s`），不标子窗/重叠/analytic；图例放波形图**右上角**（axvspan/axvline 加 label + `legend(loc='upper right', framealpha=0.9)`），**只标拾取中心时间**（如 `T-P pick 956.7s`、`T-S pick 977.6s`），不标窗口范围（`[x, y]s` 去掉）；理论虚线也进图例（`T-P th ...`）。详见 `references/pre-tp-packet-polarization.md`（含沿理论扫描最径向/最横向中心→固定 `T_TP_CENTER/T_TS_CENTER` 的分析流程与 VOI 结果）。
5. **PCA 理论到时偏振**：`pca_theory_polarization.py` → `{STA}_pca_theory_TP_TS.png`
   - 窗口：T-P/T-S 均按**计算到时 ±2s**（2026-09-02 前十站敏感性测试定稿：±2s 为最短稳定窗；T-S **不延后**，与 T-P 同口径；2026-08-26 起已删除此前模板的 `tts+9 ~ tts+15` 后移 12s 逻辑）。窗口重叠检查仍保留（`w_tp[1] > w_ts[0]` 时后延 T-S）。
   - 图总标题仅保留台站名、SNR、理论 T-P/T-S；子窗、重叠率、平滑与窗口策略写入代码注释。
   - 左下诊断面板显示 2–8 Hz Z/R/T 三分量：波形统一淡灰，包络 Z黑/R蓝/T红并浅填充；三分量共用展示时段的 global P98，保留真实相对振幅。R/T 玫瑰图参考线必须同时保留：R=0°绿色虚线、T=90°红色虚线；字母放在图内靠近对应参考线、加粗，避免与方位刻度重叠。径向刻度默认保留；若两图比较需统一刻度，不得无说明隐藏。
   - 短窗口理论偏振诊断的展示包络默认 **4 s** 平滑；不得将步1剖面图的6s或长时窗三分量图的8s机械套用到 ±2s 诊断窗口。
   - 用半透明色块标出实际 PCA 窗口：T-P 蓝色、T-S 紫色，并以同色虚线标理论中心时刻。
   - **时间坐标自检**：窗口在陆地走时坐标定义、图轴为相对发震时间时，绘制色块前必须加回 `tt_ocean`；不得直接把陆地走时数值画到绝对时间轴。
   - 窗口敏感性试验（T-P/T-S ±3s、后移 10/11/12s 等）必须输出独立命名的对比图，不能覆盖主结果；比较窗口时保持三分量波形、包络与彩色窗口同图显示。详见 `references/theory-polarization-figure-display.md`。
   - **理论窗口敏感性对比**：严格理论窗 (`T-P/T-S ±3 s`) 与用户指定平移窗须输出为不同文件，不能互相覆盖。对比台站必须从当前已验证的新版绘图模板生成，**不得复制或运行旧台站脚本**；左下角用相对发震时间显示 Z/R/T：global P98、灰色波形、分量彩色4s包络、蓝/紫实际 PCA 窗口长方形；总标题仅保留台站名、SNR、理论 T-P/T-S。逐站偏移、保留规则与滑动扫描判读见 `references/theory-window-comparison.md`。
   - **窗口扫描判读**：扫描 T-S 横向段或 T-P 径向段时，必须检查玫瑰图集中度与连续相邻窗，不能只依据中位 R-T 角或高线性度选窗；高线性度只代表线性运动，不代表目标波型。偏振最优但走时残差过大的晚到波包，不得改标为 T-P/T-S。详见 `references/theory-window-comparison-workflow.md`。
   - 窗口扫描判读、严格理论窗与后移窗的可复现对比、共同刻度玫瑰图规范 → `references/theory-window-polarization-comparison.md`。
   - **窗口敏感性**：若严格理论窗不清晰，保留严格理论窗对照后分别扫描 T-P/T-S 候选窗；不得用任意大后移量“挑好看结果”。候选必须满足走时顺序、方向集中度、线性度及物理可解释性。详见 `references/theory-window-sensitivity.md`。
   - **二维网格敏感性扫描（多参数联合验证，如窗宽 x 海洋声速）**：用户要求"做不同尺度对比验证可靠性"且涉及两个以上参数时，不要逐个一维扫描，改用网格扫描+热力图（loss-landscape 风格，x/y=两参数，色=残差，星标当前方案）。见 `references/window-velocity-grid-scan.md`。
   - **全时段偏振扫描（整段 T 波一视同仁，验证理论窗与边界拾取）**：脚本 `templates/full_scan_polarization.py` = **多台站批量版**（`python full_scan_polarization.py [STA ...]`，默认 10 高SNR站；multiprocessing 8 进程；每站输出到自己的 `{STA}/`；CSV 缓存命中则改图不重扫）。扫描 [理论T-P−100s, 理论T-S+500s]、步长 0.25s、每窗 ±3s PCA（同搜索方案口径），输出 `{STA}_polarization_fullscan.csv` + 三行图（**Z/R/T 三分量 2–8Hz 包络（共用 global P98、`Norm. Amp.`）/ dist-to-radial 折线 / 线性度折线**），竖线**只画理论到时**（蓝/红虚线 + 数值标注）——**不要画未标注的辅助线**，画了用户会问"这条实线是什么"。用途：① 理论窗处是否真有径向/横向偏振簇；② 贴 ±15s 搜索边界的拾取是否真实（BATG 实测：理论 T-P 无径向支持、真实径向簇晚 ~14.6s → 边界拾取是真实信号）；③ 展示整段波列径向/横向交替多簇结构（10 站全都在理论±60s 内找到两类簇，但偏移 −31.5s~+45.5s **各站不一致，不能统一时移**）。判读必须用**连续簇**（≥1s、res<10°、lin>0.99）而非全局最小残差窗口（孤立偶然窗会污染结论）。方法与判读 → `references/window-velocity-sensitivity-scan.md`；1D/2D 参数敏感性 → `references/window-velocity-grid-scan.md`。
6. 线性度/平面度：`polarization_attributes.py` → `{STA}_linearity_planarity_3D.png`
   - **只做 T 波时段**（T-P-50s ~ T-S+150s），滑动 1s 子窗 75% 重叠
   - 面板: (a) 线性度 vs 时间 (b) 平面度 vs 时间 (c) 3D 散点(R-T角/时间/线性度)
   - 分段统计: T-P 段(tp-20~tp+40) / T-S 段(ts-10~ts+60) 的 lin/pla/R-T角 中值

模板: `templates/pca_polarization.py` + `templates/pca_theory_polarization.py`（复制后改 EV/BASE/STA/NET 即可用）。

**独立玫瑰图（±25s / 振幅法 / 原始计数 / 地理 N-E 方向，论文 Fig2e/f 风格）**：窗中心用论文
走时（Equation 1: `ocean/1.48 + land/6.0`、`+ land/3.6`），**勿用 stations.csv TauP 到时**
（偏早 ~20-25s 导致窗罩不住真波包）；±25s 宽窗判别用径向/横向扇区占比而非 PCA 主轴中值；
每面板只画本相位理论线（T-P=径向 BAZ、T-S=横向 BAZ+90）、黑色实线**过圆心完整直径**、
刻度只留 N/E/S/W、标题极简；±25s 等窗口变体输出 `_pm25` 后缀新文件不覆盖原 ±2s 图。
完整配方/算法/判读 → `references/amplitude-rose-standalone.md`。

### 三方法偏振对比图（协方差 / 振幅加权 / 采样点等权）

当用户要看"不同偏振估计方法是否给出同一方向"时，用**3 行（方法）× 2 列（T-P 径向 / T-S 横向）**
的玫瑰对比图，每格标题写该方法的扇区占比（R 或 Tr）：

| 方法 | 每个样本的方向 | 权重 | 特性 |
|------|---------------|------|------|
| Covariance (PCA) | ±25 s 窗内 1 s 子窗 / 0.25 s 步，N/E 协方差最大特征向量方位 | 最大特征值 λ₁ | 最稳健、抗噪；假定子窗内单一偏振 |
| Amplitude | 逐采样点 `atan2(E, N)` | `N²+E²`（瞬时能量） | 时间分辨率最高、无窗参数；强尖峰会主导直方图 |
| Sample count | 逐采样点 `atan2(E, N)` | 每点 1 票 | 无偏、不看幅度；背景噪声被同等计入 → 系统性低估偏振强度 |

- 扇区占比定义：以 BAZ 为中心 ±22.5° 为径向 R，以 BAZ±90° 为中心 ±22.5° 为横向 Tr（角度用地理 N/E，0°=N 顺时针）；
  窗口 ±25 s，到时用论文式(1) 恒速。
- 判读：三者同向 → 结论稳健；等权被稀释而另两法有峰 → 偏振能量集中在波包内（波包特征）；
  PCA 与振幅法互不一致 → 窗内多段不同偏振叠加，先回看波形再下结论。
- 批量做法：单脚本循环台站列表（G1 = SA16/17/18/22/23/25/29），每台输出到自己的
  `process_flow/2_single_station_analysis/{STA}/`，并在终端逐台打印恒速到时便于核对；
  参见本目录 `scripts/` 的批量脚本模式（显式列表循环，不用 glob 花括号）。
- 🔴 **三条 r 轴刻度不可互比**（用户会指着图上的数值问"这些 12345 是什么意思"）：每根柱高 = 该方法投给该方向的**权重之和**，量纲各不相同——
  协方差 = λ₁ 之和（单位 (m/s)²，实测单票 1e-18~3e-17、扇区峰 ~1e-16）；振幅 = (N²+E²) 之和（同单位同量级）；等权 = 票数（**总和恒等于窗内采样点数**，±25 s @20 Hz = 1000，扇区峰 ~40）。
  所以只能比**方向分布形状**，绝不能用数值大小比方法强弱；可比的数字只有标题里的 R / Tr 百分比。
- 🔴 **判读前先对随机基线：径向 R = 12.5%、横向 Tr = 25%**（±22.5° 扇区 = 圆周 1/8，横向两个扇区 = 1/4）。接近基线的值不构成偏振证据（实测某台恒速窗 T-P 径向 12% ≈ 随机，换口径后才有 27%）。跨方法比较必须**同一方法跨口径比**——子窗平均让 PCA 的占比天然高于逐点法，拿 PCA 的 80% 去比等权的 19% 是错的。
- 🔴 **窗中心口径敏感，报告必须写明用的哪套**：偏振窗中心在论文恒速式与射线理论之间差 20–25 s，而窗半宽只有 ±25 s → 几乎错开一半，R/Tr 显著改变。G1 七台实测（协方差法）：T-P 径向 **7/7 台在射线理论口径下一致变强**（49→80%、45→73%、39→60%…），T-S 横向无一致趋势（3 升 4 降）。批量脚本把口径做成开关（`--tt-mode const|ray`，ray 版输出带 `_ray` 后缀不覆盖恒速版），并逐台打印所用到时便于核对。

- **汇总对比图（三方法 × 两口径）画法**：2 行（T-P 径向 R / T-S 横向 Tr）× 3 列（方法），每格 x=台站列表、两条折线（恒速=实线圆点、射线理论=虚线方点）、同台两口径之间连灰细线，并画随机基线水平灰点线（R 12.5% / Tr 25%）；数字同时存 CSV（如 `G1_polarization_compare.csv`）便于逐项核对；图例放**整图顶部**（塞进数据区会压到高值台站）。

- 🔴 **批量脚本里与台站有关的量必须逐台从 `stations.csv` 读**（BAZ、海洋/陆地距、到时），不得沿用首台或 SA18 的硬编码值——用户会专门问答式验证（“其它台站画的是对应的 BAZ 吧，不是锁死的吧”）。自证：用合成数据调用绘图函数后读回参考线角度，`np.degrees(ax.lines[0].get_xdata()) % 360` 应等于各台自己的 BAZ / BAZ+180（径向面板）、BAZ±90（横向面板）。

**方法原理说明的表达**：用户要"简要"时给 **3–5 行要点版**（每种方法一行：取什么方向 + 什么权重 + 一句特性），
然后只补一句"看三者是否同向"的判读。用户会明确说"太复杂"来要求压缩 —— 默认就写短，不要先长篇再等纠正。

**交叉验证**：偏振结果可用 t_phase_analysis 包的 `pca()` 复核（API 踩坑、R-T 角提取、窗口敏感性问题）→ `references/t-phase-analysis-pca.md`。

**PCA 诊断时序图**：`t_phase_analysis/pca.py` 自带 `plot_pca_diagnostics(result, z_trace=None, stacked=False)`，输出 5 行共享 X 轴诊断（Z 波形/特征值 λ₁λ₂λ₃/BAZ 逐窗散点/线性度-平面度-球形度/主轴 R-T 角）。对某连续窗（如 T-P_land−60 ~ T-S_land+80）`pca()` 后直接调用即得 `{STA}_pca_diagnostics.png`；注意图例 λ 下标字符（U+2081..）在 CJK 字体（Noto Sans CJK）缺字形会显示方框，可用 DejaVu Sans 渲染。入口：`run_pca_diagnostics.py` 模式 = 与 pca_polarization 相同加载 + `res=pca(st_win,...)` + `fig,axs=plot_pca_diagnostics(res, z_trace=trs['BHZ'])`。

**阈值选择与文献边界**：30% 是偏振窗口的经验性相对阈值，不是普适初至标准；初至推荐噪声统计阈值，且应做 20%/30%/50% 敏感性分析 → `references/envelope-threshold-selection.md`。

## 双波包独立拾取与 SNR 门槛（2026-08-18 更新）

T-P/T-S 偏振分析必须视为**两个独立的小波包**，不能在一个宽的 T 波段内用单次 `argmax`。单次全窗最大值会系统性偏向最强波包，导致另一个波包被漏拾或重复拾取。

### 拾取规则

1. 以理论到时分别建立 T-P、T-S 搜索窗；两个窗口独立找局部峰值。
2. 候选峰使用 `find_peaks` 的 `prominence`、峰高和与理论到时的距离联合筛选，不能只按绝对峰高排序。
3. T-P/T-S 不允许共用一个峰；两个波包的偏振窗口不得重叠。
4. 30% 相对峰值阈值只用于定义单个波包的窗口边界，不是波包存在性或严格初至判据。
5. 若某一理论窗口没有合格局部峰，标记该波包为不可用，不得强行把最大波包当作另一个波包。
6. 建议对 20%、30%、50% 边界阈值做敏感性分析；报告峰值时间、窗口边界和窗口间距。

### 单波包 SNR 与可用性

总 T 波 SNR 不能替代 T-P/T-S 波包质量控制。分别计算：

```text
SNR_TP = 10 log10(P_TP / P_noise)
SNR_TS = 10 log10(P_TS / P_noise)
```

正式双波包偏振分析要求 **SNR_TP ≥ 10 dB 且 SNR_TS ≥ 10 dB**；6–10 dB 仅作探索性结果并需偏振稳定性验证；<6 dB 不用于双波包结论。除 SNR 外，建议要求每个波包中位线性度 ≥0.6、偏振角圆形标准差 ≤20°，且两个窗口不重叠。具体质量控制表和伪代码见 `references/dual-packet-picking-and-snr.md`。

## ⚠️ 归一化与纵轴单位（2026-08-24 用户多次追问后定稿）

**去响应数据 vs 绘图显示是两回事**，用户会逐张图追问单位，必须口径统一：

1. **去响应数据** = 地面**速度（m/s）**：`tr.simulate(paz_remove=paz, remove_sensitivity=True)` 出速度；SACPZ 的 INSTGAIN/SENSITIVITY 标 `(M/S)`（STS-2 等速度计）；远震体波/T 波量级 ~1e-8–1e-7 m/s。
2. **绘图显示** = **P98 归一化**（无量纲）：Z/R/T 共用一个**全局波形 P98**（波形与包络共用同一因子）。归一化后纵轴量级是 O(1)（±2~±5），**这是归一化的正常表现，不是单位错误**。
3. **纵轴标签**：归一化图一律 `Norm. Amp.`；**不能写 `Velocity (m/s)` 或 `Counts`**（数据已归一化，物理单位已不存在）。
4. **纵轴范围**：动态 `±1.05×max(|norm|)`，**不要固定 ±2.5**——固定范围会把真实峰值截断成"假 ±2"视觉，用户会误判为归一化没生效。
   - ⚠️ **纯 max 定轴会被体波短时尖峰顶爆**（2026-09-02 AMPY 实测）：2–8Hz 带通后 P 波振幅仍达 11.7×ep98 → 旧 `1.05×max(波形,包络)` 得 y_max=±12.3，T 波(≈1×ep98)被压到轴高 ~5% 不可见（用户报"纵轴上下10太大了、T 波不明显"）。`process_image.py::plot_bandpass` 已改为**稳健分位数定轴**：`y_max = 1.05×max(P99(窗口内 (|波形|+包络)/ep98), 1.0)` → AMPY ±2.04，T 波占轴高 ~48%；P/S 顶部被截断是预期（到时线仍标注，T 波才是关注对象）。`plot_raw`/`plot_decon` 仍是 max 式——用户若要求 T 波在这些图也可读需同样处理。数值诊断法与新旧对比见 `references/robust-yaxis-scaling.md`
5. **新/旧流程不一致**（跨事件对比时注意）：旧模板（如 2014 Puerto Rico 事件 H65A 的 `BPF.png`）**不归一化**直接画去响应速度 → 纵轴显示真实量级 `1e-7`（Amplitude 标签隐含 m/s）；新模板归一化 → `Norm. Amp.` ±几。两套都能自洽，但单位标注风格不统一；用户统一方向确认后再批量改，**未确认前只改当前事件**。

其他事件模板（Southern Africa 2019/2020、Japan 2026 的 process_image.py）的 decon/bandpass 图**归一化了但标签仍写 `Velocity (m/s)`**——历史标签错误，统一时修正为 Norm. Amp.。
完整对比与代码模式 → `references/waveform-normalization-units.md`

## ⚠️ 图件自查（2026-08-24 用户要求"以后这种东西你要自己检查一边"）

每次出图后**必须自查**，不能等用户发现：
- 标题/表头文字是否完整、被截断或超出图边界（suptitle 过长 + tight_layout 会裁掉）
- 图例位置是否符合要求（如"左下角/右上角"）、是否与其他标注重叠、是否遮挡数据；图例背景用浅白半透明（`framealpha≈0.9`）
- 图内标注（震源/地震信息、单位标签、台站名）是否超出画面边缘
- 纵轴标签与实际数据是否一致（归一化图不能标 Velocity (m/s)）
- 快速数值自检（窗口/角度/线性度打印值）与图面一致
- **pyGMT 图（转换点/台站地图）标注布局**：`fig.text` 与 `fig.legend(position='jBL...')` 都默认落左下角 → 文字会被塞进图例框（`top10_obs_vs_theory_conversion.png` 案例）。文字放图例上方（`y=R[2]+5.0`）或震源机制球旁（`x=EVLO, y=EVLA+0.9, justify='CB'`）；白色描边无背景框 `pen='0.4p,white'`（pen=文字轮廓，不要 fill）。详见 `references/pygmt-annotation-layout.md`

## ⚠️ 数据质量预检（偏振分析前必做）

跑偏振分析前必须先检查三分量数据质量，坏数据会污染偏振结果。

**新事件数据准备**（从 EarthScope 下载 XK/XV 等台网数据、补 SAC 坐标、提取 PAZ 的完整流程与踩坑）→ `references/earthscope-download-workflow.md`

```python
from obspy import read
import numpy as np

for comp in ['BHZ', 'BHN', 'BHE']:
    st = read(f'{sta}..{comp}.*_resp.sac')
    tr = st[0]
    neg_pct = (tr.data < 0).sum() / len(tr.data) * 100
    print(f'{comp}: min={tr.data.min():.2e} max={tr.data.max():.2e} '
          f'mean={tr.data.mean():.2e} std={tr.data.std():.2e} neg={neg_pct:.1f}%')
```

### 异常判据

| 指标 | 正常范围 | 异常信号 |
|:---|:---|:---|
| 负值占比 | ≈ 50% | **>90% 或 <10%** → DC 偏移（如 VOET Z 分量 100% 负值） |
| Z/H 振幅比 | 0.3–3 | **>5** → Z 分量增益异常或传感器故障 |
| 分量间 std 比 | 0.3–3 | **>5** → 某分量异常放大/缩小 |

### 处置

- 仅 Z 分量异常 → 偏振分析只用水平分量（N/E），R-T 偏振不需要 Z
- 水平分量异常 → 排除该台站
- 全部异常 → 检查 SAC_PZ 是否正确，重新下载数据

## 偏振分析参数

| 参数 | 值 |
|------|-----|
| P/S 窗口 | 理论到时 ±2s |
| T-P/T-S 窗口 | 观测版（run_polarization_guided_picks.py 批拾取）：偏振引导 ±15s 搜索→中心±2s（2026-09-02 定稿）；理论版：计算到时 ±2s 不延后 |
| 包络阈值 | 30%（仅作显示/包络拾取参考，2026-08-26 起观测到时偏振改用偏振引导，不再用包络峰定窗） |
| 包络拾取 | 峰值法（旧版 v2，已被偏振引导替代） |
| 子窗口 | 1.0s，75% 重叠（Flinn 1965） |
| 方法 | 滑动协方差 PCA (`np.linalg.eigh`) |
| 平面 | **仅 R-T 平面**（v3，移除 Z-R/Z-T；T波偏振不画 Z-R 玫瑰） |
| 窗口基准 | **陆地走时**（v4：TTP_LAND = tt_tp - ocean_km/1.48，t=0=波上岸），不用绝对到时 |
| T-P 搜索窗 | 观测版：理论 ±15s（2026-08-26 定稿，替代旧 -100~+100s 宽窗）；AMPY 等理论偏差大的事件需人工核查 |
| T-S 搜索窗 | 观测版：max(TTS_LAND-15, t_obs_tp+5) ~ TTS_LAND+15（2026-08-26） |
| 配色 | P=绿, S=橙, T-P=蓝, T-S=紫；窗口方块 T-P=蓝 `#1f77b4`, T-S=紫 `#9467bd` |
| 判据 | R-T 0°=径向(P/SV), ±90°=横向(SH/Love) |

🔴 **偏振分析必须做在 2–8 Hz 上（2026-10-20 用户明确："偏振不是对 2-8Hz 的分析吗？"）**：
`pca()` 之前先对**整段**波形做 2–8 Hz 零相位带通，再切窗（只滤 4 s 窗会因窗边瞬变污染主轴角）。
不要拿去响应后的宽频位移直接算：4 s 窗内 2–8 Hz 只占位移 RMS 的 **0.05~0.5%**，宽频口径的 R-T 角
与 0.998~0.999 的线性度反映的是长周期背景（微震/去卷积分量），不是 T 波本身。实测同窗两口径差异极大
（AMPY T-P 拾取窗 +0.5° vs +101.9°；AMPY T-S 理论窗 −115° vs +13°），线性度从 ~0.999 掉到 0.54~0.87。
实现：`run_polarization_guided_picks.py` 参数区 `POLAR_BAND = (2.0, 8.0)`，在 `load()` 里整段带通，
置 None 即退回旧宽频口径（旧拾取表归档在 `polarization_guided_picks/_superseded/<日期>_rawdisp/`）。
⚠️ 2–8 Hz 下 R-T 角中值抖动大，单看中值选窗不稳（判据宜改用扇区占比）；
⚠️ `sector_pct()` 的随机基线：轴向折叠后径向/横向都是 **25%**（±22.5° 双边）；本文件他处写的
"R 12.5% / Tr 25%" 是不折叠、径向只算单瓣的口径，两条不能混着引用。

### 入射角（要 Z 分量参与）：R-T 平面偏振给不出它

🔴 **别混口径**：本 skill 的 R-T 平面 PCA 给的是**水平面方位角** `atan2(v_T, v_R)`（径向/横向判据），
**不是入射角**。用户问“前面测得的入射角是多少”时，那个角其实是方位角——拿它去答垂直入射角会被当场纠正。

要入射角就加 Z：**3D PCA**（Z,R,T 协方差主轴）取 `i = arccos( |p_Z| / |p| )`（p = 最大特征值对应主轴；
单位向量 `|p|=1` ⇒ 写成 `arccos|p_Z|` 等价，分母是防未归一化的冗余项）；自垂直方向量。
**第二独立估计器 = 水平/垂直振幅比**：`i_amp = tan⁻¹( √(RMS_R² + RMS_T²) / RMS_Z )`（= `atan(v_H/v_Z)`）。
🔴 ❌ **不是 `atan(v_Z/v_H)`** —— 那是"从水平量"的余角，会把入射角报成它的补角（同一 SA18 ±5 s 窗：
正解 57.46°、反写 32.5°）。两者差 ≤3° 才算稳。
🔴 **根号必须在分子内、逐分量平方后相加**：读成 `RMS(√(R²+T²))/RMS_Z`（先开根再取标准差）在 SA18
同一窗给 **38.24° vs 正解 57.46°**（差 19°）——`√(RMS_R²+RMS_T²)` 的平方 = `mean(R²+T²)` 是水平矢量模的
**均方根**，而 `std(√(R²+T²))` 扣掉了水平方向的"直流"分量，必然偏小。
📐 两个三角函数不是两种方法：`tan⁻¹` 就是 `atan`；`atan(H/Z) ≡ acos(Z/√(Z²+H²))` 是**恒等式**
（同一向量两种写法，逐点数值可证）；`atan2(H,Z)` 与 `atan(H/Z)` 在 H,Z≥0 时等价、只是不除零。
⚠️ 式子里写 `RMS`、代码里写 `np.std` **只在窗内均值为 0 时等价**（带通窗近似成立、实测逐位相同；
未去趋势或有偏置的窗不等价）——写进论文时统一成一个口径。
⚠️ `|p_Z|` 取绝对值 ⇒ 只给"离垂直多远"，**不分上行/下行**；要判初动方向不能这么写。
🔴 **椭圆偏振时两估计器分道扬镳，PCA 会给出错值**：相位差 90° 的椭圆运动使协方差非对角项被抹掉，
特征向量退化成 Z 轴与 R 轴本身、主轴取方差大的那个（合成 i=55.4°、R 方差更大 ⇒ PCA 报 **90.0°**，
而振幅比仍给 55.4°）。⇒ **直线度低（<0.7）时只信 `i_amp`**；报 i 前用"已知 i 的直线偏振合成波包"
跑一遍同窗同带自检（真值应原样复原）。
⚠️ **自由表面效应**使地表质点运动的视入射角大于真实射线入射角（偏陡）⇒ 由 i 反推的
v 偏小（4.5 km/s 是下界）；报 i 必须带上是哪个窗、哪条频带（≥5 Hz 直线度塌陷就弃）。

**用户贴一张公式图问"这个对吗 / 这个是 arctan 吗 / 这样吗"时的答法**（连着追问三遍的那种）：
① 先给代数判定（是不是定义式/恒等式、定义域与量级对不对）→ ② 再用**合成真值 + 真实数据各算一遍**
给数字（如真值 55.4°→复原 55.40°、SA18 恢复 README 的 57.5°）→ ③ 最后列**容易读错的地方**
（根号位置、atan/atan2、RMS/std、多余括号与绝对值号）。只答"对"或只讲道理都不合格——用户要的是
能直接粘进论文的干净式子 + 数值背书。
公式链、恒等式证明、椭圆反例、复算代码与**全子阵（多台）扩展配方** → `references/incidence-angle-from-polarization.md`。

🔴 **会不会循环验证**（用户必问）：`i` 是**水平旋转不变量**（换错 baz 也逐位不变，可数值自证），
所以“用 baz 定 R/T”不污染 i；真正要交代的是：①窗的选取——用偏振引导窗测 i 再验证偏振就是循环，
解药 = 滑动窗盲检（不引用拾取，只按直线度门筛）；②闭环 `v = sin i / u` 里的 `u` 来自 R/T **分量包络**
⇒ 依赖 baz（彻底做法 = 对原始 N/E 做矢量波束形成，或只用 Z 做 F-K）。
口径、盲检数值与闭环门禁表 → `fk-array-analysis` 的 `references/fk-parameter-selection.md` §6d/§6e。

🔴 **慢度符号统一 = u**：慢度无论写 s/km 还是 s/°，一律记 **u**（矢量，指向传播方向）；`s/°` 里的 s 是
"秒"这个**单位**、不是符号；Seizmo 代码的 `smax` / `svol` / beam 列名 `x(s/deg)` 是它的参数命名，
**不搬进公式**；本 skill 的小写 `p` 专指 3D PCA 主轴向量（`p_Z, p_R, p_T`），`p_ray` 这个符号已废除。
**改符号 / 改口径必须一次改全**：案例 README、脚本 docstring 与打印/CSV 标签、本 skill 与 references、
F-K 推导文档（`fk4d_source.md` 及其 `docs/` 副本各一份）；改完**重跑脚本核对数值逐位不变**，再
`ls -la --time-style=full-iso` 确认文件真写进去了（/mnt/e 写盘会静默失败：脚本打印 ✓ 而 mtime 不动）。
🔴 **论文的"预测慢度"是它给的主量，括号里的 km/s 只是它的换算——别把方向搞反**：原文是
`predicted value of 13.7 s/° (8.1 km/s)`，即 **13.7 s/° 是预测慢度**、8.1 km/s = `111.195/13.7` 是它的**等效视速度**
（式(1) 预设的是 1.48 / 6.0 / 3.6 km/s，根本没有 8.1）。⇒ ❌ 不要写成 `u_pred = 111.195/v, v 取 8.1`（把论文的
输出当输入 = 循环，不能算复现）；❌ 也不要给 13.7 乘 sin i。要跟它比就**在同一口径上比**：
由偏振 i 算的 `u = 111.195·sin i/v` 与论文的 13.7 是**同一类量**（都是慢度），差别只在等效视速度
（G1 七站 15.3 s/° ↔ v_app ≈ 7.3 km/s vs 论文 8.1 km/s）。要独立复现论文那个数只能回到式(1) 走时场（见 fk-array-analysis）。
▶ 复算脚本：`scripts/check_incidence_formula.py`（合成真值 + 真实窗一次跑完，给数字背书用）。

存档（用户要求“公式和结果存到文件夹”）：`process_flow/2_single_station_analysis/{STA}/incidence_angle/`
= README（公式链 + 结果表 + 可靠性结论 + 循环性分析）+ 一个可重跑脚本 + 逐窗/分频带/SNR/闭环 CSV
+ `i(t)` 与直线度双面板图（**图内文字用英文**，Arial 无中文字形会出方框；中文只留 README/日志）。

## 观测到时偏振拾取方法（v5，2026-08-26 用户定稿）

**偏振引导拾取**（替代旧包络峰+30%阈值法，`pick_by_polarization()` 见 `templates/pca_polarization.py`）：
- **在陆地走时坐标扫描**（`t_land = t - TT_OCEAN`，t=0=波上岸）
- 搜索范围：T-P 理论 ±15s；T-S `max(TTS_LAND-15, t_obs_tp+5)` ~ `TTS_LAND+15`（保证 T-S 在 T-P 之后、不重叠）
- 对每个候选中心 ±2s（含窗内所有 1s 子窗）计算 R-T 角中值作为该时窗偏振：
  - T-P 目标 = **径向 R**（0°/±180°，角度距离折叠到 [0,90]）
  - T-S 目标 = **横向 T**（±90°）
  - 同角度距离时取线性度更高者
- 最终窗口 = 选中中心 ±2s（与理论版同口径），窗口重叠检查仍保留
- 拾取窗口转回绝对坐标仅用于画图：`w_abs = w_land + TT_OCEAN`（axvspan/axvline）

⚠️ **搜索范围教训（VINA 实测 2026-08-26）**：±100s 宽窗会选到远处偶然径向/横向的时段（T-P 偏到理论+99s、T-S 拾取角 -1° 反而径向）；±15s 是定稿值，不要放宽。

⚠️ **run_pca 双重偏移 bug（VINA 实测 2026-08-26）**：`run_pca` 内部已加 `TT_OCEAN`，调用必须传**陆地走时坐标**（`w_tp0_l`），若传绝对坐标（`w_tp0`）会再偏移 +936s → 最终 R-T 角与拾取角矛盾。自检：打印的最终 R-T 角中值应 ≈ 拾取角。

⚠️ 理论走时不可靠实例（2026-08）：AMPY 理论 T-P=879s，但全 trace 包络最大峰在 ~780s（T 波实际到达）；若偏振引导在 ±15s 内找不到合格径向/横向窗，需人工核查该台站。

**窗口半宽（±N s）敏感性（2026-09-02 AMPY 实测）**：最终窗口半宽对结果稳定性有影响——
±1s 过窄时 R-T 角不稳定（AMPY T-P 角跳到 114°，偏离径向），±2~5s 稳定（~-21.5°）、
线性度随窗口增大略升（0.993→0.997）。**2026-09-02 前十站（VOI/SOLA/FOMA/AMPY/BKTA/BATG/KIRI/BAND/CPSM/VINA）敏感性对比后定稿：默认半宽 ±2s（进入稳定平台的最短窗）**，批量脚本 `window_size_top10_comparison.py`（输出根目录 `window_size_top10_comparison.csv/.png`）。若用户要求做
±1/2/3/4/5s 对比，复用 AMPY_window_size_comparison.py 模式（逐半宽跑 pca、汇总 R-T
角/线性度/平面度 + 折线图 + 玫瑰图对比；玫瑰图子图须用 polar=True 的 add_subplot，
不能放普通笛卡尔子图）。
⚠️ **对比指标必须是"偏离径向/横向角度"（用户 2026-09-02 明确要求），不是 R-T 角绝对值**：
```
dev_tp = min(abs(ang_tp % 180), abs(180 - (ang_tp % 180)))  # T-P 偏离径向(0°/±180°)
dev_ts = abs((ang_ts % 180) - 90)                            # T-S 偏离横向(±90°)
```
越小越接近理论偏振方向，图才直观（一眼看出哪个窗口偏振最接近径向/横向）；画 R-T 角绝对值
没有明确物理含义。CSV 同时输出 dev_tp/dev_ts 列。AMPY 实测：±1s T-P 偏离径向 65.8°（不稳定），
±2-5s 稳定 ~21.5°；T-S 偏离横向 33-45°（横向性一般，说明接收端波型非纯 SH）。
- **"±2s 最合适"定稿（2026-09-02 用户确认并已同步代码）**：±2s 是**进入稳定平台的最短窗**（±1 不稳、≥±2 的 T-P 相同 ~21.5°；T-S 偏横向 ±2s=33.7° 略优于 ±3s=38.4°）。已改：`run_polarization_guided_picks.py` HALF=3→2、skill 模板 `templates/pca_polarization.py` 与 `templates/pca_theory_polarization.py` 及各站相关副本统一 ±2s（理论版为纯理论中心 ±2s，T-S 不后移；观测引导版为扫描中心 ±2s），引导拾取 CSV 已按 ±2s 重跑（accepted 9 站不变），VINA 观测版单台图已重跑。2026-09-02 已补齐：10 站 `pca_polarization.py` 副本全部换成模板引导版±2 并重跑 `{STA}_pca_TP_TS.png`；`plot_polarization_map.py` / `plot_polarization_roses.py` 改 ±2 并重跑（`polarization_map_tp.png`、`polarization_roses_top10.png`）。
- **论文/汇报表述的文献边界**：Shearer (2019) 教科书第 13.5.1 节协方差偏振法只给操作流程（去均值→开窗→协方差→特征向量）与"窗内=定向信号+非相干噪声"模型假设，**未给具体窗长秒数**——引用时只能落在"窗须短到只含目标波列、长到样本足够"的原理上，不得声称书中有 ±2s/±3s 字样。美东（PR2014，USArray TA）实测点数教训在 seismic-twave-analysis `references/merged-polarization.md`（W57A：±2s 只 5-6 个测量点方位偏 7°，±5s 13-14 点偏 1°）。可用的示例措辞与数值见 `references/polarization-window-justification.md`。

**台站分布图 + T-P 偏振方向（旋转前，2026-09-02 用户要求）**：用户要求\"在台站分布图上用箭头标注偏振方向，改为旋转前的样子\"。做法：
- 每站 T-P 理论到时 ±2s 窗口内（2026-09-02 起），用 **N/E 分量（BHN/BHE）做 2D PCA（不旋转到 R/T）**，主轴方位角 = `atan2(E主分量, N主分量) mod 180`（0–180°，线性无方向性）——这就是\"旋转前\"的绝对方位角，箭头方向直接对应地图上的真实方向；线性度 = `(λ1-λ2)/(λ1+λ2)`（`np.linalg.eigh` 升序，最大特征值在 `vecs[:,1]`）
- pygmt 底图（地形/海岸线），每站画**双头箭头**（偏振主轴无方向性，双头才准确；单向箭头会误导）——pygmt 画箭头用 `fig.plot(data=[[x0,y0,x1,y1]], style='v0.45c+e+a30', pen=..., fill=...)`（4 列 data 形式；`x,y + direction` 形式报 list index out of range）
- 颜色/长度 ∝ 线性度（红>0.7 橙 0.4-0.7 蓝<0.4），长度用固定度值（MAX_LEN≈0.4°，避免密集区重叠、出界）
- 马达加斯加实测：南部台站（CPSM/AMPY/FOMA/MAHA）偏振方向指向震中 → 符合远震 P 波径向偏振；中北部方向各异（接收端波型混合，是 F-K\"哪个分量准\"的线索）
- 脚本：`plot_polarization_map.py`（输出 `polarization_map_tp.png`）

旧版（包络峰值法 + 30% 阈值 ±5s、±100s 宽窗）见 `references/envelope-polarization-analysis.md`。

## SNR 计算

Welch PSD 频谱能量比、2-8 Hz 积分 —— **⚠️ 已于 2026-10-14 废弃**，现行口径见上文「SNR 口径」与
seismic-twave-analysis `references/snr-unified-noise-baseline.md`（2–8 Hz 整段带通 → 窗内 RMS → 除以阵列噪声中位数）。
下列旧参数仅作历史对照：
- ~~数据源：去响应未带通 resp.sac~~；~~信号窗 `[tt_tp − 20, tt_ts + 100]`~~；~~噪声窗 P 前等长 `[tt_p − L, tt_p]`~~
- 仍有效的一条：时间轴基准用 npz time / `b-o`（SAC starttime 不可信）；R/T 先按与 BHZ 的 starttime 差样本对齐再 rotate

## BH1/BH2 分量映射

```python
if chn == 'BH1': chn = 'BHN'
elif chn == 'BH2': chn = 'BHE'
elif chn == 'BHX': chn = 'BHN'
elif chn == 'BHY': chn = 'BHE'
```

常见非标准分量：IU 台站 BH1/BH2，PS 台站 BHX/BHY。

## 运行

模板脚本参数写死在顶部（EVLA/EVLO/EVDP、STA/NET、BASE），**不支持 CLI 传参**。标准流程：

1. 复制模板到事件目录：`cp .../templates/polarization_analysis.py <evt>/process_flow/2_single_station_analysis/`
2. 用 python 字符串替换参数（不要手写整个脚本）：
   - `EVLA, EVLO, EVDP = 19.0428, -66.8096, 20.0` → 事件坐标（查 events.yaml）
   - `EV_ORIGIN = UTCDateTime(2019, 1, 22, 19, 1, 43)` → 事件发震时刻（文件名 `D.2019.022.190143` → 19:01:43；SAC 缺 `o` 字段时**必改**，否则时间轴偏移）
   - `STA, NET = 'K62A', 'TA'` → 目标台站/网络
   - `BASE = "/mnt/f/sharedate/case/2014-01-13-mw63-process"` → 实际事件目录（注意 `-station_low` 后缀）
3. 运行：`~/miniforge3/envs/obspy/bin/python <脚本>`；⚠️ 若报 `ModuleNotFoundError: No module named 't_phase_analysis'`，运行前加 `PYTHONPATH=/home/sy/codex:/home/sy/codex/t_phase_analysis`（polarization_analysis / polarization_attributes / plot_waveform_windows 依赖 codex 的 t_phase_analysis 包，2026-09-07 SA18 实测）
4. 重命名输出为事件目录约定（其他台站一致）并清理中间图：
   - `{STA}_p_s_polar_P*_S*.png` → `{STA}_PS_polar.png`
   - `{STA}_tp_ts_polar_*.png` → `{STA}_TP_TS_polar.png`
   - 删 `*_p_s_fullwindow*.png`、`*_p_s_onset10s_*.png`

**v3 旧七张图（legacy，需 `--all`）**（预处理图拷贝也已归入 `--all`）：
1. 复制预处理图（先做）：`sacprocess/waveform_image/{STA}/{STA}_3_spectrogram.png` → `{STA}_spectrogram.png`；`{STA}_bandpass_2-8Hz.png` → 同名
2. `templates/plot_waveform_windows.py`（T波三分量窗口，不画P/S）
3. `templates/polarization_analysis.py` → PS_polar
4. `templates/pca_polarization.py` → pca_TP_TS（观测到时偏振引导版：±15s 搜索最径向/最横向→中心±2s，2026-09-02 定稿）
5. `templates/pca_theory_polarization.py` → pca_theory_TP_TS（理论窗口，T-P/T-S 均 ±2s，T-S 不延后）
6. `templates/polarization_attributes.py` → 线性度/平面度

## 已知问题

- 🔴 **`--out-dir` 静默失效、图写进事件目录（2026-09-15，Python 3.14 实测）**：`full_scan_polarization.py` 用 `mp.Pool` 并行，子进程入口必须**显式接收** `out_root/flat/method`（已改为 `_worker(job)`）；
  不能靠“main() 里设全局变量、子进程读全局”的写法——Python 3.14 在 Linux 上默认 start method 已从 fork 改为 **forkserver**，子进程重新 import 模块，拿不到 main() 赋的全局量，输出目录会静默退回模块级 `ROOT`（= 事件目录），`图跑到 case 目录而 /tmp 里只有脚本`。
  同理 `run_single_station.py` 给**每个**模板都必须传 `--out-dir <out_dir>`（漏传则模板用脚本内 OUT → 写进 case 目录）。烟测自检：跑完 `diff` 一下事件目录的文件 mtime 应保持不变。
- **⚠️ 驱动默认只跑三张图**（2026-09-15 用户定）：`run_single_station.py --event <key> --stations <STA>` 默认只出 `PS_polar` / `pca_theory_TP_TS` / `polarization_fullscan`（+CSV）；旧图（含预处理图拷贝）需 `--all`。
- **⚠️ stations.csv 的 SNR 列（2026-10-14 现况）**：`step1_map_and_table.py` 输出的 csv 是 **19 列**：`…,snr,snr_z,snr_r,snr_t,snr_u,snr_uz,snr_ur,lon,lat,baz`；主列是 **`snr_u`**（统一噪声基准 T 分量，2026-10-14 起），`snr/snr_z/snr_r` 是自噪对照列。
  分分量重算用 `compute_rms_pair(…, t_off=)`（SAC `o` 字段口径，数据起点相对发震偏移，本项目 −600 s）。
  单台分析/转换点选站前先确认 csv 是新版（19 列），否则会选到旧口径的台站。
  详见 seismic-twave-analysis `references/snr-unified-noise-baseline.md`；旧 Welch 口径的 `references/per-component-snr.md` 已存档。

- **⚠️ 批量多台站跑模板的替换清单**（2026-08-04 实测，漏一项就错位/写错目录）：
  每个模板脚本顶部有 4 处必须全部替换，**不能只改 STA,NET**：
  1. `EVLA, EVLO, EVDP` → 事件坐标（查 events.yaml）；漏改 → 到时全错（实测 WTBG P=0s/T-P=6488s，因为模板默认是 2014 加勒比事件坐标）
  2. `EV_ORIGIN = UTCDateTime(...)` → 事件发震时刻；漏改 → 时间轴偏移（`polarization_attributes.py` 模板默认 2013 事件发震，直接复用时 T 波窗口切空 → `pca() requires exactly 3 traces`）
  3. `STA, NET` → 目标台站/网络
  4. `OUT = .../{STA}/` 硬编码台站名（pca 模板默认 `.../SOE`、polarization_attributes 模板曾默认字面 `.../STA`）→ 漏改结果写进错误目录
  注意：pca/pca_theory/polarization_attributes 的 `OUT` 行是在 `STA` 赋值**之前**出现的 f-string，直接 `re.sub` 替换 OUT 行本身即可，不要依赖 STA 变量。
  **批量脚本用显式列表循环**（`for sta in ['A','B','C']`），不要用 `glob.glob('{A,B,C}_*.py')`——Python glob 不支持花括号展开，静默返回空列表、替换全部没执行（本次实测踩坑）。
  完整清单 + 正则替换模板 + 自检 → `references/batch-template-replacement-pitfalls.md`
- **✅ 第 7 张模板已补齐**（2026-08-04）：`templates/polarization_attributes.py` 已创建（只做 T 波时段 T-P-50s~T-S+150s 的线性度/平面度 + 3D 散点）。完整 7 张图可输出。注意：`polarization_attributes.py` 也遵循"OUT 硬编码台站名"坑——批量跑时同样要替换 OUT 行。
- **⚠️ pca_polarization.py 的 run_pca 传参坐标坑**（2026-08-23 VOI 实测、2026-08-26 VINA 复现）：模板 `run_pca(w0, w1)` 内部执行 `w0_abs = w0 + TT_OCEAN`（**期望陆地走时坐标**），但模板曾把拾取窗口先转成绝对坐标（`w_tp0 = w_tp0_l + TT_OCEAN`，供 axvspan/打印）再直接传入 → 实际 PCA 窗口 = 绝对 + 929s，错位到 ~1900s，玫瑰图全是伪结果（VOI 曾输出 T-P 66.8° 不径向、T-S 48.9° 不横向；VINA 偏振引导版拾取角 98° 横向但输出 31.5° 矛盾）。修复：`run_pca(w_tp0_l, w_tp1_l)` 传**陆地走时坐标**；axvspan 仍用绝对坐标。自检：最终 R-T 角中值 ≈ 拾取角，T-P 应径向（≈0°）、T-S 应横向（≈90°）。用户手动固定窗口时（"直接调参数"）设 `T_TP_CENTER / T_TS_CENTER / HALF_WIN` 三个变量（绝对时间中心，±5s 默认），同样按"run_pca 传陆地坐标、axvspan 用绝对坐标"规则。
- **⚠️ pca 模板 OUT 硬编码台站名**（2026-08-04）：`pca_polarization.py` / `pca_theory_polarization.py` 里 `OUT = f'{BASE}/process_flow/2_single_station_analysis/SOE'` 硬编码 **SOE**。批量跑多台站只改 `STA, NET` 不够——必须同时替换 OUT 行，否则结果写进 SOE 目录并报 FileNotFoundError：
  ```python
  s = re.sub(r"OUT = f'\{BASE\}/process_flow/2_single_station_analysis/SOE'",
             f"OUT = f'{{BASE}}/process_flow/2_single_station_analysis/{STA}'", s)
  ```
  验证输出：日志行 `=> .../{STA}/{STA}_pca_TP_TS.png` 应含目标台站名。
- **⚠️ polarization_analysis.py 输出命名坑**（2026-08-04）：它同时输出 `{STA}_p_s_polar_P*_S*.png`（P/S ±2s，**正确来源**）和 `{STA}_p_s_onset10s_*.png`（P/S onset+10s）。重命名时用前者 → `{STA}_PS_polar.png`，删后者 + `*_p_s_fullwindow*.png` + `*_tp_ts_polar*.png`（⚠️ **仅 v3→v5 命名迁移时适用**，属一次性动作；日常目录清理由用户决定，别顺手删）。不要误把 onset10s 当 PS_polar。
- **⚠️ 剖面图全局 P98 被异常强台站绑架**（2026-08-04 南非 2019-01-22 实测）：`waveform_proflie.py`（seismic-twave-analysis 模板）全局 P98 归一化时，个别振幅异常台站（CRAIL 包络峰为其他台 79×、VOET 25×；去响应全局 std 大 10-50×，属响应/增益异常而非真实强 T 波）会把其余 20+ 台 T 波压成"平直噪声"——看起来像无信号，实为显示假象。诊断：① 各台算 T 波窗口 (T-P-50~T-S+100) 包络峰找 outlier（>中位数 10×）；② 对比去响应全局 std 确认异常；③ 加入 `PROFILE_EXCLUDE` 重画。信号强弱以单台峰/噪比为准，别只看剖面图视觉平直度。本事件实测该排除：SOE/POGA（用户指定）+ CRAIL/VOET（振幅异常）
- **⚠️ waveform_proflie.py 时间轴 bug**（2026-08-12 马岛 2013-07-22 实测）：旧版模板在 SAC 缺 `o` 字段时走 `else: t = np.arange(npts)*delta`，时间轴从**数据开始时刻**起算（常为发震前 600s），而 T-P/T-S 理论到时相对**发震时刻** → 整体偏移 → 剖面图 T 波"看不见"（全是平直噪声），但单台三分量图 T 波很明显。修复：模板已加 `_parse_origin_from_filename()`，缺 o 时用 `t = (starttime - origin) + np.arange(npts)*delta`。自检：T 波应在 T-P/T-S 虚线附近有纺锤形能量团。另注意强台站归一化绑架（上一项）。完整诊断/修复/频段选择 → `references/waveform-profile-time-axis-bug.md`
- SAC 头可能缺 `stla/stlo`，从 StationXML inventory 获取
- **`sac.o` 存在且合理 → 用 `b - o`**（2026-08 南非 2013-08-17 最终定稿，与顶部"时间轴铁律"一致）：数据常从发震前开始记录，`o=600, b=0` → `t_off = -600`，`starttime - EV_ORIGIN = 0` 反而是**错误**的（ObsPy starttime 被 SAC reference time 误导，reference 写成了发震时刻）。判断：SAC 有 o 字段时先试 `b - o`，看 `t[0]` 是否为负（发震前开始）；若 `t[0]=0` 但数据实际含发震前 600s 噪声 → 坐标系错了。完整复盘 → `references/sac-time-axis-b-o.md`
- **新事件先查 o 字段**：`python3 ~/.hermes/skills/seismic/single-station-analysis/scripts/check_sac_o_field.py <raw目录>` 批量检测该事件 SAC 是否缺 `o` 字段并给 t_off 建议。实测缺 o 事件：2019-01-22、2023-05-21；有 o：2013/2018/2020 各事件（与 D/M 型无必然关系，取决于下载数据中心/批次）
- **BASE 目录名后缀坑**（2026-08）：事件目录名常带 `-station_low` 等后缀（如 `2019-01-22-mww67-prince-edward-islands-station_low`），但 events.yaml 与 run_single_station.py 里的 BASE 常缺该后缀 → 替换 BASE 后先 `ls -d` 验证目录存在
- **同事件多 case 分支副本坑**（2026-08-05 实测）：同一事件可能在 `case/Southern_Africa/` 和 `case/Madagascar/` 下各有副本，**数据/process_flow 只存在于其中一个**（2013-08-17 事件只有 `Madagascar/` 分支有 stations.csv 与完整产物）。事件目录内脚本（如 Madagascar 专用台站图 `singer.py`）的 BASE 常硬编码指向无数据分支 → FileNotFoundError。处置：先 `ls -d` 确认实际数据目录，再改脚本 BASE。singer.py 修复记录（`snrs` 未定义 NameError、region 需含震源、图例左下角 jBL+w6c、全台站标注）→ `references/madagascar-singer-station-map.md`
- PMCC 要求 SAC header 含 `t0` 且 Stream 恰好 3 个 trace
- 脚本路径硬编码 `/mnt/e/` 或 `/mnt/f/`，新事件需检查
- Z 分量异常诊断 → `references/z-component-diagnostics.md`
- **台阵反演（PMCC/F-K）⏸️ 暂定不执行**（2026-08-05）：PMCC 是**水听器/OBS 台阵方法**（Cansi 1995，CTBTO IMS 标准），陆地宽频台网（~100 km 间距）对 2-8 Hz T 波空间采样严重不足（需 <1.5 km），实测 4 种方法 0 组有效。且 PMCC ≠ 频域 F-K 波束成形（`t_phase_analysis.pmcc` 是时域互相关时延反演，包内无真正 beamforming）。陆地台网用偏振 R-T 替代 → `references/fk-spacing-limits.md`；方法记录 → seismic-twave-analysis `references/fk-africa-attempts.md`
