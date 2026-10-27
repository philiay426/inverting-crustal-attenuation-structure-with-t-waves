# 出图规范（全链）：步 0-4 的标准图、版式、命名、返工点

所有图都由 `templates/` 的脚本产出，改参数不重写。默认 PNG（300 dpi），要矢量再加 PDF。
幅度只做**脚本级自查**（轴范围/单位/标签/图例是否出框/刻度是否截断），不做视觉审美迭代。

## 步 0：逐站六图（`waveform_image/{STA}/`，脚本 `preprocess_one_station.py`）

**全部 1:1：`FIG_WAVE = FIG_SPEC = (7.2, 7.2) in`，`DPI = 300` → 2160 x 2160 px。无总标题，站名进文件名。**

| 图 | 内容 | 要点 |
|---|---|---|
| `{STA}_1_raw.png` | 原始三分量 | **只标发震时刻线，不标到时**；看空道/缺口/毛刺 |
| `{STA}_2_remove_response.png` | 去响应后三分量（默认速度） | 四条到时线；量纲合理（1e-6~1e-9 m/s） |
| `{STA}_3_spectrogram.png` | **Z/R/T 三面板** STFT，0-10 Hz | 共用一条色标，标签只写 `Power (dB)`；`viridis`；`10log10(P)`、`vmin/vmax = P2/P95`；`hann`、`nperseg=12.8 s`、`noverlap=93.75%`、`nfft=2x`、`mode='psd'`；横轴与波形图同一绝对计时（`XLIM_T0` ~ `tt_ts+XLIM_TAIL`）；`ylim 0~fs/2` |
| `{STA}_4_envelope.png` | Z/R/T 三分量包络 | **只画包络**：不画波形、不镜像、只正半边；黑线 `lw=0.8`；`ylim 0 ~ 1.02x` 全局峰；ylabel `Envelope`，标题不带 `P98 norm.` |
| `{STA}_bandpass_2-8Hz.png` | 2-8 Hz 带通（按 `baz` 旋转成 Z/R/T） | 黑波形 `lw=0.4`；**面板内右上角只标 Z / R / T**，无 ylabel、无频带/单位文字 |
| `{STA}_packet.png` | **SNR 计算图**：Z 分量包络 + 信号窗/噪声窗色块 | 信号窗橙 `#ff7f0e`（`alpha=0.25`）、噪声窗红 `#d62728`（`alpha=0.18`）；标题 `SNR = … dB (Welch 2-8 Hz)`；图例在轴外上方 |

**到时线**：全项目唯一来源 `ARRIVAL_COLORS` —— P `#4daf4a`、S `#ff9848`、T-P `#4d9fff`、T-S `#ff4d4d`。
谱图的到时虚线额外加黑色描边：**彩色 `lw=0.8` + `pe.withStroke(1.5)`**（再粗就糊成一条带）。
示范包/剖面母本另用红蓝两色，属于独立口径，不与本表混用。

**期刊 rcParams**：字体 Arial→Helvetica→DejaVu Sans；刻度 9 pt、轴标签 11 pt 加粗、标题 12 pt、
图例 8.5 pt；四边全框、刻度朝内（`top=True, right=True`）、轴线 0.6 pt；`pdf.fonttype=42`；图内文字全英文。

**横轴**：`XLIM_T0 = 0`（发震时刻）起，到 `tt_ts + XLIM_TAIL`（默认 +500 s）；纵轴保持整段自动。

**SNR**：`20log10(信号 RMS / 噪声 RMS)`，Welch PSD 在 2-8 Hz 上积分；信号窗 T 分量 `[T-S±25]`、Z/R `[T-P±25]`；
噪声窗 `[tt_p-60, tt_p-10]`（`SNR_GAP = 10` 前移；贴 P 会吃进初动 → 噪声底虚高、SNR 虚低，实测可差 20 dB 量级）。
**先整段滤波再切窗**（只滤 50 s 窗时窗边瞬变使噪声虚高，实测差 23 dB）。`SNR_BAND=None` 可留不滤波对照。

一次判一台能不能用，比任何 SNR 数字都直观。图上不写标题，站名进文件名。

## 步 1：台站分布图与记录剖面

### `station_map_all_names.png`（`step1_map_and_table.py`，**默认唯一一版**）

- **全站台名都标、不突出任何单台**：`LABEL_ONLY = []`、`HIGHLIGHT_STA = []`。
  想只标前 3 台那版：`MAP_TOPN=1`（默认不出）。
- 三角按 `snr` 着色，范围 = 实际 SNR 的 `floor`~`ceil`；`0.5 cm` 三角 + 黑边 `0.3p`。
  colorbar 写 `SNR (dB, common noise baseline)`（统一噪声基准口径）或 `SNR (dB)`（自噪口径），
  **不要**写成 "colored by SNR"。
- **配色定稿**：台站 `viridis` + 地形 `geo`（−8000~5000 m、100 m 间隔）+ 山影 `+a315+nt0.8+m0`，
  海 `lightblue@60`、陆 `white@80`。`MAP_STA_CMAP` / `MAP_TERRAIN_CMAP` 可临时换（试过 batlow / gray /
  magma，仍是 viridis+geo 最稳）。岸线 gray30 `0.3p`、国界 gray50/gray70、经纬网 gray40 虚线、
  `-1000 m` 黑虚线 `0.8p`。
- **标签**：8p 粗体黑字 + `WHITE_STROKE = 0.03 cm` 四向白描边。
  ★ 描边字号必须 = 正文字号（`STROKE_FONT = LABEL_FONT.replace("black","white")`）；
  6p 描边 + 8p 正文会让白字缩进笔画里，看起来像空心字/脏白边（旧脚本常见这个错）。
- **避让**：`LABEL_POS = "auto"` 时候选顺序 = 贴符号四向（右/左/下/上）→ 半径 0.95 cm 的 8 方向 →
  半径 1.55 cm 的 8 方向，逐个试第一个不压住已放标签的位置。数十台全标时都能放下；
  只用 fixed `bottom` 会有若干标签互相压住（横排密集台站的物理限制，不是 bug）。
- 图例只写 `Station`，不带彩色说明括号；`DPI = 600`（`MAP_DPI` 可降，预览用 220 够）；`PDF_OUT=False`
  （栅格底图，PDF 约 100 MB）。
- ⚠ 已知瑕疵（未改，会改观感）：`Elevation (m)` 色标画的是 `geo` 原色阶，而海/陆被 `lightblue@60`、
  `white@80` 覆盖过 → **色标读数与图上颜色不完全对应**。修法二选一：色标加同样透明，或去掉覆盖层。

### `fig3_record_sections_3x2.png` / `..._proj.png`（`step1_record_sections_3x2.py`）

- 3 列 x 2 行：上排 = 震中距纵轴 + 相对发震时间（0-1600 s）；下排 = **陆段距离**纵轴 + 相对上岸时间
  （`t - L_ocean/1.48`，-1000~500 s）；三列 Z / R / T。
- 每道画黑包络 `lw=0.35 / alpha=0.6`，**按自身峰值归一**（母本口径；不是全局 P98，也不是逐台 max 混用）；
  振幅比例换算成等效 km（`AMP_KM_DIST` / `AMP_KM_LAND`）。
- 到时：**T-P 红、T-S 蓝**，逐台短竖虚线（**不要**通栏长竖线）。注意这套配色与步 0 的 `ARRIVAL_COLORS`
  不是一套，各按各的母本走，不混用。
- `TT_MODE='proj'`（`stations.csv` 投影走时，输出自动加 `_proj`）与 `'const'`（论文式(1) 陆段恒速 6.0/3.6）
  两套口径各出一版；两者的横向错位量（可差 20 s 量级）是判断"用哪套口径"最直观的证据。
- y 轴从数据实际下限起（取该事件的 `dist` / `land` 最小值），别从 0 起留大片废白。

## 步 2：偏振图（skill `single-station-analysis` 产出，本 skill 不画）

张图/表都由 `run_single_station.py` 出，规格以那个 skill 为准，这里只记读图要点：

| 文件 | 看什么 |
|---|---|
| `<STA>_PS_polar.png` | P/S 理论到时 ±2 s 的偏振，检查 P 段是否径向、S 段是否 45-90° |
| `<STA>_pca_theory_TP_TS.png` | T-P / T-S ±2 s 的 PCA 偏振，T 波窗内主轴是否稳 |
| `<STA>_polarization_fullscan.png` | 三面板全时段：(a) 主轴 R-T 角（左轴）+ 径向/横向扇区占比（右轴）(b) 线性度 (c) 平面度；竖线只画理论到时 |
| `<STA>_polarization_fullscan.csv` | 逐窗数值，判读用连续簇而非全局最小窗 |

读图三条：
1. 接近 45° 的偏振，`rt_angle_med_deg` 会在 0 与 90 之间乱跳 —— 看扇区占比的对比，别读单点角。
2. 判据是**连续簇**：≥1 s、`residual < 10°`、`lin > 0.99`。
3. 线性度/平面度在 2-8 Hz 整段上可能全程恒定（见 0.999 / 0.001），没有分辨力时不作判据。

**不要用步 0 的包络 npz 画偏振**：npz 里的 R/T 是加工过的分量，同一窗算出的 `frac_T` 与
从 `remove_response` 三分量现算能差一大截（实测 0.76 对 0.46）。偏振方向一律现算。

## 步 3：F-K 图组

| 图 | 内容 | 要点 |
|---|---|---|
| 图/表 | 内容 | 要点 |
|---|---|---|
| `fk_map_<TAG>.png` | 慢度图（填色/等高线 + 峰标 + 理论 BAZ 射线 + 送进 fkmap 的实时窗波形） | 峰标注写 `s/deg` 与 `deg`，**不用 @ 或 Δ 符号**；标出峰 −6 dB 内**全部**局部极大，不只 argmax；BAZ 不报小数 |
| `fk_peaks_<TAG>.csv` | 逐峰表：`dB`、`|u|`、`BAZ`、x/y 格点 | 第一行注释写清子阵标签、台数、`smax/spts/frng`；dB 相对峰 |
| `fkinput_<TAG>.txt` | 送进 fkmap 的波形（逐台一行） | 画"时窗落在哪些台上"用；同一 `TAG` 的图与数据必须同名 |
| `fig_<TAG>_subarray_map.png` | 子阵分布：本组红三角 / 另一组蓝三角 / 其余灰三角 + 沙滩球 + 震中→质心大圆虚线 | 600 dpi；-1000 m 等深线虚线；图例只写组名 |
| `fig_fk_window_<TAG>.png` | 时窗画在波形上 | 窗内波形 + 窗边界虚线；常要陆地计时与绝对计时**两种横轴各出一版** |

## 步 4：衰减四图

### 图 1 定案图：Q^-1 随频率（`plot_q_inverse.py`）

| 项 | 规格 |
|---|---|
| 尺寸 | 5.0 x 3.8 in，单轴 |
| 坐标 | 双对数；x 主刻度 2.0 / 2.5 / 3.0 / 4.0，xlabel `log frequency (Hz)`；ylabel `$Q_S^{-1}$` |
| y 刻度 | 统一科学计数法，`0.2-3.0 x 10^-3` 取整数倍，**不出现裸 10^n 标签** |
| 数据点 | 黑实心圆 `ms=5.5`、`mfc/mec='k'`、`mew=0.8` |
| 拟合线 | 黑实线 `lw=1.7` |
| 标注 | 线旁 `a = 0.47 +- 0.39`，**沿直线角度旋转**，自动避让数据点（脚本打印避让间距，要求 > 8 pt） |
| 其它 | 无标题；`tick_params(direction='in', top=True, right=True)`；网格 `alpha=0.22, lw=0.4` |
| 命名 | `<组名>_Q倒数_<口径>_<n>点_拟合线.png` |

用 `Q^-1` 而不是 `Q`：`1/Q = (1/Q0) f^-alpha` 是直线，且截距小、范围窄，低压点不挤在一起。

### 图 2 对比图（`compare_q_objects.py`）

- 同一坐标画 `1/Q = (1/Q0) f^-alpha`，只靠颜色区分（黑/红/蓝），不画点、不画误差带。
- `CASES` 每项 `(名称, Q0, alpha, sigma_alpha, 颜色, 线型, 拟合右端频率)`；论文值直接填发表值。
- 无标题，图例在图内；命名 `<组名>_Q倒数_对比_<对象简称>.png`。

### 图 3 尾波窗落点（`plot_coda_window_section.py`）

- 竖版 5.4 x 7.6 in；`2S` 口径用绝对计时 0-1500 s。
- 逐台画 2-8 Hz 包络（`waveform_t`），逐台按自身峰值归一化（只看形态）。
- 窗条画在基线上：信号窗橙 `#ff7f0e`、尾波窗红 `#d62728`，半高 0.035 道间距；
  理论到时短竖虚线 T-P 蓝 `#1f77b4`、T-S 紫 `#9467bd`。
- 纵轴按**理论 T-S 到时**排，刻度 `STA (d km)`（真实距离轴会让相距 <12 km 的台压在一起）。
- 脚本末尾自带数值自检（每台窗起止 + 与 T-S 的差）：先看数字再看图。

### 图 4 陆地计时剖面（`plot_record_section_landtime.py`）

- 13 x 9 in，每种时间轴上下两栏：(a) 全时段 (b) 放大到直达窗-尾波窗。
- **两套时间轴各出一版**：`land`（x = t - L_ocean/1.48，0 = 上岸）与 `abs`（相对发震）。
- 纵向按陆地距离升序，刻度 `STA (land km)`；`AMP_FRAC=0.45`、`lw=0.25`。
- 图例在坐标区上方（`bbox_to_anchor=(0.0, 1.03)`、`ncol=3`、`frameon=False`）。
- 陆地计时下窗与到时都要**同样扣掉海洋段**（`_off`），否则整体错位。

## 全局约定

- 图用中文名进组目录，表用英文名进 `数据/`（`q_<tag>_detail.csv`、`q_<tag>_byfreq.csv`）。
- 图里不写标题（标题进正文/文件名）；轴标题不加粗。
- 改图先 `mv` 旧图，绝不原地覆盖；旧版进 `_superseded/<日期>_清理/` + `README_去向.md`。
- `/mnt/e` 被 Windows 看图软件占用时 `Errno22` 或静默失败 → 写完 `ls -la --time-style=full-iso` 核对。
- 表头与图内注记一律纯文字：不用 `->` `@` `Δ` `°` 之类符号（`s/deg`、`deg` 直接写词）。
