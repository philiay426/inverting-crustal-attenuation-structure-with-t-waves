---
name: plotting-standards
category: seismic
description: >
  地震学出版级制图标准：rcParams 设置、配色规范、波形/频谱/玫瑰图/包络图布局、
  文件名约定、P98归一化、到时标注。画图时加载此 skill 参考。
triggers:
  - "画图标准"
  - "期刊出图"
  - "publication figure"
  - "matplotlib rcParams"
  - "波形绘制"
  - "频谱图"
  - "玫瑰图"
  - "包络图"
  - "P98归一化"
  - "配色方案"
  - "figure formatting"
  - "制图规范"
---

# 地震学出版级制图标准

> ** 预处理标准（2026-08-01 更新）**
>
> **去响应**：统一使用 `tr.simulate(paz_remove=paz_dict, remove_sensitivity=True)` 直接出速度。
> `simulate_seismometer() + np.gradient()` 已全面废弃。
>
> **多台站预处理四图**：`seismic-twave-analysis` → `templates/process_image.py`
> — 3×1 布局，per-component P98，每站四张独立图。
>
> **单台分析并排图**：`seismic-twave-analysis` → `scripts/process_single_station.py`
> — 3×2 布局（左=带通包络 ZRT，右=频谱 0-10Hz ZRT），输出到
> `process_flow/3_single_station_analysis/{STA}/`。
>
> `process_obs.py` 模板已标记 DEPRECATED。

> **优先级**：对其他图类（剖面图、偏振分析、地形图等），本文件是唯一权威标准。

---

## 〇、出图规范分两级（2026-09-28 用户明确，先看这条）

用户原话：“首先要按照那个期刊要求的标准来…有些图我只是处理过程中检查一下对错，所以不用太详细。”

- **一级 = 过程检查图**（分析中间图 / QC 图 / 自查图，**默认走这一级**）：只需满足四条——
  ① 文字不重叠、不出框（跑 `fk-array-analysis/scripts/check_fig_text_overlap.py` 与 `fig_ink_margins.py`）；
  ② 轴标签带单位、刻度可读；③ 字号不缩到看不清；④ 不覆盖既有图（新文件名，或先 `mv` 旧文件）。
  **不必**套下面的期刊全套（尺寸 / DPI / 矢量 / 文字极简 / 降饱和配色等）。
- **二级 = 期刊交付图**：按本 skill 正文全部要求做。**执行时机：用户明确说“按期刊要求检查/修改图”
  “文字要符合期刊要求”时才走二级**；他没点名就不要把过程图按期刊规格返工。
- 拿不准时先出一级图给他看；要投/要交时再走二级。
- 其余出图偏好（配色降饱和、全框、文字极简、不标测量值、面板字母…）均属二级，详见下文。

---

## 一、通用规范（所有图都适用）

### 0. 定版基线（2026-10 用户定稿，优先于下文各节的"范围值"）

下文 §1–§7 给的是**范围**（如 8–10 pt）。T 波流程一律用下面这套**固定值**，与范围冲突时以本表为准：

| 项 | 定值 |
|---|---|
| 栏宽 | 单栏 **3.5 in**、双栏 **7.1 in**（对应 A4 双栏版心 8.9 / 18.0 cm）；高度按面板数 × 2.4 in |
| 字号 | 刻度 **9**、轴标签 **11 粗体**、标题 **12**、图例 **8.5**、面板标 `(a)` **9**、图内注记 **7.5–8** |
| 分辨率 | PNG **600 dpi**（底线 300）；需要矢量才出 PDF，`pdf.fonttype=42` |
| 轴框 | 四边全框 `0.6 pt`、刻度朝内、无次刻度（对数轴若开次刻度用 0.4 pt） |
| 线宽 | 数据 1.0 pt、拟合 1.0 pt、参考线 0.6 pt 虚线；标记 4–5 pt 实心、`mew=0.4` |
| 图例 | 轴外上缘 `bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left", mode="expand", ncol=<项数>, frameon=False`；不放轴内 |
| 图内文字 | 全英文；禁 `->`、`@`、`Δ`、`°`（写 `at`、`deg`）；千分位不加逗号 |
| 轴标签 | `物理量 (单位)`，如 `Frequency (Hz)`、`Distance (km)`、`$Q^{-1}$`；无量纲不写单位 |
| 命名 | `<对象>_<内容>_<条件>.png`（条件如 `_w500`、`_proj`、`_landtime`）；不覆盖旧图 |
| 自检 | 脚本级：轴范围覆盖数据、刻度数 = 标签数、`savefig` 后核对文件时间戳 |

**可直接复制的 rcParams（与上表一一对应，不许各脚本自改）**

```python
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9,                 # 刻度
    "axes.labelsize": 11,           # 轴标签
    "axes.labelweight": "bold",
    "axes.titlesize": 12,
    "axes.titleweight": "normal",
    "legend.fontsize": 8.5,
    "axes.linewidth": 0.6,
    "axes.edgecolor": "black",
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,       # 四边全框
    "xtick.major.size": 3, "ytick.major.size": 3,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.minor.visible": False, "ytick.minor.visible": False,
    "lines.linewidth": 1.0,
    "savefig.dpi": 600,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})
```

> 若目标期刊有明确模板（如 GRL/JGR 的栏宽），**以期刊模板为准**，并回来改本表。


### 1. 尺寸

| 栏宽 | 最大宽度 | 用途 |
|------|---------|------|
| 单栏 | **9.5 cm（3.74 in）** | 小图、频谱图 |
| 双栏/通栏 | **19 cm（7.48 in）** | 大图、多面板 |

- 图高上限约 23.7 cm（9.3 in）
- **不要用屏幕好看但印刷失真的尺寸**（如随手写的 18×12 in）
- 提前换算好 `fig_w = width_cm / 2.54`

```python
FIG_W_CM = 19.0          # 双栏
FIG_W_IN = FIG_W_CM / 2.54
FIG_H_IN = FIG_W_IN * 0.62   # 高宽比按内容定
fig, ax = plt.subplots(..., figsize=(FIG_W_IN, FIG_H_IN))
```

### 2. 字体

- **无衬线**：Arial → Helvetica → DejaVu Sans（退化链）
- **字号**：坐标轴标签/图例/刻度 **8–10pt**，标题 **10–12pt**
- **矢量嵌入**：`pdf.fonttype=42, ps.fonttype=42`（禁用 Type3）
- **🔴 图内文字一律用英文**：投稿发英文期刊，标题/轴标签/图例/图内文本全用英文（如 `1s win/75% overlap`），中文只允许出现在脚本注释和 print 中。WSL/Linux 默认字体无中文字形，中文渲染成豆腐块（2026-08-04 实测）。仅个人查看用图需要中文时，装 Noto Sans CJK 并设 `rcParams['font.sans-serif']=['Noto Sans CJK SC']`

```python
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9,
    'axes.labelsize': 10, 'axes.titlesize': 10,
    'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})
```

### 3. 分辨率与格式

| 图类型 | 分辨率 | 格式 |
|--------|--------|------|
| 纯矢量（波形/散点/线条） | — | **PDF/EPS**（首选） |
| 纯栅格（照片类） | ≥300 dpi | PNG |
| 组合图（矢量线条+栅格热图/地形） | **600 dpi** | PNG + PDF（双输出） |

- **投稿优先 PDF/EPS**
- `pcolormesh` / `imshow` 热图存 PDF 时加 `rasterized=True`，避免文件过大
- 文字和线条保持矢量可编辑

```python
# 组合图：PNG 600dpi + PDF 矢量
fig.savefig(out_png, dpi=600, bbox_inches='tight')
fig.savefig(out_pdf, bbox_inches='tight')  # 矢量版
```

### 4. 线条与边框

- 线宽 **≥0.5pt**（`lw=0.5` 起步，太细印刷断线）
- 坐标轴刻度 **朝内**（`xtick.direction='in'`）
- 去掉多余装饰网格线和边框；一般只保留左、下边框

```python
plt.rcParams.update({
    'axes.linewidth': 0.6,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
})
# ⚠️ 用户偏好（2026-09-23，他直接纠正过本行）：**坐标框一律全框** —— 四条边都画
#   （`axes.linewidth=0.6`）、并开四边朝内刻度：`'xtick.top': True, 'ytick.right': True`。
#   不要半边框（半框是期刊常见做法但不是他要的）；表格同样是全框细黑线。
```

### 5. 配色

| 用途 | 色板 | 禁止 |
|------|------|------|
| 连续数据（频谱/热图） | **viridis / plasma / inferno / magma** | jet / rainbow |
| 离散类别（台站/震相） | **tab10 / Set2** | 连续色板表示离散类别 |
| 波形线条 | 黑/灰 `#555555` | 除非有物理含义 |
| 到时标记 | P=亮绿 `#4daf4a`, S=亮橙 `#ff9848`, T-P=亮蓝 `#4d9fff`, T-S=亮红 `#ff4d4d` | |

```python
CMAP = 'viridis'                # 频谱图
ARRIVAL_COLORS = {
    'P':   '#4daf4a',           # 亮绿
    'S':   '#ff9848',           # 亮橙
    'T-P': '#4d9fff',           # 亮蓝
    'T-S': '#ff4d4d',           # 亮红
}
```

### 6. 多子图排版

- 子图编号 `(a)(b)(c)` 放**左上角**（除非明确说不要）
- 共享坐标轴用 `sharex/sharey`，避免重复刻度
- **图例、色标、标题之间预留独立间距**（最常见排版事故就是叠压）
- 对齐：同一列子图的坐标轴用 `constrained_layout` 或手动 `subplots_adjust`

```python
fig, axes = plt.subplots(3, 1, figsize=(w, h), sharex=True)
plt.subplots_adjust(hspace=0.08, top=0.92, bottom=0.10, left=0.07, right=0.94)
```

### 7. 坐标轴与标注

- **所有物理量带单位**，写在轴标签里：`Frequency (Hz)`、`Time (s)`、`Velocity (m/s)`
- 不要只写变量名
- 刻度数量适中，避免过密
- 科学计数法用 `ScalarFormatter(useMathText=True)`

```python
ax.set_xlabel('Time from origin (s)', fontsize=10, fontweight='bold')
ax.set_ylabel('Frequency (Hz)', fontsize=9, fontweight='bold')
```

---

## 二、分图类规范

### A. 地形图 / 台站分布图（PyGMT / Cartopy）

> **🔴 必须输出地理地图，不是极坐标图！**
> "台站分布图" = 经度-纬度地理地图 + 地形色标。极坐标方位-距离图不符合出版标准，
> 会被直接要求重画。用户说"台站分布图"时，默认期望就是 §A 规范的地理地图。

| 要素 | 规范 |
|------|------|
| **投影** | 墨卡托 `M18c`（大陆尺度） |
| **地形网格** | `@earth_relief_15s`（约450m分辨率，已有缓存 583MB），回退 `@earth_relief_01m` |
| **地形着色** | geo 色标，独立 CPT 文件，底部水平色标 `JBC+w10c/0.3c+h`，标注 `Elevation (m)` |
| **台站着色** | SNR 着色 with viridis，右侧垂直色标 `JMR+w7c/0.3c+v`，**每 10 一档** (`xa10f1`) |
| **等深线** | -1000m 黑色虚线，图例注明 `-1000 m depth` |
| **经纬刻度** | **每 10° 标注 + 每 5° 小刻** (`xa10f5, ya10f5`)，不画满屏网格 |
| **图例** | 右下角白框：台站(colored by SNR) + 震中 + -1000m 各独立条目 |
| **台站标签** | 标注 land_km 最近 + 最远 + SNR 最高的 2-3 个台站（全站名标注时逐台），**引线灵活判断**：标签贴台站(<0.6°)不画引线，偏离远才画**细实线**，不遮挡数据点 |
| **文字描边** | **4 方向白字偏移 0.05cm + 黑字**（勿用 `fill="white"` 白底——难看；勿用 8 方向白字——offset 大了叠成白块/锯齿，实测 4 方向 0.05cm 最干净） |
| **方位角射线** | 从震中到标注台站的虚线 + 方位角度数标签 |
| **DPI** | 200（栅格地形用） |

**代码模板（权威版）：** `/home/sy/codex/scripts/plot_station_map_final.py`

```python
# 关键 pygmt 调用序列
fig = pygmt.Figure()
elev_cpt = os.path.join(tempfile.gettempdir(), "elev_relief.cpt")
pygmt.makecpt(cmap="geo", series=[-8000, 5000, 100], continuous=True, output=elev_cpt)
fig.grdimage(grid="@earth_relief_15s", region=R, projection="M18c",
             cmap=elev_cpt, shading="+a315+nt0.8+m0")
fig.coast(shorelines="0.3p,gray30", borders=["1/0.3p,gray50","2/0.2p,gray70"],
          resolution="h", water="lightblue@60", land="white@80")
fig.basemap(frame=["xa10f5", "ya10f5", "WSne"])  # 10°大标 + 5°小刻
fig.grdcontour(grid="@earth_relief_15s", region=R, levels=[-1000], pen="0.8p,black,--")

# SNR 色标
pygmt.makecpt(cmap="viridis", series=[floor(snr_min), ceil(snr_max), 1])
fig.plot(data=tmpfile, style="t0.5c", cmap=True, pen="0.3p,black")

# 双色标
fig.colorbar(frame="xa10f1+lSNR (dB)", position="JMR+w7c/0.3c+o1c/0c+v")  # 每10一档
fig.colorbar(cmap=elev_cpt, frame="xa2000f500+lElevation (m)",
             position="JBC+w10c/0.3c+h+o0c/1.6c")

# 图例
legend_spec = (
    "G 0.1c\\n"
    "S 0.3c t 0.5c gray30 0.3p,black 0.8c Station (colored by SNR)\\n"
    "G 0.15c\\n"
    "S 0.3c a 0.4c red 0.7p,black 0.8c Epicenter\\n"
    "G 0.15c\\n"
    "S 0.3c - 0.6c - 0.8p,black,-- 0.8c -1000 m depth\\n"
)
fig.legend(spec=lf.name, position="jBR+w4.3c+o0.3c/0.3c", box="+gwhite+p0.5p,black")
```

**`grdcontour` 等深线标注：4 条实测结论（2026-09-29，马达加斯加台站图）**

- **标注后线会变实线**：一开 `annotation`，GMT 把被标注的那条线改成 **annotated pen**（默认 `0.75p,black` 实线）。
  要保持虚线必须同时给两把笔：`pen=['a0.7p,black,--', 'c0.7p,black,--']`（a = 标注线 / c = 普通线）。
- **标注文字只能是纯数字**（如 `-1000`）：`+u m` 加单位、把 `annotation` 写成间隔式（`'1000'`）实测都不生效；
  要写“水深 1000 m”这类整句只能自己 `fig.text` 放文字。
- **字号只有一条可控路**：`label_placement='n2+f6p,Helvetica-Bold,black'`（`n` 模式可接 `+f`）；
  `'d4c'`（沿线每 4 cm 一个）**不能**接 `+f`（GMT 报 `not a valid number` 并静默回退），
  而 `pygmt.config(FONT_ANNOT_PRIMARY=...)` 对等深线标注**完全无效**（4p/6p/16p 暗像素像素一模一样）。
- **疏密二选一**：`d<dist>` 按距离（小图约 10 个/幅，推荐）；`n<n>` 每条**线段** n 个——
  地形碎的海岸线有上百个短线段 ⇒ 标注泛滥（实测大区域图 `n2` 出 234 个，`d4c` 只 11 个）。
- **无视觉验收的客观核验法**：`fig.savefig('x.eps')` 后数 PostScript 明文
  （`re.findall(rb'\(([^)]{0,20})\)')` 里含 `1000` 的串数）+ 看 `setdash` 确认虚线仍在；
  对比字号时先量**暗像素**（画布尺寸须把 frame 排除在 config 之外，否则被坐标轴字号干扰）。

**pygmt/GMT 文字「白描边 + 加粗」：实测 3 条（2026-09-29，震中信息文字）**

- **描边写法**：字体串后接 `=`*pen*（GMT pstext 官方：*to draw outline fonts you append =pen to the font
  specification*），如 `font='9.5p,Helvetica-Bold,black=0.6p,white'`。pygmt **无 `stroke` 参数**，且 `-F` **没有 `+o`**
  （传 `+o` 直接报 `Unrecognized modifier +o`）。
- ⚠️ **白描边会把笔画“啃瘦”，而且必须在真图尺度上量才看得出来**：同一串 9.5p 字，在 **600 dpi + 9.8 cm 宽投影**
  （成品尺度）下加 `0.6p,white` ⇒ 暗墨迹 **−62%**（字发虚）；换到 300 dpi + 12 cm 的测试画布只 **−7%**。
  ⇒ 验收字重/描边的测试画布必须用**成品的 dpi 与图幅尺寸**，否则会漏掉这个量级的偏差。
- **“白描边 + 加粗”要画两遍**（GMT 没有更粗的字重，加粗只能靠同色描边）：
  ① `font='<base>,white=1.0p,white'` 先画白衬底；② `font='<base>,black=0.25p,black'` 再压黑字加粗。
  实测（真尺度，相对无描边纯粗体）：黑 0.15p = +17%、**0.25p = +27%（推荐）**、0.5p = +48%（偏糊）；
  黑底上白圈仍在（白像素 >0）⇒ 两遍法同时保住“加粗”与“白圈”。
- 无视觉验收的量化判据：白底图上数**暗像素**（= 字重）、黑底图上数**白像素**（= 白圈）——两者分开量，
  才能分辨“真的加粗了”与“白圈被吃掉了”。

### B. 波形图（Record Section / 三分量）

#### B.1 预处理概览图（每台站 4 张）

**数据源**：到时从 `stations.csv` 读取，不要自己计算。
CSV 列：`station,snr,dist_km,land_km,ocean_km,P,S,T-P,T-S`

**共用参数**：
- 分量：Z/R/T（NE→RT 旋转）
- 横轴：`t_min=0, t_max=T-S+500s`（发震时刻为零点）
- 发震标记：灰点线 `ax.axvline(0, color='gray', lw=0.5, ls=':', alpha=0.5)`
- 分量标签：右上角 `ax.text(0.985, 0.88, label, ha='right')`

##### 图1: 原始波形（Raw）

| 参数 | 值 |
|------|-----|
| 数据 | 原始 SAC（不 demean） |
| 归一化 | 三分量统一 P98 |
| 线型 | 黑色 `lw=0.3` |
| 纵轴 | ±2.5 |
| 到时标注 | **不标**（仅发震线） |
| 单位标签 | 中间面板 `Counts` |

##### 图2: 去响应波形（Remove Response）

同图1，额外：
- 单位：`Velocity (m/s)`
- **标到时**：P/S/T-P/T-S 彩色虚线 + 图例（顶部面板）

##### 图3: 带通滤波 + 包络（Bandpass） — 2026-07-31 定稿

| 参数 | 值 |
|------|-----|
| 滤波 | 2–8 Hz Butterworth 4th |
| 波形 | 深灰 `#555555`，`lw=0.3, alpha=0.5` |
| 包络 | 红色 `#d62728`，`lw=0.8, alpha=0.8`，**不填充** |
| 包络算法 | **Hilbert → abs → 6s boxcar 平滑** |
| 归一化 | 三分量**包络**统一 P98 |
| 图例 | **左上角**，标注到时数值：`P=263s S=479s T-P=1207s T-S=1258s` |
| 虚线 | `lw=0.6, ls='--'`（比其余图细） |

```python
# Hilbert 包络
from scipy.signal import hilbert
analytic = hilbert(signal)
env_raw = np.abs(analytic)
win = int(6 * sr)
kernel = np.ones(win) / win
env = np.convolve(env_raw, kernel, mode='same')
```

##### 图4: 频谱图

见 §D。

#### B.2 记录剖面图（2026-08-04 最新）

**实现**：`seismic-twave-analysis/templates/waveform_proflie.py`（陆地距离分箱取 SNR 最高台站）

- **频段默认 0.5-2 Hz 低频**：T 波低频成分衰减慢，远处台站（480km+）仍可见；2-8Hz 高频几百公里就耗尽——这是"包络画成直线"的真正根因（归一化不是）
- **3s 平滑**（Hilbert 包络）：8s 会把峰值抹平、纺锤形变钝
- **收紧窗口**：`t_min = min(tt_tp)-30, t_max = max(tt_ts)+60`（信号充满图面）
- **包络独立归一化**：每条 trace 先扣噪声基线（t<t_min 中位数），再算全局 env P98（**不复用波形 global_p98**）
- **主图**：灰色填充包络(`#888888` alpha=0.45) + 黑色细波形叠加(`lw=0.25`)；另有 per-bin P98 辅助图 + 纯波形版
- 横轴：`Time from shore arrival (s)`（陆地走时对齐，t=0=波上岸）
- `--band` 参数切换频段（`0.5-2.0` 默认 / `2.0-8.0` 近处）

### B.3 单台分析并排图（Bandpass + Spectrogram）— 2026-08-01 新增

**单台分析标准输出**：3×2 GridSpec 布局，一张图包含带通包络和频谱。

| 参数 | 值 |
|------|-----|
| 布局 | 3行×2列：左=带通包络(ZRT)，右=频谱(ZRT) |
| 尺寸 | `fig_w=19/2.54`, `fig_h=fig_w*0.60` |
| 左列-波形 | 灰 `#555555` lw=0.25 + 红包络 `#d62728` lw=0.6, 60点平滑 |
| 左列-归一化 | 三分量统一 P98，纵轴 ±2.5 |
| 右列-频谱 | 0-10Hz, viridis, P2-P95 vmax |
| 频谱参数 | `nperseg=min(sr*12.8, max(256,n//8))`, `noverlap=nperseg*0.9375`, `nfft=nperseg*2` |
| 到时 | P/S/T-P/T-S 虚线标于左列 |
| 文件 | `{STA}_bandpass_spectrogram.png`, 300 dpi |
| 路径 | `process_flow/3_single_station_analysis/{STA}/` |

```python
from matplotlib.gridspec import GridSpec

fig = plt.figure(figsize=(fig_w, fig_h), dpi=150)
gs = GridSpec(3, 2, fig, hspace=0.10, wspace=0.06,
               top=0.90, bottom=0.12, left=0.06, right=0.92)

for i, (lbl, d) in enumerate(comps):
    ax_l = fig.add_subplot(gs[i, 0])  # 左: 带通包络
    ax_r = fig.add_subplot(gs[i, 1])  # 右: 频谱图
```

| 规则 | 原因 |
|------|------|
| `ax.set_aspect('equal')` | 否则椭圆/直线度被拉伸变形 |
| 两正交分量坐标轴范围必须一致 | N-E 平面 N和E刻度范围相同 |
| 标出时间窗口 + 方位角/入射角 | 可复现 |
| 多窗口对比用同一色标按时间渐变 | 不分开画 |

### D. 频谱图 / 时频图（Spectrogram）— 2026-07-31 定稿

#### 纯频谱图（3行×1列，Z/R/T）

| 参数 | 值 |
|------|-----|
| 色标 | **viridis** |
| vmin/vmax | 三分量共享 **p2–p95** |
| 布局 | 3行×1列，sharex |
| 时间窗 | `t_min=0, t_max=T-S+500s` |
| 横轴 | `Time from origin (s)` |
| 频率范围 | 0–10 Hz（`f_mask=f<=10`, `set_ylim(0,10)`） |
| 分量标识 | 右上角白色文字 Z/R/T |
| 频率标签 | 中间面板左侧 `Frequency (Hz)` |
| 到时虚线 | P/S/T-P/T-S，`lw=1.0, ls='--', alpha=0.85` |
| 到时图例 | 顶部面板上方水平排列，无边框 |
| 标题 | `STATION — Spectrogram  Fmin–Fmax Hz` |

```python
# 频谱参数（细颗粒）
nperseg = min(int(sr * 12.8), max(256, min_len // 8))
noverlap = int(nperseg * 0.9375)
nfft = nperseg * 2   # 最多 2x，4x 超时
```

#### 🔴 SAC 时间轴对齐：`sac.o` 偏移陷阱

SAC 头 `O` = 数据起始到发震的秒数（常为 600s）。
**未减此偏移** → 所有到时标记全错。

```python
t0_offset = tr.stats.sac.o       # 如 600s

# ✅ 正确：切片后计算发震相对时间
i1 = int((t_min + t0_offset) * sr)
i2 = int((t_max + t0_offset) * sr)
t = np.arange(i2 - i1) / sr + (i1 / sr - t0_offset)

# ❌ 常见错误：漏掉 i1/sr
# t = np.arange(i2 - i1) / sr - t0_offset
```

⚠️ **`sac.o` 可能未设置**（多事件实测：BOSA 2020 有、Prince Edward 2019 无）。使用 fallback 函数：

```python
def get_t0_offset(tr, default_origin):
    try:
        t0 = tr.stats.sac.o
        if t0 and t0 > 0: return float(t0)
    except (AttributeError, TypeError): pass
    return default_origin - tr.stats.starttime
```

详见 `seismic-twave-analysis` → `references/sac-o-fallback.md`。

### E. 数据图表 / 统计分析图

| 规则 | |
|------|-----|
| **误差棒/置信区间** | 数据支持就画，图注说明是 SD/SE/CI |
| **显著性** | p 值/星号标注方式在图注说明 |
| **图例** | 优先空白区域，其次图外，不遮挡数据 |
| **基线** | 柱状图/折线图不要截断非零基线，除非有理由并在图注说明 |
| **双纵轴** | 慎用（AGU 审稿人常要求拆成两个子图） |

---

## 三、代码模板

### 3.1 rcParams 完整模板

```python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9,
    'axes.labelsize': 10,
    'axes.titlesize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'axes.linewidth': 0.6,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})
```

### 3.2 三分量子图模板

```python
FIG_W_IN = 19.0 / 2.54
FIG_H_IN = FIG_W_IN * 0.55
fig, axes = plt.subplots(3, 1, figsize=(FIG_W_IN, FIG_H_IN), sharex=True, dpi=150)
plt.subplots_adjust(hspace=0.08, top=0.92, bottom=0.10, left=0.07, right=0.94)
```

### 3.3 到时标注模板

```python
ARRIVALS = [
    ('P',   263.1, '#4daf4a'),
    ('S',   479.2, '#ff9848'),
    ('T-P', 1207,  '#4d9fff'),
    ('T-S', 1258,  '#ff4d4d'),
]

# 虚线 + 图例（仅顶部面板）
for i, ax in enumerate(axes):
    for j, (ph, tt, clr) in enumerate(ARRIVALS):
        ax.axvline(x=tt, color=clr, lw=1.0, ls='--', alpha=0.85,
                   label=ph if i == 0 else None)
    if i == 0:
        ax.legend(loc='upper right', fontsize=7, frameon=False, ncol=4,
                  handlelength=1.2, handletextpad=0.4, columnspacing=0.8,
                  bbox_to_anchor=(0.985, 1.18))
```

### 3.4 P98 归一化

```python
# 三分量统一 P98
all_data = np.concatenate([z, r, t])
p98 = np.percentile(np.abs(all_data), 98)
data_norm = data / p98
```

### 3.5 ZRT 旋转

```python
baz = 127.5  # 从 stations.csv
br = np.radians(baz)
r = n_comp * np.cos(br) + e_comp * np.sin(br)
t_comp = -n_comp * np.sin(br) + e_comp * np.cos(br)
```

---

## 五、Excel 台站表格格式（2026-08-02 定稿）

### 5.1 台站信息表（Excel）

| 要素 | 设置 |
|:---|:---|
| 字体 | **微软雅黑** |
| 标题行 | **16号，加粗，居中** |
| 中文表头（第2行）| **14号，黑色加粗，居中** |
| 英文表头（第3行）| **14号，白色加粗，蓝底 `#4472C4`，居中** |
| 数据行 | **14号，居中** |
| 框线 | **全框细线，黑色**（含标题行） |
| SNR 高亮 | SNR≥阈值行浅绿底 `C6EFCE` |

### 列序（13列）

```
中文: 台站 | 方位角 | 震中距 | 陆地距 | 海洋距 | 信噪比 | P波 | S波 | T-P | T-S | 经度 | 纬度 | 反方位角
英文: STA  | AZ(°)  | DIST   | LAND   | OCEAN  | SNR(dB)| P(s)| S(s)|T-P(s)|T-S(s)| LON | LAT | BAZ(°)
```

### 代码模板

```python
from openpyxl.styles import Font, Alignment, Border, Side

thin_border = Border(
    left=Side(style='thin', color='000000'),
    right=Side(style='thin', color='000000'),
    top=Side(style='thin', color='000000'),
    bottom=Side(style='thin', color='000000')
)

# 标题行
info_cell.font = Font(name='微软雅黑', bold=True, size=16)
info_cell.alignment = Alignment(horizontal='center', vertical='center')

# 表头（中文 + 英文双行）
c.font = Font(name='微软雅黑', bold=True, size=14, color='000000')
c.alignment = Alignment(horizontal='center', vertical='center')

# 数据行
cell.font = Font(name='微软雅黑', size=14)
cell.alignment = Alignment(horizontal='center', vertical='center')

# 全框（含标题行）
for row in ws.iter_rows(min_row=1, max_row=n_rows, max_col=13):
    for cell in row:
        cell.border = thin_border
```

执行脚本 `step1_station_table_and_map.py` 已按此更新，统一输出：
1. `stations.xlsx` — 格式化表格
2. `stations.csv` — 数据
3. `figure_map.png` — 全站名标注版 (15s 地形)
4. `figure_map_final.png` — 期刊级版 (双色标, 图例, 方位角射线)
5. `profile_Madagascar.png` / `profile_waveform_Madagascar.png` — 马达加斯加剖面
6. `profile_Mainland.png` / `profile_waveform_Mainland.png` — 非洲南部大陆剖面

剖面按经度 43° 分界：> 43°E = 马达加斯加，≤ 43°E = 非洲大陆。

地图模板: `/home/sy/codex/scripts/plot_station_map_final.py`（权威版）。
SNR 色标每 10 一档 (`xa10f1`)，刻度每 10° 大标、5° 小刻 (`xa10f5`)。

### 台站标签防重叠（2026-09-06 密集台阵 XA 38 台定稿）

每个被标注台站尝试多个偏移方向，选第一个不与已有标签重叠的方向。
特殊台站可用字典指定固定偏移（如 `{"AMPY": (-0.8, -0.5)}` 左下）。

**密集台阵（台站经度间距 ~1°，如南非 XA 38 台）的完整方案**：
1. **偏移候选全方向螺旋**（32 个，勿按台站相对震源方位 `sx/sy`，否则密集区全部同向偏移互相挤）
2. **纵向大偏移提前**（台站横向间距 ~1° 放不下时纵向错开，最大 ±3.6°）
3. **fallback 螺旋网格搜索**（候选用尽仍重叠时：半径 0.6~5.6°、步 0.4°、方位步 40° 找更远空白，杜绝 fallback 到第一个 offset 硬放导致重叠）
4. **字号 8p + box 尺寸贴合文字**（8p 粗体 "SA74" 宽约 1.1°，box 用 1.45°，避免 box 过宽误判/漏判）
5. **画布放大 M18c→M26c**（比例尺变大，文字相对度数变小，更宽松）
6. **文字白色描边**：4 方向白字偏移 0.05cm + 黑字（`fig.text(font="8p,Helvetica-Bold,white", offset=...)` 画 (1,0)/(-1,0)/(0,1)/(0,-1) 四方向，再画黑字）。勿用 `fill="white"` 白底（难看）；勿用 8 方向白字（offset 大了叠成白块/锯齿，实测 4 方向 0.05cm 最干净）
7. **引线灵活**：`dist = hypot(dx,dy)`，`dist >= 0.6°` 才画细实线 `0.3p,gray45`，紧贴台站不画

间距阈值 `LABEL_MIN_DIST = 1.5`（度）。

### 期刊版台站分布图：只标关键台站（2026-09-24 Zhou2021 复现定稿）

投稿版不必全站标注，只标关键台（示例台/方法所用台，如 SA18），图面立刻干净、无重叠问题：
`plot_station_map(..., label_only=["SA18"])`，参数区 `LABEL_ONLY`。关键台另用**红三角**突出
（叠在 SNR 着色三角之上，0.60c），图例同步加一条 `Station SA18`（论文 Fig.1 约定）。

- **图幅换算**：`MAP_WIDTH="15.5c"` + 右侧色标 → 18.97 × 19.88 cm，满足双栏 19 cm / 限高 23.7 cm。
  ⚠️ Mercator 下本区 h/w ≈ 1.24（纬度比例 ≈ 经度比例/cos(中纬)）：地形色标若仍横排置底，
  总高会到 ~26 cm 直接超限 → 两条色标都竖排右侧（`JMR` + 不同 dy 偏移错开）。
- **标签偏移一律以 cm 定义再换算成度**（`_dpx = 经度跨度/图宽cm`，`_dpy = _dpx/cos(中纬)`）。
  写死度数会在换图幅时失效：26c 下 −0.49° 刚好，缩到 15.5c 后三角半高已 0.7°，标签会压住三角。
- **事件文字不要锚在震中右侧**：12p 全串约 4.2 cm 宽，小图幅下折合 ~9° 经度，靠右必被图框裁掉
  → `justify="BR"` 锚在沙滩球左上方，文字向左延伸。
- **图例符号必须与图上一致**：图上震中是沙滩球时，删掉「红星 Epicenter」条（GMT 图例画不出沙滩球），
  不要保留图上不存在的符号。
- **PDF 陷阱**：地形底图占 80%+ 面积的图，GMT 出 PDF 会内嵌高分辨率栅格 → 单张 ~100 MB，
  不适合投稿；直接给 600 dpi PNG（满足期刊要求）。确要 PDF 时先降 `GMT_GRAPHICS_DPU`。
- 交付前用 PIL 读 PNG 的 `info['dpi']` 反算物理尺寸（px/dpi*2.54）核对是否落在栏宽/限高内。

### 台站质量控制

`step1` 默认排除毛刺台站（`DEFAULT_EXCLUDE` 集合）。
另支持矩形排除框 `EXCLUDE_BOXES`（如 `[[32.0, 37.5, -11.5, -7.0]]` 过滤北部噪声簇）。
**不要用纬度一刀切**（会误伤好台站）。

SNR 计算：动态噪声窗 `max(20, min(300, P_time × 0.8))`，数据强制 `astype(np.float64)` 避免溢出。
回退值 0.0（非 NaN）。

> **地图改进细节**：排除框过滤、标签防重叠算法、SNR 动态噪声窗 → `references/station-map-improvements.md`

### 波形剖面图范围与图例（2026-09-07 Zhou 2021 复现定稿）

- **y 轴（距离轴）不留大片空白**：从数据实际下限附近起（`y_min = max(0, y[0] - amp_km*3)`），
  勿从 0 起——若台站 land/dist 从 600+ km 开始，0~600 全是废白。距离范围贴合数据 + 振幅余量。
- **全时段剖面**：震中距 y 也从 ~min(dist) 起（如 XA 台站 2252-2800 → y 从 2200 起）。
- **图例默认放图内右上角**（`loc="upper right", borderaxespad=0.4`，去掉 `bbox_to_anchor=(1.01,...)` 外挂）。
  图内右上只轻微压到最上方波形尾部（尾波衰减段），可接受。
- **图宽**：剖面图 18" 过宽（横长条），13" 左右合适（字号经 FONT_SCALE 等比缩放保持协调）。
- 全时段剖面标题/x轴不加粗、仅 y 轴加粗（用户偏好）；标题精简如 `T-wave profile | 日期 Mw`。

## 四、工作流

- 一张图一个独立脚本或合并版 `plot_all.py`（内部每张图独立函数，数据加载一次）
- 脚本名：`plot_1_raw.py` / `plot_2_deconvolve.py` / `plot_3_bandpass.py` / `plot_4_spectrogram.py`
- 输出放项目 `BASE` 目录

### 图件文件名用中文（2026-10-14 用户要求）

- **交付/存档的图件文件名用中文**（例：`Q随频率_南部6台去CPSM_两窗对照.png`、`S8八台_波形质控.png`），
  不要 `fig_q_s6_noCPSM_qf.png` 这类纯英文标签；图内文字仍一律英文（见 §一.2）。
- 生成脚本留参数口子，别把名字写死在 savefig 里：
  `FIG_STEM = (os.environ.get('Q_FIGNAME') or f'Q随频率_{RUN_TAG}').strip()`，
  图名 = `{FIG_STEM}.png` / `{FIG_STEM}_逐频点拟合.png`；数据文件（csv/log）仍用英文标签 `q_<TAG>_*` 便于脚本检索。
- 已有产物改中文名时：`mv` 保持 mtime（内容不变）→ **grep 同步文档里全部引用**（产物表、跑法命令块）→
  `ls -la --time-style=full-iso` 核对落盘、再 grep 一次确认没有残留旧名。

## pygmt 多行文字：`fig.text` 不认 `\n`（会把多行拼成一整行）

- 症状：写成 `fig.text(text='A\nB\nC', ...)`，图上却是**一整行**（换行被吞掉）。内容越长越糟：
  5 行约 105 字符 → 一行约 16 cm 宽，`no_clip=True` 时把**画布**从 10.9 cm 撑到 16.8 cm
  （表现为"图例出框 / 图例跑到图外"，而图例自身 anchor 完全没问题，很容易误判成定位错）。
- 实测（pygmt 0.19）：`text='A\nB\nC'` 与 `text='A'` 的墨迹范围**完全一样**（都只有一行）；
  `text=['A','B','C']`（列表）报 `GMTValueError: Invalid value` —— 两种写法都不行。
- 正确做法：多行图例/说明改用 **GMT legend 机制**（`fig.legend(spec=io.StringIO(...))`），
  每行一个 `S <dx> <symbol> <size> <fill> <pen> <text>`；多行与自动尺寸都由 GMT 管。例：
  ```python
  spec = ['G 0.06c', 'H 8.5p,Helvetica,black 0 Station groups']
  spec += [f'S 0.30c t 0.34c {color} 0.3p,black {name} ({n} sta)' for name, color, n in rows]
  fig.legend(spec=io.StringIO('\n'.join(spec)), position='jBL+o0.35c/0.35c',
             box='+gwhite@20+p0.3p,gray45')
  ```
- 验收（不靠肉眼）：**两版取差** —— 关掉该图例再渲一版，两图差集 = 图例足迹；
  再与地图边框比，即可判"在框内 / 是否与另一图例压盖 / 是否把画布撑大"。

## 六、出图前规避清单（11 类高发问题 · 2026-09-15 并入）

⚠️ 自查口径（2026-09-14 用户定）：**不做自我视觉检查** —— 只有用户明确要求看图、或任务本身是图件质检，
才调 `vision_analyze`；其余一律只做**脚本级自查**（轴范围/单位/标签/是否出框/表格截断）+ 数值逻辑核查。
下表是"**画之前就规避**"，不是画完让视觉模型挑问题。

| # | 症状 | 修法 |
|---|---|---|
| 1 | 超长 suptitle 左右被裁 | 标题 ≤ 60 字符；细节放面板标题/文件名，不进总标题 |
| 2 | `ax.table` 列宽溢出、行标签被切 | `colWidths` 控制（总和 ≤ 1）、去掉过宽的 rowLabels 列（关键信息并入第一列）、或表格通栏 |
| 3 | 每面板图例压住数据 | 改**图级共享图例**：`fig.legend(handles=proxy, loc="lower center", ncol=N, frameon=False)` |
| 4 | 色标标签挤到相邻面板 | ①共享轴标签（只最左列标 y、最下行标 x）②`wspace ≥ 0.3` ③色标 label 用短词（`dB`） |
| 5 | 相邻面板标题横向相压 | 标题拆 2–3 行、字号降到 8；每行 ≤ 30 字符（面板宽 ~5 cm 时） |
| 6 | 右侧台名/标注被右边缘裁 | 文字放轴**内**：`ha="right"` + x 取最右样本 |
| 7 | 特殊字符（★ 等）渲染成乱码 | 图内一律用英文词（`(default)`）不用符号 |
| 8 | 对数轴刻度 `3×10¹` 与 `4×10¹` 重叠 | 显式 `set_xticks([10,20,30,50,80])` + 纯数字标签，关掉 minor formatter |
| 9 | 台站地图缩到 1/3 宽后台名互压 | 按面板尺寸单独生成地图（投影宽度 = 面板宽度）+ 只标 2–3 个代表台站 |
| 10 | 面板编号贴在标题里/表格左边缘（看不出属于哪块） | 用 `panel_label(ax, "a")` 把 **(x) 加粗放面板内部左上角**，带白描边；编号不重复出现在标题/行标签里 |
| 11 | 图内自定义简写看不懂（如 `slow@BAZ`、`24.7@126`） | **不用缩写符号**：表头写全并定义顺序（`injected: slow (s/°) / BAZ (°)`），单元格用 `24.7 / 126.4`；文字叙述用 `slow 24.7 s/°, BAZ 126.4°` |

**文字宽度要量不要猜**：`fig.canvas.draw()` 后用 `t.get_window_extent(renderer)` 把串宽换成 mm，和容器宽（面板/栏宽）比，超了直接报警。
两个高频事故：① 面板标题 `set_title(pad)` 与紧贴上框的 `text(y≈1.01)` 抢同一段高度（默认 pad=6 pt 会叠成一团）—— 把标题 pad 提到 12–13 pt、实测值行贴框外 0.3 mm 即可分开；② 串宽超过容器宽（例：`measured slowness 25.13 s/deg (BAZ 133.1 deg)` 8.5 pt 实测 **65.7 mm** > 48 mm 面板）会横跨图面压住 y 轴标签 —— 缩短串或降到 7 pt（同串 41.9 mm），保证 < 容器宽。
自检用的 figure 必须显式 `dpi=` 与存盘一致，否则按默认 100 dpi 量出的宽度与 600 dpi 成品不可比（差 6 倍）。

**脚本级自检清单（画完过一遍）**：图内文字全英文、单位标签与实际口径一致（归一化后标 `Norm. Amp.` 不是 `Velocity (m/s)`）、
多子图/图例/标题/色标无叠压出界、台站标签不重叠裁切、配色合规（连续量 viridis 系，不用 jet/rainbow）、
异常台站（P98 天文数字）已剔除未绑架全局归一化、剖面振幅标尺未把波形压成直线或过度放大、
矢量用 PDF / 栅格 ≥300 dpi。
