---
name: seismic-twave-analysis
category: seismic
description: >
  T波分析全流程：预处理(去响应/滤波)、波形剖面、偏振分析、台站地图、
  衰减分析。画图标准见 plotting-standards skill。
---

# T波分析全流程

> **本流程定版模板（`templates/`）**：步 0 `preprocess_one_station.py`（六图 + 包络 npz）、
> `run_batch_preprocess.py`；步 1 `step1_map_and_table.py`（表 + 分布图）、`step1_record_sections_3x2.py`（3x2 剖面）。
> 全部走环境变量，必填缺失即报错；变量表见 `templates/README.md`。
> 目录里 `process_image.py`、`step1_station_table_and_map.py` 是旧版单图脚本，仅参照，不参与本流程。

## 🔴 执行原则：改参数跑，不重新写

**永远不要从零写代码代替现有模板。** 「能跑 5 行 matplotlib」远不如「模板 50 行 pyGMT 出版级」。

❌ **反面案例**：用户要台站分布图 → 你写 `plt.scatter(lons,lats)` 5行 → 没地形、没弧线、没T-P着色 → 用户炸了。

✅ **正确做法**：找到 `templates/` 对应脚本 → 改路径/坐标 → 直接跑。

步骤概览（"处理xxx事件"标准执行链，**目录编号 2026-09-28 用户定**）：

**编号 = 目录编号，都挂在 `process_flow/` 下**：`0` 数据处理 / `1` 台站表·台站分布图·剖面图 / `2` 单台分析 / `3` F-K / `4` 衰减（Q0/alpha）。**步 4 的口径、判据与模板脚本全在 `seismic/twave-crustal-attenuation` skill**（2026-10-23 起规范）。

**0 数据处理（`process_flow/0_data_processing/`）**

| # | 脚本 | 输出 |
|---|------|------|
| 0 | `process_image.py` | 去响应SAC (`0_data_processing/remove_response/*_resp.sac`) + 包络缓存 (`0_data_processing/envelope/*_envelope_2-8Hz.npz`) + 每站四图 (`0_data_processing/waveform_image/{STA}/`，raw/去响应/频谱/带通) |

**1 台站表 / 台站分布图 / 剖面图（`process_flow/1_station_distribution/`，🔴 全部平铺，不再分 `maps_info/`、`profiles/` 子文件夹）**

| # | 脚本 | 输出 |
|---|------|------|
| 1 | `step1_station_table_and_map.py` + 事件目录图件脚本 | stations.xlsx + stations.csv + **2 张台站图**（全台站名、Top3台站名）+ **3 张剖面图** + 记录剖面 `record_section_abs_land_3x2.png` |

**2 单台分析（`process_flow/2_single_station_analysis/{STA}/`，见 single-station-analysis skill）**

默认只做**台站分布图上标注的那 3 个台站**（最近陆 / 最远陆 / SNR 最高）；用户明确指定台站时以用户名单为准。

| # | 输出（默认每台站 3 张） |
|---|------|
| 2a | `{STA}_PS_polar.png`（P/S 偏振） |
| 2b | `{STA}_pca_theory_TP_TS.png`（T-P/T-S 理论窗偏振） |
| 2c | `{STA}_polarization_fullscan.png` + `.csv`（全时段偏振/线性度/平面度） |

推荐一键驱动：`run_single_station.py --event <KEY> --stations STA1,STA2,STA3`（旧 7 张图加 `--all`）

**3 F-K 台阵分析（反演反方位角）（`process_flow/3_fk_analysis/`，见 fk-array-analysis skill）**

进本步先把 1 号的记录剖面图 copy 进来当参照（出图脚本已自动 copy）。

> 注意：`process_single_station.py`（旧步3）已被 single-station-analysis 的流程取代，不再单独执行。

> **步1固定输出**：表格 + CSV + **5张图**，即**2张台站分布图**（全台站名、Top3台站名）和**3张波形/包络剖面图**。不得以旧版 `step1_station_table_and_map.py` 的内置图件替代；每个事件的具体文件名、筛选范围和样式以其事件目录的专用参考文件/脚本为准。

源模板位置：`/home/sy/codex/scripts/step1_station_table_and_map.py`（事件无关版，从 `config/events.yaml` 读参数）。

**改模板后不自动同步到所有事件目录**：源码库脚本更新后，各事件目录下的副本不会自动更新。改完脚本后**主动询问用户**是否需要同步到其余事件目录（用户按需决定，勿自作主张批量拷贝）。

- **预处理 4 张图**（标准见 `plotting-standards` skill → §B.1）：
  1. 图1: 原始波形 — Z/R/T，不 demean，P98归一化，黑线 lw=0.3，±2.5，不标到时
  2. 图2: 去响应波形 — Z/R/T，P98归一化，标到时+图例
  3. 图3: 2-8Hz带通+包络 — 灰波形(`#555555`, lw=0.3, alpha=0.5) + 红包络(`#d62728`, lw=0.8)，**Hilbert→abs→6s平滑**，包络P98归一化，图例左上角标到时数值
  4. 图4: 频谱图 0-10Hz — **viridis**色标（非plasma），P2-P95 vmax，详见 plotting-standards §D
- **脚本组织**：可独立4脚本或合并版 `plot_all.py`/`process_bosa.py`（数据加载一次，每图独立函数）。合并版示例见 BOSA `process_bosa.py`
- **时间窗**：动态 `t_min=0, t_max=T-S+500s`（发震时刻为零点，从 SAC 头 `sac.o` 读取偏移）
- **标签**：Z/R/T 右上角 `ax.text(0.985, 0.88, ...)`，单位标签中面板左侧 y 轴
- **NE→RT 旋转**：`baz` 从 `stations.csv` 取
- **步1图件**：固定输出两张台站图（全台站名、Top3台站名）和三张剖面图；台站范围、命名、归一化、到时标记与视觉细节由事件专用参考文件规定。不要调用旧内置剖面图逻辑覆盖既定五图。
- **波形剖面步另出记录剖面图**（`1_station_distribution/plot_record_sections_abs_land.py` → `record_section_abs_land_3x2.png`，绝对计时 + 陆地计时 3×2 六面板，Fig.3 型）：横轴窗口默认 = 上排 0–1800 s / 下排 −1000–500 s（用户 2026-09-28 定，**要改他会说，不许自行改**；`--xlim-abs/--xlim-land` 可命令行覆盖），出图后脚本自动 copy 一份到 `3_fk_analysis/`（见 §波形剖面）。
- 单台分析只做台站分布图上**显示台站名的那 3 个台站**

## 目录结构（🔴 2026-09-28 用户定稿命名，编号 = 步骤号）
```
<event>/
  raw/                          # 原始 SAC
  raw/SACPZ/                    # 响应文件
  process_flow/                 # 所有阶段的产物都挂在这里
    0_data_processing/          # 步0 数据处理 (process_image.py)
      remove_response/          #   去响应 SAC (*_resp.sac)
      envelope/                 #   与带通图一致的 Z/R/T 包络数组 (*.npz)
      waveform_image/           #   每站四图 {STA}/
    1_station_distribution/     # 步1 台站表 + 2张台站图 + 3张剖面图 + 记录剖面（🔴 全部平铺，不分子文件夹）
    2_single_station_analysis/  # 步2 单台分析 {STA}/（默认 = 台站分布图标注的 3 个台站）
    3_fk_analysis/              # 步3 F-K 台阵分析（反演反方位角）
  # 🔴 转换点定位分析：2026-09-28 用户取消（不成熟），不设目录、不进流程（旧事件里的 3_conversion_point/ 原样留存）
    4_attenuation/              # 步4 衰减 Q0/alpha（口径/判据/脚本见 seismic/twave-crustal-attenuation）
```

> ⚠️ **旧事件仍是旧名**（2026-09-28 之前处理的事件目录）：`sacprocess/`（= 新 `process_flow/0_data_processing/`）、`1_station_distribution/{maps_info,profiles}/`（新规范平铺）、`4_fk_analysis/`（= 新 `3_fk_analysis/`）。用户点名再逐个改；本 skill 的 references 里出现的旧路径描述的是**历史事件的实际路径**，不是新规范。

> ✅ **已迁移**：马达加斯加 `2013-08-17-mw62-southwest-indian-ridge` 已于 2026-09-28 完成迁移（0/1/2/3 号 + 1 号平铺，89 个脚本/文档共 117 处路径改写；证据与回滚清单在该事件 `process_flow/_migration_2026-09-28/`）。本地区其余事件（2013-07-22 等）仍是旧名。

> ✅ **脚本已按"自动探测"改造**（新事件零改动即用新名，旧事件仍读旧名）：`process_image.py`（步0 目录）、`run_single_station.py`（渲染模板时把 `sacprocess/` 改写成 `process_flow/0_data_processing/`）、`plot_record_sections_abs_land.py`（图写脚本同目录、FK 目录自动探测 3/4 号）。其余事件目录内的旧脚本副本（各事件自己的 `profiles/*.py`、`sacprocess` 常量）**仍写旧路径**，只对旧事件有效。

## 表格格式规范

`step1_station_table_and_map.py` 输出的 `stations.xlsx` 格式（2026-08-03起）：
- 字体：微软雅黑
- 标题行：16号加粗居中
- 表头（双行中文+英文）：14号加粗居中
  - 中文行：黑字无底
  - **英文行：蓝底白字**（`fill='4472C4'`, `color='FFFFFF'`）
- 数据行：14号居中
- 全框细黑线
- SNR 列位置：`stations.csv` 第 14 列 `snr_u`（统一噪声基准 T 分量）；`stations.xlsx` 主列在第 9 列，浅绿色底色标记 ≥ 阈值者
- **SNR 口径 = 统一噪声基准**（2026-10-14 起唯一口径）：20log10(信号窗 RMS / 阵列噪声中位数)，2–8 Hz 整段带通后切窗 —— 完整参数、两个坑（整体滤波、噪声窗避开 P 波）、下游接法见 `references/snr-unified-noise-baseline.md`

## T波走时计算 ⚠️

T 波在海洋 SOFAR 通道中传播，走时应基于**实际的海洋/陆地分段路径**而非震中距。

**❌ 错误做法**：`dist_km * 0.95 / 1.48`（用震中距乘粗比例，不准确）

**✅ 正确做法**：分段计算，注意路径方向取决于震源位置！
```
海洋地震: T波走时 = 震源→海岸(海洋SOFAR, 1.48 km/s) + 海岸→台站(陆地P波, ~6 km/s)
陆地地震: T波走时 = 震源→海岸(陆地P波, ~6 km/s) + 海岸→台站(海洋SOFAR, 1.48 km/s)
```

```python
from obspy.geodetics import gps2dist_azimuth

# 以海洋地震为例 (Prince Edward Islands → BOSA):
coast_la, coast_lo = -34.5, 22.0    # 南非南海岸
ocean_km = gps2dist_azimuth(evla, evlo, coast_la, coast_lo)[0]/1000  # 震源→海岸 (海洋)
land_km = gps2dist_azimuth(coast_la, coast_lo, stla, stlo)[0]/1000   # 海岸→台站 (陆地)
t_tp = ocean_km / 1.48 + land_km / 6.0   # T-P 走时 (T波→陆地P波转换)
t_ts = ocean_km / 1.48 + land_km / 3.6   # T-S 走时 (T波→陆地S波转换)
```

### 🔴 画图前必须验证 T 波到时！

不要盲信理论走时，画图前先检查数据里 T 波能量实际在哪：
```python
# 步骤1: 未滤波数据找主能量峰
pk = np.argmax(np.abs(data[600:2000])); t_pk = 600 + pk * delta
# 步骤2: 对比各频段——真实 T 波应在各频段均有能量
```

⚠️ **窄带通假峰陷阱（BOSA 案例）**：T 波能量集中 0.5-2 Hz 低频（SOFAR 通道低频衰减慢），只做 2-8 Hz 带通会漏掉远处台站的 T 波；反过来，窄带通滤波本身可能制造"假峰值"（滤波器带内噪声被放大成貌似信号的峰）。**验证方法**：对比 0.5-2 / 2-4 / 4-8 Hz 多频段，真实 T 波应在各频段均有能量且到时一致；只在一个频段出现的峰多半是滤波假象。


## 🔴 数据起点判断（所有时间轴问题的一条总规则）

**数据下载常从发震前 ~10 分钟开始**（用户习惯），因此波形数据可能比发震早 600s。所有时间轴问题的根源都是没先确认数据起点。处理数据前第一步：

```python
sac_b = getattr(tr.stats.sac, 'b', 0)
sac_o = getattr(tr.stats.sac, 'o', None)   # 读原始字段，不要只看 ObsPy starttime
if sac_o is not None and 0 < abs(sac_o) < 1e6:
    t_off = sac_b - sac_o          # o=600 → t_off=-600，数据从发震前600s开始 ✓（2013-08-17 实测）
else:
    EV_ORIGIN = UTCDateTime(年,月,日,时,分,秒)   # 从文件名 D.2019.022.190143 解析
    t_off = tr.stats.starttime - EV_ORIGIN       # 缺 o 时 fallback（2019/2023 南非）
t = t_off + np.arange(npts) * delta
```

- **`o` 有效（0<|o|<1e6）时优先 `b-o`**，它直接给出数据起点相对发震的秒数（负 = 发震前开始）
- **`o` 缺失时才用 `starttime - EV_ORIGIN`**
- ⚠️ **禁止仅凭 `starttime == EV_ORIGIN` 判断数据从发震开始**：ObsPy starttime 读的是 SAC reference time，该字段可能被写成发震时刻而实际数据从发震前 600s 开始
- 自检：`t[0]` 应为负值（发震前开始，常见 -600）或 0；若 SAC 有 o 但 t[0]=0 → 坐标系错了。窗口内应有纺锤形能量包络，若只见平坦噪声 → 先查 t_off

完整案例复盘见 `single-station-analysis/references/sac-time-axis-b-o.md`。

## 预处理

### ✅ process_image.py（2026-07-31 标准）

**所有新事件和重跑必须使用 `process_image.py`**。遵循 `plotting-standards` 规范：

- 单脚本输出：预处理（去响应 SAC） + 标准图（raw / 去响应 / 带通+包络 / 频谱 / packet 可选）
- 内建 `get_t0_offset()` 处理 SAC `o` 缺失 fallback
- 从 `stations.csv` 读取到时，不自己计算
- 按台站分子文件夹：`waveform_image/<STA>/<STA>_{1_raw,...}.png`
- 复用：复制到新事件目录，改顶部 `EVENT`/`BASE`/`UTCDateTime` 三处即可

**2026-09-15 更新（模板已同步）**：`templates/process_image.py` 换成新版（564 行，旧版留档
`templates/old/process_image_OLD_363lines_20260817.py`），新增开关：
`--only` / `--figs` / `--skip-resp` / `--packet` / `--xmax` / **`--square`** / `--env-smooth` /
`--env-mode mean|max|lp` / `--out-suffix` / `--tt-paper|--tt-both` / `--archive`。
其中 **`--square` 现在作用于五类图**（raw / decon / bandpass / packet / spectrogram；
`plot_decon` 与 `plot_bandpass` 原先 figsize 写死，已改为随 `FIG_SQUARE` 切 `(FIG_W, FIG_W)`）。
典型 1:1 出图：

```bash
python process_image.py --only SA16,SA17 --figs 1_raw,packet,2_remove_response,bandpass_2-8Hz,3_spectrogram \
       --packet --square --xmax 2000 --skip-resp
```

```bash
# 从已有事件复制模板
cp <existing_event>/process_image.py <new_event>/
# 改事件路径和发震时刻（3处）
sed -i 's|<old_event_name>|<new_event_name>|g' process_image.py
sed -i 's|UTCDateTime(2020, 8, 6, 23, 36, 44)|UTCDateTime(NEW)|' process_image.py
```

### ❌ process_obs.py（已废弃）

旧模板已从所有事件目录删除，不再使用。如需查看历史代码见 git 历史。

### 📦 代码库 (`codex/lib/processing.py`)

可复用处理库，提供标准 `parse_sacpz` / `remove_response` / `apply_bandpass` / `process_single_file` 函数。所有函数使用 **`tr.simulate(paz_remove=paz_dict)`** 方案，与 `process_image.py` 一致（2026-08-01 已对齐）。

```python
import sys; sys.path.insert(0, '/home/sy/codex/lib')
from processing import remove_response, apply_bandpass, parse_sacpz
```

### 📜 templates/ 3 脚本

| 脚本 | 步骤 | 说明 |
|------|------|------|
| `process_image.py` | 2. 预处理 | 去响应 + 标准图（raw/去响应/带通+包络/频谱/packet），`--square` 可出 1:1（2026-09-15 新版，564 行） |
| `process_single_station.py` | 3. 单台分析 | `--plot` 出 bandpass+spectrogram 并排图 |
| `waveform_proflie.py` | 剖面 | 独立波形剖面：land_km 分箱(10km)取SNR最高、分区域、0.5-2Hz低频、包络独立P98（见下「波形剖面」节） |

> ⚠️ `profile_conversion.py` 已废弃 — 剖面图功能已并入 `step1_station_table_and_map.py`（`_gen_waveform_profile()` 函数）。步1一键输出5件套。


### ⚠️ Python 环境

**统一使用 `/home/sy/miniforge3/envs/obspy/bin/python`**（pygmt 0.19.0 + obspy 1.5 + numpy/scipy/matplotlib/pandas/xarray/openpyxl/pillow 齐全）。
系统 `python3` (3.14.4) 的用户级 numpy/scipy/matplotlib 不完整或损坏。

```bash
# 直接调用（不激活）
/home/sy/miniforge3/envs/obspy/bin/python script.py
```

备用环境（勿用于绘图，pygmt 0.18 subplot 有 bug）：
- `~/codex/venv/bin/python` — 通用科学计算
- `/home/sy/envs/tphase_env/bin/python3` — T 波分析（含 future 包）

### ⚠️ ObsPy 1.5.0 API 变更

`tr.remove_response(paz_remove=...)` 在 ObsPy ≥1.5.0 不生效（报 "No response information found"）。**必须改用**：

```python
# ✅ 唯一推荐方案: tr.simulate() — BOSA 验证方案 (2026-08-01)
tr.simulate(paz_remove=paz, remove_sensitivity=True,
            water_level=water_level, pre_filt=pre_filt)
# 直接输出速度 (m/s)，无需手动 np.gradient()
```

**❌ simulate_seismometer + np.gradient 已废弃**（2026-08-01）。该方案需手动 `np.gradient(data_out, dt)` 将位移转为速度，引入额外误差和不一致。所有 skill/codex 脚本已统一为 `tr.simulate()` 方案。

### 🔴 单位决策：T 波统一用速度 (m/s)（2026-08-24 定）

**全流程（去响应 → 滤波 → 偏振 → 衰减 → 归一化）统一用速度 (m/s)**。原因：
- 去仪器响应可输出三种物理量：位移 DIS/DEF (m)、速度 VEL (m/s)、加速度 ACC (m/s²)，三者差 1/ω 因子（位移=速度/ω=加速度/ω²）
- 位移低频放大（1/ω），加速度高频放大，**速度中频信噪比最好**
- 地震计（STS-2 等）本身就是速度计（SACPZ 灵敏度单位 M/S），去响应到速度反褶积误差最小
- 2–8 Hz 内位移 vs 速度有 4 倍相对差异（2 Hz 处位移=速度/12.6，8 Hz 处=速度/50.3）→ 位移会扭曲子频带（2–3/3–4/4–5/5–6/6–8 Hz）相对振幅
- 水听器/OBS 声压类比时才用 DEF（对水听器=声压）

**⚠️ t_phase_analysis 包单位陷阱**：`t_phase_analysis/response.py` 输出 `output="DEF"`（位移 m，函数名 `remove_response_def`），但 `waveform.py` 的 `plot_waveform_panel` ylabel 硬编码 `Velocity (m/s)` —— **标签与数据不符**（数据是位移却标速度）。用该包处理时注意：去响应结果是位移；统一决策是速度，包待改 VEL（2026-08-24 未改，用户选择先改当前事件+skill）。

**图标签规范**：P98 归一化后的波形图纵轴必须标 `Norm. Amp.`（无量纲），**不能标 `Velocity (m/s)`**（数据已归一化，标签会误导）。只有显示真实速度（未归一化）的图才能标 `Velocity (m/s)`。

### 数据下载

首选 HTTP POST 批量下载，失败台站用 ObsPy fallback。按网络聚合查询避免 300s+ 超时。

### StationXML 响应文件（IRIS FDSN 2025年后标准格式）

IRIS FDSN 下载的数据使用 StationXML (`.xml`) 代替传统 SACPZ 文本文件。ObsPy 原生支持。

**去响应流程（与 SACPZ 文本版不同）**：
```python
from obspy import read_inventory

inv = read_inventory(f"{xml_dir}/IU.AFI.xml")        # 读取 XML
tr_out = tr.copy()
tr_out.detrend('demean'); tr_out.detrend('linear')
tr_out.remove_response(inventory=inv, output='VEL',   # ✅ 直接出速度
                        pre_filt=(0.001,0.005,45,50), water_level=60)
```

**⚠️ 与 SACPZ 文本版的差异**：
- 不用 `parse_sacpz()` + `tr.simulate(paz_remove=...)` — PAZ dict 不适用
- 用 `tr.remove_response(inventory=inv, output='VEL')` — ObsPy 直接读 XML 响应链
- 台站坐标从 `inv[0][0].latitude/longitude` 获取，不依赖 SAC 头
- 通道名可能有 `BH1`/`BH2` (非标准) → 映射 BH1=BHN, BH2=BHE
- 部分台站只有 BHZ 无水平分量 → 跳过多面板出图，仅输出去响应 SAC

**典型适配流程**（见 Japan 2025-08-01 案例）：
1. SACPZ 目录 = `raw/SACPZ/IU.XXX.xml`（.xml 扩展名）
2. `load_station_meta()`: 循环读所有 XML → `{net.sta: {lat, lon, inventory}}`
3. 从 meta dict 取台站坐标（SAC 头常缺 stla/stlo）
4. 检测三分量完整性：`has_3c = all(c in comps for c in ['BHZ','BHN','BHE'])`

### 数据整理前检查

**在处理前，先统计各台站的文件数量，识别异常台站：**

```bash
cd <event>
python3 -c "
import glob, os; from collections import Counter
c = Counter()
for f in sorted(glob.glob('*.SAC')):
    key = f'{f.split(\".\")[0]}.{f.split(\".\")[1]}'
    c[key] += 1
for k, v in sorted(c.items(), key=lambda x: -x[1]):
    tag = '⚠️  ' if v > 3 else '✅' if v == 3 else '⚠️  '
    print(f'  {tag} {k}: {v} files')
print(f'Normal (3-comp): {sum(1 for v in c.values() if v==3)} / {len(c)} stations')
"
```

**🔴 务必扫描非标准通道名（BH1/BH2/BHX/BHY）**，确认 SAC 和 SACPZ 都已重命名。

⚠️ **YV 网络通常全部使用 BH1/BH2**，是最容易漏改的区域台网。扫描范围必须覆盖**所有网络**（II/IU/YV/PS/AV），不只常见的 II/IU。

```bash
# 检查残留的非标准通道名（所有网络）
ls raw/*BH1*.SAC raw/*BH2*.SAC raw/*BHX*.SAC raw/*BHY*.SAC 2>/dev/null
ls raw/SACPZ/*BH1 raw/SACPZ/*BH2 2>/dev/null
# 若有输出 → 必须重命名，否则对应台站会被跳过

# 重命名后验证清零
ls raw/*BH1* raw/*BH2* raw/SACPZ/*BH1 raw/SACPZ/*BH2 2>/dev/null
# 应无输出
```

### 常见异常类型

| 异常类型 | 文件数 | 原因 | 处理方式 |
|---------|--------|------|---------|
| **多 location code** | 6–12 | 同台站有 00/10/20/60 等多个 location code（IU/II/JP 网络常见） | 默认按台站聚合，会覆盖同一 channel，需手动选择 1 个 location |
| **碎片文件** | 数十~上百 | 数据获取错误导致拆分为极短片段（如 MY.KSM 122 文件） | 直接跳过该台站 |
| **多 segment** | 6 | 同一台站有 2 个 location code，每个 location 3 分量 | 同多 location 处理 |
| **BH1/BH2 命名** | 3 | 非标准水平分量命名（II: ABPO/SUR, IU: LSZ/MIDW/TARA, YV: 全部, AV 网络） | SAC 和 SACPZ 都必须重命名 BH1→BHN、BH2→BHE。只改 SAC 漏 SACPZ → 台站被跳过 |
| **BHX/BHY 命名** | 3 | 非标准水平分量命名（部分 PS 网络，如 PS.CHUK） | 自动映射 BHX→BHN、BHY→BHE |
| **仅 BHZ 台站** | 1-3 | 缺少水平分量（如 MI.DPS, PS.CHUK 等） | 仅输出去响应 SAC，跳过四面板出图 |

### ⚠️ 非标准通道名处理 (BH1/BH2, BHX/BHY)

部分网络使用非标准水平分量命名。**在预处理前将 SAC 和 SACPZ 文件都重命名为标准名**。

🔴 **关键：SAC 和 SACPZ 必须同步重命名！** 只改 SAC 不改 SACPZ → `find_sacpz` 按 BHN/BHE 搜索找不到 BH1/BH2 命名的 SACPZ → 台站缺 SACPZ 被跳过。

**受影响的网络**: II (ABPO, SUR), IU (LSZ, MIDW, TARA), YV (全部), PS (CHUK)

```bash
cd raw/
# SAC 文件: BH1→BHN, BH2→BHE（所有网络）
for f in *BH1*.SAC; do mv "$f" "${f/BH1/BHN}"; done
for f in *BH2*.SAC; do mv "$f" "${f/BH2/BHE}"; done
# SACPZ 文件: 同样重命名！
cd SACPZ/
for f in *BH1; do mv "$f" "${f/BH1/BHN}"; done
for f in *BH2; do mv "$f" "${f/BH2/BHE}"; done
# PS 台站: BHX→BHN, BHY→BHE
for f in *BHX*.SAC; do mv "$f" "${f/BHX/BHN}"; done
for f in *BHY*.SAC; do mv "$f" "${f/BHY/BHE}"; done
```

`find_sac_triplet()` 需包含 fallback 映射（见代码库）。

### 执行预处理

`process_image.py` 一步完成预处理 + 四张标准图（遵循 `plotting-standards` 规范）。

- 按台站分子文件夹输出：`waveform_image/<STA>/<STA>_1_raw.png` 等
- 从 `stations.csv` 读取到时，内建 `get_t0_offset()` 处理时间对齐

### ⚠️ 常见陷阱

**发震时刻**：SAC 头 `nz*` 是数据开始时间，`sac.o` 常缺失。从文件名解析 `HHMMSS`：

SAC 头的 `nz*`（nzyear/nzjday/nzhour/nzmin/nzsec）是**数据开始时间**，不是发震时刻。`sac.o` 常缺失。
**正确做法**：从文件名解析 `HHMMSS`（如 `...D.2019.022.190143.SAC`）。
```python
parts = os.path.basename(fn).split('.')
time_str = parts[-2].replace('_proc','')
origin = UTCDateTime(year=int(parts[-4]), julday=int(parts[-3]),
                     hour=int(time_str[0:2]), minute=int(time_str[2:4]),
                     second=int(time_str[4:6]))
```

### ⚠️ SACPZ 文件 location code

不同台网的 SACPZ 文件 location code 可能不同：GT/BX 用 `00`，AF/ZT 用 `--`。读取时必须尝试多种：
```python
for loc in ['--', '00', '']:
    fp = f'SACPZ.{net}.{sta}.{loc}.{ch}'
    if os.path.exists(fp): ...
```

### 🔧 补下载缺失 SACPZ（IRIS FDSN）

当部分台站缺失 SACPZ 时，可尝试从 IRIS FDSN 下载：

```python
from obspy.clients.fdsn import Client
from obspy import UTCDateTime

client = Client("IRIS")  # 现在 redirect 到 EarthScope
t = UTCDateTime("2013-08-17T16:32:31")
inv = client.get_stations(network="YV", station="EURO,GLOR,...",
                          location="*", channel="BH?",
                          starttime=t-3600, endtime=t+7200,
                          level="response")
# 写入 SACPZ
for net in inv:
    for sta in net:
        for ch in sta:
            fname = f"SACPZ.{net.code}.{sta.code}.{ch.location_code}.{ch.code}"
            inv.write(fname, format="SACPZ")
```

> ⚠️ 老旧区域台网（如 YV）可能已从 IRIS 下线，返回 HTTP 204。这种情况只能跳过。

> ⚠️ **2023年及以后的波形数据**：IRIS FDSN station inventory 可能返回台站元数据，但 `get_waveforms()` 对 2023+ 的非洲区域台网（AF/BX/NR/GE）几乎全部返回 HTTP 204 "no data"。这是 IRIS 数据分发策略问题（可能涉及 embargo 或数据未提交），不是脚本 bug。**先尝试下载，若全部 204 则用已有数据直接处理**，不要反复尝试不同查询方式浪费时间。

### ⚠️ Z/R/T 分量（2026-07-30 更新）

**四张图全部使用 Z/R/T**（NE→RT 旋转后）：

- `r = n_ * np.cos(br) + e * np.sin(br)`
- `t_ = -n_ * np.sin(br) + e * np.cos(br)`
- `br = np.radians(baz)`，`baz` 从 `stations.csv` 取
- 标签：`Z` / `R` / `T`（不写 BH 前缀）

### ⚠️ 去响应质量检查

部分台站的 SACPZ 导致去响应后产生 10^12~10^18 的数值垃圾。预处理后必须检查：
```python
p98 = np.percentile(np.abs(data), 98)
if p98 > 1e6 or not np.isfinite(p98):
    # 跳过此台站
```

**🔴 预处理后必须做跨台站振幅一致性检查**（2026-08-03 FOMA 案例）：
同一事件的 `_resp.sac` 各站 P98 应在同一数量级（速度 ~1e-5 m/s）。若某站差几个数量级 → 去响应失败：
- 太大（如 FOMA P98=6e+27，比 raw counts 还大 24 个数量级）→ PAZ 参数解析错误或 simulate 反了
- 太小（如 TEZI P98=2e-10）→ 灵敏度/增益单位错误
```python
p98s = []
for f in glob.glob(f'{resp_dir}/*BHZ*_resp.sac'):
    tr = read(f)[0]
    p98s.append((os.path.basename(f).split('.')[1], np.percentile(np.abs(tr.data), 98)))
p98s.sort(key=lambda x: x[1])
# 正常: 各站 P98 都在 ~1e-5 m/s 量级；跨度 >10^4 说明有站去响应失败
```
用标准模板 `process_image.py` 重跑即可修复（旧版脚本绘图部分与模板相同，仅去响应计算可能不一致）。

**🔴 SACPZ 重复块解析 bug（FOMA 振幅爆炸根因）**：部分 SACPZ 文件（如 G.FOMA）把同一个响应块重复写 3 遍（3× `ZEROS 7` + `POLES 11`）。手工 `parse_sacpz` 若把所有块都读进去 → 21 zeros + 33 poles + 错误的 a0/sens → simulate 后振幅爆炸（P98=6e+27）。**修复：只解析第一个 ZEROS/POLES 块**，遇到第二个 `ZEROS` 即 `break`。验证：FOMA 修复后 zeros=7, poles=11, P98=2.4e-05 与其他站一致。
```python
def parse_sacpz(fp):
    zeros, poles = [], []
    const = 1.0; a0 = None; sens = None; iz = ip = False
    seen_block = False
    with open(fp) as f:
        for line in f:
            s = line.strip()
            if s.startswith('*'):
                if 'SENSITIVITY' in s.upper():
                    try: sens = float(re.findall(r'[-+]?\d+\.?\d*', s.split(':')[1])[0])
                    except: pass
                if s.upper().lstrip('* ').startswith('A0'):
                    try: a0 = float(s.split('=')[1].split()[0])
                    except: pass
                continue
            if s.startswith('ZEROS'):
                if seen_block: break  # 只取第一个块
                iz = True; ip = False; seen_block = True; continue
            if s.startswith('POLES'): ip = True; iz = False; continue
            if s.startswith('CONSTANT'):
                const = float(s.split()[1]); iz = ip = False; continue
            if iz and s:
                p = s.split()
                if len(p) >= 2: zeros.append(complex(float(p[0]), float(p[1])))
            elif ip and s:
                p = s.split()
                if len(p) >= 2: poles.append(complex(float(p[0]), float(p[1])))
    if a0 and sens: pass
    elif const and sens: a0 = const / sens
    else: a0 = const; sens = 1.0
    return zeros, poles, a0, sens
```
调试手段：直接打印 `parse_sacpz()` 返回的 `len(zeros)/len(poles)/a0/sens`，正常 3-5 zeros / 4-9 poles；若 zeros>10 或 a0>1e15 就是重复块或元数据异常。

### ⚠️ N/E 分量长度不一致

旋转前对齐：`npts = min(len(n_data), len(e_data))`

### 🔴 plot_bandpass 纵轴统一（AMPY 案例 2026-08-03 / Madagascar 2026-08-24）

`plot_bandpass` 里**必须显式设 ylim**，且三个子图共用同一个范围，否则 matplotlib 每个子图自动缩放 → Z/R/T 纵轴不一致（Z 振幅大时被压扁）。

**✅ 新标准（2026-08-24，取代固定 ±2.5）**：Z/R/T 共用全局波形 P98 归一化后，纵轴按数据自动动态设置（不固定 ±2.5，避免截断视觉假象）：
```python
all_wave = np.concatenate([np.abs(filt[l][i1:i2]) for l, _ in comps])
ep98 = np.percentile(all_wave, 98) or 1.0
y_max = 1.05 * max(max(np.max(np.abs(filt[l][i1:i2]/ep98)) for l, _ in comps),
                   max(np.max(env[l][i1:i2]/ep98) for l, _ in comps), 1.0)
ax.set_xlim(t_min, t_max); ax.set_ylim(-y_max, y_max)
```
- bandpass 归一化因子用**波形全局 P98**（不是包络 P98），波形和包络共用
- 纵轴标签：`Norm. Amp.`（归一化后不得标 Velocity (m/s)）
- 检查是否所有子图都设了 ylim（旧模板固定 ±2.5 的副本需同步更新）

---

## 台站筛选与下载（新事件启动）— 2026-08-01

### 🔴 FDSN 时区陷阱

用户报告的事件时间通常是**当地时间**（如北京时间 UTC+8）。FDSN 查询必须用 **UTC**。

- 北京时间 10:48 → UTC 02:48
- ❌ 直接用用户给的时间查 FDSN → 全部返回 204 "no data"
- ✅ 先换算 UTC，再查

### 台站筛选流程

```python
# 1. 获取台站清单（按方位角和距离过滤）
inv = client.get_stations(latitude=evla, longitude=evlo,
    minradius=0, maxradius=52, channel='BHZ', level='station',
    starttime=t1, endtime=t2)

coords = {}
for net in inv:
    for sta in net:
        dist_m, az, baz = gps2dist_azimuth(evla, evlo, sta.latitude, sta.longitude)
        dd = kilometers2degrees(dist_m/1000)
        if az_min <= az <= az_max and dd <= dist_max:
            coords[f'{net.code}.{sta.code}'] = {
                'lat': sta.latitude, 'lon': sta.longitude,
                'dist_deg': dd, 'az': az
            }
```

### 批量下载（按网络聚合）

逐个台站查询 FDSN → 91 站需 300s+ 超时。**按网络批量下载**：

```python
from collections import defaultdict
by_net = defaultdict(list)
for code in coords:
    by_net[code.split('.')[0]].append(code.split('.')[1])

for net, stas in by_net.items():
    st = client.get_waveforms(net, ','.join(stas), '*', 'BHZ', t1, t2)
```

### SNR 计算 + 方位角分选

1. 计算每站 SNR（噪声窗前 60s / 信号窗 P 波后）
2. 按方位角分组（每 3-5° 一组）
3. 每组选 SNR 最高的台站
4. 若全组 SNR < 2，选距离最近的

结果写入 `selected_stations.txt`。

### ⚠️ 台站坐标：从 inventory 获取

SAC 头常缺 `stla`/`stlo`。**始终从 FDSN station inventory API 获取坐标**，不依赖 SAC 头。

```python
# ❌ 错误：依赖 SAC 头
lat = tr.stats.sac.stla  # 可能为 0 或不存在

# ✅ 正确：从 inventory 获取
inv = client.get_stations(network=net, station=sta, level='station', ...)
lat = inv[0][0].latitude
```

## 台站分布

**期刊级地图**由 `/home/sy/codex/scripts/plot_station_map_final.py` 生成（`figure_map_final.png`）：

- **15s 地形**（约450m分辨率，GMT 缓存 583MB） + geo 色标（底部水平色标，`JBC+w10c/0.3c+h+o0c/1.6c`）
- **备选**：`@earth_relief_01m`（1.8km，global coverage）
- **GMT 缓存路径**：`/home/sy/.gmt/server/earth/earth_relief/earth_relief_15s_p/`（77 tiles）
- **注**：也缓存了 `earth_synbath_15s`（12MB，海洋高精），见 `plotting-standards/references/gmt-global-topography.md`
- **SNR 着色台站** + viridis 色标（右侧垂直色标，每10一标 `xa10f1`）
- **双色标**各自独立（elev_cpt 单独存文件避免被 SNR cpt 覆盖）
- **图例框**（右下白框黑边）：Station (colored by SNR) / Epicenter / -1000 m depth
- **-1000m 等深线** + 文字标注
- **智能标注**：land 最近 + land 最远 + SNR 最高（3个，白描边 + 虚线引线 + 方位角射线 + 角度标注）
- **标签防重叠**（2026-08-03）：每个台站从 4 个候选方向（优先远离震中）中选第一个不与已标注台站重叠的方向；`fixed_offsets` dict 可强制指定台站的固定偏移（如 `{"AMPY": (-0.8, -0.5)}` 左下角）
- **EXCLUDE_BOXES**（2026-08-03）：`[[lon_min, lon_max, lat_min, lat_max], ...]` 矩形区域内的台站不画（用于剔除噪声密集簇，如 10°S 北部）。step1 的 `DEFAULT_EXCLUDE` 同名数组也会过滤，两处需同步
- **震中信息**白描边（8方向白色偏移），14p Helvetica-Bold
- **坐标轴**：`xa10f5`/`ya10f5`（10° 标 + 5° 小刻度），图框 `WSne`（用 `fig.basemap(frame=...)` 单独调用，coast 里设 frame 可能被覆盖）
- **台站符号**：0.5c 三角形，0.3p 黑边

全站标注版（`figure_map.png`）由 `plot_station_map.py` 生成。

## 波形剖面

**实现**：`templates/waveform_proflie.py`（独立脚本，推荐新事件用，2026-08-04 最新版）。

### 核心参数（当前定版，用户还会改）

| 参数 | 值 | 说明 |
|------|-----|------|
| 分箱 | `--bin-km 10`（默认） | 按 land_km 陆地距离分箱，档内**只取 SNR 最高台站** |
| 频段 | **`--band 0.5-2.0`（默认低频）** | **T 波低频成分衰减慢，远处台站（480km+）仍可见；2-8Hz 高频几百公里就耗尽。** 近处分析可切 `--band 2.0-8.0`。这是"包络画成直线"的真正根因（2026-08-04 实测：远处台站 0.5-2Hz 峰值比 2-8Hz 高 5-25 倍） |
| 平滑 | **3s** | Hilbert 包络 `obs_env` → `uni_filt(size=int(3/delta))`。8s 会把峰值抹平、纺锤形变钝 |
| 时间窗 | `t_min = min(tt_tp) - 30, t_max = max(tt_ts) + 60` | 收紧窗口，信号充满图面 |
| 归一化 | **包络独立 global_env_p98** | 每条 trace 先扣噪声基线（t < t_min 中位数），再基于 `env_c` 算全局 env P98。**不能复用波形 global_p98**（包络是恒正平滑量，尺度不同） |
| 主图 | 灰色填充包络(alpha=0.45) + 黑色细波形叠加(lw=0.25) | 见 `plotting-standards` → 剖面图 |
| amp | `amp = min(30, max(5, bin_km*0.45))`；envelope `amp*3.0`；per-bin `amp*1.6` | |

### 输出三张图
1. `profile_env_{区域}_land.png` — 主图：包络（全局 env P98）+ 波形叠加
2. `profile_env_perbin_{区域}_land.png` — 辅助：包络 per-bin P98（远处弱信号也能看，但噪声会被放大）
3. `profile_waveform_{区域}_land.png` — 波形版

### 陆地走时对齐
横轴 = `t_land = t_abs - tt_ocean`，`tt_ocean = ocean_km / 1.48`（SOFAR），**t=0 = 波上岸**。所有到时标记同步减 `tt_ocean`。xlabel：`Time from shore arrival (s)`。

### 排除台站
`PROFILE_EXCLUDE = {'BITY', 'Z01TG', 'W12MB'}`（脚本顶部）：异常噪声/脉冲站，发现方法——检查各档 T 波时段 P98，某档 SNR 低但振幅异常大即脉冲噪声站。

### 用法
```bash
python waveform_proflie.py <event> --snr-min 3.0 --bin-km 10 [--band 0.5-2.0]
```

> **🔴 全局 P98 归一化被异常强台站绑架**（2026-08-04 南非2019实测）：单个振幅异常台站会主导 P98 把其余台站的 T 波包络压成"平直噪声"假象（CRAIL 包络峰比其他台大 79 倍、VOET ×25）。处置：算各台站 T 波窗口包络峰相对最弱台的比值（>10 倍异常）→ 把异常台加入 PROFILE_EXCLUDE 重画。自检：图上大量线平直但 SNR 不低 → 先怀疑被强台站压扁，不是信号弱。

### 记录剖面（Fig.3 型 3×2：绝对计时 + 陆地计时）— 属本步，不是独立步骤

脚本：`process_flow/1_station_distribution/plot_record_sections_abs_land.py` → `record_section_abs_land_3x2.png`（🔴 步1 平铺，图与脚本同目录；旧事件仍在 `1_station_distribution/profiles/` 下）。
上排 = 震中距 vs 时间（0 = 发震）；下排 = 陆地段距离 vs 时间（0 = 波上岸，`t − ocean_km/1.48`）；列 = Z/R/T；只画 2–8 Hz 的 6 s 平滑包络（黑色 `lw=0.25`、无灰波形、无填充），**逐台峰值归一**到绝对 km 常量；到时 T-P 红 / T-S 蓝短竖线（长度 ≈0.5×中位间距）。版式判据、留白公式、剔台名单详见 `record-section-plotting` §14。

🔴 **横轴窗口默认 = 上排 `0–1800 s`（相对发震）/ 下排 `−1000…+500 s`（0 = 上岸）**（用户 2026-09-28 定：「默认就先这个窗口，我要你改的时候会说的」）——**不许自行更改、不经商量截窄或留白**；用户点名要改窗口时才动，命令行 `--xlim-abs A,B` / `--xlim-land C,D` 可覆盖（不留参数即用默认）。

```bash
python process_flow/1_station_distribution/plot_record_sections_abs_land.py              # 默认窗口（六面板一张）
python process_flow/1_station_distribution/plot_record_sections_abs_land.py split|abs|land   # 拆成两张 3x1 / 只出一张
python process_flow/1_station_distribution/plot_record_sections_abs_land.py --xlim-land -100,400   # 用户要求改窗口时
```

出图后脚本**自动 copy 一份到 F-K 步骤目录** `process_flow/3_fk_analysis/`（旧事件 `4_fk_analysis/`；`--fk-dir none` 可关闭。F-K 用这张图圈定分析时窗 / 选子阵）。**图和代码一起走**：同步时把脚本本身也 copy 到同路径（用户 2026-09-29），1 号根本来就有一份脚本。

### 🔴 上下两排跨度统一 + 单一幅度/基准（2026-09-29 用户定）

「下排只是上排同一批波形换排序 + 换计时，不另设一套参数」——★参数区三个开关：

| 开关 | 默认 | 含义 |
|---|---|---|
| `UNIFY_SPAN` + `SPAN_S/SPAN_KM` | `True`, 1800 s / 1200 km | 两排 x 都用 SPAN_S 秒、y 都用 SPAN_KM km（各自按本排墨迹范围居中）→ 两排同 s/px、同 km/px，可直接对读 |
| `AMP_MODE` / `AMP_RATIO` | `'span_ratio'` / 1.0 | **默认**：上排 amp = AMP_RATIO × dist 中位间距；**下排 amp = 上排 amp × (下排 y 跨度 / 上排 y 跨度)** → 两排「amp / 轴跨度」同比例 ⇒ 同一条道在两排**屏幕高度一致**（对应得上）。`'ratio'`=各排按本排间距；`'fixed'`=绝对 km（`AMP_KM` / `AMP_PER_ROW`） |
| `Y_LIM_LAND` / `Y_MARGIN_BOT,TOP` / `Y_TICK` | `(0,800)` / `0.9,0.5` / `200` | 下排 y 轴由用户直接给（0 起）；上排留白 = 系数 × amp（底/顶，防贴底与切顶）；y 刻度步长两排共用 |
| `NORM_WINDOW` | `(0, 1800)` | 逐台峰值归一的**基准窗**（绝对时间，两排共用）；`None`=旧行为（各道取自己可见窗口的峰，会让远台基准变成 T 包、近台变成 P/S，包高在台间不可比） |

⚠️ 跨度统一后：`(0-data_min, data_max+amp)` 居中放进 1200 km → 因为下排陆地段距离天生只有 ~670 km，下排上下各留 ~266 km 空白；上排 ~1114 km 基本占满。

## 单台分析（步2）

只分析**台站分布图上标注的那 3 个台站**（最近陆 / 最远陆 / SNR 最高）——与步2 的默认台站选择一致。
流程见 `seismic/single-station-analysis` skill；输出落 `process_flow/2_single_station_analysis/{STA}/`。

## 衰减分析（步4，已规范化 → `seismic/twave-crustal-attenuation`）

2026-10-23 起本步有独立 skill：口径（尾波归一化 CNM、`2S` 逐台窗、rms、`sqrt(r*70)`）、
三层判据（`A_d/A_c>=1` / `|k|/sigma_k>=1` / `|alpha|/sigma>=2`）、核心引擎与四类出图脚本、
定案文档模板都在 `seismic/twave-crustal-attenuation`（含 `templates/` 可直接跑）。
本节只保留下面这条**数据层面的限制**，仍有效。

> ⚠️ **非洲南部 Q 反演限制**：实测多事件堆叠（65对）仍无法可靠约束 Q——ZT/BX/AF/GT 四网混合仪器增益不一致，单事件震源项无法消除。Zhou 2021 用统一 XA 台阵+尾波归一化才成功。建议将 Q 反演作为次要结果或不包含。

## 台阵反演（PMCC/FK）— 流程步3 = F-K（反演反方位角）

> 🔴 **2026-09-28 更新：F-K 是流程第 3 步**（目录 `process_flow/3_fk_analysis/`），**引擎实现、参数与判据全部以 `seismic/fk-array-analysis` skill 为准**（默认 method=coarray，Octave/SEIZMO fk4d）。
> 下面这段保留的是**「为什么不直接用 PMCC/水听器那套」的结论**，仍有效：稀疏陆地宽频台网做不了 PMCC，别拿它当引擎。

### 为什么不用 PMCC（水听器/OBS 方法）

- **PMCC（Progressive Multi-Channel Correlation）是水听器/OBS 台阵方法**（Cansi 1995，CTBTO IMS 标准），设计用于台间距几十~几百米到几 km 的密集水听器/OBS 阵列
- 陆地宽频台网（南非 ~100 km 间距）对 2-8 Hz T 波空间采样严重不足（需 <1.5 km）
- 实测 4 种方法全失败（0 组有效），且**三站互相关时延在 T-P 包络中心几乎同时到达（dt≈0）**，几何退化
- 注意：PMCC ≠ 频域 F-K 波束成形（Capon/Bartlett）。我们用的 `t_phase_analysis.pmcc` 是**时域互相关时延反演**，包内无真正 F-K 波束成形函数

### 方法记录（如需将来应用到 OBS/水听器数据）

| 方法 | 做法 | 适用场景 |
|:---|:---|:---|
| 传统 PMCC | 三站 2-8 Hz 波形，逐对互相关→时延→扫(baz,V) 找最小 misfit | 密集水听器/OBS 台阵 |
| 包络低通法（Zhou 2021） | 2-8Hz→Hilbert→LP 0.25Hz→PMCC | 稀疏台阵+统计堆叠（需 ≥80 站） |
| 延迟扫描 beamforming | 网格扫(baz,V)→堆叠相干最大化 | 有真实波前时 |

工具：`dependencies/t_phase_analysis/t_phase_analysis/pmcc.py`；参考 `references/fk-africa-attempts.md`、`references/fk-spacing-limits.md`（single-station skill）。

**陆地台网替代方案**：偏振 R-T 玫瑰图独立验证反方位角（single-station-analysis 流程已实现）。
## 包络缓存与剖面图对齐

`process_image.py` 在生成带通图时，同时保存步0 目录下的 `envelope/{STA}_envelope_2-8Hz.npz`（**新规范** `process_flow/0_data_processing/envelope/`；旧事件 `<event>/sacprocess/envelope/`），包括统一时间轴、Z/R/T 包络、带通波形和处理参数。`waveform_proflie.py` 优先读取该缓存，避免重复滤波/包络计算；只有旧事件没有缓存时才回退重算。

## 参考文档（`references/`，不可写成代码的知识）：
- `merged-travel-time.md` — 走时理论、SOFAR 传播、Love 波
- `merged-attenuation.md` — T波衰减、Q 值拟合
- `merged-polarization.md` — 偏振分析方法
- `merged-events.md` — 历史事件处理记录
- `merged-tools.md` — 工具/环境配置
- `merged-fk-pmcc.md` — FK/PMCC 阵列分析
- `stationxml-adaptation.md` — StationXML 适配（IRIS 2025年后 .xml 格式）
- `merged-events.md` — 历史事件处理记录（见下方各事件摘要）

## 处理事件速查表

| 事件 | 日期 | Mw | 台站 | 成功 | 缺SACPZ | 备注 |
|------|------|----|------|------|---------|------|
| Prince Edward | 2018-04-19 | 6.0 | ~ | ~ | ~ | 早前处理 |
| Prince Edward | 2023-05-21 | 6.8 | 51→6 | 6 | ~ | IRIS 2023年波形全204；BOSA仅5s不可用 |
| Southwest Indian Ridge | 2013-08-17 | 6.2 | 97 | 76 | YV×12+G×1→已恢复 | 最近站54km陆地距，SNR最高59.4dB |
| Puerto Rico | 2014-01-13 | 6.3 | ~ | ~ | ~ | 方位327°/350° |
| Mid-Atlantic Ridge | 2014-07-27 | 6.1 | ~ | ~ | ~ | 方位310° |
| Mindanao | 2026-06-07 | 7.8 | ~ | ~ | ~ | 方位180° |

预处理规则和边缘情况已写入各脚本 docstring（`process_image.py`, `process_single_station.py`）。
画图标准 → `seismic/plotting-standards` skill

## 文献地图任务

当研究需要基于论文重现台站分布图时（没有实际波形数据）：

1. **直接从论文 Figure 提取台站坐标**，不尝试 IRIS 下载（历史区域性台网数据通常不在 FDSN 库）
2. **基于标准 `map.py` 模板修改**，不要从零写脚本 — 只改顶部区域、坐标、地名
3. **海底地形**：用 `@earth_synbath_15s`（含水深数据），不用 `@earth_relief_01s`（仅陆地 SRTM）
4. **传播路径**：大圆路径虚线 + 简洁箭头，从震源指向海岸截止（不穿过岛屿）
5. **地名用白色描边标注，字体不加粗**：`'10p,Helvetica,black'`（无 `-Bold`），`WHITE_STROKE=0.02`
6. **Hawaii Koyanagi 1995 案例**：
   - 参考脚本：`/home/sy/work/Hawaii_Koyanagi1995/plot_topography.py`
   - 三个事件震源机制见 Global CMT，坐标见论文 Table 1
   - 50+ HVO 台站坐标见下方参考

## 已知陷阱 ⚠️

### 🔴 去响应异常 debug 四步法（去响应失败排查）

频谱图空白/全噪时，按顺序隔离：
1. **打印 poles/zeros 数量**（查 `ZEROS` 拼写、是否重复块）
2. **打印去响应后数据范围**（正常应 ±10^-4 量级 m/s；10^-16 则去响应失败）
3. **检查 `water_level` 和 `pre_filt`**
4. **检查 SACPZ `SENSITIVITY` 单位**（count/(m/s) vs count/m 混用）

### 🔴 STFT 频谱图时间轴偏移

`scipy.signal.spectrogram` 返回的 `t_spec` 以每个窗口的**中心**为参照，第一个 bin 偏移 `nperseg/(2*sr)` 秒。直接 `t_spec + t_off` 会让频谱图比波形滞后此偏移量。**修复**：`t_spec_rel = t_spec + t_off - nperseg/(2*sr)`（BOSA 案例：nperseg=1024 时偏移 12.8s，肉眼可见）。

### 🔴 NaN SNR

部分台站 SNR 计算返回 NaN（近台噪声窗口不足 → log(0)）。修复方案（2026-08-03）：
1. **动态噪声窗**：`actual_noise = max(20, min(noise_width, available_pre * 0.8))` — 根据 P 波时间自适应
2. **噪声中心用 `tt_p`**（非 `tt_p - 100`）— 用全部可用噪声
3. **`astype(np.float64)` 强制转换** — 避免大值溢出 float32 导致 PSD 为 inf/NaN
4. **回退值 `0.0`** — 当噪声窗口 < 20s 时返回 0.0（非 NaN 或 -999）

**验证**：FOMA 站（P=166s, land=54km）修复前 SNR=NaN，修复后 SNR=30.6 dB。

### 🔴 XK 台网毛刺（规律脉冲噪声）

XK 台网的 B*/W* 子组（2013 SW Indian Ridge）存在规律性短时高强度脉冲，伪造成高 SNR（最高 31.4 dB）。`step1` 已内置 `DEFAULT_EXCLUDE` 列表（W08KB, B04KH, W07CR, B06OR, B15MW, B07DX 等）。若新事件发现类似台站，追加到该列表。

### 🔴 SAC 头补 stla/stlo（模板前置步骤）

`step1` 模板从 SAC 头读取台站坐标。若 SAC 头缺 `stla/stlo`（常见于新下载数据），模板会返回 0 个台站。**运行模板前补 SAC 头**：从 FDSN `get_stations()` 取坐标，写回 `tr.stats.sac.stla/stlo/evla/evlo/evdp`。

### 🔴 太平洋事件：日期线跨越

`step1` 自动推算区域 `[min(lon)-pad, max(lon)+pad]`。台站跨日期线时（如日本→阿拉斯加）经度跨度 > 360°，pyGMT 报 `Map region exceeds 360 degrees`。**修复：手动设区域**（以日本 142°E 为中心）：`reg_lon_min = 120; reg_lon_max = 260`（跨日期线用 >180）。

**地形下载缓慢**：15s 地形单 tile 下载 3-5 分钟，用 01d（1度）分辨率快 10x：`fig.grdimage(grid='@earth_relief_01d', ...)`。

### 🔴 events.yaml 冲突

事件名在 `config/events.yaml` 已存在时（如 2025 和 2026 同区域事件），用不同 key 名（如 `honshu_east_2025` / `honshu_east_2026`）。不要覆盖已有条目。

### 其他小陷阱

- **路径**: 输出到 `/mnt/<盘符>/sharedate/`，用 Windows 盘符绝对路径告知用户
- **SAC 头缺失**: `stla/stlo` 常缺失（见上）；`sac.mag` 常为 -12345（从 USGS 目录或目录名取）
- **文件名时间解析**: `parts[-4]=year`, `parts[-3]=julday`, `parts[-2]` 可能带 `_proc` 需 strip
- **原始波形 DC 偏移**: 图1 保留 DC 不 demean，其他图用去响应数据（已 demean+detrend）
- **BH1/BH2 SAC 和 SACPZ 同步重命名**: 只改 SAC 不改 SACPZ → `find_sacpz('STA','BHN')` 搜索不到 → 台站被跳过。数据处理前 `ls raw/SACPZ/*BH1 raw/SACPZ/*BH2` 确认无残留
- **数据整理前全面扫描**: 检查所有网络通道命名（不仅 II/IU），YV/PS/AV 等区域台网也可能用 BH1/BH2
- **频谱图时间轴需用去响应数据**: 频谱面板用 `trs_resp`（仅去响应），不用 BP 滤波数据
- **sed 替换 f-string**: f-string 内日期需单独替换 `D.2019.022.`

> 绘图细节（ylim 统一、nperseg/nfft 分辨率、P2-P95 vmax、`spectrogram as sp` 变量覆盖、`*BHZ_resp.sac` glob 写法）已在模板代码中标注，不再在 skill 重复。

## 单台分析输出规范

**标准图**（`--plot` 标志，输出 `{STA}_bandpass_spectrogram.png`）：
- 3×2 布局：左列=带通包络(ZRT)，右列=频谱图 0-10Hz(ZRT)
- 左列：灰波形(#555555, lw=0.25) + 红包络(#d62728, lw=0.6, 60点平滑) + P98归一化
- 右列：viridis 色标，P2-P95 vmax
- 标注 P/S/T-P/T-S 到时虚线
- 输出路径：`process_flow/2_single_station_analysis/{STA}/`
- 调用：`process_single_station.py <STA> --network NET --raw-dir RAW --plot`

### 原始波形 DC 偏移

原始 SAC 数据有 DC 偏移（计数均值非零，如 BOSA BHZ 均值 ~14000 counts）。
图1（原始波形）**保留 DC 偏移不 demean**，其他图使用去响应后的数据（已 demean+detrend）。
