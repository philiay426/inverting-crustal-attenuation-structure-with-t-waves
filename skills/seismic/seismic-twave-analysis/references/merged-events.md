# 事件相关笔记

# 2014-01-13 Mw 6.3 Puerto Rico — Waveform Profile Learnings (June 2026)

## T-P Arrival: Systematic Bias in `tt_t()`

**Finding**: `tt_t()` from `t_phase_analysis.traveltime_profile` computes T-P arrivals that are systematically **~130s (median)** earlier than the actual T-wave envelope peak. All 34 stations in AZ 341–360° showed offsets of 56–224s (mean 137s, std 39s).

**Root cause**: The function uses great-circle path to -1000m isobath intersection with fixed V_OCEAN=1.48 km/s. The actual SOFAR channel propagation includes:
- Non-great-circle path deviations
- Scattering and reverberation at the coastal conversion zone
- Bathymetry-dependent velocity variations

**Fix options**:
1. Apply +130s bulk correction (simple, ~39s residual scatter)
2. Use data-driven envelope onset picking instead of theoretical T-P
3. Use observed T-P from manual pick for key stations

## Bahamas Bathymetry: Azimuth-Dependent Ocean Path Length

**Finding**: For the Puerto Rico → US East Coast path, lower azimuth stations (341–347°) have **~180 km shorter ocean paths** than higher azimuth stations (350–354°).

| Azimuth group | Ocean distance | T-P range | Cause |
|---|---|---|---|
| 341–347° | ~2165–2275 km | 1529–1580 s | Bahamas shallow banks → -1000m isobath further offshore |
| 350–354° | ~2325–2343 km | 1619–1638 s | Deep Florida offshore waters → isobath close to coast |

**Impact**: When sorting stations by land distance, T-P markers appear ~100s misaligned between azimuth groups. This is physically correct but visually confusing.

**Stations affected** (AZ<347°): PAL, M61A, L59A, M58A, M59A, K59A, G62A, G63A, G64A, I63A, H61A, J60A, H58A

## Waveform Profile Plotting Conventions

For figures showing amplitude attenuation with land distance, use these settings:

### Normalization
- **Global P98 normalization** (NOT per-trace): preserves relative amplitudes between stations
- Per-trace normalization masks attenuation — all traces look equally strong

### Amplitude clipping
- Soft-clip stations with per-trace P98 > 3× median P98
- Prevents near-coast stations (e.g., HRV at >10,000× median) from dwarfing all others
- Apply before global P98 computation for consistent scaling

### Visual style
```
waveform:  'k-', lw=0.5, alpha=0.9
envelope:  color='gray', lw=0.5, alpha=0.6
fill:      color='gray', alpha=0.08
```

### Annotations
- Remove station labels when the goal is to show attenuation only
- T-P markers: black dots (ms=7, white edge)
- T-S markers: gray dots (ms=7, white edge)
- No vertical lines, no phase text labels

### Legend
- Upper left corner, fontsize=18, single column
- Four entries: waveform, envelope, T-P (black dot), T-S (gray dot)

### Title
- Minimal: `AZ 341°–360° | 34 stations | SNR ≥ 3`
- No explanatory text, no normalization mention

## Script: profile_wide.py Features

Located at: `process_flow/2_waveform_profile/profile_wide.py`

Key parameters:
- `--az-min` / `--az-max`: azimuth filter
- `--snr-min`: SNR threshold (signal window: T-P to T-S+60s; noise: T-P-300 to T-P-60)
- `--t-min` / `--t-max`: time window
- `--amp`: waveform amplitude scaling (default 0.10)
- `--spacing`: vertical trace spacing (default 0.35)
- `--sort-by`: 'land' (land distance) or 'dist' (epicentral distance)
- `--align-tp`: shift time axis so T-P=0 for all traces
- `--exclude`: comma-separated station blacklist (e.g., "X51A,X60A")

### Path fixes required for current environment
```python
EVLA, EVLO = 19.0428, -66.8096  # NOT 18.99, -66.87
EVDP = 20.0                       # NOT 10.0
PROC_DIR = "/mnt/f/sharedate/case/2014-01-13-mw63-process/sacprocess"
OUT_DIR = "process_flow/2_waveform_profile"
sys.path.insert(0, '/home/sy/codex/t_phase_analysis')
```

# 3C Column Detection in Station Table (2026-07-07)

## Problem

Need to mark which stations have complete 3-component data (BHE or BH1)
for polarization analysis eligibility.

## Detection Code

Add to `step1_station_table_and_map.py` after station identification:

```python
parts = os.path.basename(fn).split('.')
sta = parts[1]
net = parts[0]

# Check for complete 3C (BHE/BH1)
bhe_files = glob.glob(f"{RAW_DIR}/{net}.{sta}..BHE*.SAC") + \
            glob.glob(f"{RAW_DIR}/{net}.{sta}.*.BHE*.SAC") + \
            glob.glob(f"{RAW_DIR}/{net}.{sta}..BH1*.SAC") + \
            glob.glob(f"{RAW_DIR}/{net}.{sta}.*.BH1*.SAC")
has_3c = 'Y' if bhe_files else 'N'
```

**Critical**: `RAW_DIR` must be defined! Add alongside `PROC_DIR`:
```python
PROC_DIR = f"/mnt/e/sharedate/case/{case_dir}/sacprocess"
RAW_DIR = f"/mnt/e/sharedate/case/{case_dir}/raw"  # ← must add
```

## Excel Formatting

Add 3C column as column 14:

```python
hdrs_cn = ['台站', '方位角', '震中距', '陆地距', '海洋距', 
           'P波', 'S波', 'T-P', 'T-S', '信噪比', '经度', '纬度', '反方位角', '三分量']
hdrs_en = ['STA', 'AZ(°)', 'DIST(km)', 'LAND(km)', 'OCEAN(km)',
           'P(s)', 'S(s)', 'T-P(s)', 'T-S(s)', 'SNR(dB)', 'LON', 'LAT', 'BAZ(°)', '3C']

# Data row
vals = [r['sta'], r['az'], r['dist'], r['land_km'], r['ocean_km'],
        r['tt_p'], r['tt_s'], r['tt_tp'], r['tt_ts'], r['snr'],
        r['lon'], r['lat'], r['baz'], r['has_3c']]

# Highlight
if ci == 14:  # 3C column
    if v == 'Y':
        cell.fill = PatternFill('solid', fgColor='C6EFCE')  # green
        cell.font = Font(bold=True, size=14, name='Microsoft YaHei', color='006100')
    else:
        cell.fill = PatternFill('solid', fgColor='F2F2F2')  # gray
        cell.font = Font(size=14, name='Microsoft YaHei', color='808080')
```

## CSV Format

Add `has_3c` column:
```python
csv_lines.append("sta,az,dist,land_km,ocean_km,tt_p,tt_s,tt_tp,tt_ts,snr,lon,lat,baz,has_3c\n")
# Format string: add {}, at end
"{},{:.1f},{:.0f},{:.0f},{:.0f},{:.1f},{:.1f},{:.1f},{:.1f},{:.2f},{:.4f},{:.4f},{:.1f},{}\n"
```

## Example Results

- 2015-12-04 Mw 7.1: 44 Y, 8 N (WAA stations mostly missing BHE)
- 2022-02-04 Mw 6.3: varies by network
- 2022-11-22 Mw 7.0: varies by network

# Event Parameter Verification (2026-07-22)

## Critical Pitfall: Wrong Epicenter Coordinates
Several Australian events had incorrect EVENT_LAT/EVENT_LON in process_obs.py, leading to completely wrong T-wave arrival times. The correct coordinates were in the `stations.csv` header line.

### Example
**2022-02-04 Mw6.3 Southwest Indian Ridge:**
- Wrong (in script): `-37.773, 51.268`
- Correct (in CSV): `-48.033, 99.4949`

With wrong coordinates, RKGY station T-P was computed as 3835s instead of correct 1349s.

### Verification Workflow
```python
import re, csv

# Read script's event params
with open('process_obs.py') as f:
    content = f.read()
py_lat = float(re.search(r'EVENT_LAT\s*=\s*([-\d.]+)', content).group(1))
py_lon = float(re.search(r'EVENT_LON\s*=\s*([-\d.]+)', content).group(1))

# Read table's event params
with open('process_flow/1_station_distribution/stations.csv') as f:
    header = f.readline()
m = re.search(r'([-\d.]+),?\s*([-\d.]+)', header)
csv_lat, csv_lon = float(m.group(1)), float(m.group(2))

if abs(py_lat - csv_lat) > 0.01 or abs(py_lon - csv_lon) > 0.01:
    print(f"MISMATCH! py: {py_lat},{py_lon} vs csv: {csv_lat},{csv_lon}")
```

Always run this check before doing any T-wave analysis on a new event.

# 事件清单扫描与处理状态检查

当需要快速了解工作目录下有哪些事件、各处于什么处理阶段时，使用 SAC 头扫描模式：

```python
from obspy import read
import glob

for d in sorted(glob.glob('*/raw')):
    label = d.split('/')[0]
    sac_files = glob.glob(f'{d}/*.SAC') or glob.glob(f'{d}/**/*.SAC', recursive=True)
    if not sac_files:
        print(f'{label}: 无SAC文件')
        continue
    
    # 只读头信息，不加载波形数据
    st = read(sac_files[0], headonly=True)
    tr = st[0]
    
    # ⚠️ sac 属性可能缺失，全部用 getattr 防御
    evla = getattr(tr.stats.sac, 'evla', None)
    evlo = getattr(tr.stats.sac, 'evlo', None)
    evdp = getattr(tr.stats.sac, 'evdp', None)
    mag = getattr(tr.stats.sac, 'mag', None)  # mag 常缺失或为 -12345
    
    # 台站统计
    stations = set()
    nets = set()
    for f in sac_files:
        try:
            s = read(f, headonly=True)[0].stats
            stations.add(s.station)
            nets.add(s.network)
        except:
            pass
    
    print(f'{label}: {len(stations)}台站 {sorted(nets)}')
    print(f'  震中: {evla}, {evlo}, {evdp}km')
```

**处理状态判断**（目录中有什么 = 做了什么）：
| 存在内容 | 处理阶段 |
|---------|---------|
| 仅 `raw/` | 只下载了数据 |
| `raw/` + `download_data.py` | 下载完成 |
| `sacprocess/` | 预处理完成 |
| `waveform_image/` 或 `*_raw.png` | 已出图 |
| `process_flow/` | 完整流程进行中 |

# 事件参数验证 + Python 环境

## 震中参数必须从 stations.csv 验证

每个事件的 `process_flow/1_station_distribution/stations.csv` 第一行包含正确的事件参数：
```
# Event: 2022-02-04 Mw 6.3 | -48.033N, 99.4949E | Depth 10.0km
```

**必须提取并对比 process_obs.py 中的 EVENT_LAT/LON！**

2022-02-04 事件曾因震中写错（-37.773/51.268 → 正确是 -48.033/99.4949），导致 T 波走时错误 2500s。

## 修复步骤
```bash
# 1. 提取 CSV 中的正确坐标
head -1 stations.csv | grep -oP '[\-\d.]+N.*[\-\d.]+E'

# 2. 对比脚本中的值
grep "EVENT_LAT\|EVENT_LON" process_obs.py

# 3. 如果不一致，修改后重新运行
python process_obs.py --all
```

## Python 环境问题

当 WSL 上有多个 Python 时（如 3.11 和 3.14），conda 环境可能被 ~/.local/lib/python3.14/ 下的包干扰。

**修复**：在 conda 环境内重装
```bash
conda activate obspy
pip install --force-reinstall numpy matplotlib scipy
```

或升级整个 conda 环境到同一 Python 版本：
```bash
conda install python=3.14 -y
pip install obspy geographiclib
```

# Event Processing Workflow (2026-06)

## Entry Procedure (When User Says "Process Event X")

When the user asks to process a specific event, follow this CHECKLIST before running any code:

1. **Get event info**: Query USGS FDSN for coordinates, depth, magnitude
2. **Check event directory exists**: `/mnt/f/sharedate/case/<event-label>/`
   - If not found, check for variations in naming (different date format, magnitude rounding)
   - Ask user for the exact directory path if ambiguous
3. **Check what's already done**:
   - Look at directory structure: `raw/`, `sacprocess/`, `process_flow/`
   - Check for existing Excel files (`stations.xlsx`) and output figures
   - Check the event-processing-workflow events table below
4. **Load relevant skills FIRST**:
   - `seismic-twave-analysis` (full pipeline reference)
   - Read `event-processing-workflow` (this file)
   - For single-station work: `single-station-analysis`
5. **Determine next step**: Based on what's already present:
   - Raw SAC only → Step 1 (preprocessing)
   - Has sacprocess/ → Step 2 (profile / maps)
   - Has process_flow/1_station_distribution/ → proceed to Step 2 or 3
   - Has process_flow/2_waveform_profile/ → Step 3 (single-station)
6. **Verify correct event directory**: Before running ANY script, confirm `cd` target and file paths match the event being processed. Running a previous event's code is the most common error.
7. **Then proceed with the 4-step workflow below**.

## Pipeline Overview

For a new earthquake event with 3-component broadband data:

```
raw SAC files (user downloads from IRIS FDSN)
    ↓ process_obs.py or 01_preprocess.py (codex)
sacprocess/
├── remove_response/*_resp.sac   ← Deconvolved velocity ONLY (no bandpass)
├── *_proc.sac                    ← Deconvolved + 2-8 Hz bandpass
└── waveform_image/               ← Per-station PNG previews
    ↓ 02_extract_stations.py
stations.xlsx + gmt_stations.txt
    ↓ pyGMT / 06_profile_section.py
arrival_map.png + record_section_*.png
```

## Events Processed (8 events as of July 2026)

| Event | Date | Mw | Coords | Stations | Status |
|-------|------|----|--------|----------|--------|
| **Barbados** | **2014-02-18** | **6.5** | **14.67°N, 58.93°W** | **TA+IU+US+LD=268 sta, US East Coast** | **✅ All steps done** |
| Central MAR | 2014-07-27 | 6.0 | 23.72°N, 45.58°W | TA+IU+US+LD, US East Coast | 📄 Raw SAC only |
| Northern MAR | 2015-02-13 | 7.1 | 52.65°N, 31.90°W | (raw download pending) | 📥 Pending |
| **SE Indian Ridge** | **2015-12-04** | **7.1** | **47.62°S, 85.09°E** | **AU/6C/S1/YS/II/IU/G/GE=~118 sta, Australia** | **✅ All steps done** |
| **SE Indian Ridge** | **2022-02-04** | **6.3** | **48.03°S, 99.50°E** | **AU/S1/IU/G=~53 sta, Australia** | **✅ Steps 1-2 done** |
| **Solomon Islands** | **2022-11-22** | **7.0** | **9.82°S, 159.60°E** | **AU/S1/YW=~80 sta, Australia** | **✅ Steps 1-2 done** |
| **SE Indian Ridge** | **2025-04-16** | **6.6** | **47.81°S, 99.73°E** | **AU/S1/II/IU/G=~119 sta, Australia** | **✅ Steps 1-3 done** |
| **Mindanao** | **2026-06-07** | **7.8** | **5.59°N, 125.05°E** | **AU/S1/II/IU/G, Australia** | **✅ Steps 1-3+ done** |

## Codex Workflow (4 Steps)

| Step | Name | Scripts | Output |
|:----:|------|---------|--------|
| 1 | 站台分布 | `01_preprocess.py` + `02_extract_stations.py` | `_resp.sac` + `_proc.sac` + stations.xlsx + gmt_stations.txt |
| 2 | 波形剖面 | `06_profile_section.py` / pyGMT map + `profile_wide.py` | record_section.png + arrival_map.png |
| 3 | 单台分析 | `three_component_analysis.py`, `polarization_analysis.py`, `plot_waveform_windows.py` | **7张图**: 三分量时频谱 + PS放大窗 + T波放大窗 + P/S±5s偏振 + P/S全窗偏振 + T-P/T-S偏振(Flinn协方差) + T波全窗偏振 |
| 4 | 衰减分析 | `attenuation_analysis.py` | 子频段 ln(A) vs distance + t* + Q估计 |
| 5 | F-K分析 | `07_fk_analysis.py` (PMCC) | (future) |

### Step 1: 预处理多网络数据

对于跨多个网络的台阵（如 Mindanao 2026 涉及 AU, S1, II, IU, G 五个网络），`01_preprocess.py` 会自动发现所有网络的 SAC 文件。需要注意：

- **Location code 差异**：不同网络使用不同的 location code（AU.KDU.00.BHZ vs S1.AUWSH..BHZ），预处理脚本必须支持多 pattern glob
- **采样率差异**：AU 网络 40 sps、G.CAN 20 sps、部分 AU 台站 10 sps — 2-8 Hz BPF 对 20 sps 台站会在 Nyquist 边界工作（ObsPy 自动降为高通），不影响处理流程
- **BH1/BH2 vs BHN/BHE**：II/IU 网络使用 BH1(0°)/BH2(90°) 替代 BHN/BHE，通过 `cmpaz` 判断
- SACPZ 文件搜索使用网络+台站+通道匹配，覆盖 489 个文件

### Step 1: 台站表格格式标准（2026-07 定稿）

| 项目 | 标准 |
|------|------|
| 字体 | Microsoft YaHei 14pt |
| 边框 | thin（浅灰色） |
| 表头 | 行2中文+行3英文（蓝底白字） |
| 数据起始 | 行4（无空白行） |
| 列宽 | 自适应（`max_len * 1.8 + 2`） |
| 3C列 | 第14列，Y=绿底/N=灰底 |
| SNR列 | 2位小数，≥阈值绿底高亮 |

**脚本：** `step1_station_table_and_map.py`

### Step 1: 台站分布图标准（2026-07 定稿）

**同时保留两种参考线：**
1. **方位角弧线**：黑色虚线，最近+最远海岸线，标注角度
2. **震中距圆弧**：灰色虚线，每1000km一圈，标注距离

**放大图：** 无绿色圆圈标记

**台站尺寸按SNR分档：**
- 全域图：0.12/0.22/0.35/0.50 cm
- 放大图：0.30/0.45/0.60/0.75 cm

**脚本：** `plot_station_map.py`

### Step 1: 清除无用中间文件

每次预处理运行后清理 `process_obs.log` 和其他临时输出。用户倾向干净的目录结构。

### Step 3: 偏振方法（Flinn 1965 协方差法）

偏振分析已全部改为 Flinn (1965) 协方差法。在分析窗口内滑动5s子窗（75%重叠），每个子窗计算3×3协方差矩阵并特征分解，用主方向做玫瑰图、线性度做权重。P和T-P窗口自动生成Z-R+R-T双玫瑰图，其余窗口仅R-T。

详见 `single-station-analysis` 的 `references/flinn-1965-covariance-polarization.md`。

| 分析类型 | 滤波频段 | 原因 |
|---------|---------|------|
| P/S偏振 | **0.01-50 Hz** (宽带) | 保留体波全频段偏振信息 |
| T-P/T-S偏振 | **2-8 Hz** | T波能量集中在此频段 |
| 三分量/波形窗 | **2-8 Hz** | 目视识别用 |

`polarization_analysis.py` 的 `load_station(fmin, fmax)` 和 `plot_polarization(fmin, fmax)` 均接受 filter 参数。
main() 调用时 P/S 传 `fmin=0.01, fmax=50.0`，T-P/T-S 用默认 `fmin=2.0, fmax=8.0`。

### Step 3: 自定义走时模型（非大西洋路径）

对于岛弧/浅海路径（Mindanao→Australia），`tt_t()` 的 `ref_el=-1000m` 阈值误将 Timor Sea（200-800m 深）分类为陆地，导致 T-P 系统偏晚 ~50s。

**解决方案**：对非开阔大洋事件，创建 `custom_arrivals.py` 使用经验速度模型：
```python
# V_app = 2.44 km/s (KDU 标定)
tt_TP = dist_km / V_app
tt_TS = dist_km / V_app_s     # V_app_s = 3.6 km/s
```

适用的路径类型：
- **开阔大洋**（大西洋、深太平洋）→ 继续使用 `tt_t()`（误差 <10s）
- **岛弧/浅海**（印尼群岛、菲律宾→澳大利亚）→ 使用 `custom_arrivals.py`

详见 `references/custom-arrivals-empirical-velocity.md`。

### Step 3: 偏振窗口迭代校准工作流

用户每次改窗后需要看到结果再决定下一步调整，**不允许一次性给出最终窗**：

```
(1) 先运行 plot_waveform_windows.py → 输出P_window/S_window/Twave_window 3张图
    → 用户观察包络和到时标注，判断窗口位置是否合理
(2) 用户反馈："窗太大"→收窄, "试试±5s"→精确到秒级修改
(3) 每次只保存 polarization_analysis.py 的 wps/wts 参数，不替用户运行
(4) 最终窗宽通常远小于理论±20-50s
```

### Step 3b: 视速度分析（波长传播模式判断）

通过线性回归 `T-P_obs ~ land_km` 判断T波在陆地的传播模式：

| V_land 范围 | 传播模式 |
|:-----------|:---------|
| 3-4 km/s | Love 波 |
| 4-6 km/s | S 波 |
| 7-9 km/s | P 波 |

**核心方法**（moveout，非残差法）：
```python
# 海洋路径近乎相同 → intercept 吸收海洋走时
# slope = 1/V_land（s/km）
A = np.vstack([land_km, np.ones_like(land_km)]).T
slope, intercept = np.linalg.lstsq(A, t_obs, rcond=None)[0]
v_land = 1.0 / slope
```

⚠️ **不要用 "T-P 减去海洋走时"**：陆地走时是两大数相减的小残差（30-70s vs 1620s总走时），噪声完全主导结果。

详见 `references/velocity-analysis-moveout.md`。

### Step 2: 波形剖面的到时标记

P/S/T-P/T-S 标记必须使用**固定 y 偏移量**，不能用局部振幅值：

```python
# ✅ 正确
ax.plot(tt, y0+AMP*0.6, 'o', zorder=5, mec='white', mew=0.5)
# ❌ 错误：T波包络/全局P98极小 → 标记不可见
na = s['z'][idx]/p98*AMP
ax.plot(tt, y0+na, 'o')
```

### Step 2: 波形剖面图标准（2026-07 更新）

#### 核心原则

1. **无方位角限制** — 使用全部台站，不按方位角筛选
2. **横轴：海岸转换点对齐** — `t - tt_ocean`，t=0 = 波上岸
3. **纵轴：陆地距离** — 近在下，远在上
4. **SNR 从 CSV 读取** — 使用 step1 预计算的值，不重新计算
5. **左侧无台站名标注** — 保持简洁
6. **左侧紧贴轴** — 无留白

#### 脚本

标准脚本：`profile_conversion.py`（位于每个事件的 `process_flow/2_waveform_profile/`）

```bash
~/envs/tphase_env/bin/python3 profile_conversion.py <event_key> --snr-min 3.0
```

#### 关键实现

```python
# tt_ocean 直接从CSV的ocean_km计算 (SOFAR通道速度1.48km/s)
for sta, info in arrivals.items():
    info['tt_ocean'] = info['ocean_km'] / 1.48

# 转换点对齐: t' = t - tt_ocean，t=0 表示波上岸
t_aligned = t - tt_ocean
traces.append(dict(
    tt_tp=info['tt_tp'] - tt_ocean,
    tt_ts=info['tt_ts'] - tt_ocean,
    ...
))

# 纵轴设置 (近在下，远在上)
ax.set_ylim(min_land - margin, max_land + margin)

# 画布高度 (拉长版)
fig_h = min(30, max(10, (max_land - min_land) / 30))

# 左侧紧贴
fig.subplots_adjust(left=0.02)
ax.set_xlim(args.t_min, args.t_max)  # 无额外留白
```

#### 包络线颜色

```python
ax.fill_between(tw, y0, y0+ey, color='red', alpha=0.08)      # 包络填充：淡红
ax.plot(tw, y0+ey, color='red', lw=0.5, alpha=0.6, zorder=4) # 包络线：红
```

#### 字体标准（期刊）

```python
ax.set_xlabel("Time from ocean arrival (s)", fontsize=14, fontweight="bold")
ax.set_ylabel("Land distance (km)", fontsize=14, fontweight="bold")
ax.set_title(f'{evdate} Mw {evmw}  |  {n} stations  |  SNR ≥ {args.snr_min} dB  |  Aligned by conversion point',
             fontsize=14, fontweight='bold')
```

#### 转换点对齐（核心方法）

波形剖面统一采用**转换点对齐**：每条道减去各自的海洋走时 `ocean_km/1.48`，x轴为"转换点后时间(s)"，所有台站的海岸到达对齐到 t=0。

```python
V_SOFAR = 1.48
t_coast = ocean_km / V_SOFAR          # 海洋走时
t_rel = t_abs - t_coast                # 对齐后时间轴
tp_land = tt_TP - t_coast             # T-P陆地走时
ts_land = tt_TS - t_coast             # T-S陆地走时

# 紫色参考线标注转换点
ax.axvline(0, color='purple', ls='-', lw=1.5, alpha=0.6)
```

好处：
- 消除海洋路径差异 → 纯粹反映陆地传播特征
- T-P/T-S moveout 斜率 = 1/V_land
- 同方位角群对齐效果最佳

模板脚本：`plot_az_split_aligned.py`

#### 排序与筛选

- 排序：按陆地距升序
- SNR筛选：≥8 dB（频谱能量dB）
- 方位角分群：如 >341° 和 <341° 分别画图
- 台站去重：同名不同 location code 只保留一个

#### 视觉风格

| 元素 | 规范 |
|:-----|:-----|
| 包络 | 只画上半部，从 y0 向上填充 |
| 到时标记 | 大圆点 ms=7，白描边，上半部 `y0+amp*0.6` |
| 台站标注 | 右侧 `{STA} ({land}km)`，bold |
| 转换点 | 紫色竖线 t=0 + 顶部文字 |
| x轴 | Time from coast conversion (s)，50s主/10s副 |

### Step 3: plot_waveform_windows.py 输出

输出2张独立PNG（P+S合并为一张图）：

| 文件 | 内容 | 刻度 |
|------|------|------|
| `{STA}_PS_window_{t1}-{t2}s.png` | P+S合并单图, 覆盖P到S范围, 两相峰值同轴显示 | 50s主刻度/10s副刻度 |
| `{STA}_Twave_window_{窗}s.png` | T波垂直波形 | 50s主刻度/10s副刻度 |

窗口范围使用**动态值**，自动适应各台到时：
```python
ps_min, ps_max = int(tt_P-90), int(tt_S+100)     # P+S合并窗
tw_min, tw_max = int(tt_TP-80), int(tt_TS+80)     # T波波形窗
```

### Step 3: 执行

```bash
# 通用CLI模板（推荐）：无需修改脚本，命令行传参数
cp ~/shared/codex/templates/single_station/*.py ${OUT_DIR}/${STA}/
cd ${OUT_DIR}/${STA}
python3 three_component_analysis.py ${STA} ${NET} ${EVLA} ${EVLO} ${EVDP} ${BASE}
python3 plot_waveform_windows.py   ${STA} ${NET} ${EVLA} ${EVLO} ${EVDP} ${BASE}
python3 polarization_analysis.py   ${STA} ${NET} ${EVLA} ${EVLO} ${EVDP} ${BASE}
```

### Step 3: 文件命名规则

所有输出文件名需包含时间窗范围（秒），如 `J62A_p_s_polar_P308-328_S569-589.png`。

### Step 3: 多网络单台分析

模板脚本默认使用 TA 网络。对于混合网络事件，拷贝模板后需修改网络名：

```bash
sed -i "s/STA, NET = 'K63A', 'TA'/STA, NET = '{STA}', '{NET}'/" {STA}/*.py
sed -i "s/K63A/{STA}/g" {STA}/*.py
```

同时单台分析脚本必须支持多种 SAC 文件命名模式（location code = `00` / `--` / 无）：
```python
def find_sac(base, net, sta, comp, suffix='_resp.sac'):
    for pat in [f"{base}/{net}.{sta}.*.{comp}*{suffix}",
                f"{base}/{net}.{sta}..{comp}*{suffix}",
                f"{base}/{net}.{sta}.{comp}*{suffix}"]:
        m = glob(pat)
        if m: return m[0]
```

### Step 2: 波形剖面的台站筛选

对于方位角群剖面（如 AZ351±1°），需做以下筛选：

1. **Station dedup**: 同一台站（STA 名称相同）即便有不同 location code 也只保留一个。Python: `df = df.drop_duplicates(subset='sta')`。HRV 是典型重复案例（location 00 vs 10）。
2. **SNR 过滤**: 使用频谱能量dB法（Welch PSD, 信号窗 T-P~T-S+80s, 噪声窗 P-300~P-60s），SNR≥8 dB 作为信号可靠性阈值。
3. **时间窗聚焦**: 根据 T-P 到 T-S 的范围设置窄窗口（如 1500-2000s），使用 50s 刻度而非 500s。
4. **最终检查**: 筛选后台站数量应能有效展示 T 波到时趋势（10-30 个为佳）。

### Step 2: pyGMT 地图区域自动计算

对于非大西洋事件（如 Mindanao 2026），`03_map_arrival.py` 的硬编码 `-R-85/-45/20/50` 不适用。需从台站分布计算区域：

```python
lons = [r['lon'] for r in rows] + [EVLO]
lats = [r['lat'] for r in rows] + [EVLA]
margin = 3
REGION = [min(lons)-margin, max(lons)+margin, min(lats)-margin, max(lats)+margin]
```

同时台站密度较高时（>100）使用 label dedup：`iloc[::2]` + `min_sep=2.0°`。

## Azimuth-Targeted New-Event Startup (2026-06 Pattern)

When a new event has **raw SAC only** (no sacprocess/ directory) and you only need a subset of stations by azimuth:

### Workflow

```python
# Step 1: Event info from USGS
curl -s "https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson&starttime=2014-07-27&endtime=2014-07-28&minmagnitude=6&maxmagnitude=7"
# → 23.7240°N, -45.5809°W, depth=10km, Mw 6.0

# Step 2: Scan all stations → compute azimuths → identify target range
# Uses obspy.read(headonly=True) for fast header scan
# See references/azimuth-filtered-station-identification.md

# Step 3: Write TWO scripts for the azimuth subset:
```

| Script Type | Purpose | Key Features |
|------------|---------|-------------|
| `process_{az-range}.py` | Batch preprocess only target stations | Self-contained SACPZ parsing + simulate_seismometer + bandpass. No codex dependency (can run in `~/envs/tphase_env/`). Hardcodes TARGET_STATIONS list. |
| `02_extract_stations_{az-range}.py` | Station info → Excel | Uses `lib.arrivals.compute_theoretical_arrivals()` for T-P/T-S. Imports from `/home/sy/shared/codex/lib/arrivals`. Must run from codex project or add to sys.path. |

### When to Use vs Full Pipeline

| Use Case | Full Codex Pipeline (`01_preprocess.py --all`) | Azimuth-Targeted Scripts |
|----------|-----------------------------------------------|--------------------------|
| Event with 250+ stations, all needed | ✅ | ❌ |
| Only 21 stations in 315-317° needed | ❌ (processes all 250) | ✅ (processes only target) |
| First look at a new event | ❌ | ✅ |
| Event with multiple azimuth clusters | ❌ | ✅ (one script per cluster) |

### Script Location Convention

Place azimuth-targeted scripts in the **event root directory** with the azimuth range in the filename:

```bash
case/2014-07-27-mw61/
├── process_315_317.py           # ← batch preprocess: 21 stations
├── 02_extract_stations_315.py   # ← station info: 21 stations (uses codex lib)
├── TA.J59A..BHZ.M.2014.208....  # ← raw SAC files alongside
├── SACPZ/
└── sacprocess/                   # ← created by process_315_317.py
    ├── remove_response/
    └── *_proc.sac
```

### Key Design Decisions

1. **Hardcoded TARGET_STATIONS list**: Rather than scanning raw SAC and filtering by azimuth at runtime, hardcode the 21 station names. This makes the script:
   - Self-documenting (reader sees exactly which stations)
   - Fast (no glob + read for all 250+)
   - Reproducible (no runtime sensor selection variation)

2. **Self-contained preprocessing script**: The `process_315_317.py` reimplements SACPZ parsing and does NOT depend on the codex library. This lets it run under any Python env with obspy installed (e.g. `~/envs/tphase_env/`).

3. **Import from codex lib for arrivals**: The extraction script imports `from lib.arrivals import compute_theoretical_arrivals, dist_az`. This requires:
   - Running from `/home/sy/shared/codex/` or
   - Adding `/home/sy/shared/codex/` to sys.path
   
   The `compute_theoretical_arrivals()` function internally calls `get_track()` + `tt_t()` from `t_phase_analysis`, which handles ocean/land separation.

### Example: 2014-07-27 Mw 6.1 northern MAR

- **Event**: 23.7240°N, -45.5809°W, Mw 6.0, depth 10 km
- **Data**: 250+ raw SAC files in event root (no `raw/` subdir), TA + IU + US + LD networks
- **Target**: 21 TA stations at azimuth 315-317° (New England, USA)
- **Scripts**: `process_315_317.py` + `02_extract_stations_315.py`
- **Boundary check**: K62A (az=317.06°) and J60A (az=317.13°) are just outside the range; verify with user whether to include

## Pitfalls

- **⚠️ `get_track()` 第5参数 `dist_int` 单位是度(degrees)，不是 km**：传 `dist_km`（如 2966）会只产生2个点的 track，`tt_t()` 返回 None。正确：`dist_int=0.1`（0.1°≈11km）。详见 `references/ttt-api-pitfalls-202607.md`。
- **⚠️ `tt_t()` 返回 key 不是 `tp_min`**：实际 key 为 `tt_t`（秒）、`tt_ocean`（秒）、`tt_land_p`/`tt_land_s`（秒）、`dist_ocean`/`dist_land`（**度**不是km）。T-S = `tt_ocean + tt_land_s`。距离km = `dist_ocean * 111.0`。详见 `references/ttt-api-pitfalls-202607.md`。
- **⚠️ 多源 SACPZ 整合**：部分事件（如2022-02-04）的 SACPZ 文件分散在 IRISDMC/、RESIF/、GEOFON/ 等网络子目录中，需合并到 `raw/SACPZ/`。详见 `references/australian-events-202607.md`。
- **⚠️ 根目录 SAC 副本清理**：已处理事件的根目录可能残留 SAC 副本（与 raw/ 完全相同），应删除。验证方法：比较文件大小是否一致。
- **⚠️ profile_wide.py 方位角过滤**：脚本默认方位角可能不匹配事件实际台站分布。必须先检查 stations.csv 的 az 列，再用 `--az-min`/`--az-max` 指定正确范围。
- **⚠️ plot_station_map.py 配置路径**：脚本从 `../config/events.yaml` 读取事件参数。需在 `process_flow/config/` 创建 events.yaml。
- **⚠️ 脚本同步模式**：修改脚本后必须同步到所有事件目录。使用循环：`for d in */; do cp script.py "$d/path/"; done`。config中的case_dir要正确（含Australian/前缀）。
- **⚠️ 2015-12-04事件目录名**：实际目录是`2015-12-04-mww71-southeast-indian-ridge_out`（有`_out`后缀），不是`2015-12-04-mww71-southeast-indian-ridge`。config中case_dir需要更新。
- **⚠️ analyze_station.py 路径修复**：CODE_ROOT路径从`/home/sy/shared/codex`改为`/home/sy/codex`；文件名模式从`{base}..BHZ`改为`{base}.*BHZ`（支持位置代码如00）；使用glob动态查找文件而非硬编码日期。
- **P/S偏振必须用宽带(0.01-50 Hz)**：去仪器响应后只做0.01-50 Hz带通（防假频），保留体波全频段偏振信息。T-P/T-S偏振保持2-8 Hz。拷贝旧脚本时注意 `load_station()` 和 `plot_polarization()` 必须支持 `fmin,fmax` 参数。
- **T-P/T-S 理论到时 ≠ 实际包络峰值**：对经过浅海/岛弧路径的事件（如Barbados→Massachusetts），tt_t() 可能偏晚~50s。HRV案例：理论T-P=2072s，实测包络峰值~2120s（+48s，V_app=1.56 km/s）。单台分析前先用包络检测验证理论到时，必要时手动校准窗口。
- **T-P/T-S 偏振窗必须独立对准各自理论到时**：旧模式 `('T-S', tt_TP+10, tt_TP+40, ...)` 使T-S窗从T-P结束开始而非对准T-S自身。正确：`('T-S', tt_TS-5, tt_TS+15, ...)`，各窗口独立对准各自理论到时。标准窗口：T-P±10s, T-S(-5~+15s)，于Barbados 2014会话统一修正。
- **剖面图OOM**：200+台站时 `06_profile_section.py` 的 matplotlib 画布过大会触发 OOM killer。解决：先用 SNR≥2 过滤、或缩小方位角范围、或使用自定义脚本（详见 seismic-twave-analysis 的 `references/profile-oom-mitigation.md`）。
- **plot_waveform_windows.py 输出2张PNG**：PS_window（P+S合并）和 Twave_window，不再输出独立的P/S单窗。
- **PP dominates P at 15-30° epicentral distance**: At these distances, direct P is weak (upper mantle triplication) while PP arrives 30-40s later with 10-100× more energy. Always query `phase_list=['P','PP','S','SS']` and center P-wave windows on PP.
- **自定义走时模型**：`tt_t()` 不是万能的。岛弧/浅海路径的误差可达 50s（~5.6%）。新事件首次处理时，先用一个近距离台站（如 KDU）的实测走时标定 V_app，然后应用到全阵列。
- **SAC 文件发现**：多网络事件（AU, S1, II, IU 等）的 location code 不统一，单台分析脚本必须支持多 pattern glob，否则会因找不到文件而静默失败（`KeyError` 或空白图）。
- **03_map_arrival.py 硬编码区域**：该脚本的 `-R-85/-45/20/50` 仅对大西洋/US East Coast 事件有效。新事件必须改为动态区域计算（见 Step 2）。

# Legacy Event Directory: 2014-01-13-mw63-process (Pre-Codex)

**⚠️ This event directory predates the codex unified toolkit.** Its structure, paths,
and script organization differ from the standard codex workflow (`shared/codex/`)
described in the main SKILL.md. When modifying figures for this event, consult
this reference for the correct file layout.

## Directory Location

```
/mnt/f/sharedate/case/2014-01-13-mw63-process/
```

Symlinked from: `/home/sy/shared/case/2014-01-13-mw63-process/`

## Structure v.s. Codex

| Aspect | This Dir (Legacy) | Codex Standard |
|--------|------------------|----------------|
| Script location | Inside event dir (35+ `*.py` files) | `codex/scripts/01-07` |
| PROC path | `~/shared/case/sacprocess/` (shared with other events) | Event-isolated `case/<event>/sacprocess/` |
| Raw data | `~/shared/case/raw/` (shared) | `case/<event>/raw/` |
| Library path | `sys.path.insert(0, '/mnt/f/sharedate/T_Phase_Analysis-main')` | `codex/t_phase_analysis/` |
| Output root | `waveform_contrast/` inside event dir | `oi/<topic>/` (top-level) |
| Script count | ~35 scripts (many parallel versions) | 7 numbered scripts (unified) |

## Scripts Inventory

### Record Sections (Waveform + Envelope)

| Script | Purpose | Key Params |
|--------|---------|------------|
| `record_section_waveform.py` | Pure waveform profile, CLI param | `--azimuth 327/350 --snr-min --norm p98/max` |
| `record_section_350_twave.py` | 350° zoom (1550-1850s) | Hardcoded 350°/SNR≥2.5 |
| `plot_profile.py` | Source→station elevation profile | `--sta J59A` or from SAC file |

### Amplitude Feasibility (4-panel)

| Script | Purpose |
|--------|---------|
| `figure_amplitude_feasibility.py` | 327° 4-panel (a) section (b) P98 (c) ln amp vs land (d) band attenuation |
| `figure_350_amplitude_feasibility.py` | Same, 350°, parameterized |

Both exist in two places: event dir AND `case/code/`. The `case/code/` copies
are the more recent versions (e.g., use `~/shared/case/sacprocess/` path).

### Polarization Analysis

| Script | Purpose |
|--------|---------|
| `polarization.py` | Single station Flinn (`--sta J59A`) |
| `polarization_all.py` | Batch 11 original stations |
| `run_all_polar.py` | 11 stations → JSON + GMT |
| `run_all_polar_full.py` | Full 312-station batch |
| `compute_twave_polar.py` | Add T-P/T-S to JSON + ZRT polar |
| `plot_combined_rose.py` | Top-10 combined rose (mod 180°) |

### Other

| Script | Purpose |
|--------|---------|
| `plot_hrv_full.py` | HRV full 3C ZRT waveform |
| `plot_TP_TS_separation.py` | J59A T-P/T-S separation |
| `time_frequency.py` | Single-station spectrogram |
| `run_all_traveltime.py` | 343-station travel time → GMT map |
| `plot_top10_map.py` | Top-10 station location map (pyGMT) |
| `plot_top10_rose.py` | Top-10 P/S polarization roses |
| `plot_polar_map_gmt.py` | GMT map + polarization arrows |

### Low Code Subdir

`low_code/` contains exploratory/throwaway scripts (not maintained).

## Output Organization

Output goes to `waveform_contrast/` inside the event dir:

```
waveform_contrast/
├── one/                               ← Individual station outputs
│   ├── record_section_327_waveform.png
│   ├── record_section_350_twave_zoom.png
│   ├── figure_327_amplitude_feasibility.png
│   ├── figure_350_amplitude_feasibility.png
│   ├── polar_analysis_*.png               (J59A, HRV, X57A, etc.)
│   ├── polar_*.png                        (older polarization plots)
│   ├── HRV_full_3component.png
│   └── record_section_11stations_p98.png
├── oi/twave_paper_reproduce/          ← Publication-quality figures
│   ├── figure1_map_labels.png
│   ├── figure2a_record_section_v3.png   (two-column, az split)
│   ├── figure2a_record_section_v3.pdf
│   └── t_wave_detections_final.json
├── oi/all_polar_analysis/             ← (currently empty, see top-level oi/)
├── (top-level PNGs hosted here too)
│   ├── record_section_az327.png
│   ├── record_section_6stations.png
│   ├── record_section_3000s.png
│   └── K60A_3C_comparison.png, K60A_nearby_3000s.png
```

The top-level `oi/` (`/mnt/f/sharedate/oi/`) also holds:
- `twave_arrival_map.pdf` — original GMT arrival map
- `twave_record_section.pdf` — original record section
- `all_polar_analysis/` — earlier polarization batch output

## Path Conventions (Legacy)

All scripts in this directory reference shared paths, NOT event-isolated paths:

```python
PROC_DIR = "/home/sy/shared/case/sacprocess"   # Shared across ALL events
RAW_DIR  = "/home/sy/shared/case/raw"           # Shared raw data
OUT_DIR  = "/home/sy/shared/case/waveform_contrast"  # Or waveform_contrast/one/
```

This means the scripts in this dir read from the **same** `sacprocess/` that also contains
data from other events processed through the shared `case/raw/`. Any new preprocessing
for this event MUST go through the shared `case/sacprocess/`, NOT create a local
`sacprocess/` inside the event dir.

## SAC File Pattern

```python
# BHZ proc files
f'{PROC_DIR}/TA.{sta}..BHZ.M.2014.013.035104_proc.sac'

# BHN/BHE (for 3C)
f'{PROC_DIR}/TA.{sta}..BHN.M.2014.013.035104_proc.sac'
f'{PROC_DIR}/TA.{sta}..BHE.M.2014.013.035104_proc.sac'

# IU network
f'{PROC_DIR}/IU.HRV.00.BHZ.M.2014.013.035104_proc.sac'

# Raw files (counts)
f'{RAW_DIR}/TA.{sta}..BHZ.M.2014.013.035104.SAC'
```

**Date code**: `2014.013.035104` = Year 2014, Julian day 013, 03:51:04 UTC.

## Key Publication Figures

| Figure | Script | Description |
|--------|--------|-------------|
| Figure 1 | `plot_twave_map_labels.py` (at `oi/`) | Station distribution map with labels |
| Figure 2a | `plot_fig2a.py` / `plot_fig2a_all.py` (at `case/code/`) | Two-column record section: AZ<341° left, AZ≥341° right, 1200-2160s window |

plot_fig2a.py was the precursor to `codex/scripts/06_profile_section.py`. The figure2a workflow was:
1. Read all proc.sac files from shared `sacprocess/`
2. Split by azimuth (341° threshold)
3. Left panel: ~232 stations (near-side azimuth)
4. Right panel: ~118 stations (far-side azimuth, 350° group)
5. Time window: 1200-2160s (20-36 minutes), x-axis in seconds
6. Hilbert envelope + P98 norm, every 4th station subsampled

# T-Wave Research Planning: Event & Station Selection

## USGS FDSN Event API

Query format for finding suitable T-wave source events:

```
curl -s "https://earthquake.usgs.gov/fdsnws/event/1/query?\
starttime=YYYY-MM-DD&endtime=YYYY-MM-DD&\
minlatitude=-60&maxlatitude=-30&\
minlongitude=0&maxlongitude=70&\
minmagnitude=5.5&format=text&orderby=magnitude"
```

### Key Parameters for T-Wave Events
- **Magnitude**: ≥5.5 for observable T waves at regional distances
- **Depth**: ≤15 km preferred (shallow events couple better to SOFAR channel)
- **Location**: Must be oceanic (mid-ocean ridges, transforms, oceanic intraplate)
- **Path**: Ocean path to coast required for T-P/T-S conversion

## IRIS FDSN Station API

Query format for finding available stations:

```
# Networks in a region
curl -s "https://service.iris.edu/fdsnws/station/1/query?\
minlatitude=-38&maxlatitude=-15&minlongitude=10&maxlongitude=40&\
format=text&level=network"

# Stations in a network
curl -s "https://service.iris.edu/fdsnws/station/1/query?\
network=AF&format=text&level=station"

# Stations in a region (all networks)
curl -s "https://service.iris.edu/fdsnws/station/1/query?\
minlatitude=-38&maxlatitude=-22&minlongitude=12&maxlongitude=38&\
format=text&level=station"
```

## Southern Africa Networks (as of 2026)

| Code | Name | Stations | Period | Notes |
|------|------|----------|--------|-------|
| AF | AfricaArray | 67 | 2004-present | Permanent, South Africa & neighbors |
| BX | Botswana Seismological Network | 21 | 2001-present | Permanent, Botswana |
| 6H | Mozambique Rift Tomography | 38 | 2011-2013 | Temporary, east coast |
| 8A | AfricaArray Namibia | 19 | 2015-2018 | Temporary, southwest coast |
| ZT | Bushveld Complex | 24 | 2015-2020 | Temporary, northern SA |
| MZ | National Seismic Network of Mozambique | 3 | 2022-present | Permanent |
| II | GSN | SUR | 1990-present | Sutherland, South Africa |
| XA | Kaapvaal Craton | 82 | 1997-1999 | Zhou et al. (2021) used this |

## T-Wave Conversion Location Requirements

Converted T waves (T-P, T-S) are typically observed:
- **Distance from coast**: ≤4° inland (Buehler & Shearer, 2015)
- **Conversion zone**: Continental slope/shelf where ocean acoustic waves couple to solid earth
- **Coast orientation matters**: East coast receives Indian Ocean T waves, west/southwest coast receives Atlantic T waves

## Event Selection Criteria for Validation Studies

When selecting events to validate Zhou et al. (2021) methodology:
1. **Different event**: Not the same Mw 6.6 Prince Edward Islands event (1997-04-28)
2. **Different stations**: Use permanent AF network or other temporary arrays, not XA
3. **Similar path characteristics**: Ensure ocean-to-continent conversion path exists
4. **Signal strength**: Larger magnitude → better SNR for identification
5. **Multiple events**: Consider using 2-3 events for robustness

## Regional Event Catalog (M≥6.0, shallow, 2000-2026)

### Prince Edward Islands Region (Indian Ocean)
- 2023-05-21: Mw 6.8, -43.4°N, 39.3°E, 28 km
- 2019-01-22: Mw 6.7, -43.1°N, 42.4°E, 13 km
- 2020-08-06: Mw 6.3, -41.7°N, 42.4°E, 10 km
- 2024-05-31: Mw 6.2, -43.3°N, 41.7°E, 10 km

### Bouvet Island Region (Atlantic/SW Indian Ridge)
- 2014-04-15: Mw 6.8, -53.5°N, 8.7°E, 11 km
- 2017-10-10: Mw 6.7, -54.3°N, 8.6°E, 9 km
- 2018-01-28: Mw 6.6, -53.1°N, 9.7°E, 10 km
- 2003-04-17: Mw 6.5, -54.6°N, 1.4°E, 10 km

### South of Africa (SW Indian Ridge)
- 2024-07-10: Mw 6.6, -53.3°N, 25.3°E, 4 km
- 2004-09-28: Mw 6.4, -52.5°N, 28.0°E, 10 km
- 2023-12-22: Mw 6.1, -52.1°N, 28.0°E, 10 km
- 2024-01-20: Mw 6.2, -39.9°N, 46.2°E, 10 km (Southwest Indian Ridge)

## References

- Zhou et al. (2021): Determining Crustal Attenuation With Seismic T Waves in Southern Africa. GRL.
- Buehler & Shearer (2015): T phase observations in global seismogram stacks. GRL.
- de Groot-Hedlin (2020): Seismic T-wave observations at dense seismic networks. SRL.

# US East Coast T-Wave Candidate Events (2014–2016)

Catalog of M6+ earthquakes in the Atlantic/Caribbean region suitable for US East Coast
T-wave crustal attenuation analysis. Each event's back-azimuth (BAZ) relative to the US
East Coast is computed to assess path diversity.

## Reference BAZ Values

BAZ computed via `obspy.geodetics.gps2dist_azimuth(stla, stlo, evla, evlo)` to four
representative TA stations on the US East Coast:

| Station | Lat | Lon | Region |
|---------|-----|-----|--------|
| J59A | 43.46°N | 74.50°W | New York |
| H57A | 41.81°N | 77.65°W | Pennsylvania |
| N63A | 44.54°N | 71.56°W | New Hampshire |
| Y59A | 39.04°N | 76.82°W | Maryland / DC |

## Candidate Events (sorted by BAZ diversity)

| Event | Date | Mw | Lat | Lon | Depth | BAZ(J59A) | BAZ(H57A) | BAZ(N63A) | BAZ(Y59A) | Description |
|-------|------|----|-----|-----|-------|-----------|-----------|-----------|-----------|-------------|
| PR M6.3 | 2014-01-13 | 6.4 | 19.04°N | 66.81°W | 20km | **163°** | **155°** | **170°** | **154°** | Puerto Rico *(already processed)* |
| MAR M6.0 | 2014-07-27 | 6.0 | 23.72°N | 45.58°W | 10km | **120°** | **114°** | **126°** | **111°** | Central MAR *(already processed)* |
| Barbados M6.5 | 2014-02-18 | 6.5 | 14.67°N | 58.93°W | 15km | 150° | 144° | 156° | 142° | Caribbean |
| MAR-N M7.1 | 2015-02-13 | **7.1** | 52.65°N | 31.90°W | 17km | **57°** | **55°** | **58°** | **51°** | Northern MAR ★ best |
| Barbados M6.5 | 2015-07-16 | 6.5 | 13.87°N | 58.55°W | 20km | 150° | 144° | 156° | 142° | Caribbean (same path as Feb 2014) |
| MAR M6.0 | 2015-09-18 | 6.0 | 15.28°N | 45.99°W | 10km | 130° | 125° | 135° | 123° | Central MAR |
| Antigua M6.0 | 2016-03-19 | 6.0 | 18.00°N | 60.70°W | 26km | 151° | 143° | 158° | 142° | Caribbean |
| MAR M6.1 | 2016-06-21 | 6.1 | 22.66°N | 45.14°W | 10km | 121° | 115° | 126° | 112° | Central MAR (similar to 2014 MAR) |

## Top Recommendations by Azimuth Group

### Group A: NE Approach (BAZ 50-60°) — **Most Different from PR**
**2015-02-13 M7.1 Northern Mid-Atlantic Ridge** (52.65°N, 31.90°W)
- BAZ ~55° → wave arrives from NE, enters North America via Newfoundland/Quebec
- **Land path**: crosses Canadian Shield → New England, completely different crustal
  structure from the Puerto Rico event's Appalachian/southern path
- **M7.1** → strongest T-wave generation of any candidate
- **Why useful**: orthogonal azimuth (~100° difference from PR), tests crustal attenuation
  across the Canadian Shield vs Appalachian province

### Group B: ESE Approach (BAZ 110-130°)
**2014-07-27 M6.0 Central MAR** (already have this)
**2016-06-21 M6.1 Central MAR** (22.66°N, 45.14°W)
- BAZ ~115-125° → wave arrives from the east, hits the coast near Carolinas
- **Land path**: crosses Piedmont → Southern Appalachians

### Group C: SSE Approach (BAZ 140-170°)
**Puerto Rico 2014-01-13 M6.3** (already have)
**Barbados/Caribbean events** (BAZ ~142-158°, similar paths to PR)

## Key Selection Criteria

1. **Azimuth diversity**: events in Groups A & B provide the most independent paths
   from the already-processed Puerto Rico event (Group C)
2. **Magnitude**: M6.5+ preferred for reliable T-wave detection across 3000+ km
3. **Data availability**: IRIS FDSN (service.iris.edu) has these events; NA network
   (TA/US/IU/LD) stations along the US East Coast provide SACPZ files

## For Crustal Attenuation Analysis

Compare the ln(A) vs inland-distance slope between:

- **PR event (Group C)**: path through Southern Appalachians / Atlantic Coastal Plain
- **MAR-N 2015 M7.1 (Group A)**: path through Canadian Shield / New England
- **MAR Central (Group B)**: path through Piedmont / Coastal Plain

Each path crosses different crustal provinces with distinct Q values, providing
natural experiments for T-wave attenuation tomography.
