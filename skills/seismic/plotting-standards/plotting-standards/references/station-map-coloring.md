# 台站地图着色方法

---
name: station-map-coloring
description: Station distribution map coloring methods — tt_t, wavefront model, dual conversion point, STA/LTA detection, and pitfalls
version: 2.0
tags: [seismic, twave, map, coloring]
---

# 台站分布图着色方法

## 检测方法对比

| 方法 | 脚本 | 原理 | 站数 | 范围 | 可靠性 |
|------|------|------|------|------|--------|
| 包络峰值 (arrivals.py) | 02_extract_stations | np.argmax(env) ±120s | 350 | - | ❌ 噪声峰 |
| STA/LTA | 02_stalta_detect.py | STA/LTA on env, 5/500s | 290 | 22.4-31.8 | ✅ DGH 一致 |
| 包络初动 20% | 03_onset_detect.py (v1) | 20% onset | 316 | 18.4-29.3 | ⚠️ celerity 敏感 |
| **峰值+空间平滑** | **03_onset_detect.py (v2)** | **peak + 局部中值 3°半径** | **336** | **19.2-30.9** | **✅ 推荐** |

**⚠️ 用户原则**：不要绕过问题找替代方案。先证明计算正确性，再调整可视化。

## 当前步1脚本清单 (2026-06-29 清理后)

| 脚本 | 用途 |
|---|---|
| `station_map_theoretical.py` | **主画图脚本** — 读JSON(STALTA检测) → 双面板分布图, hot色标, 全部站标注 |
| `02_stalta_detect.py` | STA/LTA检测 → 生成JSON + Fig.1/2a/2b |
| `arrivals.py` | 到时计算函数库 (compute_theoretical_arrivals等) |
| `traveltime_profile.py` | get_track() + tt_t() |

项目根还有 `02_extract_stations.py` (理论到时提取, 输出Excel+gmt_stations)。

## 六种方法

| 方法 | 脚本 | T-P 来源 | 范围 | 差异 |
|------|------|---------|------|------|
| 转换点CSV | station_map_theoretical.py | step4 CSV tp_s | 16-29 min | 基准 |
| tt_t实时 | station_map_tt_t.py | tt_t() 逐站 | 16-29 min | <2s |
| 匀速公式 | station_map_simple.py | ocean/1.48+land/6.0 | 17-29 min | <2s |
| 双转换点 | station_map_dual_conv.py | 327+350分组 | 20-29 min | 中等 |
| 波前模型 | test.py (DGH风格) | 方位角函数+Love 3.6 | 21-31 min | 显著 |
| STA/LTA | plot_twave_observed_v2.py | 1-5Hz STA/LTA | 22-32 min | 最大 |

## 陷阱

- 前三种方法视觉几乎相同(差<2s)，非缓存bug——MD5验证过
- gmt_stations.txt存的是实测(t_obs)，不可靠。画理论图从Excel T理论列或CSV读
- 用户不接受换陆地距着色绕过问题——必须先验证计算正确
- 跨方位角T-P排序: ocean_time占94%，312→FL(1536km)<336→NC(1943km)<350→NE(2352km)。内陆FL站可能比岸边NC站更早——物理正确非bug
- **GMT 15s地形数据损坏**: S60E000等tile可能是HTML文件而非JP2（下载失败）。症状：`ERROR 4: ...jp2 not recognized as being in a supported file format`。修复：删除损坏文件后重新下载 `rm -f ~/.gmt/server/earth/earth_relief/earth_relief_15s_p/S60E000*`，然后重新运行脚本会自动下载
- **方位角线不要画太多**: 用户明确要求只画2条（最近/最远海岸），不要每30°一条。这是视觉偏好，不是技术限制
- **震源机制解留空**: 如果没有震源机制解数据，不要跳过，留个位置以后补上

## 波前模型 (test.py, D:\\uu)

DGH 2020方法。海洋距=方位角线性函数:
- VO=1.48, VL=3.60 (Love wave速度)
- az<327: Dc=1800km
- az>350: Dc=2350km
- 327-350: Dc=1800+(az-327)/23*550
- T-P=Dc/VO+land/VL (分钟)

使用pyproj.Geod(非geographiclib)。输出gmt_stations_front.txt。
参数经用户调整(原2000/2500偏高2-3min，下压至接近W61A实际海岸线)。
327°转换点~32.46N/76.74W, 350°~39.44N/71.92W。

## STA/LTA检测 (02_stalta_detect.py, 源自D:\\uu\\plot_twave_observed_v2.py)

原始观测到时方法。比arrivals.py的包络峰值法更可靠(使用STA/LTA而非np.argmax)。
流程:
原始SAC→1-5Hz BP→Hilbert包络→STA/LTA(5s/500s)
→15-38min窗口找STA/LTA峰值→Celerity<1.8km/s→空间离群过滤(5近邻中值±3min)
→输出t_wave_detections_final.json

结果: 343站→311检测→290保留。到时22.4-31.8min(与DGH 2020 Fig.2a 24-32min吻合)。
同时生成Fig.1地图(含327/350转换边界弧+2min等时线)、Fig.2a记录剖面、Fig.2b振幅图。
JSON格式: {n_total, n_detected, n_removed, stations[{net,lon,lat,dist,az,tt_min,celerity,max_amp}]}

## 包络峰值+空间平滑检测 (03_onset_detect.py v3)

在STA/LTA和包络初动之间——峰值检测 + 局部中值空间平滑(3°半径, >2min调整):
- 336站, 19.2-30.9min, 仅调整10站
- 输出 `t_wave_peak_smooth.json`
- 画图: `station_map_theoretical.py` 改JSON路径

**空间平滑用局部中值，不要用RBF**。RBF(gaussian)在边缘外推产生负值(-34.8min)和极大值(102min)。

## 制图风格偏好 (2026-07-23 南非事件最终定稿)

**用户对台站分布图的明确偏好**:

1. **台站名标注**: SNR>6dB 才标注。6pt Helvetica-Bold，白色底框（`fill='white@70'`）+灰色边框（`pen='0.2p,gray50'`）。固定偏移 `offset='0.1c/-0.1c'`（右上方）
2. **方位角参考线**: **只画2条**——最近海岸台站和最远海岸台站的方位角。两条都画**深灰色**（`gray30`），不区分颜色。延伸到画布边缘，在外端标注角度值
3. **区域padding**: **3度**（不是5度或10度），pad=3足够
4. **地形**: 15s（@earth_relief_15s），色标 `geo`（蓝绿海洋+棕色陆地）。全球数据已缓存（~/.gmt/server/earth/earth_relief/earth_relief_15s_p/，2.7GB）
5. **海洋覆盖**: 80%透明度白色（`water='white@80'`）
6. **水深1000m等深线**: 黑色虚线（`pen='0.8p,black,-'`），**必须放在coast之后**绘制，才能显示在白色覆盖层上面
7. **震中距圆弧**: 每500km一条灰色虚线（`pen='0.5p,gray40,--'`），东侧标注距离值
8. **SNR色标**: viridis，`math.floor(min)` 到 `math.ceil(max)` 整数边界
9. **震源机制解**: 有数据画沙滩球，无数据**留空**不画
10. **地震信息标注**: 左上角白底黑字框（`fill='white@70'`），含：发震时间UTC、震级、深度、震中经纬度
11. **坐标轴标签**: 必须有 "Longitude" 和 "Latitude"
12. **通用脚本**: `~/codex/scripts/plot_station_map.py` 从 `events.yaml` 读取参数，`figure_map.png` 输出。用法: `python3 plot_station_map.py <event>`

**示例代码片段**:

```python
# 台站名偏移
def label_offset(lon, lat, az):
    if az < 45:       return (0.5, 0.35)
    elif az < 90:     return (0.35, 0.3)
    elif az < 135:    return (0.3, 0.0)
    elif az < 180:    return (0.3, -0.3)
    elif az < 225:    return (-0.6, -0.3)
    elif az < 270:    return (-0.6, -0.2)
    elif az < 315:    return (-0.6, 0.2)
    else:             return (-0.6, 0.35)

# 方位角参考线: 只画2条（最近/最远海岸），不要每30°一条
sorted_by_land = sorted(stations, key=lambda s: s['land'])
min_land_sta = sorted_by_land[0]  # 最近海岸
max_land_sta = sorted_by_land[-1] # 最远海岸

def plot_azimuth_line(fig, az, pen_style, label):
    """绘制方位角线并标注"""
    line_lons, line_lats = [], []
    for dkm in range(0, 6000, 50):
        g = geod.Direct(EVLA, EVLO, az, dkm * 1000)
        lo, la = g['lon2'], g['lat2']
        if REGION[0]-2 <= lo <= REGION[1]+2 and REGION[2]-2 <= la <= REGION[3]+2:
            line_lons.append(lo)
            line_lats.append(la)
        else:
            break
    if len(line_lons) > 2:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tf:
            for lo, la in zip(line_lons, line_lats):
                tf.write(f'{lo} {la}\n')
            tn = tf.name
        fig.plot(data=tn, pen=pen_style)
        os.unlink(tn)
        # 标注方位角（在圆弧外端）
        label_dist = 4800  # km
        g = geod.Direct(EVLA, EVLO, az, label_dist * 1000)
        if REGION[0] <= g['lon2'] <= REGION[1] and REGION[2] <= g['lat2'] <= REGION[3]:
            fig.text(x=[g['lon2']], y=[g['lat2']], text=f'{az:.0f}°',
                     font='8p,Helvetica-Bold,black', fill='white@60', pen='0.3p')

# 最近/最远海岸方位角（深灰色，两条同色）
plot_azimuth_line(fig, min_land_sta['az'], '1.2p,gray30,-', 
                  f"{min_land_sta['sta']} ({min_land_sta['land']:.0f}km)")
plot_azimuth_line(fig, max_land_sta['az'], '1.2p,gray30,-', 
                  f"{max_land_sta['sta']} ({max_land_sta['land']:.0f}km)")

# 震中距圆弧（每500km一条）
for rr in [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000]:
    radius_deg = rr / 111.12
    fig.plot(x=[EVLO], y=[EVLA], style=f'c{radius_deg:.1f}c', pen='0.5p,gray40,--')
    # 标注距离（放在右侧）
    g = geod.Direct(EVLA, EVLO, 90, rr * 1000)
    if REGION[0] <= g['lon2'] <= REGION[1] and REGION[2] <= g['lat2'] <= REGION[3]:
        fig.text(x=[g['lon2']], y=[g['lat2']], text=f'{rr}',
                 font='7p,Helvetica,gray30', angle=0, justify='LM')

# 地震信息标注（左下角）
info_text = (
    f"{EVENT_DATE} {EVENT_TIME}\\n"
    f"Mw {EVENT_MW}  Depth {EVDP:.0f} km\\n"
    f"{abs(EVLA):.2f}°{'S' if EVLA<0 else 'N'} {abs(EVLO):.2f}°{'E' if EVLO>0 else 'W'}"
)
fig.text(x=[REGION[0]+1], y=[REGION[3]-1.5], text=info_text,
         font='10p,Helvetica-Bold,black', fill='white@70', pen='0.5p,black', justify='TL')

# 震源机制解（如果有数据）
# fig.meca(spec={'strike1': ..., 'dip1': ..., 'rake1': ..., ...},
#          longitude=EVLO, latitude=EVLA, depth=10,
#          scale='0.5c', compression_fill='red', extension_fill='white', outline='0.5p,black')
# 如果没有数据，留空不画
```

## 彩色地形图

使用 pyGMT + ETOPO1 绘制带等深线/等高线的彩色地形台站分布图。见 `references/colored-topographic-map.md`。

### 增强对比度参数

彩色地形在强阴影下等高线可能被掩盖，需要调整 shading 强度：
- 标准：`shading='+a45+nt0.5'`（灰色地形够用）
- 增强：`shading='+a30+nt1.5'`（彩色地形需要，nt增大3倍）

等高线在彩色地形上也需要加粗加亮：
- 标准：`pen='0.2p,blue'`（灰色底图）
- 增强：`pen='0.4p,cyan'`（彩色底图，青色/橙色比蓝色/棕色更醒目）

## 脚本整合模式

当一个目录有多个相似的地图脚本（如 station_map_tlove.py、station_map_tp.py、station_map_theoretical.py），
应整合为一个脚本，通过命令行参数选择模式：

```python
# station_map.py — 整合版本
MODES = {
    'tlove_topo': {'values': t_love, 'label': 'T-Love travel time (min)', 'topo': 'etopo1', 'contours': True, ...},
    'tp':         {'values': t_tp,   'label': 'T-P arrival time (min)',   'topo': 'gray',   'contours': False, ...},
    'theoretical':{'values': dist,   'label': 'Time after event (min)',   'topo': 'gray',   'contours': False, ...},
}

# 用法: python station_map.py tlove_topo | tp | theoretical | all
```

关键差异通过配置字典管理：地形着色(cmap)、色标范围(series)、标签(label)、等高线(contours)、
shading参数、transparency。其他元素（海岸线、台站、震中、方位角线）完全不变。

**优势**：3个脚本→1个，减少代码重复，统一维护。

## SNR 着色标准（2026-07-23 南非事件最终定稿）

台站色标按 **SNR** 着色（不用 T-P 走时），边界取整数：

```python
import math
snr_values = [s['snr'] for s in stations]
snr_cpt_min = math.floor(min(snr_values))   # 下限取左边整数（如0.0→0）
snr_cpt_max = math.ceil(max(snr_values))     # 上限取右边整数（如19.1→20）
pygmt.makecpt(cmap='viridis', series=[snr_cpt_min, snr_cpt_max, 1])

# 色标标签
fig.colorbar(frame='xa2f1+l"SNR (dB)"',
             position='jBL+w6c/0.25c+h+o1c/0.6c+ml')
```

**标签阈值：** 仅显示 SNR>6dB 的台站名，6pt Helvetica-Bold，白色底框（`fill='white@70'`）+灰色边框（`pen='0.2p,gray50'`），固定偏移 `offset='0.1c/-0.1c'`。

**可调参数脚本模式（推荐）：** 将台站大小、字体、透明度等放在脚本顶部作为变量，方便批量修改：

```python
# 可调参数
STATION_SIZE  = 0.45        # 台站三角形大小 (cm)
LABEL_FONT    = "8p,Helvetica-Bold,black"
LAND_ALPHA    = "@80"       # 陆地白色透明度
MAP_WIDTH     = "16c"       # 地图宽度
DPI           = 200
```

**与旧方法对比：**
| 项目 | 旧方法 (T-P) | 新方法 (SNR) |
|------|-------------|-------------|
| 着色值 | `tt_tp / 60` (min) | `snr` (dB) |
| 色标 | `hot` reversed | `viridis` |
| 边界 | 自动范围 | floor/ceil 整数 |
| 含义 | 走时早晚 | 信号质量高低 |


# 彩色地形图台站分布图 (2026-07-09, 2026-07-23 更新)

## 方法

使用 pyGMT + ETOPO1 地形数据绘制彩色地形图，包含等深线/等高线。

### 关键代码模式

```python
import pygmt

fig = pygmt.Figure()

# 地形底图（15s高精度，已下载全球数据）
fig.grdimage(
    grid="@earth_relief_15s",  # 或 @earth_relief_01m 用于快速迭代
    region=region,
    projection="M15c",
    cmap="geo",              # 或 'etopo1'（更彩色）
    shading="+a45+nt1.0",    # 中等增强阴影
)

# 海洋覆盖（可选透明度）
fig.coast(
    shorelines="0.5p,black",
    borders=["1/0.3p,gray50", "2/0.2p,gray70"],
    resolution="h",
    water="white@80",        # 80%透明白色覆盖
)

# 水深等深线（用负值！）
fig.grdcontour(
    grid="@earth_relief_15s",
    levels=[-1000],          # 水深1000m，用负值！
    pen="0.8p,gray50,-",     # 灰色虚线
)
```

### Colormap 选择 (2026-07-23 更新)

| Colormap | 海洋 | 陆地 | 风格 |
|----------|------|------|------|
| `etopo1` | 蓝色（浅→深） | 绿/棕/黄 | 彩色地形，视觉丰富 |
| `geo` | 蓝绿色 | 棕色/褐色 | 经典地图风格，对比清晰 |

**用户偏好**: `geo` 更专业干净，推荐用于出版图。

### 关键参数

#### 海洋覆盖透明度
```python
water="white@5"    # 5%透明，几乎不覆盖
water="white@80"   # 80%透明，轻覆盖（用户偏好）
water="lightblue@10"  # 浅蓝色覆盖
# 不设置water参数：原始地形颜色
```

#### 等深线（海洋深度用负值！）
```python
# 正确 — 负值表示海洋深度
fig.grdcontour(grid="@earth_relief_15s", levels=[-1000], pen="0.8p,gray50,-")
fig.grdcontour(grid="@earth_relief_15s", levels=[-2000, -4000], pen="0.3p,cyan")

# 错误 — 正值是陆地海拔
fig.grdcontour(grid="@earth_relief_15s", levels=[1000], ...)  # 这是1000m海拔，不是水深！
```

#### 阴影增强
```python
shading="+a45+nt0.5"   # 标准（灰色地形够用）
shading="+a30+nt1.5"   # 增强（彩色地形需要，nt增大3倍）
shading="+a45+nt1.0"   # 中等增强（用户偏好）
```

### 台站着色（按距离）

```python
# 写临时文件（lon, lat, color_value）
tmp_xy = "/tmp/stations_xy.txt"
df[['lon', 'lat', 'dist']].to_csv(tmp_xy, sep=' ', header=False, index=False)

pygmt.makecpt(cmap="rainbow", series=[f"{df['dist'].min():.0f}/{df['dist'].max():.0f}/100"])

fig.plot(
    data=tmp_xy,
    style="t0.12c",
    cmap=True,
    pen="0.2p,black",
)
```

### ⚠️ pyGMT v0.18 陷阱

**不能用 `color=` 参数**，必须用 `data=` + 临时文件 + `cmap=True`：

```python
# ❌ 错误（v0.18会报错 GMTInvalidInput）
fig.plot(x=lon, y=lat, color=dist, style="t0.1c", cmap=True)

# ✅ 正确
df[['lon', 'lat', 'dist']].to_csv(tmp_xy, sep=' ', header=False, index=False)
fig.plot(data=tmp_xy, style="t0.1c", cmap=True)
```

### 等深线标注陷阱

标注值不能为负数：

```python
# ❌ 错误（annotation不能是负数）
fig.grdcontour(levels=[-200, -1000, ...], annotation="-2000+f8p,blue")

# ✅ 正确（用正数，GMT自动处理海洋负值）
fig.grdcontour(levels=[200, 1000, 2000, ...], annotation="2000+f6p,blue")
```

### 增强对比度（彩色地形必须）

彩色地形（etopo1）的阴影和颜色会掩盖等高线，需要增强参数：

```python
# 1. shading 强度增大3倍
fig.grdimage(..., cmap='etopo1', shading='+a30+nt1.5')  # nt从0.5→1.5

# 2. 等高线加粗+高对比度颜色
fig.grdcontour(..., pen='0.4p,cyan')    # 海洋：cyan比blue更醒目
fig.grdcontour(..., pen='0.4p,orange')  # 陆地：orange比brown更醒目

# 3. 海洋区域加白色半透明遮罩
fig.coast(..., water='white@30')
```

| 场景 | shading | 等深线 | 等高线 | transparency |
|------|---------|--------|--------|-------------|
| 灰色地形 | `+a45+nt0.5` | `0.2p,blue` | `0.2p,brown` | 50 |
| **彩色地形** | **`+a30+nt1.5`** | **`0.4p,cyan`** | **`0.4p,orange`** | **0** |

### 脚本整合模式

多个相似地图脚本应整合为一个，通过参数选择模式：

```python
MODES = {
    'tlove_topo': {'values': t_love, 'topo': 'etopo1', 'contours': True, 'shading': '+a30+nt1.5', ...},
    'tp':         {'values': t_tp,   'topo': 'gray',   'contours': False, 'shading': '+a45+nt0.5', ...},
}
# python station_map.py tlove_topo | tp | all
```

**优势**：3个脚本→1个，减少代码重复，统一维护。差异通过配置字典管理。

### 输出

- 文件：`figure_map_color_topo.png`
- DPI：300
- 大小：~5-8 MB（取决于区域）

## 快速改造：灰色地形→彩色地形+等高线

当原脚本已用 `cmap='gray'` 画地形时，只需两步改造：

```python
# 1. 把 cmap='gray' 改成 cmap='etopo1'
fig.grdimage(grid='@earth_relief_15s', ..., cmap='etopo1', shading='+a30+nt1.5')

# 2. 在 grdimage 之后、coast 之前加等深线/等高线
# 等深线（海洋）— 青色（彩色底图上更醒目）
fig.grdcontour(grid='@earth_relief_15s', levels=[2000, 4000], pen='0.4p,cyan',
               annotation='4000+f6p,cyan')
# 等高线（陆地）— 橙色
fig.grdcontour(grid='@earth_relief_15s', levels=[500, 1000], pen='0.4p,orange',
               annotation='1000+f5p,orange')

# 3. coast 加 water='white@20' 提高海洋区域可读性
fig.coast(..., water='white@20')
```

**其他元素（台站着色、色标、震中等）完全不变。**

**注意**：WSL环境下需要设置代理才能下载GMT数据：
```bash
export http_proxy=socks5://127.0.0.1:10808
export https_proxy=socks5://127.0.0.1:10808
```

示例脚本：`plot_station_map_final.py`（南非Mw6.3事件，完整版含方位角线+震中距圆弧+地震信息）

## 台站名字标注 (2026-07-23)

```python
# 带白色底框，过滤低SNR
for sta in stations:
    if sta['snr'] > 6:  # 只显示高SNR台站
        fig.text(
            x=[sta['lon']], y=[sta['lat']],
            text=sta['sta'],
            font='6p,Helvetica-Bold,black',
            justify='TL', offset='0.1c/-0.1c',
            fill='white@70',    # 白色底框
            pen='0.2p,gray50'   # 细边框
        )
```

## SNR色标 — 整数边界 (2026-07-23)

```python
import math
snr_cpt_min = math.floor(min(snr_vals))  # 0.0 → 0
snr_cpt_max = math.ceil(max(snr_vals))   # 19.1 → 20
pygmt.makecpt(cmap='viridis', series=[snr_cpt_min, snr_cpt_max, 1])
```

## 常见错误

1. **等深线用正值**: 正值是陆地海拔，海洋深度必须用负值
2. **忘记proxy**: 15s数据下载需要代理，否则超时
3. **S60E000损坏**: 这个tile可能是HTML文件，删除后重新下载
4. **geo色标灰色**: 如果geo显示灰色，可能是tile损坏，不是色标问题
5. **water透明度**: `white@80` 是80%透明度（即20%不透明），不是80%不透明

# GMT 地形数据下载与代理设置 (2026-07-23)

## 代理设置

WSL环境下下载GMT数据需要代理。v2rayN混合端口默认10808：

```bash
export http_proxy=socks5://127.0.0.1:10808
export https_proxy=socks5://127.0.0.1:10808
```

测试代理是否通：
```bash
curl -I https://oceania.generic-mapping-tools.org
```

## 全球15s地形数据下载

GMT 15s数据按10°×10°分tile，全球约648个tile，总计约3GB。

### 下载策略

分区域下载，避免单次请求过大：

```python
import pygmt

regions = [
    ([-180, 180, -90, -30], 'S90-S30'),
    ([-180, 180, -30, 30], 'S30-N30'),
    ([-180, 180, 30, 90], 'N30-N90'),
]

for region, name in regions:
    fig = pygmt.Figure()
    fig.grdimage(grid='@earth_relief_15s', region=region, projection='Q15c', cmap='etopo1')
    print(f'{name} 完成')
```

### 缓存位置

```
~/.gmt/server/earth/earth_relief/earth_relief_15s_p/
```

GMT会自动使用本地缓存，无需额外配置。

## GMT 15s数据损坏修复

### 症状

```
ERROR 4: `/home/sy/.gmt/server/earth/earth_relief/earth_relief_15s_p/S60E000.earth_relief_15s_p.jp2' not recognized as being in a supported file format.
```

### 原因

某些tile下载失败，保存的是HTML错误页面而非JP2文件。

### 诊断

```bash
file ~/.gmt/server/earth/earth_relief/earth_relief_15s_p/S60E000*
# 如果输出 "HTML document, ASCII text" 而非 "JPEG 2000" 则是损坏文件
```

### 修复

```bash
# 1. 删除损坏文件
rm -f ~/.gmt/server/earth/earth_relief/earth_relief_15s_p/S60E000*

# 2. 重新运行脚本会自动下载
python3 your_script.py
```

## 15s vs 01m 对比

| 精度 | 数据源 | 文件大小 | 适用场景 |
|------|--------|----------|----------|
| 15s (~460m) | @earth_relief_15s | ~3GB (全球) | 出版级地图 |
| 01m (~1.8km) | @earth_relief_01m | ~300MB (全球) | 快速迭代 |
| 02m (~3.7km) | @earth_relief_02m | ~100MB | 概览图 |

**建议**：开发时用01m，最终出图用15s。

# STA/LTA T-wave Detection Method

From `plot_twave_observed_v2.py` (D:\uu). More robust than simple envelope peak detection.

## Detection Pipeline

```
raw SAC → 1-5 Hz bandpass (4th order, zerophase)
        → Hilbert envelope
        → STA/LTA (STA=5s, LTA=500s)
        → search window 15-38 min
        → find max STA/LTA peak
        → celerity check: dist / (t_peak * 60) < 1.8 km/s
        → spatial outlier removal: compare to 5 nearest neighbors
          if |tt - median(nn5)| > 3 min → reject
```

## Results (PR2014)

- Input: 343 stations
- Raw detection: 311
- After spatial filtering: 290
- Arrival time: 22.4-31.8 min
- 327° profile: 34 stations, 350° profile: 32 stations
- DGH 2020 reference: Figure 2a shows 24-32 min range — match!

## Output (t_wave_detections_final.json)

```json
{
  "n_total": 343,
  "n_detected": 290,
  "n_removed": 21,
  "stations": [
    {
      "net": "IU.CCM",
      "lon": -91.2446,
      "lat": 38.0557,
      "dist": 3169.6,
      "az": 317.1,
      "tt_min": 29.47,
      "celerity": 1.793,
      "max_amp": 167.6
    }
  ]
}
```

## Key Differences vs arrivals.py

| | arrivals.py | plot_twave_observed_v2.py |
|---|---|---|
| Detector | `np.argmax(envelope)` | STA/LTA peak on envelope |
| Bandpass | None | 1-5 Hz |
| Window | theory_s ± 120s | Fixed 15-38 min |
| Spatial filter | None | 5-NN median ± 3 min |
| Failures | L63A: theory 1624s, obs 1510s (113s error) | L63A likely rejected as outlier |

## Script Location

- Original: `D:\uu\plot_twave_observed_v2.py`
- Copy in step 1: `process_flow/1_station_distribution/02_stalta_detect.py`

Paths to modify before running:
- `RAW_DIR` → raw SAC directory
- `OUT_DIR` → step 1 output directory
