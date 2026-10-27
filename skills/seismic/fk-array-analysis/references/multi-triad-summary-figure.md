# 多三角形 F-K 结果图：热力图 + 单张综合地图

马达加斯加 2013-08-17 Mw6.2 实战定稿（2026-08-27）。用户对 F-K 结果表达的两种
偏好形态，脚本模板在事件目录 `process_flow/4_fk_analysis/`：
`fk5_multi_triads.py`（热力图）、`fk6_summary_map.py`（综合地图，`--comp Z|T`）。

## 1. 多三角热力图（fk5_multi_triads.py）

- 选 6 个三角形覆盖南→北（Delaunay 候选 + 区域分散），Z 分量、2-8 Hz、直接窗 F-K
- matplotlib `3×2` subplots（不要 pygmt subplot——有 grdimage bug）
- 每面板 imshow 慢度图：`extent=[-35,35,-35,35]` s/deg、`origin='lower'`、
  `cmap='turbo'`、`vmin=-3 vmax=0`（dB，`10*log10(P/Pmax)`）、`aspect='equal'`
- 标注：红★=T-P 峰、橙★=T-S 峰（`annotate` 带 BAZ/v 数值）、黑虚线=理论 BAZ、
  绿圆=理论慢度（T-P 18.5 s/deg=6 km/s、T-S 30.9 s/deg=3.6 km/s）、白点线=等视速度圆
- 每面板右下角信息框（白色半透明）：`BAZ_obs/theo`、T-P/T-S 走时窗、v_TP/v_TS
- 右侧竖排共享 colorbar（`fig.add_axes([0.92, 0.12, 0.012, 0.72])`）+ 总标题

## 2. 单张综合地图（fk6_summary_map.py）——用户定稿形态（v4 期刊版，2026-08-27 最终）

一张图包含全部信息。**注意：下面每项都是用户逐条修正后的最终版。
v3（带走时标注、BAZ 数值、黄光晕、黑方块中心）已被否决——用户最终要求
"无需标注到时和反方位角，图要符合期刊标准"。**

- **底图**：白底 + 大陆轮廓线（用户明确要求，不要地形）：
  ```python
  # ⚠️ region 必须包含震中！马达加斯加事件震中 (-34.89S) 在 [43,60,-28,-11] 画布外，
  # 用户指出"原始震源没有啊，注意经纬度画布大小"。最终 [42,61,-38,-10]。
  fig.basemap(region=[42,61,-38,-10], projection='M17c', frame=['xa5f5','ya5f5','WSne'])
  fig.coast(land='white', water='white', shorelines='0.5p,gray30',
            borders=['1/0.4p,gray50'], resolution='h')
  ```
- **-1000m 等深线**（灰色虚线，T 波海陆转换参考）：
  `fig.grdcontour(grid='@earth_relief_15s', region=region, levels=[-1000], pen='0.8p,gray40,--')`
- **震源**：干净红星 `a0.7c red`（**无黄光晕**——期刊版）+ `fig.meca` 机制球 `scale='0.5c'` +
  简短斜体 `Epicenter`（8p,Helvetica-Italic，**不写坐标/深度**）
- **台站**：全部用**三角**（`style='t...'`，不要圆形）；参与三角形的台站红色
  `t0.32c`、其余灰色小三角 `t0.22c`；**不画三角形彩色连线**（用户："无需彩色线链接"）；
  **不标站名**（期刊版）
- **反演 BAZ 虚线**：从三角中心向**震源方向**延伸（**无黑方块中心标记**）
  ```python
  az = np.radians(baz)          # ⚠️ 指向震源！(baz+180) 是传播方位=远离震源，用户指出方向反了
  ```
  T-P 黑色虚线、T-S 灰色虚线（`gray45`），**延伸到 -1000m 等深线交点**，交点画灰圆点
  `c0.16c gray40`；**不标 BAZ 数值**（期刊版；数值在正文表格/汇报里给）
- **不加走时标注**（v3 的每站 `站名 拾取/理论 s` 已删——用户："无需标注到时"）
- 图例（简洁版，jBL）：Epicenter / Used station / Other station / Inverted BAZ (T-P) /
  Inverted BAZ (T-S) / BAZ x -1000m crossing / -1000 m isobath；
  `box='+gwhite+p0.4p,gray50'`
- **dpi=300**（期刊要求；fig.savefig(dpi=300)）
- 运行（2026-08-31 起为**分窗形态**）：`python fk6_summary_map.py --comp Z|R|T --win tp|ts`
  输出 `fk6_{comp}_{win}.png`（Z/R/T × T-P/T-S 窗共 6 图，每图只画对应窗的 BAZ 虚线、
  图例单条 `Inverted BAZ (T-P|T-S)`）。⚠️ 连跑多分量会 OSError 12 内存不足，必须分开
  执行（每分量一次）。旧综合双窗版（`_Z.png/_T.png` 各含 TP+TS 两条虚线）已被分窗版替代。

## 3. 陷阱

- **3 站 F-K 峰对 T-S 窗边界敏感**：fk5 用 `[min-15, max+50]`、fk6 用 `[min-15, max+45]`，
  同一三角 TS_BAZ 可差 10-30°。跨脚本/跨图必须锁定窗口参数（定稿：
  T-P `[min(tt)-15, max(tt)+45]`、T-S `[min(tt)-15, max(tt)+45]`）
- **pygmt legend spec 不要用 `L` 行**：`L 0.3c 7p,Helvetica-Oblique,black "..."` 触发
  `text [ERROR]: Record N had bad justification info`——用多条 `S` 行 + 图上 fig.text 代替
- 北部稀疏三角（BARY-ZAKA-SOLA、SOLA-ANTS-BERG）反演 BAZ 常落假峰（v>10 km/s 网格边缘）：
  如实标注，不强行挑选；南部→中北部 4 个三角偏差 6-8° 是可靠结果
- 白底综合图非白像素仅 ~5%（海岸线+标注），属正常，不要误判为空白图

## 4. 走时标注口径（用户曾要求；v4 期刊版综合图已不用，保留供 fk5 热力图/正文表格）

"走时不采用峰值，而采用之前关于偏振的搜索方案拾取的最早径向到时，结合海陆的计算方案"
= 观测走时用 `run_polarization_guided_picks.py` 偏振引导拾取（T-P 选最径向窗、
T-S 选最横向窗，中心 ±3s，CSV 的 `tp_pick_center`）；理论走时用海陆分段
（`tt_tp/tt_ts`，stations.csv 已有）。拾取 CSV 是绝对时间（相对发震秒），
海陆理论也是相对发震——两者可直接对比。

## 5. Z/T 分量选择（用户问"BAZ 用哪个分量最合适？R？"后的实测结论）

综合地图 `fk6_summary_map.py` 支持 `--comp Z|R|T`（默认 Z，分窗输出 `fk6_{comp}_{win}.png`）。
T/R 分量 = NE→RT 旋转后的切向/径向分量，**每子阵用自己中心的理论 BAZ 做旋转基准**
（不是全局统一 baz_ref——各三角中心理论 BAZ 不同，146-166°）。

马达加斯加 Mw6.2 六个三角的 Z vs T vs R 反演 BAZ（2-8 Hz 直接窗，定稿窗口参数）：

| 三角形 | 理论 | Z(TP) | T(TP) | R(TP) | 优选 |
|---|---|---|---|---|---|
| South FOMA-BKTA-CPSM | 146.4° | 138.3° ✅ | 112.1° | 130.8° | Z |
| S-Central BATG-BAND-BITY | 156.8° | 150.9° ✅ | 126.5° | 194.7° | Z |
| Central BITY-BAND-VOI | 156.0° | 179.0° ❌ | 157.2° ✅(1.2°) | 115.8° | T |
| C-North ZOBE-BATG-ABPO | 159.7° | 151.5° ✅ | 151.5° ✅ | 180.0° | 一致 |
| North BARY-ZAKA-SOLA | 164.2° | 130.6° ❌ | 152.6° ✅(11.6°) | 143.1° | T |
| NE SOLA-ANTS-BERG | 165.6° | 205.4° ❌ | 173.2° ⚠️ | 122.4° | 均不佳 |

结论：
- **单一分量无法全覆盖**——波型随台站-转换点几何而异；Z/T 双版本互证是论文级做法
- **R 分量整体最差**（无一三角偏差<10°）：`snr_r` 普遍最低（VOI 5.2 dB vs Z 39.8）、
  依赖旋转基准（真实 BAZ 偏差大时能量串扰到 T）；用户仍要求出 R 版作完整对比
- **BAZ 反演首选 T-S 窗**（用户确认"应该是S吧"）：T-P 窗混入 T-S 能量（tt_ts-tt_tp 短时
  尤其严重，FOMA 仅 6s），25 站对齐残差 F-K 的 T-S 窗残差 s≈0 最自洽；T-P 窗残差 450-740
  ms/km。单台偏振则相反（P 波偏振轴直接指 BAZ）——区分 F-K 波前方法与偏振方法
- Z 分量：SNR 最高（stations.csv snr 列）、标准 F-K 惯例、无旋转依赖 → 综合图默认
- T 分量：S 波切向能量（用户要求表述说 **S 波/T-S**，不写 "Love 型"）

## 6. BAZ×-1000m 交点计算：性能优化（避免 grdtrack 超时）

`fk6_summary_map.py` 早期版对每条 BAZ 虚线调 `pygmt.grdtrack`（12 次），每次重新加载
`@earth_relief_15s` → 400s 超时。定稿方案：**预加载地形一次 + 内存插值**：

```python
from pygmt.datasets import load_earth_relief   # ⚠️ 在 pygmt.datasets，不是 pygmt.load_earth_relief
elev_da = load_earth_relief(region=region, resolution='15s')

def elev_at(lons, lats):
    import xarray as xr
    return elev_da.interp(lon=xr.DataArray(lons), lat=xr.DataArray(lats),
                          method='linear').values

def baz_1000m_cross(lat0, lon0, baz):
    az = np.radians(baz)                        # BAZ 方向 = 指向震中
    steps = np.arange(0.02, 10.0, 0.02)         # 度, 从三角中心向外
    lats = lat0 + steps * np.cos(az)
    lons = lon0 + steps * np.sin(az) / np.cos(np.radians(lat0))
    elevs = elev_at(lons, lats)
    idx = np.argmax(elevs <= -1000) if np.any(elevs <= -1000) else -1
    return (lons[idx], lats[idx]) if idx >= 0 else None
```
