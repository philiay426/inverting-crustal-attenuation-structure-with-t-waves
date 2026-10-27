# 偏振角度 → 观测反方位角 → 海岸线转换点图（2026-08-23 马达加斯加定版）

偏振引导拾取（`polarization-guided-picking.md`）之后，用拾取窗口的偏振角反演观测反方位角，
并延长到海岸线求转换点。用户明确要求："只延长靠近震源那一侧，另一侧不延长"，且"在前十台站都做，
不必在意信噪比"（被拒站也画，但其方向不可靠）。

## 1. 反方位角反演：`obs_baz = baz + tp_rt_angle`

`t_phase_analysis.pca(back_azimuth=baz)` 把 ZNE 旋转到 ZRT，`principal_vector` 列为 `[Z, R, T]`，
即 `v[:,1]=R`（指向 back_azimuth，即台站→震源方向）、`v[:,2]=T`。

```python
a = np.degrees(np.arctan2(v[:, 2], v[:, 1]))   # 偏振轴相对 R 的角（0=径向）
# 对径向 P/SV 波，偏振轴沿传播方向 → 观测反方位角 = 理论 baz + 该 R-T 角
obs_baz = (baz + median_angle) % 360
```

- `median_angle` 就是拾取 CSV 的 `tp_rt_angle`（T-P 用）字段，直接读即可；
- 只取 `obs_baz` 这一侧（震源侧），**不画 `obs_baz+180` 分支**（用户要求）；
- 方向自检：接受站的 `obs_baz` 应集中在理论 baz 附近（本例 135–156°），被拒站（如 SOLA 径向残差 58°）
  的 `obs_baz` 会偏到另一侧（222.9°，落到西海岸），属预期不可靠结果。

（与 `conversion-point-joint-constraint.md` 的 N-E 特征向量反演不同：这里直接用 pca() R-T 角，
更简单；两者都是"理论 baz 仅作旋转参考、不做独立反演"这一原则的实例。）

## 2. 拾取玫瑰图（前十站 R-T 角分布）

对每个拾取窗口（`pick_center ±3 s`）做滑动 PCA（1 s 子窗 0.25 s 步长），收集全部 `a % 180`
（轴向折叠），画半圆直方图（`projection='polar'`，`set_theta_zero_location('N')`、
`set_theta_direction(-1)`，10° 一档）。布局 = 10 站 × (T-P 蓝 / T-S 紫) 两列。
判读：T-P 玫瑰应集中 0/180°（径向），T-S 集中 90°（横向）。

## 3. 射线延长到 -1000 m 等深线求转换点

```python
for d in np.arange(0, MAX_KM + 1e-6, STEP_KM):      # 大圆 Direct 步进
    p = GEOD.Direct(sla, slo, obs_baz, d * 1000); pts.append((p['lon2'], p['lat2']))
tr = pygmt.grdtrack(points=pname, grid=grid)         # 采样地形
elev = tr.iloc[:, 2].to_numpy(float)                 # 第3列=海拔
k = np.where(elev <= -1000)[0][0]                    # 首次进入深海的点 = 转换点
```

### 坑（均 2026-08-23 实测）

- **`MAX_KM` 要给足**：内陆台站到海岸可能 600+ km（VINA 沿 155° 到东南海岸 ~620 km 才穿过 -1000 m），
  默认 `MAX_KM=1000`；太短会误报 "no coastline crossing"。
- **pygmt 0.18 `grdtrack` 返回整数列名** `[0,1,2]`（lon, lat, 海拔），`newcolname=` 参数无效；
  用 `tr.iloc[:, 2]` 取海拔，不要 `tr['elev']`（KeyError）。
- **`@earth_relief_01m` 会尝试联网下载**（甚至请求不相干的 S90W180 瓦片）导致超时；改用本地缓存瓦片拼区域网格：
  ```bash
  gmt grdblend ~/.gmt/server/earth/earth_relief/earth_relief_01m_g/S60E030.earth_relief_01m_g.nc \
               ~/.gmt/server/earth/earth_relief/earth_relief_01m_g/S30E030.earth_relief_01m_g.nc \
               -R42/57/-38/-18 -G/tmp/elev_01m.nc -I1m
  ```
  pygmt 0.18 无 `pygmt.grdblend`（AttributeError），用 subprocess 调 `gmt grdblend`。
- 无 -1000 m 交点时回退到 0 m 海岸线再试一次。

## 4. run_polarization_guided_picks.py 运行前置

- `from t_phase_analysis.pca import pca` 会 `ModuleNotFoundError`：先
  `sys.path.insert(0, '/home/sy/codex/t_phase_analysis')`（包在 codex 下，不在 site-packages）。
- 批量脚本用 `TOP_N = 10` + 按 `snr` 降序取前 N，读取**当前** `stations.csv`（新 SNR 口径），
  覆盖旧 SNR 结果前先 `cp` 备份旧 CSV/summary 到 `backup_YYYYMMDD/`。

输出：`polarization_guided_picks/polarization_roses_top10.png`（玫瑰）、
`top10_conversion_topography.png` + `top10_observed_conversion_points.csv`（转换点）。
改图基于旧脚本（如 `plot_ampy_voi_conversion_topography.py`）时输出新文件名，**不覆盖原图**。
