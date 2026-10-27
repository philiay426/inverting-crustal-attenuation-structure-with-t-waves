# Python F-K 波束形成实现（Bartlett 频域慢度扫描）

2026-08-26 Madagascar 2013-08-17 Mw6.2 事件实战验证（合成+真实数据）。
脚本在事件目录 `process_flow/4_fk_analysis/`：
- `fk_analysis.py` — 引擎（fk_bartlett / peak_info / make_grid）+ 全阵列对齐残差 F-K + 子阵直接 F-K
- `fk3_3comp.py` — 三角形子阵三分量（Z/R/T）F-K，NE→RT 用统一 baz_ref 旋转
- `fk4_figure.py` — 论文 Figure 4 样式 2×2 图（pygmt）
- `plot_triad_map.py` — 三角形子阵选择图

## 引擎要点

- Bartlett 频域波束：`P(s) = Σ_f |Σ_i X_i(f)·e^{i2πf s·r_i}|²`，对慢度网格扫描，峰值→BAZ/视速度
- ⚠️ **steering 符号陷阱**：台站延迟 τ=s·r 使 FFT 相位为 e^{-i2πfτ}，波束补偿必须用
  **正号 e^{+i2πf s·r}**。⚠️ 但用正号时峰值落在 **s = −s_true**（不是 +s_true）：注入
  τ_i = −(s_true·r_i) 后波束在 s = −s_true 处对齐 ⇒ **本 Python 引擎的 BAZ = `atan2(sx,sy) + 180°`**
  （实测 9/9 例，偏差 ≤5.9°）。**fk4d 相反**：它的 `s4d.x/s4d.y` 字段直接给来波 BAZ、**不加 180**
  —— 两套引擎的峰方位约定差 180°，任何跨引擎核对/画图代码必须带各自偏移（漏掉会得到
  "存盘 peak ≠ 结果表"的假报警、且真值诊断点取到镜像方位）。
- ⚠️ **白化数组必须 dtype=complex**：`np.zeros(shape)` 默认 float，复数频谱赋值丢虚部
  （ComplexWarning），相位信息全丢，结果完全错误
- 相位白化（每站每窗 FFT 后 X/|X|）：稀疏台阵必需，防强台站主导
- ⚠️ **频率轴合成口径要固定并声明**：本引擎用**原始 npts（不填充）+ 线性功率和 `Σ_f P_f`**
  （按频点能量加权，高信噪频点主导）。fk4d 相反：**零填充到 2ⁿ + 逐频点转 dB 后等权平均**（见
  `references/octave-seizmo-fk4d.md`）。两者在"带内大部分是噪声"时峰位可差 **2-3 s/°**，
  而把带收到信号真正相干的段后两者**逐值一致**。⇒ ① 跨引擎对比前先统一频带与合成规则；
  ② 报峰位必须同时报频带；③ 引擎间差异先按这条规则去解释，别当成"哪套代码算错了"。
- ⚠️ **峰位对台站顺序无关**（波束 = `Σ_i`，可交换；打乱道顺序峰位不变，实测三种顺序逐位一致）：
  能把结果搞错的是别的事——① 各道必须**同一 b + 同一 npts**（否则 fk4d 直接报 equal B 错）；
  ② **坐标数组必须与数据行同步重排**（只打乱波形不打乱坐标才真会错）；③ 画图时按 τ/距离排序
  只是**显示**顺序，与送进 F-K 的读入顺序无关（用户会问这个）。
- ⚠️ **引擎向量化后必须先做回归校验再拿去下结论**：把逐格双循环改写成"逐 sx 行 + einsum"后，
  与原实现逐点比 `max|Δbeam|/max`（实测 ~1e-15）**且峰位索引相同**，才算同一引擎；
  差在 1e-13 以上就要查 exp/einsum 的维度顺序（顺序错会让整次标定白做）。
- 网格：SMAX=0.7 s/km（覆盖 v≥1.43 km/s），SPTS=101~121；s/deg = s/km × 111.19
  （v=3.6 km/s ↔ 30.9 s/deg；v=6 ↔ 18.5 s/deg；v=1.48 ↔ 75.1 s/deg）
- ⚠️ **扫描范围必须用 SMAX≈0.7 s/km（≈77.8 s/deg）量级**：fk6/fk7 批量脚本曾误用
  SMAX_DEG=35（只覆盖 v≥3.18 km/s），而真实 T 波视速度 1.5–1.7 km/s（慢度 67–73 s/deg）
  在网格外 → 慢度图整片弥散无峰（max/mean≈1.1，真峰应>3）。改 77.8 s/deg 后视速度进入
  合理量级。慢度图弥散先查扫描范围是否覆盖 SOFAR 速度（1.48 km/s ↔ 75.1 s/deg）
- BAZ 换算（**2026-09-11 合成定标修正，旧文写 BAZ=az+180 是错的**）：峰 (sx,sy) 的方位角
  az = atan2(sx, sy)（东/北，从北顺时针）= **来波 BAZ（指向震源），不要 +180**。定标方法与实测：
  注入来波 BAZ=127° 的平面波（物理延迟 τ_i = −(r_i·û_source)/v，û_source=(sin BAZ, cos BAZ)），
  用同一引擎读回 az=143°、慢度 24.5 s/deg（注入 24.7）→ az≈BAZ ✓；若按 baZ=az+180 画理论线，
  理论点会落到实测峰的反侧、差 180°。传播方位 = az + 180（只在画"远离震源"的传播箭头时用）。
  扇区掩膜必须用 atan2(东,北) 顺序。

## 合成验证（改引擎后必做！）

- 用**真实台站几何** + 宽带带通（2-8Hz）白噪声 + 已知 BAZ/v 平面波延迟（FFT 相移 e^{-i2πf·delay}）合成
- ❌ 单频正弦不行：白化后其余频点全是随机相位，噪声主导波束图
- 几何必须满足最大延迟 < 窗长（否则循环移位时间混叠）
- 马达加斯加实测：6 站+120s 窗恢复 (153.0°, 3.6 km/s)→(153.4°, 3.99 km/s) ✓；
  4 站栅瓣淹没真峰；8 站孔径过大假峰反超
- 校验注入延迟用**两两互相关 + 抛物线内插**，别用采样点 argmax：包络顶部平坦时噪声会让 argmax
  跳 1-2 个采样（2 Hz 采样 → 1 s 误差），而互相关内插可到 0.07 s 量级（判据取 ±0.15 s）。
- 9 组 case 矩阵（慢度扫描 × 方位扫描）、双引擎共用同一 `cases.csv`、逐步 PASS/FAIL 核对报告的
  完整配方与实测回收基线见 `references/forward-modelling-calibration.md`。

## 大孔径稀疏台阵（关键物理限制）

- 波前穿过阵列时间 = 孔径/v。25 站孔径 1341 km，3.6 km/s 波需 ~380 s
- 40-80 s 窗直接 F-K 全阵列 → 台站间延迟远超窗长 → FFT 循环移位时间混叠 → 每窗 BAZ 随机漂移
- 解法：
  1. **对齐残差 F-K**：每站窗对齐到各自理论 T-P/T-S 到时（[tt-10, tt+45]s）。
     模型自洽 → 峰≈原点（s≈0）；残差非零 → 方向/速度系统偏差。
     马达加斯加 T-S 对齐残差≈0 → 转换点模型独立验证 ✓
  2. **子阵直接 F-K**：选相邻台站子阵（孔径 ~300 km），窗长 ≥ 2×最大延迟
- 三角形子阵选择：Delaunay 三角剖分 + 质量筛选（边长 60-400 km、面积>2000 km²、
  最小角>20°、质量=面积/周长²、SNR 均值）。FOMA-BKTA-CPSM（159/159/194 km，近等边，
  SNR 均值 30.6）为马达加斯加南部最优

## 结果解读

- 3 站统计量小：BAZ 波动 ±10-18° 正常；**方向比慢度大小可靠**（稀疏几何下直接 F-K
  视速度偏低 1.5-2.1 vs 理论 3.6-6 km/s）
- 三分量：T 分量（Love 型）BAZ 最准；T-P 窗会混入 T-S 能量（当 tt_ts-tt_tp 小时，
  如 FOMA 仅 6s）→ T-P 残差反映混合波前
- 马达加斯加实测：南子阵 T-P 峰 BAZ 138°/北 151°（理论 146.4/156.8°）

## 论文 Figure 4 样式（用户指定表达方式）

2×2：地图（子阵红三角 + 全部台站黑三角 + 震中红星 + 理论 BAZ 蓝虚线）+
F-K 慢度图（s/deg、dB 色标 10log10(P/Pmax) -1~0）。
- 慢度图标注：T-P_local / T-S_local 峰（× + 文字）、黑实线=理论方位角、
  绿十字=理论慢度（T-P 18.5 s/deg、T-S 30.9 s/deg）
- pygmt subplot 2×2（autolabel 自动 (a)(b)(c)(d)），投影传面板内第一个绘图命令
- 慢度网格用 xarray.DataArray 直接喂 grdimage（xyz2grd 要求网格点数精确匹配，易踩坑）

## matplotlib 慢度面板：论文 Fig4(b)(d) 复刻（2026-09-07 用户逐条修正）

pygmt 2×2 适合多子阵组合图；**单张慢度图要像论文图 (b)/(d)——白底上只有几团"能量团"**：

1. **阈值遮罩（核心差异）**：论文图背景纯白、只有 T-P/T-S 两团 >-0.6 dB 的彩色晕圈。用
   `np.ma.masked_where(beam_db < vmin_db, beam_db)` 再 `pcolormesh`——不是把整个网格铺成连续云图。
   `vmin_db=-0.6`（窄动态范围，比默认 -6 敏感得多，晕圈小而集中）。
   ⚠️ **遮罩只适用于"有清晰峰"的结果**：复现出来是弥散/无峰时**绝不能遮罩**——会留大片空白、
   被用户否（"背景不对，即使是假象也不至于这样"）。此时改**全盘连续填色**，色标范围取
   `[峰, 下限]`（基准 = 背景中值），让能量分布如实显示。完整规则与实测见
   `references/slowness-figure-paper-style.md`。
2. **极坐标风格网格（手动画）**：同心圆 `Circle((0,0), r, fill=False, ls='--', color='0.72')`
   （r=10/20/30）+ 放射线（12 根，`ax.plot([0, rmax*sin(a)], [0, rmax*cos(a)])`）——
   不是方形虚线 `ax.grid()`。
3. **色标窄 + turbo**：`cmap='turbo'`（彩虹晕圈更接近原图），colorbar 0~-0.6。
4. **理论慢度十字位置（⚠️ 曾写死出错）**：必须按 BAZ 算坐标，禁止 `ax.plot(0, -slow, '+')`
   （= 硬编码理论方向为正南）：
   ```python
   def baz_to_xy(baz_deg, slow_sdeg):
       ang = np.radians(baz_deg)
       return slow_sdeg*np.sin(ang), slow_sdeg*np.cos(ang)
   tx, ty = baz_to_xy(baz_ref, slow_theory)   # baz_ref 与黑色 BAZ 方向线同一角度
   ```
   黑色 BAZ 方向线也画同一角度：`[0, smax*sin(baz)], [0, smax*cos(baz)]`，保证绿十字落在线上。
   BAZ = 子阵中心→震中反方位角（南非 G1=126°、G2=134°）。
5. **十字尺寸收敛**：`ms=8, mew=1.5`（不要 ms=15 占图 1/8），文字 `fs=8`。
6. theory 各相位点若共用 BAZ 则同线不同半径；若各自有 BAZ 则分别算（不要假设共用）。
7. **F-K 不要求三角台阵**（用户澄清）：任意已知相对位置的台站群（≥3-4 站）即可波束形成；
   Delaunay 三角是马达加斯加大孔径稀疏台阵的经验做法，不是 F-K 前提。



用户定稿规范（画布必须含震中、BAZ 虚线指向震源、台站三角红色、-1000m 交点灰点、
走时标偏振拾取 vs 海陆理论）完整配方见 `references/multi-triad-summary-figure.md`
（fk5_multi_triads.py 热力图 + fk6_summary_map.py 综合地图，均为用户逐条修正后的最终版）。

## 地形交点求法（-1000m）

- ❌ 对每条 BAZ 线分别 `pygmt.grdtrack` → 每次重新加载 15s 地形网格，12 条线直接超时
- ✅ 预加载一次：`from pygmt.datasets import load_earth_relief`（0.19 中**不是**
  `pygmt.load_earth_relief`，后者 AttributeError）；`elev_da = load_earth_relief(region, resolution='15s')`
  返回 xarray DataArray，之后 `elev_da.interp(lon=xr.DataArray(lons), lat=xr.DataArray(lats))` 内存插值
- 沿 BAZ 方向从中心每 0.02° 步进采样，取第一个 elev ≤ -1000 的点为交点

## pygmt 版本/API（2026-08 实测）

| 环境 | pygmt | 情况 |
|---|---|---|
| `miniforge3/envs/obspy` | 0.19.0 | subplot 正常（绘图首选） |
| `~/codex/venv` | 0.18.0 | subplot 机制 bug（"No subplot information file"） |

- pygmt≥0.18：`fig.subplot(nrows, ncols, figsize, autolabel=True)` 初始化后，
  `with fig.set_panel(panel=(r, c)):` 选面板；`irow/icol` 已废弃（Unrecognized parameter）
- set_panel 不接受 projection——投影传给面板内第一个绘图命令（grdimage/basemap）
- `fig.text` 无 `box` 参数（0.18/0.19 均无）
- 慢度网格喂 grdimage 用 xarray.DataArray（xyz2grd 要求点数与 -I 网格精确匹配，报
  "More than N records, only M was expected"）
