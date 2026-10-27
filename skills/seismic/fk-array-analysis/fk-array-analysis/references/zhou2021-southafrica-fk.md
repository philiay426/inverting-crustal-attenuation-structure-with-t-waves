# Zhou 2021 南非 XA F-K 复现排查（2026-09-07，绝对慢度直接 F-K 无峰）

目标：对 1997-04-28 Mw6.6、XA 台网 southern/northern 子阵，复现论文的绝对慢度 F-K
（论文 T-P ~10 s/deg、T-S ~24.7 s/deg）。

## XA 台网身份卡（写汇报/论文/图注要引「全称」时查这里）

论文正文只给缩写（"seismic array (XA) … Southern Africa Seismic Experiment"），汇报/摘要里用户要全称，
别只抄代码：

- **Southern Africa Seismic Experiment（SASE，南非洲地震实验）**，FDSN 网络代码 **XA**（数据 DOI
  `10.7914/SN/XA_1997`），是 **Kaapvaal Project**（Kaapvaal 克拉通多学科计划）的子项目。
- 运行 1997-04 – 1999-07；**55 台宽频带仪器（REFTEK/STS-2）在 82 个站点轮换**（1998-04/05 约半数迁址）
  ——Zhou2021 按 "82 broadband stations" 计；几何为 NNE–SSW 长条，约 1800 km × 600 km、台间距 ~100 km
  （论文把该间距与 f<0.25 Hz 防空间混叠写在同一段）；覆盖南非、津巴布韦、博茨瓦纳。
- 数据取自 IRIS DMC；引用该实验用 James et al. (2001, GRL) 或 James (2003, *The Leading Edge* 22(3),
  238–249, doi:10.1190/1.1564529)。论文 Data Availability 只给数据 Zenodo DOI（`10.5281/zenodo.4776417`）。
- 查全称/规模的渠道：`fdsn.org/networks/detail/<code>`、项目数据 DOI 页、上述 James 文献；
  先查开放状态可用 Semantic Scholar / Unpaywall API（`openAccessPdf.status`），别硬撞 paywall。

## 事件身份卡（写汇报/摘要引事件名与震级时查这里）

- 规范地名 **Prince Edward Islands region**（论文写单数 "Prince Edward Island region"），中文
  "爱德华王子群岛地区"；震中位于西南印度洋脊（SWIR）爱德华王子断裂带附近 —— 洋中脊地震，
  T 波激发效率高，这正是论文选它研究"T 波→陆地转换"的理由（写"选取事件"一栏可用）。
- 目录口径**不一致，引用前先定调**（1997-04-28 12:07:38 UTC）：

  | 目录 | 编号 | 震级 | 深度 | 位置 |
  |---|---|---|---|---|
  | Zhou et al. (2021) | — | **Mw 6.6** | 10 km | −42.4134°N, 42.6807°E |
  | USGS PDE | `usp00081ec` | **M 6.8** (mwc) | 10 km | 42.504°S, 42.686°E |
  | ISC | `iscgem1024295` | 6.72 | 18.9 km | — |
  | GCMT | `C042897C` | **Mw 6.70**（M₀=1.5e19 N·m） | 15 km | — |

  论文的 6.6 与各目录（6.7-6.8）对不上、坐标差 ~10 km ⇒ 汇报写"Mw 6.6（Zhou et al. 2021）"或区间
  "Mw 6.6-6.8"，被问时要能说清是目录口径差；事件编号可作溯源句。
- 查法（无 GUI，一次请求即可拼出上表）：`https://earthquake.usgs.gov/fdsnws/event/1/query?eventid=<ID>&format=geojson`，
  `features[0].properties.products` 同时含 `origin`（US/ISC/GCMT 各一条）与 `moment-tensor`（标量矩/
  节面/深度）；`properties.place` 就是规范地名。

## 数据与响应来源（可核查 —— 被问"数据从哪下载 / 去响应用什么文件"时按此答）

- **波形**：`<case>/download_xa_events.py`（+ `download_xa.log`）→ **IRIS/EarthScope FDSN**
  （`Client("https://service.iris.edu")`、`get_waveforms_bulk`），网络 **XA**（SASE）× 通道 `BHZ/BHN/BHE`，
  窗 = **发震前 600 s ~ 后 3600 s**；本事件 42 台站 → **117 个 SAC**（失败批 0）。
  文件名 `XA.<STA>..<CH>.Q.<年>.<儒略日>.<HHMMSS>.SAC`（起点时刻），SAC 头另写 `stla/stlo/stel` +
  `o=600`（发震在数组中的位置）+ `evla/evlo/evdp/mag`。
- **仪器响应**：`raw/SACPZ/SACPZ.XA.<STA>..<CH>`（**120 个 = 40 台 × 3 分量**），由 obspy
  `Client.get_stations(network='XA', station=…, channel=…, level='response')` → `inv.write(format='SACPZ')`
  **逐台逐分量**抓取（每个文件头 `CREATED` 时间戳都不同、跨十余分钟，正是"循环里一台一请求"的指纹）。
  **证明来源的做法**：同一条命令重抓一台 → `diff` 忽略 `CREATED` 行 → 逐行一致（实测一致，灵敏度/A0/CONSTANT 全对）。
  ⚠️ **头格式即指纹**：本目录是 obspy 写出格式（含 `CREATED/START/END/DESCRIPTION/DIP (SEED)` + `CONSTANT`）；
  Madagascar 案例的 `raw/SACPZ/` 是 **IRIS SACPZ web 服务**短头（`* XK B01KR BHE SAC PZ` + `* LATITUDE :`）
  ⇒ 写数据说明/方法时这两批数据**不是同一来源**，别合并成一句。
- **去响应链**（`process_flow/process_image.py`）：`parse_sacpz()` 读该文件 →
  `paz={zeros,poles,gain=A0,sensitivity=SENSITIVITY}` → `detrend('demean')`/`('linear')` → `taper(0.05)` →
  `tr.simulate(paz_remove=paz, remove_sensitivity=True, water_level=60.0, pre_filt=(0.001,0.005,45.0,50.0))`
  → `sacprocess/remove_response/*_resp.sac`（速度 m/s、20 Hz、`o=600` 保留 ⇒ `t_off = b − o = −600`）；
  台站名单走 `process_flow/1_station_distribution/stations.csv`，SAC 头缺 `stla/stlo` 时从 SACPZ 兜底。
  ⚠️ **只取第一个 ZEROS/POLES 块**：XA 元数据里有多于一个 PolesZeros 级（obspy 抓取时会警告
  `More than one PolesZerosResponseStage encountered`），全读会得到 21 zeros + 33 poles → `simulate` 振幅爆炸。

## 几何与参数核查（窗长是首要嫌疑，但不是唯一病灶）

- southern 6 站（SA16/17/18/22/23/25）最大站间距 312 km；northern 8 站 386 km
- T-S 站间理论到时差（论文 Equation 1）59/52 s；smax=40 s/deg(0.36 s/km) 波束扫描最大延迟
  56/69 s → 原 ±25 s(50 s) 窗不足（须 ≥ 到时差 且 ≥2×扫描延迟，实测需 ~150 s 窗）
- T-P 与 T-S 仅隔 ~70 s：T-P 不能用大窗，改小 smax=20（扫描延迟 28 s < 50 s 窗）
- 修正为 T-S 150 s 窗 + SMAX=35、T-P ±25 s + SMAX=20 后，Bartlett 显著性仍 1.7–2.6（<3），
  slow 值跳变不合理 → 不是窗长/网格参数问题

## 决定性检验：envelope 台间延迟 vs 单一平面波几何预测

用 F-K 实际输入（2-8Hz→hilbert→0.25Hz 低通→2Hz envelope_lp）逐台与参考台互相关：
- southern 以 SA18 为参考，T-S land 窗 [−40,120]s：各台实测延迟 ≈ −2~+3 s
- 而按"单一平面波 4.5 km/s 沿震源方位穿过台阵"的几何预测应为 +15~+51 s
- → 子阵内 T-S 波包几乎同时到达，**不是同一平面波前**（各台转换/散射点不同，或波包被
  平滑/背景主导）→ 直接 F-K 扫不出相干峰是数据特性，与马达加斯加岛内"直接 F-K 平谱、
  仅对齐残差 F-K 成立"一致

## 最终突破：全体台站对齐残差 F-K 成功（2026-09-07 晚，用户指出"F-K 不需要三角，台站群即可"）

- 用 **全体 SNR>3 的 36 台**（而非 6-8 站小簇）做对齐残差 F-K（每站 [tt-10, tt+45] 对齐各自理论
  T-S/T-P 到时，数据用全体台统一的 envelope_lp_all 缓存）：
  - T-P：峰 slow=0.0 s/deg（正好原点），显著性 **3.63** ✓
  - T-S：峰 slow=0.0 s/deg（原点），显著性 **3.01** ✓
- 而 7-8 站小簇（自选 southern/northern、及用户指定组 SA16/17/18/22/23/25/29 与
  SA31/32/37/38/39/40/45/46）对齐残差峰都不在原点、显著性 1.7-2.7 → **相干性对台站数量/
  组合高度敏感**：样本太少波束统计弱；丢最北端台站反而不相干。
- 教训：① 对齐残差 F-K 验证"到时模型自洽"时优先用**全体可用台站**（孔径大无妨，因为
  逐台对齐后不存在延迟/平面波入射问题），样本多显著度才达标；② 引擎是否有问题的反证：
  "全体 36 站对齐残差能出峰≈0 且显著"说明代码与理论到时都正确——若 F-K 代码有 bug 这一步
  不可能成功。
- 绝对慢度（统一绝对时间窗）36 站仍弥散（sig 1.5-1.7）：各台绝对到时差不满足单一平面波，
  与"按分段路径理论（ocean/1.48 + land/3.6）对齐后完全相干"并不矛盾——分段路径模型精确、
  但绝对到时空间梯度≠单平面波。论文能测 24.7 s/deg 说明其数据/子阵存在可识别的平面波成分。

## 复现性总结（诚实上报口径）

- 已具备：全体 36 站对齐残差 F-K 验证 T-P/T-S 到时模型自洽（峰≈原点、显著性>3）
- **已复现（2026-09 定案，本行旧结论「未复现」已作废）**：用她的引擎 SEIZMO `fk/fkmap.m`（coarray）+ **包络**
  + `frng=[1/1000 1/25]` + 子阵 `sf([2:9,11])` / `sf([9:16,24])`（名单由地图 `.fig` 散点仿射反查定案）
  ⇒ **25.13 s/°@133.1° / 25.90@163.1°**，与她的 `fk_result_array{1,2}-20191226.fig` **逐点一致**
  （rms ≤0.009 dB、r = 1.0000）。此前「绝对慢度扫不出峰」只对**错误口径 / 小台簇**成立
  （缺包络、或 `frng` 上界写成 1/33.3 就掉——差一项整张图变样），**不要据此断言论文结果不可复现**。
  配方与试错矩阵见 `legacy-bfa-beamforming-artifacts.md` §9；产物 `legacy_disk_params_rerun/08_reproduce_her_fkmap/`
- Text S1 物理解释：T-P_local 转换区大（南非以南大范围）、T-S_local 转换区小（东南局部）→
  转换区大则各台波前来自不同转换点，更难构成平面波（T-P 尤甚）

## 结论与诚实上报

> ⚠️ **本节是当时（6–8 站小簇 + 绝对时间窗）的结论，已被前文「复现性总结」的更正取代**：换用她的引擎 + 包络
> + `frng=[1/1000 1/25]` + 定案子阵后，绝对慢度能精确复现她的图（`legacy-bfa-beamforming-artifacts.md` §9）。
> 下面两条仅适用于「小簇直接绝对窗 F-K」这一口径。

- 参数全修正后仍平谱时，用互相关延迟 vs 平面波预测做硬判据（比单纯看 max/mean 更直接）
- 对该类数据：不要声称绝对 BAZ/慢度峰；验证模型用对齐残差 F-K（峰≈原点）；与论文对比需
  先解决"台间波包不同步"或改用更大尺度/不同子阵
- 附：每台 envelope 峰相对理论 T-S 的偏移离散 ±20~50 s（无系统性 moveout）——若理论到时
  本身对每台有随机误差，"按理论到时对齐"后波包仍不共面，也是 F-K 弥散来源之一，排查时
  一并看各台峰相对理论到时的偏移分布

## 预处理链原理（用户疑问："2-8Hz 滤波后再低通怎么会有能量？"）

链路：去响应速度(BHZ) → 2-8Hz bandpass(butter4 zerophase) → |Hilbert| envelope →
butter(4, 0.25Hz) lowpass(filtfilt) → 降采样 2Hz → 送 F-K。
中间隔了 **envelope（幅度解调）**：带通滤出 2-8Hz 载波 → 包络提取幅度慢变（低频、非负）
→ 0.25Hz 低通保留慢包络。类比 AM 收音机（载波 vs 调制信号）。实测喂给 F-K 的文件：
2Hz、能量 100% 在 <0.3Hz、基本非负 = 链路正确。解释不清时画 SA18 三步对照图
（2-8Hz 波形 / 包络 / 0.25Hz 低通包络 + 频谱能量搬移），或对照本句复述。
论文 Fig2d 的 envelope 来自 2-8Hz bandpass（Fig2c→d）；本事件 2-8 与 1-5Hz 两口径
fk4d 结果几乎无差别。

⇒ **这条三步链的原始出处 = Matsumoto et al. (2016)**（IMS 水听器三联体：6 阶高通 1 Hz → 包络 →
6 阶低通 0.375 Hz），连同无混叠判据原文、与 Zhou 的逐项对应、三处口径差异，见
`fk-literature-chain.md`（含三篇原文在本机的位置与"按内容 grep 找 PDF"的做法）。

## 论文 Fig4 慢度图绘制逐点规范（2026-09-07 用户逐点纠正）

- 背景纯白、只留能量团：`beam_db=10log10(P/Pmax)` 后
  `np.ma.masked_where(beam_db < vmin_db, beam_db)`（目标图 vmin≈-0.6dB 窄动态范围）
- colormap **turbo**（暖色晕圈更接近论文）；色标 0~-0.6
- 网格：**极坐标风格同心圆+放射线手动画**（matplotlib Circle(ls='--') + 从原点放射线），
  不是方形虚线网格
- 理论慢度十字**必须按 BAZ 算坐标，勿写死 (0,-slow)**：
  `baz_to_xy(baz,slow)=(slow*sin(rad(baz)), slow*cos(rad(baz)))`；与黑色 BAZ 方向线用同一
  角度（子阵中心→震中反方位角），十字才落在线上
- 十字要精巧小号：ms=8、mew=1.5、label fs=8（勿 ms=14 占图 1/8）
- 单相位图每图一个理论点即可（共用 baz_ref）；双团同框时两点同方向不同半径共用 BAZ
- ⚠️ matplotlib 此环境 `plt.subplots(..., sharex=[False,True,True,False])`（按行 list）
  会渲染成 34 万像素宽的爆图 → 改用 `fig.add_subplot(412, sharex=ax0)` 手动逐行共享

## 作者原始代码（yzhou_disk_fk_export）：哪一步在、口径是什么（2026-09-21 核对）

来源 `/mnt/e/sharedate/yzhou_disk_fk_export/`（作者工作站整盘导出 = `home_yzhou/` + `rootfs/`，6 GB，**只读核对、别改**）。
拿到作者代码后**先 grep 确认"哪一步在、哪一步不在"，再核常数与 md5** —— 不要假定"盘里有 = 论文就是这么做的"：

```bash
grep -ril --include='*.m' --include='*.sh' --include='*.pl' --include='*.py' \
     -e fk4d -e '<事件名或网络码>' <导出根> | grep -vi -e '/raw/' -e '\.SAC'
```

- **引擎一致**：`home_yzhou/Documents/MATLAB/Add-ons/seizmo-master/fk/fk4d.m` 与本机那份 **md5 全等**
  ⇒ "论文引擎版本不同导致复现偏差"这条排除。但**导出里没有任何脚本是在本事件上跑 fk4d 的**
  ⇒ 论文 F-K 那一步的实现仍缺（复现仍是"按引擎 + 自己定参数"）。
- **T 波到时口径（`home_yzhou/bin/Ttime`，作者 2017-06-22 自写）**：`vp=6.0, vs=3.6, vt=1.5`，
  `τ_TP = L_land/6.0 + L_ocean/1.5`、`τ_TS = L_land/3.6 + L_ocean/1.5`（写成 SAC 头的 t 标记）。
  ⚠️ **本项目沿途用的是 1.48**：L_ocean≈1630 km 时 1.48 给 1101.4 s、1.5 给 1086.7 s（差 ~15 s），
  逐台差按 `1.48/1.5` 缩放 ⇒ 陆地计时平移、理论到时线、τ 判据、(τmax−τmin)/T 都会跟着变。
  **要统一口径前先问用户**（现行图/表/文档全是 1.48 那套；改就是全链重算）。
- **L_land 的两套作者口径，别混**：
  ① 振幅/衰减链用 `ext_contdis.sh`：`ext_topo`（GMT4.5.16 `project`+`grdtrack`，地形格网 **ETOPO2**，
     0.1 km 采样）取大圆剖面 → `gawk '$4>0'` 数**高程 > 0 的点数 × 0.1 km** = 陆段长（= 海平面 0 m 岸线切分）；
  ② 交互式 `Ttime` 走 `topo2dis.py <profile> <depth>`：用**等深线**（depth 是批参数，例 −1000 m）定转换点，
     `L_land = dist − 沿线到转换点的距离`。
  ⇒ 复核本地 `stations.csv` 的 `land_km` 前先确认它属于哪一套（差几 km 量级，但直接进拟合横轴）。
- **振幅链（`ext_amp.sh`，同一事件目录）**：P 窗 `cut o 250 400`、T 窗 `cut o 1100 1600`（o 相对秒）
  取绝对峰 `depmax` → `amp_P_*.txt`/`amp_T_*.txt`；`plot_amp.m` 画 P、T、**T/P 比值 vs 距离 / 方位 /
  陆段距离**（列：dist, az, P amp, T amp, T/P, P/T）⇒ 与本项目衰减拟合的横轴组织可逐列对照。
- **独立 F-K 引擎（第三方交叉验证可用）**：`home_yzhou/opt/Beamforming/` = 张涵 2015 的 C+Fortran 套件
  （`bfa` 同时测慢度+反方位角、`vespagram` 给定 BAZ 测慢度；Perl 驱动默认 **SMAX 10 s/°、步长 0.05**，
  单位支持 s/° 与 s/km）；`Documents/zy_fk/`（刘巧霞 `fkone_circle.c`）。
  ⚠️ `rootfs/opt/fk/`（L. Zhu v3.2）是**反射率法合成格林函数**工具，不是台阵 F-K，别拿它当对照。
- **论文期 F-K 的真实参数记录在产物表头里（比找脚本更硬）**：作者那套 F-K 是 `Beamforming.pl + bfa`，
  不是 fk4d；`home_yzhou/Documents/T_waves/obs_data/1997-04-28_mw66_Prince_Edward_Islands_Region/fk_signal.out`
  的头逐行写着当时生效的参数：`unit_in_deg 1`（s/°）、`slowness range 10.0`、
  **`corrected frq_lo&hi 1.997 / 7.998 Hz`（即 2–8 Hz）**、`apparent veloc / slowness 1.153 / 0.867`、
  `back azimuth 313.42°`、`quality 4`、`rflo/rfla 26.62 / −24.63`；配套 `fk_signal.dat` 的头是逐台 SAC 清单
  （`delta 0.05` = **20 Hz 原始去响应波形**、`length 4001` ≈ 200 s 窗，单台 SA18、BHZ）。
  ⇒ 三条对照结论必须写进复现笔记：① 作者是**直接对 20 Hz 的 2–8 Hz 波形做 F-K**，不是本项目的
  「2–8 Hz → 包络 → 0.25 Hz 低通 → 2 Hz」链；② 网格/smax/入射角（10 s/°、grid 400、incidence 90°）
  与 `Beamforming.pl` 默认值一致；③ 该产物只有这一个事件、且是**单台**（SA18）测试，
  **不能当论文 Fig.4 的子阵结果**（缺子阵名单这一步仍无解）。
  这两份产物 + `beamforming.ps` 已拷进
  `yzhou_disk_fk_export/home_yzhou/Documents/T_waves_code/obs_data/1997-04-28_.../`。
  通用做法：问"当时用的参数设置是多少"时**先翻结果/输出文件的表头与注释**（`.out`/`.dat` 头、处理日志），
  比先找脚本再读默认值更快也更准（脚本里的默认值可能已被命令行覆盖）。
- **盘里写死过一条 2021 备份盘路径，值得向用户确认**：`VespagramR.pl` 顶部
  `$VESHOME = /run/media/yzhou/Backup_YZhou/T_waves_Research_backup-yzhou-21090414/Codes_ZhangHan/beamforming_v1.0/Beamforming/`
  ⇒ 作者 2021-09 做过一块名为 **Backup_YZhou** 的备份盘，含 `T_waves_Research` 全量备份与
  `Codes_ZhangHan/beamforming_v1.0`（正是论文收尾期）——比 2019 年的工作站盘更可能藏有论文期脚本/
  子阵名单 ⇒ 复现缺的那一步优先去那里找，先问用户这块盘还在不在。
- **导出里没有本事件波形**（作者按要求删了数据）⇒ 只能拿脚本、拿不到当年中间数。
  `home_yzhou/Documents/T_waves/` 是 SEM2D 平面波源生成器 + XA 事件清单（`eqs`、`earthquakes_XA_*.png`），
  与 T 波分析链无关，别混进来。
