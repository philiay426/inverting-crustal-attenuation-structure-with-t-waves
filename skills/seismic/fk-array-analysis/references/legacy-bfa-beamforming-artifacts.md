# 旧产物 / 第三方 F-K 产物（bfa + Beamforming.pl）的解析与判读

用在这几件事上：拿到别人（旧工作站、合作者）留下的 F-K 产物（`fk_signal.out` / `fk_signal.dat` /
`beamforming.ps` / `vespagram.xyz`）；用户问"这份旧产物是不是论文里的那个结果"；要把旧读数与
自己的 fk4d 结果并列。

## 1. 判"是不是论文的结果"——四项核对，任一不同就不能当论文结果

| 核对面 | 怎么查 | 不同就意味着 |
|---|---|---|
| 时间戳 | `ls -la --time-style=full-iso` 产物与引擎目录 | 早于/晚于论文期，只是某次测试 |
| 引擎 | 产物命名 + 表头字段 + 生成脚本（`bfa`/`Beamforming.pl` vs SEIZMO `fk4d`） | 换引擎 ≠ 论文口径 |
| 输入口径 | 表头写的频带（`corrected frq_lo&hi`）+ 是否包络（论文是 **0.25 Hz 低通包络** 不是原始带通） | 原始 2–8 Hz 直接 F-K ≠ 包络 F-K |
| 台站子集 + 单位 | 解析 `fk_signal.dat` 的 `filename` 行取台站；单位见 §3 | 全阵 ≠ 论文的南/北子阵 |
| **是否同一批次** | 表头 `rflo and rfla`（= 台站经纬度均值）+ `grid_pts` + `corrected frq_lo&hi` | 同目录的 dat 与 out 未必同一次运行 |

⚠️ **`fk_signal.dat` 与 `fk_signal.out` 放在同一目录 ≠ 同一次运行**：本项目顶层那份 dat 只有 3 台
（`s16/17/18.sac`，一次 3 台小测试），而同目录的 out 是 **39 台**运行（输入在 `sac_file/`）。
判"是不是同一批次"只看表头 `rflo and rfla`（台站经纬度均值，本项目 39 台 = 26.6192/−24.6287、
3 台 = 23.2422/−28.8384）与 `grid_pts`，**不要看文件名相邻**。

结论句式固定成"这是 <日期> 用 <引擎> 在 <口径> 上的一次测试，与论文 Fig.<n> 口径不同"，
**不要把旧测试写成"论文期参数"**（用户会立刻纠正）。论文自己发布的数据包也可能**只有波形、没有
F-K 产物/脚本** ⇒ 论文的 F-K 结果只能从 PDF 图里读或按口径重算，别指望在旧盘上找到。

## 2. 格式

**`fk_signal.out`**（bfa 写出）：

- 14 行 `!*` 表头：`unit_in_deg` / `apparent veloc` / `slowness` / `back azimuth` / `slowness range` /
  `rel max power` / `abs. power` / `quality` / `rflo and rfla` / `Sx & Sy of max power` /
  `corrected frq_lo&hi` / `grid_pts`
- 一行 `====Sx========Sy========Power====`，随后 `grid_pts²` 行三列 `(Sx, Sy, Power)`
- 行序 = **Sy 外层、Sx 内层**；坐标 = `SMAX − (grid_pts−k+1)·step`，`step = 2·SMAX/grid_pts`
  ⇒ 从 `−SMAX` 到 `+SMAX−step`（不是到 `+SMAX`）
- 写入值 = 源码里的 `work_spc/power` ⇒ **按其自身峰值归一，max 恒 = 1.000000**
- 表头 `grid_pts` 常是 201（脚本默认 `$grid=400`，**按产物表头为准**，别照默认值断言）

**`fk_signal.dat`**（bfa 输入）：每台一条记录 = `filename :` + `stlo:` / `stla:` / `stdp:` / `delta:` /
`length:` / `start:` + ASCII 波形采样值。台站数 = `filename` 行数；台站码用 regex 从
`*.<NET>.<STA>..<CH>.M.SAC` 取（**按 `.` 固定下标取会拿到小时字段**）。

⚠️ 这类文件常属旧机器上的 uid（如 1001），普通用户读不了 → 用 root 会话读/拷（WSL 下 `wsl.exe -u root`）。

解析脚本：本类任务的产物 `parse_fkout3.py`（网格尺寸/步长/归一化、网格 argmax 与表头峰的一致性、
输入台站清单、单位核对一次打印完），落盘在导出目录里，可直接改路径复用。

## 3. 两个必须报出的坑

- ⚠️ **表头报的峰 ≠ 网格 argmax**：表头 `Sx & Sy of max power` 来自引擎内部（`slbbfk`）的峰值搜索，
  与写出的网格 argmax 可以差 0.6 个网格步长、7°；表头峰位置的网格值不是 1.0（实例 0.776）。
  ⇒ 报"当时的峰"用**表头值**；自己从网格取峰时必须注明"这是网格 argmax，与表头报的峰不同"。
- ⚠️ **单位口径**：表头 `apparent veloc = 1/|S|` 严格成立，但 `unit_in_deg=1` 声称 s/deg 时
  |S| = 0.867 对应 128 km/s（非物理），按 **s/km** 读才是 1.15 km/s（T 波量级）。`bfa.f` 里
  `deg_to_km = 111.19`、`app_veloc/deg_to_km`、`sx_max*deg_to_km` 的换算方向与论文的
  `v = 111.195/u`（u 记 s/°）**差一个 111.19 因子**。引用旧读数前先选定口径并写明用的是哪个。

## 4. 引擎与参数（源在 `opt/Beamforming/`，来自 `Codes_ZhangHan/beamforming_v1.0`）

- `Beamforming.pl SACfiles -F frq_lo/frq_hi -T begin/duration [-A incidence] [-K 单位换 s/km]
  [-S slowness[/step]] [-O outfile] [-P psname] [-C] [-V]`
  默认：`grid=400`、`slowness_max=10.0 s/deg`、`step=0.05 s/deg`、`theta=90°`，输出
  `fk_signal.dat` / `fk_signal.out` / `beamforming.ps`（GMT 图）；要求各 SAC 的 `o` 一致
- `VespagramR.pl`（vespagram 聚束）：`N=1.0`（N 次根叠加）、`incidence=100°`、Hilbert 变换开
- 代码分工：`bfa.f` 读 stdin 参数、写表头与网格；`slbbfk` / `bbfk.f` 做白化 + 慢度扫描 + 峰值搜索；
  `beam_del.f` / `slfkms.f` / `rec2pol.f` 为同族函数
- ⚠️ **脚本里常有写死的绝对路径**（`$VESHOME=/run/media/<user>/<盘标>/…`）：既能反推"这份代码抄自谁的包"，
  也是"旧盘上资料不全时该去找哪块备份盘"的线索——见到就报给用户。

## 5. 与 SEIZMO fk4d 的可比性

| 项 | bfa / Beamforming.pl | SEIZMO fk4d |
|---|---|---|
| 慢度单位 | 文件自称 s/deg，数值上更像 s/km | s/°（= s/deg），`v = 111.195/u` |
| 网格 | 方阵 ±SMAX，`grid_pts` 写在表头 | 方阵 ±SMAX，`SPTS×SPTS` |
| 时间维 | 单窗一帧（`-T begin/duration`） | `width/overlap` 滑窗（`width<100%` 受"帧长/2"上限） |
| 输出 | 三列 `(Sx,Sy,Power)` + 表头读数 | `s4d.beam`（**已是 dB**，`[SPTS×SPTS×nfreq]`） |
| 频率 | 表头 `corrected frq_lo&hi` | 带内频点由 `Δf = 1/(nsp·dt)` 阶梯定 |

⇒ 只有"同事件 + 同窗 + 同输入口径"的旧产物才可与自己的 fk4d 结果并列；否则只能作"历史上另一次尝试"的
背景信息。本项目实例：旧盘上那次 bfa 运行 = 39 台 BHZ / 2–8 Hz **原始波形** / 全阵 / `grid_pts 201`
（表头 313.42°、0.867），论文 Fig.4 = `fk4d` + **0.25 Hz 低通包络** + 南/北子阵（南子阵 10.0 / 24.7 s/°）——
两者不可互证，别把前者写进复现对照表。

⚠️ **表头那个方位（313.42°）不是"这次运行的方位"**：它的峰 `|u| = 0.867 s/°` 是近原点针尖，
方位本身无意义（见 SKILL.md"峰落在原点时 BAZ 无意义"）；而该台阵对这次事件的真实 BAZ 是 ~137°，
旧产物里的 `az` 头段（315.43）与 313.42 是**另一端口径**——拿它当基准会把指向震源的峰误判成反向。

## 6. 把这两套老引擎拉到本机跑（沙盒 + stdin 接口 + 成本）

沙盒：`bin/` 放 `bfa` / `readsac`（另拷 `FK_analysis/*.f`、`plotbfa.sh` 备查）；`lib/` 放老运行时库
`libgfortran.so.3`（两条来源都验证可用：① 旧盘 root 卷 `<旧盘 root>/usr/lib64/libgfortran.so.3.0.0`；
② 无旧盘时从 Ubuntu 老发行版包取，**不需要 root**：
`curl -sO http://archive.ubuntu.com/ubuntu/pool/main/g/gcc-5/libgfortran3_5.4.0-6ubuntu1~16.04.12_amd64.deb &&
dpkg-deb -x *.deb ./x` → `./x/usr/lib/x86_64-linux-gnu/libgfortran.so.3`）；运行一律
`env LD_LIBRARY_PATH=<sandbox>/lib <bin>/bfa`（先 `ldd` 确认解析到旧库，不要用 libgfortran.so.5 硬顶）。

**`bin/` 必须同时在 PATH 上**：`Beamforming.pl` / `VespagramR.pl` 先查 `-x readsac && -x bfa`，不满足时
`$BFAHOME=''` 退化为按 PATH 找；PATH 里也没有就报 `Can't exec "readsac"`（产物只剩一个空 dat/表头）。

⚠️ **`-O` 选项必须贴紧写**（`-Ofk_signal.out`，不是 `-O fk_signal.out`）：带空格时那个文件名会被当成
SAC 输入 ⇒ readsac 报 `Unable to open <out>`，接着 bfa 拿到空的 outfile 名 ⇒
`Fortran runtime error: Cannot open file ''`、产物为空或只有表头。见到这两个错先查命令行选项写法，
不要怀疑引擎或输入数据。

**bfa 的 stdin 是固定字段格式读入，宽度必须照 `Beamforming.pl` 的 printf 来**（`bfa.f`：
`read(*,100) frq_lo,frq_hi` / `100 format(f10.6,f10.6)`、`read(*,'(f10.6)')` 单值）：

1. `out_file` 2. `chan_num` 3. `"%10.6f %10.6f"` → `frq_lo frq_hi` 4. `"%f"` → `slowness_max`
5. `grid_pts`（`-S slowness/step` ⇒ `grid = 1 + int(2·slowness/step)`；±10/0.1 → 201）
6. `verbosity` 7. `unit_in_deg` 8. `"%f"` → `theta`。数据文件名**硬编码** `fk_signal.dat`。

⚠️ **随手写 `2 8` / `10` 会被静默读成 `0.000028` / `0.000010`**（字段宽度不足按定点数解读）：产物照样生成、
表头照样打印，只是整张网格全 0。见到 `corrected frq_lo&hi = 0.000028 0.000000`、`slowness range 0.000010`、
网格 Power 全 0，**先查 stdin 格式**，不要怀疑数据或引擎。

**readsac 的 stdin**：第 1 行台站数，之后每台两行 `SAC 路径` + `"begin duration 0.0"`（`begin` 相对该 SAC 的 `b`）；
输出 `fk_signal.dat`（每台 7 行头 + `length` 个 ASCII 采样）。**每条记录按自身最大值归一**（`max|d| = 1`）
——旧产物里的 1.000000 不是单位问题。台站坐标取自 SAC 头段 `stlo/stla`，`stdp` 缺失写 `-12345` 即可（引擎不报错）。

**fkone_circle（刘巧霞）的 stdin**：`t1` / `nsamp`(2 的幂) / `nwin` / `"中心频点序号 ±邻点数"` /
`"smin smax sinc"` / `"参考 lat lon"` / 之后逐行 SAC 直到 EOF。`t1` 在 **SAC 自身时间轴**上（首样点 = `b`）；
`t1 + nsamp·dt` 不能越过文件数据末端，否则报 `Endtime oustide of time range` 并整套终止
（200 s 窗用 nsamp 2048/4096、500 s 窗用 8192）。

**成本**：bfa 单次 ≈ `npts → 2^⌈log₂⌉` × 带内频点数（后者 ∝ 窗长）；7-8 台 200 s 窗约 30 s，
500 s 窗要 20–30 min ⇒ **多组合写成独立进程并行**（各自工作目录 + `start_new_session=True`），别串行等。
⚠️ 删掉**正在运行**的 bfa 工作目录会让它的输出 fd 指向已删 inode、结果写丢 ⇒ 只清理未启动的组合。

⚠️ **obspy 写 SAC 时 `b` 会被覆写**：`trim()` 之后再赋 `tr.stats.sac.b = 0` 不生效（写出 `starttime − 参考时刻`）；
要么**不切文件**、把 `begin/duration` 交给 readsac（更贴近原脚本），要么直接改 SAC 头字节
（`b` 偏移 20、`e` 偏移 24，float32，按文件的字节序 `pack`）。给 fkone 这种"`t1` 全局唯一"的引擎做逐台窗时，
只能靠头字节归零后再喂。

## 7. 反解旧产物的真实窗口（零样本法 + 文件名扫描表）

旧产物的取窗口径常与脚本/注释不符，用两条独立证据定死：

1. **逐条数开头零样本数**：`nz_i = argmax(|d|>0)`、`窗起点_i = 数据起点 − nz_i·dt`。各台 `nz` 不同
   ⇒ 窗是**逐台**给的（本项目 1499 递减到 500）。
2. **读切片文件头段**：`b`（数据起点相对参考时刻）、`o`（发震时刻相对参考时刻）。本项目 `o = 324.965`
   ⇒ 发震时刻 = 文件起点 + 324.965 s，与论文正文给出的 UTC 互证，确认时间基准；
   ⚠️ 这个 `o` 也是"`-T` 被当成相对 `o`、实际按 `b` 执行"的那个差值来源。

两证据联立可判出：旧运行实际窗 = **每台切片起点 + 250 s**（= `sac_file/250_200_2_8.pdf` 那次），
即 `origin − 75 s … origin + 175 s`、逐台按 `tt_p` 平移 ≈50 s —— 逐台平移本身即"重定义阵列几何"，
再叠加 2–8 Hz 在 ~100 km 间距上完全混叠 ⇒ **那次运行的弱针尖（`quality 4`、`rel max power 0.038`）
是噪声，不是可判读的 T 相结果**。

⚠️ 窗口扫描产物名 = `{begin}_{duration}_{band}.pdf`（`0_100_2_8.pdf`、`250_200_2_8.pdf` …）：
**这批文件名就是当年扫过的窗参数表**，比任何脚本注释可靠；用它给"导出产物对应哪一次运行"定案。

3. **拿文件名表当候选、用原二进制逐组合跑，与记录表头对表**（最硬的一条，已实测）：
   把上面的 begin/duration 直接当候选（本项目 22 组），跑完比
   `|S|` / `BAZ` / `rel max power` / `abs. power` / `quality` / `corrected frq_lo&hi` / `grid_pts`，
   完全一致即定案。已验证的两套：
   - 顶层 `fk_signal.out` = 39 台 `cut_*.SAC` / 2–8 Hz / **begin 250、duration 200** / SMAX 10 / step 0.1
     → grid 201：复跑**表头逐行相同**、网格逐点一致（Sx/Sy 全同，Power 最大差 1e-6、argmax 同行同值）。
   - `sac_file/fk_signal.out` = 39 台 `cut_*.SAC` / 2–8 Hz / **begin 700、duration 100** / SMAX 10 / step 0.05
     → grid 400：复跑网格逐点一致（差 1e-6），表头只差 BAZ 末位（193.546158 vs 193.546127）。
   ⇒ `corrected frq_lo&hi` 是**窗长的指纹**（200 s → `1.997070/7.998047`；100 s → `1.992188/7.998047`）：
   先用它把候选砍到同一档，再比 BAZ/|S|，比盲扫省一大半组合。
   ⚠️ `-T begin/duration` 的 `begin` 是**相对该 SAC 自己的 `b`**（`readsac` 里 `nstart=(t1−b)/delta`），
   所以同一个 `begin` 在不同切片文件上落到不同绝对时刻——定案后要再用 §7.2 的头段把绝对窗换算出来。

## 7b. vespagram 的复跑（同一套老二进制）

- `VespagramR.pl SACfiles -Bbaz -Htn -Smin/max[/step] -Tbegin/duration [-A inc] [-N Nroot] [-O out] [-P] [-R rflo/rfla/rfdp]`；
  `duration` = 窗长，输出时间轴 = `begin … begin+duration`（本项目 `-T-100/400` → −100…+300 s、8001 点 @0.05 s），
  横轴**相对 `-H` 指定的头段时刻**（`-Ht3` ⇒ 相对 t3）；`-H` 的值逐个从 SAC 头段取（缺值会
  `Error: tN in SACfile … had bad value`）。
- ⚠️ **`vespagram.c` 把 SAC 文件名读进 50 字节缓冲**：绝对路径或长文件名（原始名 52 字符）直接堆溢出
  （`malloc(): corrupted top size`、产物只剩 220 字节表头）。**把输入拷成短本地名**（`SAxx.sac`）即通过；
  这是老代码限制，不是参数或数据错误。
- 成本：`2001 慢度点 × 台数 × 窗长样本` ⇒ 34 台 / 500 s 窗约 **15 min**、输出 467 MB 文本，别在交互里干等。
- 复跑自证：与旧 `vespagram.xyz` 比 **16,002,000 点最大绝对差 1e-6**（仅打印精度）。本机无 gmt4 ⇒ 不加 `-P`
  （只出数据文件，出图另做）。本项目峰 1.000 @ **t3+81.25 s**、`u = 1.47 s/°`（视速度 ≈76 km/s）
  ⇒ 与前两套 bfa 同族：近垂直/共模，不是 T 相。

## 8. 上新数据时的判读纪律

- 先跑完 §6 的**逐位复现自证**再上新数据；新数据用同一套参数，只换输入与窗。
- 每套结果都要给：带内**实际**频点（`corrected frq_lo&hi`）、`slowness range`、`grid_pts`、
  `rel max power` / `quality`、峰 `|u|`/`BAZ`/峰−中位 dB，并注明**峰是否压在网格边界**
  （`|u| ≥ 0.95·SMAX` 只能读作"≥ SMAX"，直射 P 的理论慢度（~13.4 s/°）在 ±10 网格上就是这种）。
- 子阵间距 ~100 km 时 2–8 Hz 全混叠（无混叠上限 `v/(2·d_min)` ≈ 0.022–0.03 Hz）：
  **单频引擎（fkone_circle 取带中心 5 Hz）在这类子阵上方位乱跳属预期**、不可判读；
  bfa 是整带平均，可以给方位，但必须用**平面波注入检验 + 多分量/双子阵一致性**背书。
- ±10 s/° 的网格**装不下 24.7 s/°**（T-S_local）⇒ 该口径不可能复现论文 Fig.4 的 T-S 能团；
  要复刻必须换口径（子阵 + `f < 0.25 Hz` + 更大网格），交付时把两套口径**分开报**。
- **四种输入口径的对照读数**（A/B 是"参数组合不对"的失败跑；正确口径见 §9）：
  A = 她的 `fk_map.m` 原样（2–8 Hz **原始波形**、9 台、`frng=[1/1000 1/25]`）、
  B = `fk_map_bak.m` 原样（+`envelope`/`getspline`、8 台、`[1/1000 1/33.3]`）、
  C/D = **只把输入换成论文的 0.25 Hz 低通包络**（2 sps，其余参数照 A/B）。全图峰分别在
  217.9° / 18.9° / 47.6° / 31.9°（理论 BAZ 137°）。两条教训：
  ① **差一个参数就整张图变样**——A 缺 `envelope`、B 的 `frng` 上界是 1/33.3 而非 1/25，任一处错位峰就离开源向；
     别把一次失败跑读成"这套引擎/这批数据出不了峰"（§9 的试错矩阵把每项的贡献单列了）。
  ② **论文文字写的口径 ≠ 她图里用的口径**：换成论文的 0.25 Hz 低通包络后，源向扇区最强只比中位高
     0.3–0.4 dB、全图峰在 47.6°/31.9°，而按 §9 口径能精确复现她的图 ⇒ 报读数必须写明用的是哪一套。
  读数口径见 `peak-reading-and-method-comparison.md` §5（扇区局部极大 + 中位）。
  产物：`legacy_disk_params_rerun/06_her_engine_fkmap_4ways/`（四份 beam txt + 脚本 + `源向扇区核对.txt`，
  图 `figures/fig_fkmap_her_engine_4ways.png`）。
- **bfa/vespagram 这一批只给近零慢度共模，没有源向峰**——bfa A `|S|=0.867 s/°`（BAZ 313.42°，指向理论反方位）、
  B `|S|=0.0027 s/°`（网格 argmax 恰在原点）、vespagram `u=1.47 s/°`。⇒ 报告里只能写成
  "历史上另一次口径不同的尝试"。⚠️ **别把这条推广到 SEIZMO 侧**：换成 §9 的正确口径（`envelope` + `frng=[1/1000 1/25]`）
  能精确复现她的 array1/array2 图（25.13 / 25.90 s/°）。
- ⚠️ **多组合复跑要先报价再开跑**：组合数 × 单次耗时（bfa 39 台 200 s 窗 ≈1.3 min/组、vespagram ≈15 min/次）
  先讲清总耗时与结束方式（后台跑 + 完成通知），并在用户开口前给进度——本项目 22 组 bfa 扫描 + 15 min
  vespagram 让用户直接问"什么情况，怎么这么久"。

## 9. 已实测：精确复现她的 `fkmap` 图（2019-12-26 array1 / array2）

这是**复现旧盘 f-k 图的唯一已验证口径**，照抄即复现；缺下面任一项，峰就离开源向。

⚠️ **引擎身份先说清（用户反复问过，别每次都从头查）**：她用的是 `fk/fkmap.m` —— 频域 **co-array 波束成形**
（白化 → 互谱 → 相位斜坡 → 求和 → `abs(real)` → dB），**峰位即估计值、不给残差**。既**不是** `fk4d`
（那是壳，逐帧调 `fkvolume`，而 `fkvolume` 的波束段与 `fkmap` 逐字相同），也**不是平面波最小二乘拟合**
（`fk/` 目录里没有 fit/lsq/inverse/estimate 类函数）⇒ 写方法时不能说成"平面波拟合"，也不必去找拟合函数。

| 项 | 值 |
|---|---|
| 引擎 | SEIZMO `fk/fkmap.m`（coarray 默认） |
| 输入 | 2–8 Hz 带通 → **`abs(hilbert)` 包络 → `getspline(0.1)` 平滑**（**不是**原始波形） |
| 采样 | 20 sps；`synchronize` → `interpolate(20,'spline',1200,2200)` |
| 预处理 | `normalize` → `removemean` → `taper`（顺序照她脚本） |
| 参数 | `smax=30 s/°`、`spts=300`、`frng=[1/1000 1/25]` Hz（上界 0.04 Hz） |
| 子阵 | array1 = 台名排序 `sf([2:9,11])`；array2 = `sf([9:16,24])` |
| 时间基 | 公共起点 = **最晚那台**的 starttime（`synchronize` 语义），窗再 +1200…+2200 s |

| 目标（从她 `.fig` 解出的网格） | 她 | 复跑 | 逐点对比 |
|---|---|---|---|
| array1（9 台，300×300） | 25.13 s/°@133.1° | **25.13@133.1°** | rms 0.008 dB、r=1.0000 |
| array2（9 台，300×300） | 25.90 s/°@163.1° | **25.90@163.1°** | rms 0.005 dB、r=1.0000 |
| 8 台版（`fkmap_0/1000.png` 那一族） | — | 25.38@133.1° | — |

**试错矩阵（哪一项变了就掉）**：原始波形（无包络）→ 16.15@202.5°；她自己的 `transed_2_8`（无包络）→ 13.88@318.5°；
把 `frng` 换成 `[1/1000 1/33.3]` → 11.63@18.9°；换成论文的 0.25 Hz 低通包络（2 sps）→ 10.24@136.2°（方向对、慢度减半）。
⇒ **决定性的是「包络」与「`frng` 上界 1/25」两项**；`spts` 200/300 不影响峰位，"用 obspy 等价带通还是她 SAC `bp n 4 p 2 c 2 8` 的产物"也不影响
（两套滤波的包络结果逐位相同 ⇒ 用 `butter(4,[2 8])+filtfilt` 外置替代 Octave 缺的 `iirfilter/dfilt` 是安全的）。

⚠️ 复现自证要**逐点**（rms dB + 相关系数 + 峰位），只比峰位会漏掉"峰对上但整图不同口径"的情况。
产物：`legacy_disk_params_rerun/08_reproduce_her_fkmap/`（V1–V8 beam txt + `repro_*.m` + 对照图 + `复现说明.md` + `逐点对比.log`）。

**出图对齐（用户会要求"色度尺对应上"）**——对照图的色标逐项照她的 `.fig`，不要自己挑：

| 项 | 值 |
|---|---|
| `CLim` | array1/array2 = **`[-0.7 0]` dB**；8 台族（`fk_8_*`） = `[-0.15 0]`（`fk_map_bak.m` 那句 `caxis([-0.15 0])` 只对应这一族） |
| colormap | 从 `.fig` 导出的 64×3（= SEIZMO `plotfkmap` 的 `flipud(fire)`）→ 存 `her_cmap.txt` 直接复用；**白 = 低、黑 = 峰** |
| 渲染 | `pcolormesh(vmin=-0.7, vmax=0, shading='auto')`（她用的是 `pcolor`）+ `set_aspect('equal')` + 轴限 `±30.1`，色标刻度**每 0.2**（`0 / −0.2 / −0.4 / −0.6`）|
| 脚本 | 定稿 `compare_repro5.py`（参数与交付约定见本节末「叠加元素 / 交付目录」） |

⇒ 换成她的 colormap / `CLim` 不改数值（逐点 rms、相关系数、峰位不变），但**"数值对上"不等于"图上对上"**：
她 `.fig` 里的显示属性同样是要交付的口径，报完数值要把图也对齐。
轴属性怎么批量读（`CLim`/`Colormap`/`YDir`/`DataAspectRatio`）见 `seismic-paper-reproduction/references/legacy-artifact-excavation.md` §1。

**极坐标网格与叠加元素（用户按这套认「对」）**：她的 f-k 图网格 = 等慢度虚线圆 10/20/30 s/° **+ 360° 每 30° 一条虚线辐条**
（`.fig` 里就是 12 条 2 点线）。用户说「**BAZ 的十字虚线**」指的就是这套辐条网格，**不要自行加垂线/十字准星**
（加了他会回「不是这个十字」）。要“对应上”时先看他指的是哪张图的哪个元素，再动手。

- 理论 BAZ 线（她 `.fig` 里那条 2 点线，长 30 s/°，方向不是 30° 整数倍）可直接读出：**array1 138.16°、array2 144.14°**；
  实测峰相对它 **ΔBAZ = −5.1° / +19.0°**（北子阵偏 19° 是真读数，要报出来）。
  这个方位不是"随手一条参考线"：她 `fk_map.m` 里就是 `baz = azimuth(子阵中心 → 事件位置)`，事件位置取自 `*SA18*BHZ*.SAC`
  的 `evlo/evla` —— **普通 `.SAC` 头段是「真震中」（本项目 42.6807/−42.4134），而 `.transed` 头段被改写成 T 波转换点
  （28.83/−32.75）**。用子阵中心复算：array1 → 震中 = **138.17°**、array2 → 震中 = **144.15°**（与她图上 138.16°/144.14°
  差 0.01°）；"指向转换点"的假设被 array2 直接证伪（中心 → 转换点 = 159.1° ≠ 144.14°）。⇒ 画/报这条线一律用**震中方位**，
  别拿 `.transed` 的头段去算；两子阵差 ~6° 只因中心位置不同。
- ⚠️ **读 `.fig` 网格时 `YData` 常是降序**（第 1 行 = 最大慢度）：先 `np.unique` 再 reshape 会把行序与 y 解耦，峰 BAZ 镜像成 `180−x`
  （133.1° ↔ 46.9°）。按**文件行序**取 y（`ys_raw = d[::nxs,1]`），降序则把矩阵上下翻转。
- ⚠️ **配色别照搬原饱和度**：她的 `flipud(fire)` 原色偏艳 ⇒ 保留她的色彩**顺序**，RGB 向自身亮度混
  （`muted = lum + SAT*(c−lum)`，脚本里一个 `SAT` 变量控制）。实测用户口味：`1.0` 嫌"好鲜艳"、`0.5` 又嫌"太暗"，
  **定稿 `SAT = 0.8`**（0.3 近灰阶）；要调就一次给两三个候选值让他挑，别一个个来回试。
- **叠加元素定稿参数**：等慢度虚线圆 10/20/30 s/° + **360° 每 30° 一条虚线辐条**；黑实线 = 理论 BAZ；
  **绿半透明十字 = 理论慢度位置（不是实测峰！）**：T-S_local **24.7 s/°**、T-P_local **13.7 s/°**
  （这两个是论文**预测**值；`10.0` 是**实测**值，别混），沿各自理论 BAZ 用 `baz_to_xy(baz,slow)` 定位；
  实测峰不画任何标记；`pcolormesh` + 等比例 + 轴限 ±30.1；`dpi 300`；十字样式 `ms 16 / mew 2.4 / alpha 0.50`。
- ⚠️ **图上不加文字标注**（用户定稿）：不画读数框、不写面板/总标题、十字旁不写字 ⇒ 实测峰
  （`|u|` / BAZ / ΔBAZ）只报在 md/日志里，别贴在图上。",
- 与论文图拼成对照图时的版式（左右**等框**、右栏**不带表头**、行名只在最右标一次、字号对齐论文图、
  **源图取本次渲染输出而非交付目录旧副本**）见 `seismic-paper-reproduction/references/reproduction-work-conventions.md` §8。
- **交付目录只留最后一版、代码与图同路径**：定稿图放目录根（不留 `figures/` 子目录）、脚本与图并列；
  中间版本与旧脚本**移出工程**（`/tmp/<case>/archived_figs/`，不删）并在说明 md 里记下归档去向。
- 报告图件更新时**只给路径 + 改动清单，不要往回复里嵌图**（用户：「不用贴图，我自己去看就行了」）。
- ⚠️ 往 `/mnt/e` 覆盖 PNG 报 `invalid argument` 或静默不更新时，先 `mv` 旧文件走再 `cp`，并 `md5sum` 核对源与目标。

## 10. 从旧盘 `.fig` / `.png` 里挖当年的网格与参数

她的读数真值来源（MATLAB `.fig` 里的 CData、PNG 图内文字块）与"子阵名单反查"的方法，
见 `seismic-paper-reproduction/references/legacy-artifact-excavation.md`——**反查参数前先读它**，
比拿脚本注释猜参数可靠（本项目正是靠它把子阵、窗、频带一次定死）。
