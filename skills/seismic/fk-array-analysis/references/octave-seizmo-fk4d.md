# Octave + Seizmo fk4d 复现（MATLAB F-K 原始脚本）

## Seizmo 来源

GitHub 官方仓库 `https://github.com/g2e/seizmo`（作者 Garrett Euler，MATLAB/Octave 均支持）；
可用 `git clone https://github.com/g2e/seizmo.git` 或 tags tarball。本地这份是 tag
`0.6.16-Otgontenger` 的 tarball 解压（目录名=tag 名，形如
`https://github.com/g2e/seizmo/archive/refs/tags/0.6.16-Otgontenger.tar.gz`），
无 .git 元数据。项目主页 http://epsc.wustl.edu/~ggeuler/codes/m/seizmo/。

2026-09-07 Zhou 2021 南非 XA 台阵复现验证。用户要求"MATLAB 原始脚本跑一遍 + Python 跑一遍对比"。
MATLAB 原始脚本 = Seizmo（用户 git 下载的脚本包）的 `fk/fk4d.m`。

## 两条路线的边界（勿接成 pipeline）⚠️

- **路线 A（本文件）**：Octave 只做"读 SAC → `synchronize` → `cut` → 调 `fk4d` → 取峰"，慢度扫描与
  波束形成全在 Seizmo `fkvolume`/`fk4d` 内部，Octave 侧没有自己的扫描代码。
- **路线 B**：Python（numpy Bartlett）自建网格、相位白化、扫描、找峰，与 Seizmo 无任何调用关系。
- 两者**各自独立可跑**、对同一批输入 SAC 交叉验证（本项目 G1 T-S 两条路线都给出 ≈29 s/deg），
  不是"Python 算网格 → 喂 Octave 做波束形成"的流水线（用户曾问"你现在是通过 python 扫描、
  f-k4d 给 octave 跑"，答案是否）。

## 环境（无 sudo 用 conda）

- `conda create -n octave -c conda-forge octave`（Octave 10.3.0，无需 root）
- ⚠️ **必须 `conda run -n octave octave-cli --eval "..."`** 激活环境。直接 `export PATH=.../envs/octave/bin` 时 octave 的 path 只有 `.`，连 `num2str` 都 undefined。
- Seizmo 用 `startup_seizmo`（cd 到 seizmo 目录后 `run`），**勿用 `install_seizmo`**（交互式问 TauP Y/N 会挂死 + 会 rename 顶层目录，路径全乱）。
- MATLAB 与 Octave 跑**同一份 `.m` 源码**数值等价（差异只在商业工具箱、GUI/绘图、执行速度、
  个别语法细节如 string/table/arguments 块）；用户问"Octave 和 MATLAB 有区别吗"时答：语言/算法
  层面基本等价、fk4d 这类纯 .m 实现结果一致（浮点末位量级），别暗示 Octave 结果不可靠。

## 在 Linux 终端运行（用户常直接问"哪个是 octave 的、怎么在终端跑"）

- 分工按扩展名答，不要逐个脚本解释：**`.m` = Octave/MATLAB（F-K 计算，全部在 `fk_code_octave/`）**；
  **`.py` = Python（预处理 / 检验 / 绘图）**。第三方 Seizmo 工具箱不在事件目录里，别混入清单。
- 标准三条命令（脚本内部自 cd 到 Seizmo 并 startup，所以任意 cwd 都能跑）：
  ```bash
  conda activate octave                        # 交互式终端：~/.bashrc 已有 conda initialize 块，直接可用
  cd <事件目录>/process_flow/4_fk_analysis
  conda run -n octave octave-cli --eval "run('fk_code_octave/export_Tcomp_TS_beams.m')"
  ```
- ⚠️ **别再说"每个新终端要先 `export PATH=/home/sy/miniforge3/bin:$PATH`"**——用户交互式终端本来就配好了
  （`bash -lc` 实测 `conda activate octave` / `conda activate obspy` 直接可用）。**只有非交互 shell**
  （agent 的工具调用、脚本里的 `subprocess`）不加载 `~/.bashrc`，那里裸 `conda` 可能根本不在 PATH
  （实测报 `conda: 未找到命令`）⇒ **用绝对路径 `/home/sy/miniforge3/bin/conda run -n <env>`**，
  或先 `source /home/sy/miniforge3/etc/profile.d/conda.sh`。答环境问题先分清这两类 shell，
  别把工具的限制说成用户的限制（用户会当场反问"我的依赖库里没法直接切环境吗"）。
- ❌ 直接调 `.../envs/octave/bin/octave-cli`（不经 conda）**连 Octave 自带的 m-file 都找不到**：
  实测报 `'fileparts' undefined`（比 `fk4d` undefined 更早失败——环境变量未设、Octave 库路径不全）；
  只 `export PATH=.../envs/octave/bin` 同样不行。⇒ 一律走 `conda run -n octave` / `conda activate octave`。
- ⚠️ **别为省事把两套环境合并**：用户嫌"出图还要换环境"时提议"把 octave 装进出图环境"——不要照做。
  obspy 环境 187 包/2.0 GB 是预处理+绘图+文档的主力环境；octave 环境 1.2 GB 的依赖树含 qt/hdf5/fltk，
  混装会让 conda 求解时降级/替换现有包。正解 = 留一个**一键入口脚本**，内部自动切环境（用户仍在
  obspy 环境里敲一条命令）：
  ```python
  CONDA = "/home/sy/miniforge3/bin/conda"
  def run_fk(tags):        # Octave 侧：算波束，写 fk_code_octave/output_data/
      expr = f"TAGS={{'{tags}'}}; run('{FK}/fk_code_octave/export_Tcomp_TS_beams.m')"
      subprocess.run([CONDA, "run", "-n", "octave", "octave-cli", "--eval", expr], check=True)
  def run_plot(method, tags):   # Python 侧：当前解释器就是 obspy 环境
      subprocess.run([sys.executable, f"{FK}/fk_code_octave/plots/plot_Tcomp_TS.py", "BOTH", method, tags], check=True)
  ```
  用法 `python run_fk.py --window G`（`--only fk|plot` 只算/只画、`--window all` 全跑）；
  参数仍只在 `.m` 顶部参数区，wrapper 只做调度、不藏参数。脚本放在 `4_fk_analysis/run_fk.py`。
- 交互式：`conda run -n octave octave-cli` 进去后 `run('<绝对路径脚本>')`。
- Python 侧固定 `/home/sy/miniforge3/envs/obspy/bin/python <脚本>.py`。

## 默认 method = coarray（2026-09-13 用户改定，原为 center）

- 用户要求 **F-K 一律先用/默认 `coarray` 跑**（原来定的 center 已作废，别自行改回去）。落实方式（批量、可复核）：
  ```bash
  sed -i "s/methods={'center','coarray','full'}/methods={'coarray','center','full'}/" *.m
  sed -i "s/,false,'center')/,false,'coarray')/" export_*.m      # 导出/出图脚本直接 coarray
  grep -h "methods=" *.m | sort -u                                # 复核：只剩一种写法
  ```
- `center` 仍一并导出作对照（结果 csv/日志里排在 coarray 之后），**但出图/报告口径用 coarray**。
- 同一合成数据上两者都无偏、精度相当（|Δslow| 2.60 vs 2.68 s/°、bias ≤0.52 s/°）——默认选谁不是
  "谁更准" 的问题；但**换 method 会改变峰位**（三者都不是"选对就出峰"的开关 → 见下节）。

## 读 SAC + 预处理

- `readheader(files)` + `readdata(data)` 读 SAC。`d.head` 是**扁平数值数组**（非 struct），字段用 `getheader(d,'stla'/'stlo'/'delta'/'npts'/'o'/'b')` 访问。
- FDSN 下载各站 `starttime` 有 ms 级差异 → 必须 `synchronize(data)` 对齐参考时间，否则 fk4d 报 `No windows contained 2 or more records!`（interpolate 产生 nan 后整组被丢）。
- ⚠️ **写 envelope 输入 SAC 时各站 starttime/npts 必须完全一致**：obspy `trim/decimate` 的亚秒对齐不可靠（starttime 仍差 ms 级、decimate 后起点漂移）——用 numpy 对统一起点做 `[::decim]` 重采样 + 强制所有站相同 starttime/相同长度，并写 `o=发震−starttime`、`b=0`（fk4d 的 cut 依赖）。
- `cut(data, t1, t2)`：t1/t2 是**相对 reftime 的秒**，两个标量参数（`cut(data,[t1 t2])` 会当成绝对年份报 "Pre-1972 date"）。
  - ⚠️ **`synchronize` 会重定基准（改写 `b`/`o`），此后裸偏移 `cut` 的语义跟着变**：实测同一文件
    `b=−100, o=0` → `synchronize` 后 `b=0, o=100`（参考时刻被挪到记录起点）。裸偏移相对重定后的参考时刻算
    ⇒ **凡输入 `b ≠ 0` 的脚本，`cut(data, win)` 实际切到的事件段比预期早 `|b|`**（陆地计时输入 `b=−100`
    ⇒ 早 100 s；症状 = 峰塌到近零慢度、方位反向，而注入标定照样"通过"）。
    通用写法 = **`cut(data, o+win)`，`o` 在 synchronize 之后用 `getheader(data(1),'o')` 取**（绝对计时族
    一直这么写，所以只有裸偏移脚本受害）。自检：跑完打印 `getheader(d,'b')/('o')`，并用段和/指纹确认
    切到的位置（完整案例见 `references/land-time-fk.md` §10）。
  - ⚠️ **不要用 `changeheader` 单独改参考时刻字段来"换时间基准"**：`nz` 必须与 `kztime/kzdate` 自洽，
    单独改 `nz` 报 `Reference time fields … stored as noninteger values!` 或 `NZ* not in defined range!`
    且 `exit=1`（stdout 上的读数行只打了一半，别只看读数就当跑完）。要换基准就用 obspy 重写 SAC
    （`starttime` 与 `sac.b` 一次写全），或在内存里同时改 `b/o` 并保持窗选择等价。
- ⚠️ **改波形要写 `data(i).dep`，不是 `data(i).data`**（2026-09-17 合成注入踩到）：Seizmo 记录的数据字段叫
  `dep`；写 `.data` **不报错**、只是新建一个没人读的字段 ⇒ fk4d 照旧读原 `.dep`，表现是"注入/平移没生效
  但也不报错"（据此得出错误的注入标定结论）。自检：换一组参数重跑，两份结果逐位相同 = 没生效。
- ⚠️ **Octave 结构体数组要求所有元素字段完全一致**：给配置数组追加一个带新字段的元素
  （`C(3)=struct(...,'abs',true)` 而 `C(1)`、`C(2)` 没有）→ 报 `incompatible fields in struct assignment`，
  **整个脚本 exit=1**；而 stdout 上读数行只打印到一半（错误在 stderr）——**别只看读数行就以为跑完了**，
  要查退出码。规则：新增字段时给**所有元素**都补上（`'abs',false`）。
- ⚠️ **别在 `.m` 脚本末尾定义子函数当工具用**（如 `function s = compsfx(c)`）：Octave 下不可靠
  （本次报同样的 struct 错误、清掉函数才通）。用映射表代替：
  `CSFX = struct('Z','','R','_R','T','_T'); dir = [base CSFX.(cc) '/...']`。
- ⚠️ **同一份 `.m` 连跑两次、读数必须逐位一致**（本次实测两次 `diff` 完全相同 ✓）——有非确定性就说明
  读到了刚被改写的中间产物；要报稳态结论/剔除"两次不一致"的疑虑，先跑两遍再写进文档。

## 跑"老 MATLAB 脚本"（`iirfilter` / `fkmap` 这类）：Octave 侧的替代与包依赖

- **`iirfilter` 在 Octave 下跑不通**：内部用 MATLAB 的 `dfilt` ⇒ 报 `'dfilt' undefined near line 547`
  （补/不补 order 参数都一样）；省略 order 时报的是 `Corner(s) at/above Nyquist!` —— **误导性错误**，
  别当"参数写错"去调（真实原因仍是 dfilt 缺失后的分支）。
  替代：**滤波在 Python/obspy 侧做等价实现**（`butter(4,[f1 f2],fs=1/dt)` + `filtfilt` = 4 极 2 通带零相位，
  与 `'np',4`+`'p',2` 同规格），其余步骤（`synchronize` / `interpolate(rate,method,win)` / `normalize` /
  `envelope` / `getspline` / `removemean` / `taper` / `fkmap`）**全部走 SEIZMO 原函数、参数照抄**；
  在交付说明里把这一处替代写成"偏差"（用户接受等替代，但不接受不声明）。
- 缺包按报错名装进 octave 环境（conda-forge 有）：`conda install -n octave -c conda-forge octave-signal`
  （`butter` 等）、`octave-splines`（`getspline` 依赖 `csaps`，缺它报 `'csaps' undefined`；
  该步**可跳过**——本项目加与不加输出完全一致，但仍要声明）。
- `fkmap` 可用（与 `fk4d` 同工具箱）：`fkmap(data,smax,spts,frng)`；`frng` 单位 **Hz**（惯用周期倒数写法，
  如 `[1/51 1/49]` = 约 50 s 周期）；内部 `nspts=2^nextpow2(npts)`、`f=(0:nspts/2)/(delta*nspts)`、
  `fidx=find(f>=frng(1)&f<=frng(2))`（**带内实际参与频点数先算出来再报**），`frng>=1/(2·delta)` 直接报错。
  要求各道 `b utc` 相同，否则 `error: Records in DATA must have equal B (UTC)!` ⇒ 老脚本里的
  `synchronize` 是必需的（FDSN 各站 starttime 本就差 ms–s 级）。
- **无头环境画不出图**：`plotfkmap` 需要图形后端 ⇒ 出图改在 Python（读同一份导出的 beam 文本）。
- Octave 自身两个常撞的小坑：
  ① **匿名函数体只能是单表达式**：`@(s) (fprintf(...); fflush(stdout))` 报
     `anonymous function bodies must be single expressions` ⇒ 进度输出写成
     `prog = @(s) fprintf(stderr,'%s\n',s);`（stderr 无缓冲，重定向日志也能实时看到）。
  ② 被 `sed -i` 改过的 `.m` 再想用 write_file 整篇重写会被拒（detected external change）⇒ 用 `patch` 逐处改，
     或另存新文件名。

## fk4d 用法

- `s4d = fk4d(data, width, overlap, smax, spts, frng, polar, method)`
- `smax` 单位 s/deg；`frng` 单位 Hz；`method`：`'coarray'`(默认)/`'center'`/`'full'`

### fk4d 滑窗：首帧静默丢失 ⚠️（2026-10-14 实测，6 台已严格对齐的输入）

- 每帧起点 = `b + (i-1)*tstep`、且 `b = min(b)`（fk4d.m:149,170），所以 i=1 的起点与本批记录的
  第一个样点时刻**重合**；而 fk4d 内部取帧用 `interpolate(data, 1/delta, [], start, finish, nan)`
  （fk4d.m:172，nan 是窗外填充值）——实测已对齐的记录会有 5/6 台的首样点被判为窗外 -> 含 NaN ->
  第 175 行整道丢弃 -> 第 178 行"剩不到 2 道"**整帧跳过**。结果是 `numel(s4d) = nframes - 1`，
  返回帧的 butc 从 `b + tstep` 起（**静默、不报错**）。WIDTH=40/OVERLAP=50 与 40/52 都一样丢首帧。
- 规则：**fk4d 滑窗跑完先核对 `numel(s4d) == fix((1-w*o)/(w-w*o))`**；不等就是丢了首帧，
  别把返回数组的第 k 个当"第 k 帧"（窗位置整体差一个步长，读数标签会全错）。
- 绕开办法：**逐窗自己切 + `fkmap`**（= 现役单窗链：`synchronize` ->
  `interpolate(d, SR, 'spline', off1, off2)`（无 nan 填充、不丢道）->
  normalize/envelope/getspline/removemean/taper -> `fkmap(d, smax, spts, [1/T, min(C/smax,1/25)])`），
  窗口精确、不丢帧、与现役单窗读数同引擎。给数据加前导保护边**没用**（帧起点仍是 min(b)）。

### fk4d 三种 method 的语义（读 fkvolume.m 确认，用户会问）

- **相位白化对三种 method 都做，且在 method 分支之前的公共代码里**：`cs=data(fidx,:).'; cs=cs./abs(cs);`
  （源码注释即 "cs is normalized by the auto spectra"、"center/user … do need normalization"）。
  所以 center 测不出峰**不是"没白化"**，别这样解释。
- 白化之后才分叉：`coarray` = 台站两两配对交叉谱 `cs(row).*conj(cs(col))`（N(N-1)/2 对）；
  `full` = 全部配对；`center` = **不配对**，直接用各记录相对阵中心的白化谱。
  center 计算最省，但统计量只有 N 个"记录-中心"相位、没有台站对之间的互相关信息 → 抗噪最弱。
- ⚠️ **method 不是"选对就出峰"的开关**：同一子阵、同一 2-8 Hz 输入，三种 method 给出的峰互不相同
  （XA G1 T-S：coarray 4.1 / center 12.0 / full 12.0 s/deg；G2 T-S：3.5 / 37.9 / 36.2），
  且**没有一个给出与理论慢度(24.7)和理论 BAZ 同时一致的清晰峰**（峰值-中值仅 2.3-3.0 dB）。
  结论应表述为"方法不是决定性因素"，不要拿某一种 method 的单个数字当复现成功。
- ⚠️ **比较两脚本的 fk4d 结果前，先比输入（grep 实际读的目录/参数），再怀疑语言语义**：
  曾出现两脚本"同一 method 差 6 倍慢度"，追查下来是 `sed` 改 `ENVDIR` **没生效**（目标行行尾是
  `'` 而替换模式写的是 `"`）→ 一个脚本读 1-5 Hz 目录、另一个读 2-8 Hz，被误判为"cell vs 字符串
  传参 bug"，白折腾一轮。规则：`sed -i "s/A/B/"` 之后**必须 `grep -n` 复核实际值已变**（sed 静默
  不匹配、exit 0 不代表改到）；两脚本数值矛盾时先 `grep ENVDIR/参数` 对输入，一致后再看代码。
- Octave 里 `for meth={'a','b'}` 的循环变量**取到的是字符串**，脚本内既有 `m=meth{1}` 是正确写法；
  不要再画蛇添足改成 `m{1}` 索引（会报 `sq_string cannot be indexed with {`）。改动别人的可跑代码前
  先确认它本来就错在哪一步。
- `s4d.beam` 是 spts×spts×nfreq（**已是 dB**：逐频点转 dB 后整体减 normdb，使最大值 = 0 dB）；
  ⚠️ **不要拿 `mean(beam,3)` 当慢度图出图**：dB 域等权平均会把弱峰频点的票数抬上来、峰位偏到低慢度那一档
  （见下「频带合并方式决定哪个峰当第一」）；要一张图就用振幅线性平均或直接 `fkmap`。第一维=南北 y，第二维=东西 x（`beam(i,j)` → `s4d.x(j)` 东西、`s4d.y(i)` 南北，且 `s4d.y` 是 fliplr 后从 +smax 到 -smax）。
  - ⚠️ **因为 beam 已经是对数（dB）量**：算"−3 dB 覆盖占比"必须用 `beam >= max−3`；沿用线性口径
    `max*10^-0.3` 会得到 0.00%。Python 端存盘的 beam 若是线性功率，占比则用 `beam >= max*10^-0.3`
    ——**两端口径不同，公式不能互拄**。
  - ⚠️ **频率轴内部处理（与 Python 引擎差异的根源）**：`fkvolume.m` 先 `nspts=2^nextpow2(npts)`
    把每道**零填充**再做 FFT，故同一 `FRNG` 内的频点数比"原始 npts 不填充"的 Python 多近一倍；
    且 **先逐频点转 dB、再对频率等权平均** ⇒ 噪声频点与信号频点同权。用它取峰时，若`FRNG` 内
    大部分是噪声频点（例：σ=8s 高斯包络 96% 功率在 0.03 Hz 以下而带取 [0.02,0.25]），峰会被噪声
    拉动 2-3 s/°；把带收到信号真正相干的段（如 [0.02,0.05]）后，两种做法给出**完全相同**的峰。
    核对方法：在 Python 里复刻这两条规则（填充到 2ⁿ + dB 域平均）应与 fk4d 存盘**逐值一致**。
    ⇒ **报 fk4d 峰位必须同时报 `FRNG`**，别用全带结果判"数据里有没有相干波前"。
  - ⚠️ 峰位对**台站顺序无关**（波束是对台站求和，coarray 是对所有 i<j 台对求和，都可交换）：
    打乱道顺序峰位不变（实测三种顺序逐位一致）。真正会错的是别的两件事：各道必须同一 `b`+ 同一
    `npts`，且**坐标数组必须与数据行同步重排**（只打乱波形不打乱坐标才会错）。

## fkmap / fkvolume / fk4d 的关系（用户会问"区别是什么"）

**同一估计量、三层封装**：`fkvolume` 是内核；`fkmap` = 内核 + 频带平均；`fk4d` = 内核 + 沿时间滑窗。
参数（`smax/spts/frng/polar/method`）三者含义完全相同。

| | `fkmap` | `fkvolume` | `fk4d` |
|---|---|---|---|
| 自述 | beamformer **map** in fk space | beamformer **volume** in fk space | beamformer volumes in fk-**time** space |
| 调用 | `fkmap(data,smax,spts,frng,polar,method,w)` | `fkvolume(data,...)` 同签名 | `fk4d(data,width,overlap,smax,spts,frng,polar,method)` |
| 输出 | 慢度 2 维（spts×spts） | 慢度 + 频率 | 慢度 + 频率 + 时间帧 |
| 频率维 | `frng` 内**平均掉** | 保留 | 保留 |
| 时间维 | 无（整段一次算完） | 无 | `width`/`overlap` 滑窗，`nframes=fix((1-w*o)/(w-w*o))` |
| 实现 | 自带一套与 fkvolume 相同的波束循环（**不调用**它） | 内核 | **每帧调一次 fkvolume**（源码注释：each column = a time range, each row = a frequency range） |
| 用途 | 单窗一张慢度图（论文 Fig.4 就这个） | 看慢度随频率 | 看慢度随时间 / 分开先后到达的相位 |

- `WIDTH=100` ⇒ 只取一帧 ⇒ **fk4d ≈ fkvolume，没有时间分辨**（本项目默认）；要时间分辨就调小 `width`
  （受"帧长/2 ≥ 台间最大延迟"限制）。
- 计算量 `fk4d ≈ nframes × fkvolume`；`fkmap ≈ fkvolume` + 一次频带平均（最省）。
- `fkvol2map` 把体转回图；交互看体 `fkframeslide / fkfreqslide / geofkslowslide`；出图 `plotfkmap / plotfkvolume`。
- ⚠️ 白化在 method 分支**之前**完成，`coarray` 的台对互谱幅度恒为 1 ⇒ **调用前对每道做整体缩放（`normalize`）
  不改变 coarray 结果**；`taper` / `removemean` 改谱形与低频泄漏，仍有影响。别把"加不加 normalize 结果一样"
  当成脚本写错。
- ⚠️ 报 `fkmap` 峰同样要带 `frng` 与 `nspts`（内部零填充到 2ⁿ ⇒ 带内实际频点数不等于 npts 直算），
  与上文 fk4d 的"频率轴内部处理"同一条规则。
- **为配图导出"送检波形"**（用户会要"画一下送检 f-k 的波形"）：在复现脚本里加一个导出块，把引擎调用前的那一叠
  记录（同窗、同采样、`dep` 字段）原样写文本 —— 头行给 `nsta npts delta 窗起点`，其后每行 `台名 + npts 个样点`；
  配图直接读这份 ⇒ 图上与引擎输入是**同一份数据**，不用另算近似。**导出块必须证明是纯附加**：重跑后引擎产物
  md5 逐字节不变（实测加块前后波束 md5 一致）。

### 频带合并方式决定"哪个峰当第一"（跨频带合并三种口径，读数完全不同）⚠️

- 源码实证：`fkmap` 在 `frng` 内做**复数（相干）叠加**、最后才取模
  （`beam(:)=beam(:)+exp(...)*cs(:,b)` 累加的是复数，收尾 `10*log10(|Re(beam)|/nfreq)`）；
  `fkvolume` 是**逐频点各自** `10*log10(|Re(A_b)|)`，跨频点不做任何合并。
  ⇒ `mean(体积, 3)` 只能得到**非相干平均**，永远不会等于 `fkmap`，两者甚至会指向**不同的相位**。
- 同一份数据 / 同窗 / 同参数（`smax=30, spts=300, frng=[1/1000 1/25]`, coarray）实测四口径：

  | 跨频带合并方式 | 峰 \|u\| (s/°) | BAZ (°) | 峰−中值 |
  |---|---|---|---|
  | `fkmap`（复数相干叠加） | 25.13 | 133.1 | 6.72 dB |
  | `fk4d` 单帧(`width=100`) + **振幅线性平均**（旧叫法"功率平均"）`10*log10(mean(10.^(B/10),3))` | 24.00 | 133.0 | 3.70 dB |
  | `fk4d` 单帧 + **dB 平均** `mean(B,3)` | 11.50 | 144.2 | 3.77 dB |
  | `fk4d` 单帧 + 取单个频点 | 35.60 | 302.6 | 5.75 dB |

  ⇒ 规则：**不要用 dB 平均把体积压成图**（等权投票让"数量更多的弱峰频点"胜出，峰位偏到低慢度那一档）；
  出一张图要么直接 `fkmap`（相干、prom 最高），要么做振幅线性平均。**报读数必须写明用的是哪一种合并**。

⚠️ **命名口径（不要再写"功率平均"）**：`10*log10(mean(10.^(B/10),3))` 且 `B = 10*log10|Re A_b|`，
展开就是 **`mean(|Re A_b|)` = 振幅（线性）平均**，既不是能量平均（那要 `mean(|Re A_b|^2)`），
也不是功率谱的平均。三个层次要分清：

| 层次 | 式子 | 性质 |
|---|---|---|
| 相干叠加（波场层，A_b 带相位） | `\|Re Σ_b A_b\|` | 模在求和号**外**；相位失配可相消，prom 最高 |
| 振幅线性平均（非相干） | `mean(\|Re A_b\|)` | 模在求和号**内**；恒非负、不可相消，固定丢 2–3 dB 对比度 |
| dB（对数）平均 | `mean(10log10\|Re A_b\|)` | 对数域等权投票，权重最偏，会把峰位颠倒到低慢度档 |

一句话判据：**模在求和号外 = 相干，模在求和号内 = 非相干**。两个层次之比就是相位相干性
（semblance / coherency）`C = \|Σ_b A_b\|² / (N·Σ_b \|A_b\|²)`（完全相干 1、完全不相干 1/N）。
脚本/数据文件名里的 `pwravg` 是**历史标签**，指的就是振幅线性平均（名字不改，脚本与文档都硬引用它）。
- **相识别靠方位，不靠慢度**：同一张 `fkmap` 图上可以同时有两个包 —— 实测 T-S 25.13 s/°@133.1°(0 dB)
  与 T-P 9.89 s/°@143.2°(−0.19 dB) 并存，分别对上论文 T-S_local 24.7 与**实测** T-P_local 10.0。
  判据用论文原话：**T-P_local 的方位 > T-S_local 的方位**（实测 143.2 > 133.1 ✓）。
  慢度单独分不开两相：单频点主瓣半宽 `0.64·111.195/(f·D)` 在 0.0067–0.0092 Hz、D≈313 km 时是 **27–37 s/°**，
  远大于两相 15 s/° 的间距；整带相干叠加把有效半宽压到带内最高频对应的 ~7 s/°，两个峰才分得开。
- 「一张图同时扫到两相」的三个条件：(1) `SMAX >= u_max + 主瓣余量`（u_max=24.7 时 30 是刚好、35 更稳）；
  (2) 频带里两相都有相干能量（实测低频 bin 0.0067–0.0092 Hz 搬 T-S、中频 bin 0.010–0.013 Hz 搬 T-P）；
  (3) 用相干或功率合并（dB 平均会把名次颠倒）。
- 逐频点判读（定性"哪一段频带在搬哪个相"）：对体积逐频点取峰位/电平/prom
  （`for b=1:nf; M=B(:,:,b); [~,i]=max(M(:)); ...`）。落在理论位附近的频点**必须同时看 prom** ——
  周期 > 窗长/3 的最低几个 bin 主瓣极宽（实测 prom 只有 0.2–1.5 dB），它们"位置对"没有说服力；
  `f > 111.195/(2·d_min·SMAX)` 以上的 bin 已进混叠区（判据推导见 `references/fk-parameter-selection.md`），
  也不能当证据。

### 非相干合并为什么出一大片"糊图"：低频子带的功率占比（归因顺序）⚠️

`fk4d` 体积做振幅线性平均时，**图的形状由包络功率最低的那个子带决定**。同一子阵/同窗实测：包络谱随频率急降，
把带内 64 个频点按每 8 个分 8 个子带，最低子带（0.0012–0.0055 Hz）一个就占平均功率的 **91%（南部）/ 94%（北部）**，
而该子带的主瓣半宽 `0.61·111.195/(f·D)` 在 f=0.007 Hz、D≈362 km 时是 **27 s/°**，比 ±30 s/° 的图幅还大 ——
它本来就"没有分辨能力"。振幅线性平均是逐频点取**非负**的 `|Re A_b|` 再在频率维平均，**无法像相干叠加那样让相位不对齐的频点
互相抵消**，于是这层宽平台被完整保留 ⇒ 图上就是一大团。
判据（可直接算）：**全带 fk4d 图的 −3 dB 等效半径 ≈ 最低子带自己的 −3 dB 等效半径**时，说明这张图其实在画那个子带。
实测北部全带 28.0 s/° vs 南部最低子带 28.1 s/°（几乎相同）；南部全带 19.3 s/° 则介于最低(28.1)与最高(5.2)子带之间 —— 
同一个机制在北部被放大。

**遇到"fk4d 图糊成一片"按这个顺序归因**（① 先排除几何，别先怀疑数据）：

1. **几何**：用台站坐标直接算平面波点扩散函数 PSF（`A = mean(exp(-2πi·f·τ))`，`τ = r·u/111.195`，r 用本地 km 平面坐标、
   u 用 s/°），比两个子阵在同一 f 下的 −3 dB 面积。实测南北几乎相同（0.02 Hz：169 vs 153 (s/°)²；d_min 94.7 vs 95.2 km、
   孔径 362 vs 364 km、方位覆盖都完整）⇒ **几何不是原因时，不要再拿孔径/台距说事**。
2. **合并方式造成的固定损失**：prom 与 −3 dB 面积在两种合并下的变化。实测 prom 南部 6.72→3.70 dB、北部 4.55→2.42 dB；
   −3 dB 面积南部 726→1172、北部 1310→2465 (s/°)²（等效半径 20.4→28.0）。**非相干平均固定削掉 2–3 dB**，
   本来就低对比度的子阵被削完就没有动态范围了。
3. **信号质量**（与几何无关的真实原因）：逐台 Welch 2–8 Hz SNR 对比子阵。实测北部中位 19.9 dB（14.1–24.5）
   vs 南部 24.8 dB（18.4–28.6），SNR_T 中位 22.9 vs 29.8 ⇒ "北部更糊"是信号相干性差，不是台阵坏。
4. **各子带功率占比**：对**包络**做 periodogram、按子带求和、归一化后算百分比。
   ⚠️ 直接拿 2–8 Hz 带通波形做谱再取 <0.04 Hz 的"低频子带"是**错的**（那个频段里根本没有信号，会得到 −100 dB 级假数字）：
   低频子带来自包络，必须先 `|hilbert|` 再算谱。
5. **混叠票数**：±SMAX 网格无混叠要求 `f <= C/u`（`C = 111.195/(2·d_min)`）。实测 64 个频点里只有 31 个满足
   `f <= C/SMAX`；对 T-S(26 s/°) 只有 36 个，而 T-P(13.8 s/°) 全带 64 个都合格 ⇒ **这种非相干平均天然给 T-P 更多票**，
   弱子阵上会直接把全图峰从 T-S 翻成 T-P（实测北部 25.90 s/°@163.1° → **13.81 s/°@168.7°**）。
   ⇒ 出慢度图优先 `fkmap`；非要用 fk4d 就**只保留 `f <= C/u` 的频点**再做振幅线性平均，且读数只报到能分辨的量级
   （BAZ 误差 ~Δu/|u|，本例几十度，别报小数）。
6. **子阵独立性先查再报**：老盘脚本的台号切片可能让两个子阵共用台站（本例 `[2:9,11]` 与 `[9:16,24]` 共用 SA31、SA37），
   "南北两个独立子阵"这类话必须先比对台号名单再说。

## 工具箱 helper 清单（用之前先读 header）⚠️

- `fk/fk4d.m` = 窗包装层；`fk/fkvolume.m` = 波束引擎。两者源码的中文注释版 = 事件目录
  `fk_code_octave/fk4d_source.md`（入口 + 被调核心 + 依赖函数清单 + 数值约定 + 关键词/公式推导 + 行号速查）。
- `fk/select_fk_peaks.m`（231 行）**不是自动拾峰**：先 `plotfkmap(map)` 画图，再让用户**鼠标框选**，
  每个框返回 `.xrng/.xtype/.yrng/.ytype/.db/.backazi/.horzslow/.meddb`；需要可交互 axes ⇒
  **无头批处理用不了**。讲到它时必须说明这点，别把"源码自带自动判据"写进笔记/论文。
- `fk/fkdbinfo.m`（181 行）= 真正可借用的一层：算整图或**限定范围**（`frng`/`bazrng`/`srng`）的
  **max/median/min dB**，max/min 附带 `.horzslow/.backazi/.freq`。⇒ 我们的"峰−中值 dB"显著度与它
  同口径（= `maxdb.db − meddb.db`），自己用 `max()/median()` 实现即可，不需要画图/鼠标。
  若要把"零慢度针尖"排除后再统计，照它的 `srng` 限带思路做慢度环带中值（本项目实现见
  `old_error/fk_code_octave/scan_frng_smax.m`，默认仍是更省事的全图统计）。
- 其它可借用：`slowbaz2kxy.m / kxy2slowbaz.m / slowness2deg.m`（慢度↔方位/单位换算）；
  `fkfreqslide / fkframeslide / fkmapslide`（自带出图，版式不可控 ⇒ 本项目一律 Python 出图）。

## Zhou 2021 关键参数

- 预处理：2-8Hz 带通 → Hilbert 包络 → 低通 0.25Hz → 降采样 2Hz（写 SAC 给 Seizmo 读，需写 stla/stlo/o/b 字段）
- 到时：论文 Equation 1 = `L_ocean/1.48 + L_land/6.0`(P) 或 `3.6`(S)。**勿用 TauP 陆地走时**（偏早 ~20s → 窗错位 → F-K 弥散；逐台量化、T-S 包络峰校准标尺与切换方式见 `references/theory-arrival-convention.md`）。
- 子阵：southern/northern 各 6-8 站，孔径<400km
- 结果：T-S 慢度 26-34 s/deg（论文 24.7，复现成功）；T-P 慢度测不到（T-P 在 Z 分量弱，论文也承认 T-P 只弱现于 Z/R）

## 合成验证（改引擎/流程后必做）

合成平面波（真实台站几何 + 已知慢度 10 s/deg + 2-8Hz 带通白噪声延迟）→ 完整 envelope 处理链 → F-K 应恢复 ~10 s/deg 且显著性 max/mean>3。验证通过 = 引擎+处理链对，问题在真实数据；失败 = 处理链有 bug。

## 配套衰减反演坑（coda normalization，Zhou2021 Q 反演）

- coda 窗**必须用固定 lapse time**（相对发震的统一时间窗，如最晚 tt_ts 之后 [1400,1460]s），不能取"每站 tt_ts+固定滞后"——tt_ts 随 land_km 增大时 coda 窗后移、coda 也衰减，比值 ln(A/Ac) 衰减被抵消 → 斜率变正 → **Q 反演成负值**（负 Q 第一排查项）
- 几何扩散 `+γ·ln(r)` 方向易反（r 600-980km 时 ln(r) 变化 0.49 ≫ 衰减变化 ~0.2，直接加会盖过衰减）；不确定就先不加
- Zhou2021 实测：Qs=2398 f^0.55（论文 685 f^0.53，**α_s 复现成功**），Q0 偏高 3.5×（几何扩散未修 + 论文 Text S2 精确窗未知）；Q 的 α 比 Q0 可靠

## envelope npz 本地修改（去某段突起）必须沿用原 time/t_off ⚠️

给 `sacprocess/envelope/{STA}_envelope_2-8Hz.npz` 做本地修改（如 taper 掉 1085~1115s 干扰段）时：

- **必须读原 npz 的 `time`/`t_off` 作基准**（process_image 存 time 相对发震，t_off=+600=发震相对数据起点，time[0]=-600）。
- ❌ 不要自己用 `sac.b - sac.o` 重算 t_off（0−600=−600，符号反）——会导致整条 time 轴错位 +1200s、waveform 相对时间全错，剖面图里该台波形全乱（2026-09-07 SA30 实测翻车，用户一眼看出）。
- 正确流程：读 `.bak_npz` 原 npz 的 time → 用 `np.searchsorted(time, 窗)` 定位目标段 → 对 `waveform_z/r/t` 乘 cosine taper（边缘 ~4s 平滑到 0）→ 用 taper 后 waveform 重算 `envelope_*`（|hilbert| + 6s 平滑）→ 存回，元数据（time/delta/t_off/band/smooth）沿用原版。改前 `cp npz npz.bak_npz`。
- 修改后验证：time 与原版 allclose、目标段波形≈0、窗外波形与原版零差异。
- 数据源修正后要**重跑全部用该 npz 的剖面图**（step1_profiles 3 张 + plot_unbinned 全站/snr6 各 8 张 + 定制全时段图），否则旧图仍带坏数据。
- ⚠️ 诊断"毛刺/突起"先逐样本看原始速度：看似尖峰可能不是孤立单点，而是 ~0.1Hz 大幅低频运动在 2-8Hz 拐点的泄漏（SA30 1100s 案例：raw 平滑、3 点中值残差 P98=0 因数据量化、真正原因是低频运动）——先弄清形态再决定 taper 段，不要当"毛刺"单点剔除。

## 本项目的三个引擎脚本（`fk_code_octave/`，形态统一：参数区在顶 + 逐段中文注释 + `cd`/`startup_seizmo`）

| 脚本 | 用途 | 关键参数（都在文件顶部参数区） | 产物 |
|---|---|---|---|
| `export_Tcomp_TS_beams.m` | 真实数据主链：相位(TS/TP) × 分量(Z/R/T) × 窗(A/B/D/G/R/R2) × method(center/coarray) | `SMAX=35`（= T-S 理论 24.7 × 1.4 留余量）、`SPTS=201`、`FRNG=[0.02 0.25]`、`WIDTH=100`（整窗一帧）、`OVERLAP=0`、`POLAR=false`、`BAZREF=126`（只用于打印偏差）；`TAGS={'G'}` 只跑某窗 | `fk_code_octave/output_data/beam_{tsc,tpc}_{Z,R,T}_<窗>[_coarray].txt` |
| `export_beams_cb_orig.m` | 「相干带 vs 原始带」对照：验证 `FRNG` 里 0.08–0.25 Hz 的噪声段是否在稀释峰 | 两组 `sets`：cb = `[0.005 0.08]`/smax=50、orig = `[0.02 0.25]`/smax=35；G1/G2 台站名单与窗写在脚本内 | `output_data/beam_{cb,orig}_{G1,G2}_{center,coarray}.txt` |
| `run_fk4d_bandscan_G1.m` | 分频带扫描：检验全带看到的峰是不是**空间混叠折叠像**（真平面波 ⇒ 三段同向；只在一段出现/方向跳变 ⇒ 混叠） | `SMAX=78`（刻意比常规大一倍以暴露折叠像）、`SPTS=441`（步长 0.354 s/°，与常规 0.35 一致才可与全带结果比）、`METHOD='coarray'`、`BANDS/FRNGS` 三段 b1/b2/b3 = `[0.005 0.03]/[0.03 0.08]/[0.08 0.25]` Hz | `fk_code_octave/bandscan/beam_G1_{b1,b2,b3}_coarray.txt` + `summary_G1_bandscan.csv`（逐段 nsig/slow/BAZ/峰−中值/−3dB 占比） |

- 三者都靠脚本顶部 `cd <seizmo 目录>; startup_seizmo;` 启动 ⇒ 任意 cwd 下 `conda run -n octave octave-cli --eval "run('<绝对路径>')"` 都能跑；cwd 不参与定位。
- ⚠️ 出图脚本（`plots/*.py`）按固定列序 `x(s/deg) y(s/deg) beam_dB` 解析这些 txt，**改列序会静默画错图**（不报错）。
- ⚠️ 提参数/加注释后的逐字节回归要**先确认 `/tmp` 备份里确实有该文件**再判差异（`diff -q` 对缺失文件同样返回非零）——`_coarray` 那几组可能此前从未跑过、根本没有基线，别把"备份不存在"报成"改动引入差异"。
