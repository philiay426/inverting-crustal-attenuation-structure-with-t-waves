# 用 SEIZMO 自带函数跑 F-K（"纯源码"基线）+ Octave 兼容补丁

用户会要求「用纯源码跑一遍，不要自己写」—— 他要的是可复现的干净基线，用来排除
我们自写脚本引入的偏差。照这个顺序调用 SEIZMO **自己的**函数，不自写波束/绘图逻辑：

```matlab
cd('<seizmo 发行版目录>'); startup_seizmo;
addpath('<补丁目录>');                       % 见下面「Octave 兼容补丁」
data = readheader(files); data = readdata(data); data = synchronize(data);   % 各台 b 必须相同
d    = cut(data, o + win0, o + win1);        % 绝对计时：o = getheader(data(1), 'o')
s4d  = fk4d(d, WIDTH, OVERLAP, SMAX, SPTS, FRNG, false, method);
smap = fkvol2map(s4d, FRNG);                 % 体 → 面平均
fh = figure('visible', 'off'); ax = axes('parent', fh);
plotfkmap(smap, [-0.7 0], 'max', 'k', 'w', ax);   % (map, dblim, zerodb, fgcolor, bgcolor, h)
print(fh, '-dpng', '-r200', outfile);
```

- `fk4d(record, width, overlap, smax, spts, frng, polar, method)`；`method ∈ {center, coarray, full}`。
- `plotfkmap(map, dblim, zerodb, …)` 自带色标与色标尺：`dblim = [-0.7 0]` + `zerodb = 'max'`
  就是论文那种 0…−0.7 dB 的显示口径，不用自写色标。
- **必须把 axes 句柄传进 `plotfkmap`**：不传（自己 `figure()` 后直接 `print(fh, …)`）会报
  `no axes object in figure to print` —— 它会另开一个 figure。
- 把 `smap.beam` 原样导出成 txt（列序与我们既有 beam 文件一致）后，就能用同一个多峰表脚本
  （`scripts/peaks_table.py`）读它，与自写链路的读数逐项对比。

## 结果：与自写脚本逐位一致（G1）

| 配置 | 源码链路网格峰 | 自写脚本 |
|---|---|---|
| T-P 整带 0.008–0.0322、窗 1173–1292 s | 1.26 s/° @146° | 1.26 @146 ✓ |
| T-S 单频点 0.0078 Hz、窗 1243–1386 s | 35.31 s/° @138° | 35.31 @138 ✓ |

⇒ 自写 `export_*.m` 只是 SEIZMO 的壳。被问「是不是你们脚本写错了」时，跑一遍源码链路即可定案。

## Octave 兼容补丁（两处，少一就跑不起来）

1. **`interp1q` 未实现**（`plotfkmap.m` 用到）：放一个 `shims/interp1q.m`
   （`yi = interp1(x, y, xi, 'linear', 'extrap')`，语义同 MATLAB 的快速线性插值），`addpath` 后即可。
2. **`colorbar` 的位置参数顺序**：SEIZMO 写 `colorbar('eastoutside', 'peer', ax, 'fontweight', …)`。
   MATLAB 允许 LOC 在最前，**Octave 要求 LOC 必须是最后一个参数** ⇒ 改成
   `colorbar('peer', ax, 'fontweight', …, 'eastoutside')`。改前把原文件备份成 `.orig_matlab_order`。
   ⚠️ 只把 LOC 挪到第二位**不够**：后面还有 name-value 对时照样报
   `colorbar: LOC specification must occur as final argument`。

## Octave 里跑 SEIZMO 缺的东西（缺失就装，别改链路）

| 报错 | 缺什么 | 处理 |
|---|---|---|
| `'butter' undefined`（`iirfilter` 内） | Octave `signal` 包 | `conda install -n octave -c conda-forge octave-signal -y`；脚本里 `pkg load signal;` |
| `'csaps' undefined`（`getspline` 内） | Octave `splines` 包 | 同上装 `octave-splines` + `pkg load splines;`（本例装上前后输出逐项一致 ⇒ 该平滑步对结果无影响，但仍应装上以保持链路完整） |
| `'dfilt' undefined`（`iirfilter` 内） | **MATLAB 专有对象，Octave 无替代** | `iirfilter` 本体在 Octave 里**跑不通**，不要为此改 SEIZMO 源码。改在**预处理侧**做等价滤波：`butter(4,[f1 f2],fs)` + `filtfilt`（= 4 极 2 通带零相位，与原脚本 `'c',[2 8],'p',2,'np',4` 同规格），其余步骤（`synchronize/interpolate/normalize/envelope/removemean/taper/fkmap`）仍走 SEIZMO 原函数 |

- `pkg load signal/splines` 要在调用这些函数之前执行（放在 `startup_seizmo` 前后都行；不 load 就是 `undefined`）。
- 报「按原脚本参数复跑」时把这条替代写进交付，并注明"只换缺失实现、不改任何参数"。
- 用的是**单窗版入口 `fkmap`**（老脚本调的就是它；`fk4d` 只多时间维，单窗慢度图不需要）；
  `frng` 是**绝对频率 Hz**，参与平均的 bin 由 `Δf = 1/(2^⌈log₂ npts⌉·delta)` 决定（带→bin 换算见 `fk-parameter-selection.md` §0b）。
- ⚠️ 无图形环境下 `plotfkmap` 的 `print` 可能**静默不落盘**（脚本不报错但 PNG 不存在）⇒ 交付/读数图一律用 matplotlib 从导出的 beam 文本重画，数值（峰值/占比）从数据算，与图分开验证。

## 出图质量：源码出图只作对照，不作交付

本机 Octave 只装了 **gnuplot** 后端（`available_graphics_toolkits`），`plotfkmap` 出的是
约 1556×1167 的 PNG，线宽与字体都比 matplotlib 版粗。**出版级图仍用 matplotlib 脚本**
（色板/表头/±30 显示范围等口径见 `slowness-figure-paper-style.md`）；源码图只用于对照引擎行为
（色标口径、峰位、多峰结构）。

---

## 附：另一版记录（原 `references/seizmo-source-run-octave.md`，2026-09-28 整理时逐字并入）

# 用 SEIZMO 自己的源码跑 F-K（Octave 下）+ 与自写壳逐位对照

用户会说"**用纯源码跑一遍，不要自己写**"。意思是：**只用 SEIZMO 自带函数链**做完一遍，
顺带证明我们自写的 `.m` 只是壳、没有偷偷改数据。做法与必须的兼容补丁如下。

## 1. 源码链（全部是 SEIZMO 的函数）

```matlab
cd('<seizmo 目录>'); startup_seizmo;
addpath('<case>/shims');                          % 见 §2.1
files = {'<dir>/SA16_env_lp.sac', ...};
data  = readheader(files); data = readdata(data); data = synchronize(data);   % 各台必须同起点
o     = getheader(data(1), 'o');
d     = cut(data, o + WIN(1), o + WIN(2));        % 绝对计时窗；陆地计时输入 o=0 → cut(data, WIN)
s4d   = fk4d(d, 100, 0, SMAX, SPTS, FRNG, false, 'center');   % (data,WIDTH,OVERLAP,SMAX,SPTS,FRNG,POLAR,METHOD)
smap  = fkvol2map(s4d, FRNG);                     % 体 → 面（SEIZMO 自己的频率维平均）
plotfkmap(smap, [-0.7 0], 'max', 'k', 'w', ax);   % (map, dblim, zerodb, fgcolor, bgcolor, h)
print(gcf, '-dpng', '-r200', '<out>.png');
```

- `plotfkmap` 的 `dblim = [-0.7 0]` + `zerodb = 'max'` 就是项目要的色标口径；`fg='k', bg='w'` 给白底。
- driver 脚本 `pure_seizmo_G1.m` 是最薄的壳（读文件、切窗、调函数、落盘），**不含任何自写分析/绘图**。
- 跑法同既有约定：`conda run -n octave octave-cli --eval "run('pure_seizmo_G1.m')" > log.txt 2>&1`，
  在 log 里 grep 关键词取读数（**不要 `| tail` 管 stdout，会被 SIGPIPE 杀**）。

## 2. Octave 兼容补丁（不加就跑不起来，都不改 SEIZMO 的语义）

1. **`interp1q` 未实现**（`fk/plotfkmap.m` 用到）→ 在 `shims/interp1q.m` 写等价实现并 `addpath`：

   ```matlab
   function yi = interp1q(x, y, xi)
     yi = interp1(x(:), y(:), xi(:), 'linear', 'extrap');   % 含线性外推，同 MATLAB 语义
   end
   ```

2. **`colorbar` 的 LOC 位置**：MATLAB 允许 `colorbar('eastoutside','peer',ax,'fontweight',...)`，
   Octave 要求 **LOC 必须是最后一个参数** → 改成
   `colorbar('peer',ax,'fontweight','bold','xcolor',fc,'ycolor',fc,'eastoutside')`。
   `plotfkmap.m` 里有两处（约 218 / 382 行），改前 `cp` 一份 `.orig_matlab_order` 备份。
   报错原文是 `colorbar: LOC specification must occur as final argument`。
3. 本机 Octave 只装了 **gnuplot** 后端（无 qt/fltk）⇒ 源码出图是 gnuplot 渲染（约 1556×1167 px，
   白底黑字），比项目 matplotlib 脚本（300 dpi、论文实测色板）粗糙。**源码图只用于对照引擎行为**，
   出版图仍走 matplotlib 脚本。

其他 Octave 小坑：`string()` 不存在（用 `num2str`）；`figure('visible','off')` 后自己建 axes 再
把句柄传给 `plotfkmap`，否则 `print(fh,...)` 会报 `no axes object in figure to print`。

## 3. 源码跑的唯一目的：逐位对照

源码链与自写壳的读数必须**完全一致**，否则说明壳动了数据（切窗/同步/分量/时间基准）：

| 配置 | 源码（fk4d→fkvol2map 的网格峰） | 自写壳 |
|---|---|---|
| G1 T-P 1173–1292 s / 17.1 / 99 / 0.008–0.0322 | **1.26 s/° @146°，prom 10.49 dB** | 1.26 @146（一致） |
| G1 T-S 1243–1386 s / 42.2 / 243 / 0.006–0.009 | **35.31 s/° @138°，prom 9.90 dB** | 35.31 @138（一致） |

- 把 SEIZMO 的 map 原样导出成文本（列序与我们既有 beam txt 相同：`x y beam_dB`）就能直接用
  `peaks_table.py` 数多峰 → 见 `peak-reading-and-method-comparison.md`。
- 顺序：**先跑源码链，再拿自写壳对照**；只有两者一致时，才可以用自写壳的批量产物报结论。

