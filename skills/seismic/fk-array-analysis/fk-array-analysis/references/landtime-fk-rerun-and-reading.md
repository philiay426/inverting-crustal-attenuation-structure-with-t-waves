# 陆地计时 F-K：口径复制到其他分量、复算与判读

适用场景：`<case>/process_flow/4_fk_analysis/<子阵>/04_landtime_fk/` 这类「陆地计时（0 = 上岸）」目录。
入口 `export_landtime_fk.m`（★配置区 `C(1)…C(n)` + `RUN_TAGS` 过滤）× 出图 `plot_landtime_map.py`
（`PANELS` 单格 / `WIDE_SETS`+`LAYOUT=wide3` 并排）× 独立方向检验 `xcorr_direction_check.py`。

---

## 1. 「某一个分量的口径，做到其他分量上」

1. **复制原配置行，只改 `comp`，tag 保持不变**，并在 ★配置区尾部注明「口径与 X 分量一字不动」：
   ```matlab
   C(31) = struct('tag','ts120_270','phase','T-S','comp','R','FRNG',[0.005 0.013],'win',[120 270],'SMAX',45.0,'SPTS',257);
   C(32) = struct('tag','ts120_270','phase','T-S','comp','Z','FRNG',[0.005 0.013],'win',[120 270],'SMAX',45.0,'SPTS',257);
   ```
   波束文件名 = `{tag}_{comp}[_coarray|_full].txt` ⇒ **同 tag 不同 comp 不会互相覆盖**；
   会互相覆盖的是「同 tag 同 comp 的两条配置」（例如同一 tag 下先后写 T-S T 与 T-P T，后者盖前者）。
2. `plot_landtime_map.py`：`PANELS` 加对应单格；`WIDE_SETS` 加一组「同口径 N 分量并排」
   （suptitle 从集合的 `title` 键取，别把标题硬编码在 `fig.suptitle` 里，否则新组会顶着旧标题）。
3. `MEAS_BAZ` 补 `(tag, comp)` 的实测方向，否则该格回退画 F-K argmax 黑实线 —— 同一张图里各格方向线来源不同源。
4. 复算（实测 54 s：3 分量 × 3 method，SPTS 257）：
   ```bash
   cd .../<子阵>/04_landtime_fk
   cp -p summary_landtime_fk.txt /tmp/summary_landtime_fk.before.txt          # ★ 先备份，见 §2
   conda run -n octave octave-cli --eval "RUN_TAGS={'ts120_270'}; run('export_landtime_fk.m')"
   /home/sy/miniforge3/envs/obspy/bin/python xcorr_direction_check.py         # CONFIGS 里补同窗同带的 R/Z 行
   /home/sy/miniforge3/envs/obspy/bin/python plot_landtime_map.py             # 单格（放最后跑，见 §4）
   LAYOUT=wide3 WSET=ts_rtz /home/sy/miniforge3/envs/obspy/bin/python plot_landtime_map.py
   ```
   （非交互 shell 里 `conda` 不在 PATH ⇒ 写全路径 `/home/sy/miniforge3/bin/conda run -n octave octave-cli`。）
5. 收尾：把本次运行日志按「只留读数行」过滤后写进目录（`log_<入口>_<主题>.txt`）；README 追加一节
   （新配置号、复算命令、读数表、产物清单），并同步「最终保留集」表里的图名与 `C(1)…C(n)` 上界。

---

## 2. `RUN_TAGS` 会整文件重写 summary —— 备份 + 按键合并

入口脚本用 `fopen(summary_landtime_fk.txt,'w')` ⇒ 只跑 3 组时，原有 18 行被洗成 9 行。
合并键必须**正则解析**，不要按空格切 token：

```python
import re
def key(l):
    comp = re.search(r"\b([RZT])\s+\d", l).group(1)        # 分量紧邻 FRNG 数字
    meth = re.search(r"\b(coarray|center|full)\b", l).group(1)
    return (l.split()[0], comp, meth)                        # tag, comp, method
```

**为什么不能按 token 位置取**：phase 列含空格（`T-S 两能团`、`总 -T-P`），而 method 前一个 token 是 `SPTS`，
所以 `t[m-1]` 拿到的是 SPTS 而不是 comp ⇒ 同 tag 同 SPTS 的行被误判为同一键，去重后每个 method 只剩最后一行。
备份在 `/tmp`，跑完按 key 合并回写（新的优先），并核对行数 = 旧 + 新。

---

## 3. 判读：argmax 一张图不够，要过方向锥 + 独立方向检验

**① 方向锥再查。** `argmax` 常落在近零慢度的「针尖」或非震源方向；报读数前在
「理论 BAZ ±15°、慢度 ≈ 目标区间」的锥内再取一次最强值，并给出它相对全图峰的 dB 差：

```python
U = np.hypot(XX, YY); AZ = np.degrees(np.arctan2(XX, YY)) % 360
daz = np.abs((AZ - BAZ_TH + 180) % 360 - 180)
cone = (daz <= 15) & (U >= u_lo) & (U <= u_hi)
```

这样能把「该分量确实有源向峰，只是比另一峰低 x dB」与「根本没有源向峰」区分开 —— 前者不要写成「没信号」。
同一张图上再对镜像方向（BAZ+180）的锥取一次最强值：镜像弱 2–5 dB 说明该引擎的 coarray 波束**实测不严格 ±u 对称**，
方位因此可读；镜像等强时方位不可读，只能靠滞后拟合定方向。

**② 独立方向检验（`xcorr_direction_check.py`，台对互相关滞后最小二乘）。**
判据：拟合 rms < 15 s 才用它画方向线，否则 `MEAS_BAZ[(tag,comp)] = None`（不画实线），并看
镜像候选 rms / 拟合 rms 的比值（≫1 才说明方向唯一）。
**同一个分量换窗就会从「废拟合」变「可用」**：T-S 120–270 s 窗 R/Z 的 rms 7–12 s（比值 5–9×，可用），
而 60–260 s 长窗与 T-P 窗的 R/Z rms 30–49 s（比值 ≈1，方向不可定）。⇒ 说「某分量定不了方向」时必须带窗和带。

**③ 角度用 %180、BAZ 与理论 137.3° 比残差**；报告里同时给出慢度（与陆段理论 111.195/v 比，T-P 18.53、T-S 30.89 s/°）。

---

## 4. 出图脚本的 summary 只反映最后一次运行

`plot_landtime_map.py` 在两种布局末尾都会以 `'w'` 写同一份 `summary_landtime_maps.txt`：
`LAYOUT=single` 写单格那几行，`LAYOUT=wide3` 只写并排的 3 行。若单格集是并排集的超集，**最后跑 `LAYOUT=single`**，
汇总才完整（并排的读数与对应单格逐字相同，不会丢信息）。

产物自查（脚本已打印，逐条核）：色标 `pc.x1 <= 1.0`（在画布内）、`pc.x0 >= 各主图 x1`（不与主图重叠）、
刻度 = 该格 SMAX 的 10 s/° 倍数；PNG 尺寸单格 1920×1770（6.4×5.9 in @300 dpi）、并排 4050×1620（13.5×5.4 in @300 dpi）。
写 `/mnt/e` 共享盘后核对 mtime（Windows 侧占用会静默失败）。

---

## 5. 速度/合法性自检（换窗换带必过）

- 合法频点：窗长 T ⇒ `nsp = 2^ceil(log2(T/dt))`、`Δf = 1/(nsp·dt)`，带内只留 `f ≥ 1/T` 的栅格点。
- 抗混叠上限：`SMAX ≤ 0.550231/f_hi`（`f_hi` = 带内最高参与频点）。
- 例：150–200 s 窗、dt 0.5 s ⇒ nsp 512、Δf 0.00390625 ⇒ 带 0.005–0.013 Hz 只剩 **0.0078125 / 0.0117188** 两档；
  带收到 0.005–0.0095 或 0.006–0.0090 就只剩 **0.0078125** 一档。
