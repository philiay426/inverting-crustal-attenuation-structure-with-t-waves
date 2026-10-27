# Zhou2021 Fig.3 式记录剖面（3×2：三分量 × 两种距离轴）

复刻 Zhou et al. (2021) Figure 3：**3 列 = Vertical / Radial / Tangential，2 行 = 两种纵/横轴组合**。
脚本：`process_flow/1_station_distribution/profiles/plot_fig3_record_sections.py`（事件目录内）。
数据只读 `sacprocess/envelope/{STA}_envelope_2-8Hz.npz` 的 **`envelope_z/r/t`**（2–8 Hz 带通 →
|Hilbert| → 6 s 滑动平均；用户要求**只画包络、不叠波形**，无填充阴影）
+ `stations.csv`（`dist/land_km/ocean_km`），**不重算预处理**。

## 面板定义（照论文抄，不要自己发明）

| 面板 | 纵轴 | 横轴 | 到时短竖线 |
|---|---|---|---|
| 上排 (a)(b)(c) | Epicentral distance (km)，按 `dist` 升序排台 | Time (s)，**0 = 发震**，窗 **0–1600** | 红 T-P_local / 蓝 T-S_local 短竖线，位置随 `TT_MODE` 口径（默认投影走时 = `stations.csv` 的 `tt_tp`/`tt_ts`） |
| 下排 (d)(e)(f) | **Continental path distance (km)**，按 `land_km` 升序排台 | Time (s)，**0 = T 波到达大陆坡**：`x = time − ocean_km/1.48`，窗 −1000…+500 | 同口径自动换算：`tt − ocean/1.48`（= 纯陆地段走时；投影口径下 SA18 为 83.8 / 149.6 s） |

口径开关（两套到时差 20–25 s，见 `theory-arrival-convention.md`）：

```python
TT_MODE = "proj"    # 'proj' = 投影走时（stations.csv 的 tt_tp/tt_ts，陆地段 TauP）—— 当前默认
                    # 'const' = 论文式(1) 恒速 6.0/3.6
if TT_MODE == "proj":
    OUT = OUT.replace(".png", "_proj.png")   # 两口径产物分名，互不覆盖
```

下排两轴**不要写死速度**，用 `x_tp, x_ts = d["tt_tp"] − d["tocean"], d["tt_ts"] − d["tocean"]`，
两排随口径自动一致。用户可能要求把两套画在一张图上（对照用），做法见 SKILL.md 的 `--tt-both` 条。

- 论文 Fig.3 的配色与 Fig.2 一致：**红 = T-P_local，蓝 = T-S_local**（本项目剖面其它图用的
  蓝 T-P / 紫 T-S 配色在这里不适用，别照搬）。
- 布局：列标题 Vertical/Radial/Tangential **居中置顶**；面板标签 `(a)…(f)` 左上角；
  **每一列**（不是只最右列）都标 `T-P_local` / `T-S_local`；只最左列给 y 轴标签、每个面板都给
  `Time (s)`。
- ⚠️ 列标题与相位标签各占一行：轴坐标里标题置 `y≈1.11`、相位标签 `y≈1.00`，否则挤成一团
  （改版式后必须重跑 bbox 自查）。
- 画布用论文近似的横向比例（≈1.1 宽高比）；画完 `PIL` 核 `size` 比值。

## 🔴 逐台峰值归一化（不是 P98）

> ⚠️ **本章的线型/幅度/时间窗数值都是“某一轮”的定稿，不是固定规范** —— 用户会逐轮要求
> “变短/变淡/变细”、改幅度、拉长时间轴。复刻新案例时**取最近一次定版**（参数集中在脚本顶部
> “★参数区”），别照抄本文的历史值，也别把旧值当“论文原样”。
> 最新一次（马达加斯加 2013-08-17，abs+land 3×2）：到时短线 `SEG_FRAC=0.50`、
> 包络 `lw=0.25, alpha=0.65`、上排横轴拉到 1800 s（下排陆地计时窗不动）、
> `AMP_KM_DIST=50` / `AMP_KM_LAND=20`。

上排窗口 0–1500 s 必然含强 P/S。实测该窗口内**波形峰值 / P98 = 5.2（中位）～7.8（最大）**
⇒ 用 P98 归一化时 37/37 台的波形全都溢出相邻台站、叠成一团。
论文 Fig.3 用**逐台峰值归一化**，波形峰值正好等于其显示高度：

```python
amp  = AMP_KM_DIST        # 上排半高（km）；下排用 AMP_KM_LAND
norm = np.max(w[m]) or 1.0                             # 窗口内峰值（包络非负，不用 abs）
ax.plot(x[m], y0 + w[m] / norm * amp, color="black", lw=0.35, alpha=0.60, rasterized=True)
```

- **幅度用绝对 km**（用户会直接报数，不说"间距的几倍"）：`AMP_KM_DIST`（上排，本轮定稿 50）、
  `AMP_KM_LAND`（下排，25）——允许重叠（用户明说"允许重叠"，当前 ≈4× 中位间距）。
- 纵轴留白必须跟着幅度走：`pad = max(0.6 × 间距, 1.05 × 幅度)`，否则放大后边缘台站被裁。
- 线型定稿（用户逐轮提的"变短/变淡/变细"）：到时短线 `SEG_FRAC=1.05`（线长 = 1.05×间距）+ `lw=0.8, alpha=0.45`；
  包络 `lw=0.35, alpha=0.60`；列间距 `wspace=0.22`（用户："每一列图多空一点距离"）。
- 改成其它归一化（P98 / 全局）前先重算一次"峰值/基准"，确认 ≤1 才不会被裁。

## ⚠️ 台站分布不均 ⇒ 部分波形必然重叠（如实报告，别默默换布局）

论文台阵沿剖面基本均匀（间距 ~15–20 km），图面整齐；XA 37 站中位间距 dist 12.5 km / land 6.0 km，
但**最小间距只有 0–4 km**（2–3 台挤在几 km 内）。按 `0.40 × 中位间距` 定高时，
**dist 轴 ~15 对、land 轴 ~18 对**相邻台站波形会接触/重叠 —— 这是台站分布的性质，不是绘图 bug。

出图前/后跑一次一次性检查并**在交付说明里给出数字**：

```python
for name, v in (("dist", dist_sorted), ("land", land_sorted)):
    g = np.diff(v); amp = 0.40 * np.median(g)
    print(name, "min %.1f p10 %.1f median %.1f km" % (g.min(), np.percentile(g, 10), np.median(g)),
          "| 接触对数(间距<2*amp):", int((g < 2 * amp).sum()))
```

要"条条分明"时可切**等间距布局**（排序后固定间隔，刻度仍标真实距离值，记录剖面常规画法），
代价是纵轴不再是真实距离 —— **换不换由用户定，别擅自换**（用户偏好真实距离轴/不等间距纵轴）。

## 出图后自查（视觉模型关闭时靠这个）

1. `fig.canvas.draw()` → 取 `fig._suptitle` / 各 `ax.get_window_extent(renderer)` 的 bbox，
   判重叠与越界（矢量级 bbox 自查，见 SKILL.md `scripts/check_figure_whitespace.py` 条目）。
2. 核归一化后每条波形峰值 = `amp`（防再次溢出重叠）。
3. `ls -la --time-style=full-iso` 核产物 mtime（`/mnt/e` 被 Windows 看图软件占用会静默不更新）。
