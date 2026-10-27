# T-P 前峰包偏振检查（pre-TP packet polarization）

事件：Madagascar 2013-08-17 Mw6.2，台站 VOI（2026-08-23 实测）。
触发场景：用户在 `{STA}_pca_theory_TP_TS.png` 或包络图上指出 T-P 拾取窗口前有一个峰包，要求看该峰包的偏振。

## 步骤

1. 载入三分量（同 `pca_theory_polarization.py`：2–8 Hz bandpass，zerophase, corners=4；包络 = |hilbert|，8s 平滑）。
2. 定位峰包：在 **T-P 窗口起点之前**的一段（VOI: 940~956s，即拾取窗起点 957.15 前）取包络 argmax → `tpk`。不要在 T-P 窗口内或全窗找。
3. 对 `tpk` 做 R-T PCA（`t_phase_analysis.pca`，1s 子窗 75% 重叠，analytic，baz 旋转，仅 R-T 平面）：
   - **±1s 峰值窗**（5 子窗）：表征峰包核心偏振；
   - **±3s 搜索窗**（21 子窗）：表征含两侧能量的混合偏振。
4. 判读：±1s 角集中（p25/p75 窄）→ 单相；±3s 角跨大半圆 → 峰核干净、两侧混合。R-T≈90°=横向（Love/SH），≈0°=径向（P/SV）。
5. 输出图 `{STA}_pca_preTP_peak.png`：2×2（三分量波形 / 包络+窗口 / ±1s 玫瑰 / ±3s 玫瑰）。玫瑰用半圆轴向折叠（0–180°，0/180=径向 R，90=横向 T），角度取 `arctan2(v[:,2], v[:,1])`（R-T 平面）。

## VOI 实测结果

峰包 `tpk=952.70s`（T-P 拾取 960.15 前 7.45s，理论 964.4 前 11.7s）：

| 窗口 | R-T 角中值 | 轴向判读 | lin | p25/p75 |
|------|-----------|---------|-----|---------|
| ±1s | 91.0° | 横向（T） | 0.995 | 80.5 / 95.0° |
| ±3s | −18.6° | 分散 | 0.999 | −103.2 / 80.5° |

- 与 T-P 自身的 −0.9°（径向）**相反** → 该峰包是独立横向（Love/SH 型）到达，不是 T-P 前驱；
- 距 T-S 拾取 977.6s 有 24.9s，与 T-S 无关；
- 科学后续：跨台站普查各站 T-P 前 ~10s 内是否普遍存在横向峰包，可判断是否为区域性 Love/SH 转换到达。

## 调试/绘图坑

- **去响应后包络量级 ~1e-8**（速度单位 m/s）：`print(f'{env:.1f}')` 显示 0.0 会误判"无峰"。打印峰值/包络先归一化（`envs/envs.max()`）或用科学计数法。
- **matplotlib 极坐标子图**：`fig.add_subplot((2,2,3), projection='polar')`（tuple 形式）在新版 matplotlib 报 `TypeError: Bbox.from_bounds() missing height`（tuple 被当作 Bbox）；必须分开传 `fig.add_subplot(2,2,3, projection='polar')`。玫瑰图统一用后者（`rose_rt(fig, nrows, ncols, idx, ...)` 模式）。

## 延伸：沿理论到时扫描最径向/最横向 → 固定窗口分析（2026-08-23 VOI）

用户在看完峰包偏振后问"沿理论到时搜索，最符合径向/横向的是多少s"——直接在终端 heredoc 里对理论中心 ±15s 以 0.25s 步长扫 ±3s 窗口（同 `run_polarization_guided_picks.py` 的 candidates 逻辑，评分函数可只算轴向残差）：

```python
def axial_rad(a): return abs(((a+90)%180)-90)     # 0°=径向
def axial_trans(a): return abs(axial_rad(a)-90)   # 90°=横向
# 对每个中心 c: pca(±3s 窗) → med R-T 角 → 径向/横向轴向残差
```

VOI 结果（窗口 ±3s、2-8Hz、lin 均 ≥0.998、SNR 38dB+）：

| 范围 | 最径向 | 最横向 |
|------|--------|--------|
| 理论 T-P 964.4 ±15s | **956.65s**（R-T −0.9°，残差 0.9°） | 975.40s（88.4°） |
| 理论 T-S 991.6 ±15s | 981.10s（10.0°，不纯） | **977.60s**（R-T 89.0°，残差 1.0°） |

- 搜索方案选的 T-P=960.15 与 956.65 残差同为 0.9°（960.15 是联合排序结果）；波形序列：~952s 横向峰包 → 956–960s 径向 T-P → 975–978s 横向 T-S。
- **用户要求"把时窗调为这两个时间做一个分析"**：在 `pca_polarization.py` 里只改 `T_TP_CENTER=956.65, T_TS_CENTER=977.60`（绝对时间中心），重跑 → T-P −0.9°（径向）/ T-S 89.0°（横向），lin 0.999/0.998，完美命中。
- ⚠️ 改脚本前先重读文件：用户可能在外部手动改过参数（本次实测 HALF_WIN 已被用户从 5.0 改成 3.0、中心改成 930/966）；`patch` 的模糊匹配会替换成我写的 old_string 语义，必须重读确认再动手。
- 窗口参数速查：当前脚本 ±3s（HALF_WIN=3.0）；`t_phase_analysis.pca` 包默认 window_length_s=10s/window_step_s=2s，**脚本均显式覆盖为 1.0/0.25**（75% 重叠）；搜索方案 SEARCH_HALF_S=15、STEP=0.25、HALF=3。

## 图标注偏好（用户 2026-08-23，pca 类图通用）

- **表头（suptitle）只标理论到时与实际选择时间**，不标子窗/重叠率/analytic 等叠加参数（子窗细节写代码注释即可）：`GE.VOI PCA  T-P th 964.4s → pick 956.7s  |  T-S th 991.6s → pick 977.6s`
- **理论/实际时间作为图例放左下角波形图**（2×2 布局的 ax3 包络面板）：`axvspan(..., label=f'T-P pick [{w0:.1f},{w1:.1f}]s')` + `axvline(ttp, ..., label=f'T-P th {ttp:.1f}s')` + `ax3.legend(loc='upper left', fontsize=9, framealpha=0.9)`（loc 用 upper left 避开 950–980s 峰区；色块图例直接对应时间轴，信息量最大）。
- 推荐做法：表头总览 + 波形图图例细节，两者都做；若只能留一处，留波形图图例。
