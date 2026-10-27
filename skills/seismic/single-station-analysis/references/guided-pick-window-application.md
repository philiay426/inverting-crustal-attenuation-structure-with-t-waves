# 偏振引导搜索窗口回填单台主图（guided-pick → pca_theory_polarization.py）

事件：Madagascar 2013-08-17 Mw6.2（2026-08-23 实测）。
搜索方案脚本：`process_flow/2_single_station_analysis/run_polarization_guided_picks.py`。
结果 CSV：`.../polarization_guided_picks/polarization_guided_tp_ts_picks.csv`。

## 搜索方案参数（脚本内硬编码）

| 参数 | 值 |
|------|-----|
| 理论先验 | T-P/T-S 各自理论中心 ±15 s |
| 滑窗步长 | 0.25 s |
| 窗口半宽 | ±3 s |
| PCA 子窗 | 1.0 s / 75% 重叠 / analytic / ZNE→ZRT(baz) / 仅 R-T |
| 门槛 | 线性度 ≥0.8；波包 SNR ≥10 dB；角残差 ≤20°；spread p75 小 |
| 物理顺序 | T-S ≥ T-P + 3 s（违反则 order_fail，禁止反序对） |
| 评分 | score = max(0,1-res/90)*3 + lin*2 + min(max(snr,0),30)/30 − spread/90 |
| 联合选对 | 按 (res_tp+res_ts, spread_tp+spread_ts, −lin 和, −snr 和) 字典序最小 |

波包 SNR：2–8 Hz Welch 能量比，信号窗 center±3s，噪声窗 −300~−60s（相对发震）。

## CSV 窗口是绝对时间，模板内部是陆地走时

`pca_theory_polarization.py` 模板：
- 内部窗口变量 w_tp/w_ts 是**陆地走时坐标**（t=0 = 波上岸），slice 时 `w_abs = w + TT_OCEAN`；
- 图 ax3 的 x 轴是**相对发震绝对时间**（`t = off + np.arange(npts)*delta`）。

回填搜索窗规则：
1. 定义 `PICK_TP_ABS / PICK_TS_ABS`（CSV 的 window_start/end，绝对时间）；
2. PCA slice 用 `w = PICK_ABS - TT_OCEAN`（保持模板 run_pca 内部 +TT_OCEAN 不变）；
3. axvspan/axvline 用 `PICK_ABS` 直接画——⚠️ 原模板曾把陆地走时值（~35/63s）画到 900–1000s 信号轴上，色块完全错位，必须用绝对时间。

## VOI 实例（2026-08-23）

stations.csv: tt_tp=964.4, tt_ts=991.6, ocean_km=1375 → TT_OCEAN=929.05s。

| 震相 | 理论 | 拾取中心 | 偏移 | 窗口(绝对, s) | 陆地走时窗口(s) | R-T 角 | lin |
|------|------|---------|------|--------------|---------------|--------|-----|
| T-P | 964.4 | 960.15 | −4.25s | [957.15, 963.15] | [28.1, 34.1] | −0.9°(径向) | 0.998 |
| T-S | 991.6 | 977.60 | −14.0s | [974.60, 980.60] | [45.6, 51.6] | 89.0°(横向) | 0.998 |

结论：VOI T-S 拾取比理论早 14s，理论后移12s 窗（≈1000.6–1006.6s 绝对）完全落在实际横向波包之后——验证"不得把 +12s 机械套用"规则；凡搜索方案给出大偏移，主图必须用搜索窗重画。脚本输出 R-T 角/线性度与 CSV 逐位一致即视为自检通过。
