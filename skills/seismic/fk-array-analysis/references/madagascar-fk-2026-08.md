# 马达加斯加 Mw6.2 F-K 分析实录（2026-08-26）

## 数据与几何

- 事件: 2013-08-17 Mw6.2 SW Indian Ridge, evla=-34.887, evlo=54.128
- 25 站去响应速度 SAC（Z 分量）: `<event>/sacprocess/remove_response/*BHZ*_resp.sac`
- 孔径: 南北 1341 km × 东西 726 km；两两间距 min 36 / median 482 / max 1427 km
- 理论 BAZ: 137.8-169.2°（25 站阵列中心 159.8°；南部 6 站子阵 151.4°）
- T 波能量 2-8 Hz（SNR 4-6.5 dB）；低频 0.5-1 Hz SNR 为负（岛内台，高频转换波主导）
- 网络构成: XV×21（MACOMO 岛内）、GE×2、G×1（FOMA 塞舌尔）、II×1（ABPO）

## stations.csv 列索引（易错）

表头: `sta,az,dist,land_km,ocean_km,tt_p,tt_s,tt_tp,tt_ts,snr,lon,lat,baz`
→ **lon=idx10, lat=idx11, baz=idx12**（不是 idx11/12/10！首次写成 lon=11/lat=12/baz=10 直接 ValueError 纬度越界）

## 调试时间线

1. v1 全 25 站 + 40s 窗直接 F-K：300s 超时 → 优化（预计算 G、限 T 波窗 [min(tt_tp)-30, max(tt_ts)+60]、SPTS 101）
2. 结果乱（BAZ 286-307° vs 理论 160°）：三个问题叠加
   - peak_info 扇区掩膜 `atan2(SY,SX)` 参数顺序错（从东逆时针而非从北顺时针）
   - 全阵列波前穿越 ~380s ≫ 40s 窗 → 时间混叠（物理限制，非代码 bug）
   - 早期窗多数台站未到 T-P，白化后噪声相位主导
3. 合成验证（单频正弦）假失败：白化后其他 239 个频点全是随机相位 → 换宽带合成（带通白噪声 + 平面波延迟）
4. 合成验证（宽带）：4 站仍失败（栅瓣淹没真峰）；**6 站 + 120s 窗成功恢复 (153°, 3.6 km/s)**；8 站孔径拉大到 930km 假峰反超
5. 最终两路分析（对齐残差 + 子阵直接）通过

## 最终结果

| 分析 | 结果 |
|---|---|
| T-S 对齐残差 F-K（25 站） | 峰值≈原点 s≈0 → 转换点模型自洽 |
| T-P 对齐残差 F-K（25 站） | s=450-740 ms/km（T-P 窗混入 T-S 能量，FOMA 的 tt_ts-tt_tp 仅 6s） |
| 子阵 2-4Hz 堆叠 | BAZ=141.1° vs 理论 151.4°（偏差 10°） |
| 子阵 4-8Hz 堆叠 | BAZ=160.0° vs 理论 151.4°（偏差 9°） |
| 子阵初至窗 t=880s | **BAZ=153.4°**（偏差 2°），v=1.88 km/s |
| 子阵各窗 2-8Hz | BAZ 波动 125-174°，v 1.1-7.0 km/s（初至段最可信） |

## 输出文件

`<event>/process_flow/4_fk_analysis/`
- `fk2_aligned_tp_*.png` / `fk2_aligned_ts_*.png`（对齐残差慢度图，2-4/4-8/2-8 三频段）
- `fk2_subarray_*.png`（子阵堆叠慢度图）+ `fk2_subarray_baztime_*.png`（BAZ/v 演化）
- `fk2_*.csv`（峰值汇总/每窗峰值）

## 一句话结论

F-K 确认 T 波来自震中方向（初至窗 BAZ 153.4° vs 理论 151.4°），T-S 转换点模型在 25 站间自洽；视速度偏低是稀疏台阵 + 波前畸变的固有局限，F-K 对方向约束强于对速度大小的约束。
