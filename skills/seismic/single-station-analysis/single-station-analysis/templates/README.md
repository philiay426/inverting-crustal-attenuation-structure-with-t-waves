# templates/ —— 步 2 单台分析（本流程定版两张图）

> 环境变量进、路径不进代码；必填缺失 `SystemExit`。

| 脚本 | 作用 | 必填 | 可选 |
|---|---|---|---|
| `polarization_rose_geo.py` | 两个地理方位玫瑰（T-P / T-S）+ Z/R/T 包络时窗；PCA 滑窗 | `SRC`、`STA` | `CSV`、`NET`、`OUT`、`HALF=25`、`PWIN=1`、`PSTEP=0.25`、`FMIN/FMAX`、`EVENT_TIME`、`PAD_PRE/PAD_POST`、`SMOOTH_N` |
| `polarization_fullscan_3panel.py` | 竖排三栏（R 占比 / T 占比 / 线性度+平面度）全时段扫描 | `SRC`、`STA` | `CSV`、`NET`、`OUT`、`PSTEP=0.25`、`FMIN/FMAX`、`HALF` |

目录里另有 `polarization_analysis.py`、`pca_polarization.py`、`full_scan_polarization.py` 等**旧版**，
本流程只认上表两个（PCA 滑窗 = `--method cov` 口径）。

**输出契约**：可用台清单 + 每台剔除理由；`rt_angle_med_deg` 有 180° 模糊，只报 `%180` 轴向角。
