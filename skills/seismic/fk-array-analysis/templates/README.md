# templates/ —— 步 3 子阵 F-K（本流程定版）

> 引擎是 Octave/SEIZMO（`fk4d`）；本目录只放“算几何 + 出图”。必填缺失 `SystemExit` / `error`。

| 脚本 | 作用 | 必填 | 可选 |
|---|---|---|---|
| `fk_geometry.py` | **先跑这个**：台站坐标 + `V_T` → `d_min`/孔径/`C`/`u_pred`/`smax`/`spts`/分辨率地板/理论 BAZ，打印并写 `fk_geometry.csv` | `FK_STAS`、`FK_VT`、`T_CSV` 或 `FK_SAC`（其一） | `FK_SMAX`、`FK_SPTS=301`、`FK_FMAX=0.04`、`EVENT_LAT/LON`、`T_OUTDIR` |
| `step3_fk_driver.m` | Octave 驱动：SAC 目录 → `beam_<TAG>.txt` + `fkinput_<TAG>.txt` | `FK_SEIZMO`、`FK_DATA`、`FK_TWIN="W0,W1"`、`FK_VT` 或 `FK_SMAX` | `FK_C`、`FK_FMAX`、`FK_TAG`、`T_OUTDIR`、`FK_SPTS`、`FK_FRNG`、`FK_SPS=20` |
| `step3_fk_demo.py` | 读波束网格出慢度图 + `fk_peaks_<TAG>.csv`，列全部局部极大并打印几何判据 | `T_BASE`、`FK_STAS` | `FK_TAG`、`T_OUTDIR`、`EVENT_LAT/LON`、`FK_WIN`、`FK_VT` 或 `FK_SMAX` |
| `step3_fk_subarray_map.py` | 子阵分布图（子阵 + 沙滩球 + 大圆射线，600 dpi） | `T_BASE`、`SUB1` | `SUB2`、`T_OUTDIR`、`EVENT_*`、`CMT_*`、`REGION` |

**参数铁律**：`smax = 1.4·u_pred`、`u_pred = 111.195/V_T`、`f_hi = min(F_MAX, C/smax)`、`C = 111.195/(2·d_min)`。
换地区只改 `FK_VT`（视速度）与 `FK_TWIN`（时窗），不许手填 `smax` 与 `frng`。
