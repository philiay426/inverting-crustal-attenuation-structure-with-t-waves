# templates/ —— 步 4 衰减反演（本 skill 只保留这一步）

步 0/1 的模板在 `seismic-twave-analysis/templates/`，步 2 在 `single-station-analysis/templates/`，
步 3 在 `fk-array-analysis/templates/`。本目录只有衰减计算与它的图。

| 脚本 | 作用 | 必填 | 可选 |
|---|---|---|---|
| `twave_q_inversion.py` | **核心引擎**：逐台逐频点振幅 → 逐频点斜率 `k` → `Q(f)` → `Q0/α`（1σ） | `Q_CASE` | `Q_STAS`、`Q_CODA_MODES`、`Q_AMP=rms`、`Q_HM=70`、`Q_FREQ_EXCL`、`Q_TAG`、`Q_PCOMP=z` |
| `plot_q_inverse.py` | 定案图：`Q^-1` 对 f + 拟合线 + `α` 注记 | `Q_CSV` | `Q_DETAIL`、`Q_OUTPNG`、`Q_FREQ_EXCL` |
| `plot_coda_window_section.py` | 尾波窗落点（绝对计时剖面） | `Q_STAS`、`Q_BASE` | `Q_CODA_MODE`、`Q_ENV_DIR`、`Q_FIG_DIR` |
| `plot_record_section_landtime.py` | 陆地计时 / 绝对计时剖面（含窗与理论到时） | `Q_STAS`、`Q_BASE` | `Q_ENV_DIR`、`Q_FIG_DIR` |
| `compare_q_objects.py` | 多对象 `Q^-1` 对比（自己的 / 别的台阵 / 文献） | `Q_OUTPNG` | 曲线参数填进脚本顶部 `CASES` |
| `select_station_set.py` | 台组子集枚举筛选（纯算术 + 引擎产物） | `Q_BASE`、`Q_POOL` | `Q_SCAN`、`Q_ENGINE` |

**输入契约**：去响应后的**速度** SAC（`remove_response*/` 或 `.transed_<f>` 成品）+ `stations.csv`（含 `tt_tp/tt_ts`）。
`Q_BASE` = 事件的 `process_flow/` 目录（图脚本读其下 `1_station_distribution/stations.csv`）；`Q_DETAIL` 不给时由 `Q_CSV` 的 `_byfreq` 自动换成 `_detail`。
**输出契约**：逐台明细 + 逐频点 `f_hz, Q, sigma_Q, k, se_k, n_used`。

两个硬约束：取窗索引 `i0 = (t − b)/dt`（扣 `b`）；多事件直接 pooled OLS，不扣事件常数、不剔点。
