# templates/ —— 步 0 预处理 + 步 1 台站表/剖面（本流程定版）

> 全部走环境变量，**不改代码里的路径**；必填缺失直接 `SystemExit`。解释器 `/home/sy/miniforge3/envs/obspy/bin/python`。
> 目录里另有 `process_image.py`、`step1_station_table_and_map.py` 等**旧版**单图脚本，仅供参照；本流程只认下表这 4 个。

| 脚本 | 作用 | 必填 | 可选 |
|---|---|---|---|
| `preprocess_one_station.py` | SAC → 去响应（速度）→ 2-8 Hz 带通 → 包络；出六张 1:1 图 + npz | `SRC`、`STA` | `OUT`、`FIG_OUT`、`ENV_OUT`、`FIG_STA`、`OUTPUT=VEL\|DISP` |
| `run_batch_preprocess.py` | 逐台 `subprocess` 批量调用，日志 `batch_log.txt` | `SRC`、`STA_LIST=A,B` | `OUT` |
| `step1_map_and_table.py` | `stations.csv` + `.xlsx` + 台站分布图（pygmt，600 dpi） | `T_BASE` | `T_OUTDIR`、`T_MAPDIR`、`T_CSV`、`EVENT_*`、`EXCLUDE_STA`、`LON_RANGE`、`CMT_*`、`MAP_*` |
| `step1_record_sections_3x2.py` | 三分量 3x2 记录剖面（上排震中距 / 下排陆段） | `T_BASE` | `T_OUT`、`T_CSV`、`T_ENV`、`TT_MODE=proj\|const` |

**交接物（下游契约）**

- `envelope/<STA>_envelope_2-8Hz.npz`：`time`、`waveform_z/r/t`、`envelope_z/r/t`
- `stations.csv`：`sta, lat, lon, dist_km, land_km, ocean_km, tt_p, tt_s, tt_tp, tt_ts, snr_db, baz`

**口径**：时间轴 `t = b - o`；速度 m/s；`land_km` 是投影陆段（不是震中距）；SNR 口径全表只准一种。
