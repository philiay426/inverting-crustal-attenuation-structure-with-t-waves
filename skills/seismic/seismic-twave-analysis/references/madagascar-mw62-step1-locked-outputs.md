# 2013-08-17 Mw6.2 马达加斯加事件：步1图件锁定配置

## 适用范围
仅适用于：
`/mnt/e/sharedate/case/Madagascar/2013-08-17-mw62-southwest-indian-ridge/`

这是一轮数据筛选后的当前输出，不应外推为其他马达加斯加事件或通用事件流程。

## 当前输出目录
`process_flow/1_station_distribution/`

## 步1两个脚本（2026-08-23 起，旧脚本已删除）

1. `step1_map_and_table.py` — 台站表格（`stations.xlsx` 18 列 + `stations.csv` 19 列）+ 两张台站分布图。
   参数在文件顶部：事件信息、SNR 阈值（6.0）、排除台站/排除框、经度范围（43–60E）、地图参数、CMT 解，
   以及 ★SNR 口径参数区（`SNR_HALF` / `SNR_NOISE_LEN` / `SNR_NOISE_GAP` / `SNR_BAND` / `SNR_THRESHOLD`）。
   用法：`python step1_map_and_table.py`（或 `--map-only` 只重画图）。
   SNR 口径（2026-10-14 起）= **统一噪声基准**：20log10( 信号窗 RMS / 阵列噪声中位数 )，
   2–8 Hz 整段带通后切窗；T 分量信号窗 T-S±25 s、Z/R 用 T-P±25 s，噪声窗 [tt_p-60, tt_p-10]。
   完整参数、两个坑、下游接法 → `references/snr-unified-noise-baseline.md`（含被取代的 Welch 频段积分口径）。
   当前结果：≥6 dB → 统一 T 19 / Z 17 / R 15 台（对照自噪口径 T 17 台）；top3 = MAHA / BITY / ZAKA。

2. ~~`step1_profiles.py`~~ **已删除**（当年生成的三张剖面图即下方锁定产物，仍保留归档）。
   现在在用的剖面脚本：`plot_unbinned_waveform_profiles.py`（未分箱，`--snr-min` 按统一噪声基准筛台）
   与 `plot_record_sections_abs_land.py`。勿再按 `step1_profiles.py` 的参数名（SNR_MIN/BIN_KM/AMP_BASE…）找文件。

## 台站图（两张）
由 `step1_map_and_table.py` 输出：

- `station_map_all_names.png`：标注当前入选的所有台站名；标签防重叠；只保留 land_km 最近和最远台站的两条灰色虚线；不标方位角数字。
- `station_map_top3.png`：同底图，只标注当前 SNR 最高的三个台站名。

## 剖面图（三张）
由 `step1_profiles.py` 输出，当前锁定使用：
- `profile_1_envelope_globalP98.png`
- `profile_2_waveform_global_normalized.png`
- `profile_3_waveform_full_origin.png`

它们与 `archive_non_madagascar/process_flow/1_station_distribution_all74/` 中指定的三个源图字节一致。改标准图后必须用 `cp` 同步归档源图（`profile_envelope_globalP98_black_Madagascar_land.png` / `profile_global_norm_waveform.png` / `profile_full_origin_distance_waveform.png`），否则下次复制链路会回退旧版。

## 数据筛选说明
“当前25个岛内台站”只是本次数据筛选结果，不是该事件固定要求，也不是通用规则。后续台站范围以用户最新筛选结果为准。
