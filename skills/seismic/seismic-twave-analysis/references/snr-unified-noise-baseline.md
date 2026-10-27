# SNR 口径：统一噪声基准（2026-10-14 用户定，现行口径）

## 一句话

SNR_dB = 20·log10( **信号窗 RMS / 该分量的阵列噪声中位数** )；
波形先做 **2–8 Hz 整段带通**，再按窗口切出来算 RMS。

## 为什么换掉「自噪口径」

旧口径的分母是**该台自身**噪声窗 RMS → SNR 排序 = 台站安静程度，不是信号强弱。
2013-08-17 Mw6.2 马达加斯加实测：各台自身噪声底相差 **384 倍**
（AMPY 7.88e-11 m ~ BITY 3.03e-08 m，2–8 Hz 去响应位移），
于是最安静的 AMPY/VOI 虚占 SNR 前二，而它们的信号振幅只排第 13/11 名；
信号最强的 MAHA/ZAKA/BITY 反被压到第 5/8/22 名 —— 用户指出"voi 和 ampy 还是异常，
这不符合传播的规律"。换共同分母后：SNR 排序 = 绝对信号振幅排序，
log10(信号 RMS) 与陆地距 r=−0.71（约 −4.2 dB/100 km），符合 T 波沿陆段衰减。

**判据**：任务是"哪台信号强 / 按信号强弱筛台"→ 必须用统一基准；
"这台自身有多干净"是另一个问题（保留对照列，别删）。

## 完整参数（`1_station_distribution/step1_map_and_table.py` 参数区）

| 项 | 值 | 说明 |
|---|---|---|
| 带通 | 2–8 Hz（`SNR_BAND`） | 必须**整段滤完再切窗** |
| 信号窗 | T 分量 `T-S±25 s`；Z/R 分量 `T-P±25 s`（`SNR_HALF=25`） | T 分量以 T-S 为中心；三分量窗长统一 50 s |
| 噪声窗 | `[tt_p-60 s, tt_p-10 s]`（`SNR_NOISE_LEN=50`, `SNR_NOISE_GAP=10`） | 右端比直达 P 早 10 s，避开 P 波 |
| 分母 | 各分量「噪声窗 RMS」在台阵上的**中位数** | 本事件 T 1.612e-09 m / Z 1.132e-09 m / R 1.301e-09 m |
| 阈值 | `SNR_THRESHOLD=6.0`，作用在 `snr_u` 上 | ≥6 dB 等价于 信号振幅 ≥ 2× 阵列噪声中位数 |

## 两个坑（都实测过）

1. **整段滤波再切窗**：只对 50 s 窗单独带通，窗内长周期分量在窗边产生滤波瞬变、
   噪声窗 RMS 虚高。FOMA 实测「切窗带通 5.5 dB vs 整段带通 28.7 dB」，差 **23 dB**，排名被带偏。
2. **噪声窗不能盖到 P 波**：原 `[P-50, P]` 吃到 P 波初动尾部，AMPY 因此虚低 **24 dB**
   （改 `[P-60, P-10]` 后 +23.9 dB）。噪声窗还要避开数据起点
   （本项目数据从发震前 600 s 起；ZE 数据集只有 300 s 前导，不能用同一窗）。

## 产物列结构（`1_station_distribution/`）

- `stations.csv` **19 列**：
  `sta,az,dist,land_km,ocean_km,tt_p,tt_s,tt_tp,tt_ts,snr,snr_z,snr_r,snr_t,snr_u,snr_uz,snr_ur,lon,lat,baz`；
  头 3 行注释写入统一噪声基准值（下游脚本据此解析分母，勿删）。
  主列 **`snr_u`**（统一基准 T 分量）；`snr/snr_z/snr_r` 是**自噪对照列**（T/Z/R）。
- `stations.xlsx` **18 列**：新增「统一信噪比(T分量)/统一Z信噪比/统一R信噪比」，
  标题同时给两口径的计数。
- 台站图按 `snr_u` 着色，colorbar 标签 `SNR (dB, common noise baseline)`。
- 本事件结果：≥6 dB → 统一 T **19** 台 / Z **17** 台 / R **15** 台（自噪口径 T 17 台）；
  top3 = MAHA(40.2) / BITY(28.8) / ZAKA(26.2)。

## 下游脚本怎么跟着用（2026-10-14 已改）

- `plot_unbinned_waveform_profiles.py`：新增 `compute_snr_unified()` ——
  **直接 import `step1_map_and_table` 调 `compute_rms_pair` / `snr_db`**，避免两处口径漂移；
  `noise_baseline()` 从 csv 头注释解析分母（失败退回常量）。
  `--snr-min`（comp=None）读 `snr_u` 列；`--comp z/r/t` 与南非台表走 `_resp_snr_components()`
  （resp.sac 原始波形 → 同一实现）。旧 `compute_snr_welch_time()` 保留但标注【已废弃】。
  验证：`--print-snr` 重算 25 台的 Z/R/T 与 csv 的 `snr_uz/snr_ur/snr_u` 最大差 **0.1 dB**
  （仅四舍五入），≥6 dB 计数完全一致（19/17/15）。
- 其它：`plot_island_waveform_qc.py` 的 QC 列表改打印 `snr_u`。

## 被本口径取代的旧写法（勿再引用）

- Welch PSD 2–8 Hz 频段积分 (`compute_snr_welch*`, `references/snr-computation-fix(es).md`,
  `record-section-plotting/references/per-component-snr.md`, `single-station-analysis/references/snr-source-unbandpassed.md`)
- 信号窗 `[T-P−20, T-S+100]`、噪声窗 P 前动态窗 `max(20, min(300, 0.8·tt_p))`
- "SNR 必须用去响应未带通波形" 的铁律 —— 现口径是"2–8 Hz 带通后算 RMS"；
  「未带通」只保留为**禁止双重带通**的注脚（npz 里 `waveform_*` 已是 2–8 Hz 缓存，不得再滤一遍）。
