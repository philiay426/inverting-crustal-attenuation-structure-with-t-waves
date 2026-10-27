# ZNE→ZRT 旋转约定统一到 obspy（2026-09-02 全项目）

## 背景

排查"带通图 P 波被截断"时顺带发现：马达加斯加事件整个项目（process_image.py、
lib/rotation.py、t_phase_analysis fallback、各台站单台脚本、F-K 脚本、templates）
的 ZNE→ZRT 旋转都用**手写公式**，它与 obspy 标准约定差 **180° 反相**。

## 两套约定

| | 手写公式（旧，全项目） | obspy `rotate_ne_rt`（标准） |
|---|---|---|
| R | `R = +N·cos(baz) + E·sin(baz)` | `R = -N·cos(baz) - E·sin(baz)` |
| T | `T = -N·sin(baz) + E·cos(baz)` | `T = -E·cos(baz) + N·sin(baz)` |
| R 正方向 | 指向震源 (toward source) | 背离震源 (away from source) |

关系：`项目(R,T) = −obspy(R,T)`（整体反相 180°）。

## 影响矩阵

| 分析类型 | 是否受影响 |
|---|---|
| 能量/SNR/包络剖面/频谱 | **无影响**（\|振幅\|² 不变） |
| R/T 波形显示 | 形态相同，整体反号 |
| 偏振 R-T 角 atan2(T,R) | 轴向 %180 处理 → 无影响；带符号 %360 → 差 180° |
| 玫瑰图/方向角（带符号） | 差 180° |

T 波偏振引导拾取/BAZ 反演用 %180 轴向判别，故**历史结果不受影响**——
当时能得出正确结果正是因为这个约定。

## obspy 陷阱

- **obspy 没有 `rotate_zne_rt`**（只有 rotate_ne_rt / rotate_rt_ne / rotate2zne /
  rotate_zne_lqt）。`from obspy.signal.rotate import rotate_zne_rt` 必然 ImportError
  → 静默走手写 fallback，掩盖约定不一致。
- 统一做法：**直接调用 `rotate_ne_rt(n, e, baz)`**，不要手写公式。

## 逆变换（pca.py `_horizontal_back_azimuth`）

obspy 约定下由 (u_r, u_t) 反算地理 (u_n, u_e)：

```
n = -r·cos(baz) + t·sin(baz)
e = -r·sin(baz) - t·cos(baz)
```

注意正负号必须与 rotate_ne_rt 匹配；用 %180 轴向处理可消除特征向量 ± 歧义。

## 验证方法（可复用）

1. 合成纯径向/纯切向信号验证分离：纯径向（沿 baz 方向）→ R 保留、T≈0。
2. 正变换对比 `np.max(|project_R - obspy_R|)` 应 = 0。
3. 逆变换自洽：纯 ±R 经 %180 后都还原 baz_true。
4. 缓存与直接旋转相关 = 1.0（注意三分量长度对齐）。

## 已修改文件

- `/home/sy/codex/lib/rotation.py` — `rotate_zne_to_zrt` 改调 rotate_ne_rt
- `/home/sy/codex/t_phase_analysis/t_phase_analysis/pca.py` —
  `rotate_stream_to_zrt` 改调 rotate_ne_rt；`_horizontal_back_azimuth` 逆变换修正
- 事件目录：`process_image.py`（get_zrt/plot_bandpass/plot_spectrogram 三处）、
  `step1_map_and_table.py`、`plot_spectrogram_sections.py`、
  `plot_unbinned_waveform_profiles.py`、`fk3_3comp.py`、`fk6_summary_map.py`
- 模板：`/home/sy/codex/templates/single_station/*`（3 个）——防新事件带回旧约定
- 备份：`lib/rotation.py.bak_pre_obspy`、`pca.py.bak_pre_obspy`

重跑清单：process_image 全 25 站 → 批量偏振拾取 → 单台偏振图 → R/T 剖面图/
频谱剖面图。

---

## stations.csv 解析陷阱（勿用固定列索引）

stations.csv 自 2026-08 起为 **15 列**：
`sta,az,dist,land_km,ocean_km,tt_p,tt_s,tt_tp,tt_ts,snr,snr_r,snr_t,lon,lat,baz`

⚠️ 旧 13 列格式 `p[12]` 是 baz；15 列格式 `p[12]` 是 **lon**（baz 在 p[14]）。
任何按固定索引解析的旧代码会把经度当 baz → R/T 旋转全错（曾实际发生，
process_image `load_stations` 中 p[12] 取到 lon=46.79 而非 baz=151.6）。

**必须用 `csv.DictReader` 按列名解析**（跳过 `#` 注释行），参考：

```python
import csv
with open(csv_path) as f:
    for r in csv.DictReader(line for line in f if not line.startswith('#')):
        baz = float(r['baz'])   # 不是 p[12]！
```

---

## process_image bandpass 纵轴 P 波截断修复

`plot_bandpass` 纵轴原用"窗口内 99% 分位"定标，但 P 波初动常达
6–12×ep98（远大于 1.5×ep98 的 99% 分位）→ **P 波被完全截断**。

修复：`y_max = 1.05 × max(99% 分位, P/S 窗内峰值, 1.0)`
（P/S 窗 = [tt_p−10, tt_s+60]，三分量取最大）。代价是 T 波被相对压小
（占 15–40%），可接受；如要 P/T 同轴兼顾可用断轴/分窗。
