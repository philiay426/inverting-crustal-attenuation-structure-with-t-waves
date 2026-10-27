# 台站分布图最新改进 (2026-08-03)

## 排除框过滤

大尺度地图常有大面积低质量台站簇。用经纬度矩形框排除：

```python
EXCLUDE_BOXES = [
    [32.0, 37.5, -11.5, -7.0],  # 10°S 北部密集噪声簇
]

def _excluded(lon, lat):
    for lo, hi, la, ha in EXCLUDE_BOXES:
        if lo <= lon <= hi and la <= lat <= ha:
            return True
    return False
```

排除后的台站不进入 CSV/xlsx/剖面/地图，保持输出一致。
剔除的 raw 和 sacprocess 数据也应清理。

## 标签防重叠算法

每个被标注台站尝试 4 个方向（优先远离震中），选第一个不与已标注标签重叠的：
- 台站在震中左侧 → 优先向右偏移
- 台站在震中右侧 → 优先向左偏移
- 台站在震中上方 → 优先向下偏移

特殊台站可固定偏移：
```python
fixed_offsets = {"AMPY": (-0.8, -0.5)}  # 左下角
```

`LABEL_MIN_DIST = 1.5`（度）防止标签重叠。

## SNR 动态噪声窗 —— ⚠️ 已于 2026-10-14 废弃

> 现行口径 = 「统一噪声基准 SNR」：2–8 Hz **整段带通**后取窗内 RMS，
> 除以**阵列噪声中位数**（T 1.612e-09 / Z 1.132e-09 / R 1.301e-09 m，本事件）；
> 噪声窗 `[tt_p-60, tt_p-10]`（等长 50 s），信号窗 T 分量 T-S±25 s、Z/R 用 T-P±25 s。
> 参数与两个坑见 seismic-twave-analysis `references/snr-unified-noise-baseline.md`；
> 台站图按 `snr_u` 着色、colorbar 标 `SNR (dB, common noise baseline)`。
> 下面这段（动态噪声窗 + Welch）仅作历史对照——固定/动态窗 + 单台自噪分母会把
> SNR 排序变成"台站安静程度"排序（AMPY/VOI 虚高）。

~~近台 P 波到得早，固定 300s 噪声窗可能超出数据范围；改用动态窗口 `max(20, min(300, P_time * 0.8))`~~

~~`compute_snr_welch` 参数：~~
~~- `signal_center = tt_ts`（T-S 理论到时）~~
~~- `noise_center = tt_p`（P 波理论到时，不用 tt_p-100）~~
~~- 噪声窗 = P 之前 80% 可用数据~~

以上已被 `compute_rms_pair`（整段 2–8 Hz 滤波 + 窗内 RMS）取代。

> 仍有效的教训：去响应后数据量级很大，算 PSD/RMS 前必须 `astype(np.float64)`；
> 回退值用 0.0（非 NaN），避免 Excel 表格空白。

## 剖面图分区域 + 到时标记

波形剖面按经度 43° 分两组（lon_range 参数）：
- > 43°E = 马达加斯加岛台站（profile_Madagascar.png）
- ≤ 43°E = 非洲南部大陆台站（profile_Mainland.png）

时间轴对齐到发震时刻：`t = sac.b + arange(n)*delta`，`t_min=0`。

每条波形道标记四组到时（仅 waveform 模式）：
- P=#4daf4a(绿), S=#ff9848(橙), T-P=#4d9fff(蓝), T-S=#ff4d4d(红)
- 圆点 ms=3, markeredgewidth=0.5

`tr['sac_b']` 必须存入 traces 字典，用于到时标记点时计算采样索引。

图例三列布局：`ax.legend(ncol=3, fontsize=10)`

## 地形缓存

| 分辨率 | 路径 | 大小 |
|:---|:---|:---|
| 15s (~450m) | `~/.gmt/server/earth/earth_relief/earth_relief_15s_p/` | 583 MB, 77 tiles |
| 01m (~1.8km) | `~/.gmt/server/earth/earth_relief/earth_relief_01m_g/` | fallback |
| GMT 总缓存 | `~/.gmt/server/earth/` | 8.4 GB |

pygmt 用 `@earth_relief_15s` 自动处理缓存，无需手动 grdcut。
