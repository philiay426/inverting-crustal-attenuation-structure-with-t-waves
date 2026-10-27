# 波形剖面图参考（2026-08-04 最新）

## 实现

| 实现 | 脚本 | 适用 |
|------|------|------|
| **独立脚本（推荐）** | `templates/waveform_proflie.py` | 新事件标准，land_km 分箱取 SNR 最高 |

> ⚠️ `profile_conversion.py` 已废弃。step1 内置 `_gen_waveform_profile()` 仅作备份。

## 参数速查

```bash
python waveform_proflie.py <event> --snr-min 3.0 --bin-km 10 \
    --regions "10:43:Mainland,43:60:Madagascar" [--band 0.5-2.0]
```

| 参数 | 默认 | 说明 |
|------|------|------|
| `--bin-km` | 10 | 陆地距离分箱宽度；档内只取 SNR 最高台站（不叠加） |
| `--regions` | "10:43:Mainland,43:60:Madagascar" | lon:lon:label 逗号分隔 |
| `--band` | **0.5-2.0** | **低频：T 波低频衰减慢，远处台站(480km+)仍可见；2-8Hz 几百公里就耗尽。近处可切 2.0-8.0** |
| 纵轴 | land_km（陆地距离） | dist 版已弃用 |
| 时间窗 | `min(tt_tp)-30 ~ max(tt_ts)+60` | 收紧窗口，信号充满图面 |
| 滤波 | `butter(4)` band + filtfilt | 必须 cast float64 防溢出 |
| 包络 | Hilbert → **3s 平滑** | 8s 会把峰值抹平、纺锤形变钝 |
| 归一化 | **包络独立 global_env_p98**（先扣噪声基线） | **不复用波形 global_p98**（尺度不同） |
| 主图 | 灰色填充包络(alpha=0.45) + 黑色细波形叠加(lw=0.25) | 输出 `_env` 版 |
| 辅助图 | per-bin P98（`_env_perbin` 版） | 远处弱信号也能看，但噪声会被放大 |
| 波形版 | 纯波形（`_waveform` 版） | |
| 排除台站 | `PROFILE_EXCLUDE` 集合 | 异常噪声站（SNR 低但振幅异常大） |

## 🔴 时间轴对齐（2026-08-04 关键修复，08-17 修正规则）

数据下载常从发震前 ~10 分钟开始（`o=600` → 数据从发震前 600s 起）。统一：
```python
sac_o = getattr(tr.stats.sac, 'o', None)
if sac_o is not None and 0 < abs(sac_o) < 1e6:
    t_off = tr.stats.sac.b - sac_o        # o=600 → -600，数据从发震前开始（2013-08-17 实测正确）
else:
    EV_ORIGIN = UTCDateTime(年, 月, 日, 时, 分, 秒)   # 从文件名 D.2019.022.190143 解析
    t_off = tr.stats.starttime - EV_ORIGIN            # 缺 o 时 fallback
t = t_off + np.arange(len(data)) * delta
```
⚠️ 勿仅凭 `starttime == EV_ORIGIN` 判断从发震开始（ObsPy reference time 可能被写成发震时刻）。
自检：`t[0]` 应为负（发震前开始，常见 -600）或 0；窗口内应有纺锤形能量包络，若只见平坦噪声 → 先查 t_off。
