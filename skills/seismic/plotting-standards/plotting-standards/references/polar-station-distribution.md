# 极坐标台站分布图（方位-距离）

不画地理底图时，用极坐标展示台站在方位-距离空间的分布。适合快速筛选台站覆盖范围。

## 标准参数

| 参数 | 值 |
|------|-----|
| 投影 | `subplot_kw={'projection': 'polar'}` |
| 尺寸 | 双栏 19cm，`figsize=(19/2.54, 19/2.54*0.95)` |
| 距离圈 | 每 10° 灰细线 `lw=0.3, alpha=0.4`，外端标注 `{r}°` |
| 震中 | 红色星号 `*`, `s=120, c='#d62728'`, 黑边 `lw=0.6` |
| 台站 | 灰色小圆 `#aaaaaa`(s=5)，精选台站红色大圆 `#d62728`(s=45) + 黑边 `lw=0.4` |
| 方位范围 | `set_thetamin(az_min) / set_thetamax(az_max)` |
| 方向 | `set_theta_zero_location('N')`, `set_theta_direction(-1)` |
| 标注 | 精选台站名：`ax.annotate()` + arrow, fontsize=7 bold |
| 图例 | `loc='lower left'`, framealpha=0.8 |
| 标题 | `EVENT — Station Distribution  (az=…°, dist≤…°)` |
| rcParams | 按 §一 通用规范 |

## 代码模板

```python
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9,
    'axes.labelsize': 10, 'axes.titlesize': 10,
    'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'axes.linewidth': 0.6,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})

fig_w = 19.0 / 2.54
fig, ax = plt.subplots(figsize=(fig_w, fig_w * 0.95),
                        subplot_kw={'projection': 'polar'}, dpi=150)

selected = {'STA1', 'STA2', ...}

for s in all_stations:
    is_sel = s['code'] in selected
    ax.scatter(np.radians(s['az']), s['dist'],
               c='#d62728' if is_sel else '#aaaaaa',
               s=45 if is_sel else 5,
               alpha=0.9 if is_sel else 0.4,
               edgecolors='black' if is_sel else 'none',
               linewidth=0.4, zorder=10 if is_sel else 1)

ax.scatter(0, 0, c='#d62728', s=120, marker='*',
           edgecolors='black', linewidth=0.6, zorder=20)

for r in [10, 20, 30, 40, 50]:
    az_arr = np.linspace(np.radians(52), np.radians(180), 100)
    ax.plot(az_arr, [r]*100, 'gray', lw=0.3, alpha=0.4, zorder=0)
    ax.text(np.radians(180), r + 0.5, f'{r}°', fontsize=7, ha='left', color='gray')

ax.set_ylim(0, 55)
ax.set_thetamin(52); ax.set_thetamax(180)
ax.set_theta_zero_location('N'); ax.set_theta_direction(-1)
ax.tick_params(labelsize=8, length=3)

for s in all_stations:
    if s['code'] in selected:
        name = s['code'].split('.')[1]
        ax.annotate(name, xy=(np.radians(s['az']), s['dist']),
                    xytext=(np.radians(s['az']), s['dist'] + 2.5),
                    fontsize=7, fontweight='bold', ha='center',
                    arrowprops=dict(arrowstyle='->', color='black', lw=0.5))

ax.set_title(f'EVENT — Station Distribution\n'
             f'{len(selected)} selected / {len(all_stations)} total',
             fontsize=10, fontweight='bold', pad=22)

from matplotlib.lines import Line2D
ax.legend(handles=[
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#d62728', markersize=7,
           markeredgecolor='black', markeredgewidth=0.4, label=f'Selected'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#aaaaaa', markersize=5,
           label=f'Available'),
    Line2D([0],[0], marker='*', color='w', markerfacecolor='#d62728', markersize=10,
           markeredgecolor='black', markeredgewidth=0.4, label='Epicenter'),
], loc='lower left', fontsize=8, framealpha=0.8)

fig.savefig('station_map.png', dpi=300, bbox_inches='tight')
```
