# StationXML 响应文件适配指南（2026-08-01）

IRIS FDSN 2025年后下载的数据用 StationXML (`.xml`) 代替传统 SACPZ 文本。

## 关键差异

| | SACPZ 文本 | StationXML |
|---|---|---|
| 去响应 | `tr.simulate(paz_remove=paz)` | `tr.remove_response(inventory=inv, output='VEL')` |
| 台站坐标 | SAC 头 stla/stlo | `inv[0][0].latitude/longitude` |
| 通道名 | BHZ/BHN/BHE | 可能有 BH1/BH2/BHX/BHY |

## 适配模板

从 `process_image.py` 修改以下部分：

### 1. 读取响应
```python
from obspy import read_inventory

def load_station_meta():
    meta = {}
    for xml_f in glob.glob(os.path.join(XML_DIR, '*.xml')):
        inv = read_inventory(xml_f)
        net, sta = inv[0].code, inv[0][0].code
        meta[f"{net}.{sta}"] = {
            'lat': inv[0][0].latitude,
            'lon': inv[0][0].longitude,
            'inventory': inv
        }
    return meta
```

### 2. 去响应
```python
def remove_response_xml(tr, inv):
    tr_out = tr.copy()
    tr_out.detrend('demean'); tr_out.detrend('linear')
    tr_out.taper(max_percentage=0.05)
    tr_out.remove_response(inventory=inv, output='VEL',
                            pre_filt=(0.001,0.005,45,50), water_level=60)
    return tr_out
```

### 3. 分量映射
```python
comp_map = {'BHZ':'BHZ', 'BHN':'BHN', 'BHE':'BHE',
            'BH1':'BHN', 'BH2':'BHE', 'BHX':'BHN', 'BHY':'BHE'}
```

### 4. 仅BHZ台站跳过出图
```python
has_3c = all(c in comps for c in ['BHZ','BHN','BHE'])
if not has_3c:
    print(f"BHZ only — skip 4-panel plot")
    continue
```

## Japan 2025-08-01 案例

事件: lat=41.30°N, lon=142.00°E, depth=80km (Honshu East, 陆地地震)
台站: 19 太平洋岛台站 (IU + JP + MI + PS), StationXML 响应
结果: 11 台站全四图, 3 台站仅 BHZ, 4 台站缺 XML 元数据

## pyGMT 台站分布图

台站分布图必须用 pyGMT — 出版级要求：地形底图、方位角弧线、T-P 着色、色标。

```python
import pygmt
from geographiclib.geodesic import Geodesic

# 跨日期线：经度跨度 >300° 时转为 0-360
if max(lons)-min(lons) > 300:
    lons = [(lo+360) if lo<0 else lo for lo in lons]
    evlo360 = EVLO if EVLO>0 else EVLO+360

fig = pygmt.Figure()
fig.grdimage(grid='@earth_relief_01d', region=reg, projection=f'M{w}c',
             shading='+a45+nt0.3', cmap='geo')
fig.coast(shorelines='0.5p,gray30', borders='1/0.5p,gray40')

# 方位角弧线（每 15° 虚线）
for azi in range(0, 360, 15):
    # ... Geodesic 大圆投影 ...
    fig.plot(x=cl, y=ca, pen='0.3p,gray30,--')

# 台站 — T-P 走时 hot 着色
pygmt.makecpt(cmap='hot', series=[t_min, t_max, 2], reverse=True)
fig.plot(x=lons, y=lats, fill=t_tp_min, style='t0.25c', cmap=True, pen='0.15p,black')

# SNR top 5 标注
fig.text(x=lo+0.3, y=la+0.15, text=key, font='5p,Helvetica,black', fill='white@70')

# 震中
fig.plot(x=[evlo360], y=[EVLA], style='a0.4c', fill='red', pen='0.5p,black')
fig.colorbar(frame='xa10f5+l"T-P (min)"', position='jBL+w8c/0.3c+h+o1c/0.8c+ml')
```

**反面案例**：`plt.scatter(lons, lats)` — 无地形、无弧线、无着色 → 不合格。❌
