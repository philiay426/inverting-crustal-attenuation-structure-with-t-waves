# 文献地图标注技巧（基于 pyGMT map.py 模板）

## 传播路径标签分散放置

当多条传播路径汇聚到同一个小区域时，沿路径的标签容易挤在一起。

**策略**：每条路径用不同的 `frac`（路径长度比例），分散到画布不同象限：

```python
# idx 0 = 最靠海岸的路径（frac 高）
# idx 1 = 最远的路径（frac 中低）  
# idx 2 = 中间路径（frac 中高）
positions = [0.90, 0.55, 0.70]  # 按 idx 排列
frac = positions[idx]
mi = int(npts * frac)
# 向后搜索确保在画布区域内
for i in range(mi, npts):
    if R[0] <= plon[i] <= R[1] and R[2] <= plat[i] <= R[3]:
        mx, my = plon[i], plat[i]; break
```

## 路径终点选择

箭头应终止于海岸附近（论文中 T 波转换区），而非岛内部台站。

## 海底地形

使用 `@earth_synbath_15s` 而非 `@earth_relief_*`，后者仅含陆地高程。

## 图例格式

pyGMT legend spec 中虚线格式：`"S 0.3c - 0.8c - 1.2p,#d62728,- 0.8c Label"`
- 第 3 字段 `0.8c` = 短线长度
- 第 4 字段 `-` = 占位（无符号）
- 第 5 字段 = 线样式（`1.2p,#d62728,-`）
