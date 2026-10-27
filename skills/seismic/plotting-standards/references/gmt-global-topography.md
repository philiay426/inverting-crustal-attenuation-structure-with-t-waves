# GMT 全球地形缓存

## 路径

| 分辨率 | 路径 | 大小 | 说明 |
|:---|:---|:---|:---|
| 15s (~450m) | `/home/sy/.gmt/server/earth/earth_relief/earth_relief_15s_p/` | 583 MB, 77 tiles | **主用** |
| 15s 海洋 | `/home/sy/.gmt/server/earth/earth_synbath/earth_synbath_15s_p/` | 12 MB, 2 tiles | 海洋高精 |
| 01m (~1.8km) | `/home/sy/.gmt/server/earth/earth_relief/earth_relief_01m_g/` | — | 回退 |
| GMT 总缓存 | `/home/sy/.gmt/server/earth/` | 8.4 GB | 全部 |

## pygmt 用法

```python
# 直接引用，GMT 自动从缓存加载，缓存未命中时自动下载
fig.grdimage(grid="@earth_relief_15s", ...)
```

## 不再手动 grdcut

GMT 6.x + pygmt 支持 `@remote_dataset` 语法自动处理缓存，
不再需要手动 `gmt grdcut` 或保存临时文件。
