# 区域概览地形图（多台阵 + 克拉通 + 等深线）

用途：画"研究区地形图"——南部非洲 + 马达加斯加 + 西南印度洋脊这一类的多台阵总览图。
参考实现：`case/Madagascar/2013-08-17-mw62-southwest-indian-ridge/process_flow/1_station_distribution/plot_region_overview_map.py`

## 1. 克拉通/太古宙地块边界数据源（两个都验证可用）

**首选：Hasterok et al. (2022) Earth-Science Reviews「global_tectonics」cratons.shp**
- 定义：新太古代基底区（>2500 Ma），多边形主要来自 global_gprv.shp；属性 `prov_name` / `continent` / `reworked`（是否后太古宙改造）/ `area`。
- 下载（jsDelivr 在国内秒下；raw.githubusercontent 极慢会超时）：
  `https://cdn.jsdelivr.net/gh/dhasterok/global_tectonics@main/plates%26provinces/shp/cratons.<ext>`（ext = shp/shx/dbf）
- 读：`shapefile.Reader(path, encoding="latin-1")`（dbf 含拉丁字符，默认 utf-8 会炸）。pyshp 若缺：装进 obspy 环境即可。
- 图幅 [10,66,-44,-8] 内 12 个：Kaapvaal / Zimbabwe Craton、Limpopo Belt、Angolan Shield、Kasai Craton、Cuvette Centrale、Kibaran Belt、Tanzania Craton、West Congo and Kimezian Belts、马达加斯加 Antenanarivo(原文拼写，实为 Antananarivo) / Antongil / Itremo Block。

**备选：EarthByte Craton_Boundaries v1.1**（Zenodo record 15003747，DOI 10.5281/zenodo.15003746，GPL-3.0，Shirmard et al. 2025）
- 定义：地球物理（岩石圈厚度/热流反演）划出的克拉通根，与地质定义的边界不完全重合。
- Zenodo 直下快，含 Craton_Boundaries.shp / 沉积盆地 / 矿床等。

引用写法：图注里写 "Craton boundaries after Hasterok et al. (2022), Earth-Science Reviews" 并给 Zenodo DOI 10.5281/zenodo.5093930。

## 2. 文字标注自动避让（多台阵图上人工摆必压台站）

102 个台站符号 + 12 条注记时，手工摆位置必出错。做法：
- `text_box(lon,lat,text,fs)`：字宽按 `0.55 * 字号` 估，`PT_PER_DEG = 图宽cm/2.54*72 / 经度跨度`；纬向高再除 `cos(lat)`（Mercator）。
- `place_label()`：在给定搜索窗内按 0.5° 步长枚举，打分 = (框内障碍物数, 是否撞已放框, 距偏好点距离)，取最优。障碍物要把**台站与震源星标**一起算进去（只算台站会漏掉压住星标的情形）。
- 图例框也作为一个已占框喂进去（`LEGEND_BOX`），否则注记会压图例。
- 检查输出现成：本脚本会打印每条注记的落点与"框内障碍=N"，全部为 0 才算过。

## 3. 出图后不靠眼睛的核验

按像素统计代替目视：按每个台阵的填充色做容差匹配计数，可反推三角符号个数（含图例 1 个）；星标同理分主/次两种颜色；克拉通线用 `gray20` 单独计。渲染静默失败（图例语法错、颜色名错）时脚本会报错或计数为 0，能第一时间发现。

## 4. pygmt/GMT 两个顺序陷阱（各踩过一次）

- `-R` 的纬度是 **south/north** 顺序：`R=[west, east, -40.0, -11.0]`。写成北在前（`[-11, -40]`）会报 `Option -R parsing failure`，而不是静默出错，但容易误以为是数值格式问题。
- 色标 `-W` 不论横竖都是 **<长>/<粗>**：竖排写 `JMR+w9c/0.45c+v`。反过来写 `+w0.45c/9c+v` 不会报错，但会自动把 9 cm 当成横向尺寸，整张图白白宽出十几厘米（控制输出像素尺寸可当场发现）。
- 横版图（PPT 半页用）的做法：图内区域先做到 ~2:1（经纬范围 + Mercator 自动定高），色标改放右侧竖排，避免底部色标把整体高度顶起来。

## 5. GMT/pygmt 排不了中文（实测结论）

- 试过 6 种字体写法（`SimHei` / `STSong-Light` / 字体文件路径 `/mnt/c/Windows/Fonts/msyh.ttc` / 系统字体名 `Microsoft YaHei` / `Noto Sans SC` / 默认），输出与默认字体**逐像素完全一致** → GMT 直接忽略这些指定，中文变成空白或方框。判据：同一串中文用两种字体渲染后比较像素差异，0 差异＝字体没生效（这套差分法同样适用于验证别的“字体到底用上没”）。
- 做法改成两段式：pygmt 只画底图（地形/等深线/克拉通线/台站/震源/色标条与刻度数字）存 `base.png`；matplotlib 在 `base.png` 上叠中文标注、图例、色标标题，再存 PNG/PDF（PDF 里底图是栅格、文字是矢量）。
- 经纬度 → 像素换算不需要额外的坐标输出：底图里**地图图框是最长的横/竖黑线**（占图宽/图高 75% 以上即命中；色标条只有 9 cm 高、图例框更短，不会误判），检出的图框像素范围 + 球面 Mercator 公式就能换算。自检办法：检出图框的宽高比应等于图幅的 Mercator 宽高比（实测 1.954 vs 理论 1.963）。像素 y 向下，纬度用 Mercator 内插。
- 中文描边（白边样式）：`path_effects=[pe.withStroke(linewidth=2*HALO_PT, foreground="white")]`（描边总宽取 2 倍对应向外 HALO_PT）。中英混排的框宽估算：中日韩字符按 1.0 em、拉丁/数字按 0.55 em（布局避让搜索要用）。
- 可用中文字体文件：`/mnt/c/Windows/Fonts/msyhbd.ttc`（微软雅黑 Bold，与台站表/交接文档字体一致）。matplotlib 用 `FontProperties(fname=...)` 直接吃 .ttc。
- 竖排色标的刻度数字占满右侧，标题别硬塞右边——放色标**正下方**横排。
