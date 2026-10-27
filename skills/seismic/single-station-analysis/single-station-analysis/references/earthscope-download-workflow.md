# EarthScope 台网数据下载流程（2026-08-04 实测）

给单台分析/新事件准备数据时，从 EarthScope FDSN 下载台网（如 XK/XV）三分量数据的完整流程与踩坑。

## 1. 确认台网运行期（先查元数据）

```python
from obspy.clients.fdsn import Client
c = Client('EARTHSCOPE')   # IRIS 已改名 EarthScope，用 'EARTHSCOPE' 短名（'IRIS' 也可但会告警）
inv = c.get_stations(network='XK', level='network')
for net in inv:
    print(net.code, net.start_date, net.end_date, net.description)
```

实测（2026-08）：
- **XK = SAFARI**（Seismic Arrays for African Rift Initiation），2012-01 ~ 2014-12
- **XV = MACOMO**（Madagascar/Comores/Mozambique），2011-01 ~ 2013-12

⚠️ **台网元数据声称运行 ≠ 实际有数据**：XK 声称运行到 2014-12，但 2014-11-17 事件查询返回 204 no data（实际已停摆）。先试下载，204 就跳过，不要反复重试。

## 2. 下载波形

逐台站×分量下载慢（90 站需 10+ 分钟），按网络批量查询更快：
```python
st = c.get_waveforms(network=net, station=','.join(stas), location='*',
                     channel='BH?', starttime=t0-PRE, endtime=t0+POST)
```
时间窗：发震前 600s ~ 发震后 7200s（T 波需要长记录）。

## 3. 🔴 SAC 头缺 stla/stlo（必踩坑）

FDSN 下载转 SAC 后**通常没有台站坐标头**（stla/stlo 缺失），step1/模板按坐标过滤会返回 0 台站。
必须先从台网元数据取坐标补写：
```python
inv = c.get_stations(network=net, starttime=t0-3600, endtime=t0+7200, channel='BH?')
coords = {f'{n.code}.{s.code}': (s.latitude, s.longitude, s.elevation)
          for n in inv for s in n.stations}
for f in glob.glob(f'{raw}/*.SAC'):
    net, sta = os.path.basename(f).split('.')[0], os.path.basename(f).split('.')[1]
    if f'{net}.{sta}' in coords:
        tr = read(f)[0]
        tr.stats.sac.stla, tr.stats.sac.stlo, tr.stats.sac.stel = coords[f'{net}.{sta}']
        tr.write(f, format='SAC')
```

## 4. 🔴 get_sacpz() 不可用 → 从 StationXML 手动提取 PAZ

`client.get_sacpz()`（EarthScope）实测返回 0 文件。改用 `get_stations(level='response')` 拿 StationXML，
手动提取 PolesZeros stage 写 SAC PZ 文本（process_image.py 的 parse_sacpz 读取格式）：

```python
inv = c.get_stations(network=net, starttime=t0-3600, endtime=t0+7200,
                     channel='BH?', level='response')
for st in inv[0].stations:
    for ch in st.channels:
        resp = ch.response
        pz = None
        for sg in resp.response_stages:
            if hasattr(sg, 'poles') and hasattr(sg, 'zeros'):
                pz = sg; break
        if pz is None: continue
        norm = pz.normalization_factor
        # ⚠️ obspy 新版本 PoleZero stage 无 .sensitivity 属性！
        sens = resp.instrument_sensitivity.value   # 从 instrument_sensitivity 取
        # 写 SACPZ.{net}.{sta}..{ch} 文本: 头注释 + ZEROS n + 复数行 + POLES n + CONSTANT norm*sens
```

⚠️ **obspy 版本差异**：`response.get_paz()` 在某些版本报 `'PolesZerosResponseStage' object has no attribute 'sensitivity'`。
手动遍历 `response_stages` 找 `hasattr(sg,'poles')` 的 stage 最稳，sensitivity 统一从 `resp.instrument_sensitivity.value` 取。

## 5. 生成 stations.csv 的捷径

SAC 头补齐坐标后，可直接复用同区域旧事件的简化脚本
`<旧事件>/process_flow/step1_compute_arrivals.py`（sed 改 EVLA/EVLO/EVDP + date + mw + BASE 4 处）。
⚠️ sed 替换顺序坑：`6.2→6.3` 会把 `mw62` 变成 `mw63`，导致后面匹配 `...mw62-southwest...` 的替换失效——BASE 路径用独立 replace 或最后统一检查。

## 6. 台站来源判断（事件归属）

判断事件属于哪个区域/文件夹（如 Madagascar vs Southern_Africa），按 raw SAC 文件的台网前缀：
```bash
ls <event>/raw/*.SAC | sed 's|.*/||' | cut -d. -f1 | sort -u
```
实测：含马达加斯加 MACOMO 台站的事件用 **XK/XV** 台网；纯非洲大陆事件用 ZT/BX/AF/GT。
移动事件目录后，events.yaml 里该事件的 `case_dir` 必须同步改为新路径，否则脚本找不到数据。
