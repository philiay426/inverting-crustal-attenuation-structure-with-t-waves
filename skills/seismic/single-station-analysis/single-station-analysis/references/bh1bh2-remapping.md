# BH1/BH2 → BHN/BHE 映射

## 背景

US 网络等台站使用 BH1/BH2 命名水平分量（而非 BHN/BHE）。分析脚本的 `load_station()` 默认搜索 `BHN/BHE`，需要添加分量名映射。

## 分量映射

```python
for ch in ['BHZ','BHN','BHE','BH1','BH2']:
    for pat in [f"{resp_dir}/{net}.{sta}.*.{ch}*_resp.sac",
                f"{resp_dir}/{net}.{sta}..{ch}*_resp.sac"]:
        m = glob.glob(pat)
        if m: trs[ch] = read(m[0])[0]; break

# BH1/BH2 → BHN/BHE 映射（读取后、trim前）
if 'BHN' not in trs and 'BH1' in trs: trs['BHN'] = trs['BH1']
if 'BHE' not in trs and 'BH2' in trs: trs['BHE'] = trs['BH2']
```

## SACPZ 下载修复

BH1/BH2 的 SACPZ 文件可能在首次下载时失败（`400 BAD_REQUEST: level=resp + format=text` 不兼容）。

### 正确下载方式

```python
import requests
url = f'https://service.earthscope.org/fdsnws/station/1/query?' \
      f'net=US&sta={STA}&loc=00&cha={CH}&level=response&format=xml' \
      f'&starttime=2014-01-13&endtime=2014-01-14'
r = requests.get(url)  # XML格式
# 从XML中提取PolesZeros和InstrumentSensitivity
```

### 快捷方案（STS-2传感器）

US 网络的 LBNH 等台站使用 STS-2 传感器，**BH1/BH2 与 BHZ 具有完全相同的仪器响应**（相同 NormalizationFactor、Poles/Zeros、Sensitivity）。可直接复制：

```bash
cp raw/SACPZ/SACPZ.US.LBNH.00.BHZ raw/SACPZ/SACPZ.US.LBNH.00.BH1
cp raw/SACPZ/SACPZ.US.LBNH.00.BHZ raw/SACPZ/SACPZ.US.LBNH.00.BH2
```

然后重新运行预处理（去仪器响应）。
