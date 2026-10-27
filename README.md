# T波地震数据处理工具包

## 3 步事件处理流程（批量统计）

```
加新事件: 编辑 config/events.yaml → 按顺序跑 3 个脚本
```

| 步骤 | 脚本 | 功能 | 输出 |
|:----|:-----|:-----|:-----|
| 1 | `step1_station_table_and_map.py` | 台站信息+到时+SNR+台站分布图 | stations.xlsx, stations.csv, figure_map_tp.png |
| 2 | `profile_wide.py` | 波形剖面图(全局P98+软裁剪+海洋对齐) | profile_az*.png |
| 3 | `find_conversion_point.py` | 转换点计算+走时残差验证 | conversion_points.csv, residuals.csv |

### 用法

```bash
# 前提: 激活虚拟环境
source /home/sy/codex/venv/bin/activate

# 步骤1: 台站分布与到时表
python3 scripts/step1_station_table_and_map.py puerto_rico
python3 scripts/step1_station_table_and_map.py mindanao_2026 --snr-min 5.0

# 步骤2: 波形剖面图
python3 scripts/profile_wide.py puerto_rico --snr-min 3 --y-land --align-tp
python3 scripts/profile_wide.py puerto_rico --snr-min 3 --y-land --align-tp --exclude X51A,X60A

# 步骤3: 转换点计算与走时残差
python3 scripts/find_conversion_point.py puerto_rico
```

### 添加新事件

编辑 `config/events.yaml`，copy 一段：

```yaml
new_event:
  label: "new_event"
  date: "2015-01-01"
  lat: 18.0
  lon: -65.0
  depth_km: 10.0
  mw: 6.5
  data_dir: "/mnt/e/sharedate/new_event_raw"
  sac_date: "2015.001.000000"
  case_dir: "2015-01-01-mw65-new_event"
  snr_threshold: 3.0
  v_ocean: 1.48    # 可选, 默认1.48
  v_love: 3.6      # 可选, 默认3.6
```

然后把原始 SAC + SACPZ 放到 `data_dir`，跑预处理(`lib/processing.py`)后开始 3 步流程。

## 目录结构

```
codex/
├── README.md
├── config/events.yaml          ← 事件参数
├── lib/                        ← 核心函数库
│   ├── processing.py           → 去响应、滤波、SACPZ
│   ├── arrivals.py             → 到时计算+SNR(Welch PSD)
│   └── rotation.py             → ZNE→ZRT、能量比
├── scripts/                    ← 3 步主脚本 (通用, 参数化)
│   ├── step1_station_table_and_map.py
│   ├── profile_wide.py
│   └── find_conversion_point.py
├── scripts/old_scripts/        ← 旧版脚本备份
├── t_phase_analysis/           → T波走时+PMCC+时频谱
└── case/                       → 各事件工作目录
```

## 关键参数

| 参数 | 值 | 说明 |
|:----|:----|:-----|
| 频带 | 2-8 Hz | T波特征频带 |
| water_level | 60 dB | 去仪器响应 |
| pre_filt | 0.001-50 Hz | 去仪器响应 |
| celerity_max | 1.8 km/s | T波与众波区分 |
| V_ocean | 1.48 km/s | SOFAR声速 |
| V_land | 3.6 km/s | Love波速 |
| SNR方法 | Welch PSD 2-8Hz | 频域功率谱, 非包络峰值 |

## Hermes skills（T 波衰减分析工作流）

`skills/seismic/` 下是可直接装入 Hermes Agent 的 skill（`~/.hermes/skills/seismic/`）：

| skill | 作用 |
|---|---|
| `twave-crustal-attenuation` | **编排**：五步顺序、步间契约、Q(f) 反演（步 4） |
| `seismic-twave-analysis` | 步 0 预处理 + 步 1 台站表/走时/剖面 |
| `single-station-analysis` | 步 2 单台偏振与四类证据 |
| `fk-array-analysis` | 步 3 子阵 F-K（几何→引擎→慢度图） |
| `plotting-standards` | 制图定版基线（栏宽/字号/线宽/图例/命名） |

全流程按五步走；只做其中一段时单独加载对应模块。详见 `skills/README.md`。

## 运行环境（WSL / Linux 优于 Windows 盘）

流程是"几千个小文件"型 I/O，实测同一批 SAC：`/mnt/e`（9p）30.2 ms/个、ext4 1.7 ms/个（约 20 倍）。
所以：

| 内容 | 放哪 |
|---|---|
| 代码 / 仓库 / skill | Linux 本地盘（`~/codex`、`~/.hermes/skills/`） |
| 原始数据、成品图 | 共享盘 `/mnt/e/sharedate`（容量大，顺序读够用） |
| 桌面上的 skill 包 | 只作分享镜像，**不要在里面跑脚本** |

两个辅助脚本：

```bash
scripts/stage_event.sh <事件目录名> [子目录...]        # 数据缓存到 ~/data_cache/<事件>/（默认 raw）
scripts/stage_event.sh --back <事件目录名> <子目录>    # 本地产出回写共享盘
scripts/sync_skills.sh "提交信息"                      # ~/.hermes/skills/seismic -> 仓库 -> 双推
```
