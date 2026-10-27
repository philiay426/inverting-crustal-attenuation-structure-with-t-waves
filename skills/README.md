# T 波衰减分析工作流（Hermes skills）

一套可拆可合的 T 波分析流程：**做完整衰减反演时按五步走；只做其中一段时单独用对应模块**。
这些目录是 Hermes Agent 的 skill（`SKILL.md` + `templates/` + `references/`），拷进
`~/.hermes/skills/seismic/` 即可生效，不需要改代码里的路径——所有模板都吃环境变量。

## 分层

| 层 | skill | 管什么 | 什么时候加载 |
|---|---|---|---|
| **编排** | `twave-crustal-attenuation` | 五步顺序、步间契约、跨步铁律、图件总清单、定案文档、Q 反演（步 4） | 做完整 Q 反演时 |
| **模块** | `seismic-twave-analysis` | 步 0 预处理（去响应/滤波/包络/六图）、步 1 台站表与走时/剖面 | 只要数据、只要表 |
| | `single-station-analysis` | 步 2 单台四类证据、PCA 偏振、拾取 | 只判断单台可用性 |
| | `fk-array-analysis` | 步 3 台阵 F-K（几何→引擎→慢度图→判据） | 只求来向/慢度 |
| **共用** | `plotting-standards` | 制图定版基线（栏宽/字号/线宽/图例/命名/rcParams） | 画任何图 |

## 五步

```
0 预处理      原始 SAC -> 去响应(速度) -> 2-8 Hz 带通 -> 包络 + 六张 QC 图
1 台站表      走时(TauP)、陆/海段(-1000 m 等深线)、SNR、分布图、记录剖面
2 单台分析    波形/SNR/偏振/质点运动 四类证据 -> 可用台清单 + 剔除理由
3 子阵 F-K    fk_geometry.py 现算 u_pred/smax -> Octave `fkmap` 引擎 -> 慢度图 + 峰表
4 衰减反演    尾波归一化 A_d/A_c -> 逐频点斜率 k -> Q(f)=pi f/(V|k|) -> Q0, alpha
```

步 3 的预测慢度**由算法现算**（`u_pred = 111.195 / V_T`，`smax = 1.4·u_pred`，
`f_hi = min(F_MAX, C/smax)`），换地区只改视速度与时窗，不手填网格。

## 安装

```bash
git clone git@github.com:philiay426/twave-toolkit.git
mkdir -p ~/.hermes/skills/seismic
cp -r twave-toolkit/skills/seismic/* ~/.hermes/skills/seismic/
```

依赖：Python（ObsPy / NumPy / Matplotlib / PyGMT）、Octave + SEIZMO（仅步 3 引擎）。
模板解释器统一 `python`，Octave 用 `conda run -n <env> octave-cli`。

## 模板清单

| 步 | 脚本 | 必填环境变量 |
|---|---|---|
| 0 | `preprocess_one_station.py` | `SRC`、`STA` |
| 0 | `run_batch_preprocess.py` | `SRC`、`STA_LIST` |
| 1 | `step1_map_and_table.py` | `T_BASE` |
| 1 | `step1_record_sections_3x2.py` | `T_BASE` |
| 2 | `polarization_rose_geo.py` | `SRC`、`STA` |
| 2 | `polarization_fullscan_3panel.py` | `SRC`、`STA` |
| 3 | `fk_geometry.py` | `FK_STAS`、`FK_VT`、`T_CSV`(或 `FK_SAC`) |
| 3 | `step3_fk_driver.m` | `FK_SEIZMO`、`FK_DATA`、`FK_TWIN`、`FK_VT`(或 `FK_SMAX`) |
| 3 | `step3_fk_demo.py` | `T_BASE`、`FK_STAS` |
| 4 | `twave_q_inversion.py` | `Q_CASE` |

缺必填直接 `SystemExit`，不会静默用某个案例的默认值；每个目录的 `templates/README.md`
有完整的可选变量表。

## 说明

- 本仓库只放**方法与模板**：不含具体事件的路径、台站名单、测值或结论。
- 判据写在 skill 里（如 `A_d/A_c >= 1`、`|k|/sigma_k >= 1`、主次峰差 > 1 dB、BAZ 差 < 5 deg），
  执行时按判据出结论，不按结果挑台。
