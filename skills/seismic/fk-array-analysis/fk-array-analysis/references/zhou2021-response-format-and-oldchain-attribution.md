# Zhou2021: 响应格式差异 与 旧链 Qs 偏差归因 (2026-10-12 实测)

补充 `zhou2021-ecnm-exact-repro.md`。数字全部来自脚本可复算的结果, 不是推理。

## 1. 响应格式: SACPZ vs RESP
- 导师: RESP + evalresp (`bin/z_trans` 第 36 行 `transfer from evalresp ... to vel freqlim 0.01 0.02 10 20`) -> 速度, 2 个零点。
- 旧链: FDSN `raw/SACPZ/SACPZ.<net>.<sta>..<ch>` + 手写 `parse_sacpz` + obspy `simulate(paz_remove=..., remove_sensitivity=True)` (`process_flow/process_image.py:102,135`) -> 3 个零点, `INPUT UNIT: M` (位移型)。
- 实测 (1997 事件 39 台 BHZ, 2.0/2.5/3.0/3.5/4.0 Hz): 两者振幅响应之比 = 1.2511 / 1.4997 / 1.7537 / 2.0017 = **f/2**, 39 台完全一致 (一个干净的 2*pi*f 因子)。
- 含义: 格式差是真的, 但在窄带 (+-0.1 Hz) 内对信号窗与尾波窗是同一个常数 -> 在 A_sig/A_coda 里约掉 -> **不解释旧链 Qs 差一倍**。它对 Z 与 T 同等影响, 而旧链当年 Qp 是对的, 这本身就排除了它。
- 灵敏度: 个别台 SACPZ 的 SENSITIVITY 错 (SA62 = 3.07e6, 其他约 4.2e8, 差 136 倍) -> 常数因子, 比值里同样约掉, 但会让绝对振幅/SNR 闸门/出图错。

## 2. NE->RT 旋转口径 (实测, 用于纠正 "csv 的 baz 不可信" 的旧说法)
- 用 `stations.csv` 的 `baz` 列 + obspy `rotate_ne_rt(n, e, baz)`: 与导师 `.transed` 的 BHT 比 signed corr **+0.9998** (35/37 台 |corr|>0.99)。
- 用真实大圆 BAZ (`gps2dist_azimuth(ev, sta)[1]`): signed corr **-0.9880** (180 度反相)。
- 手写 `T = -n*sin(b) + e*cos(b)` 配 csv baz: **-0.9998** (与 obspy 差 180 度)。
- 事实: 导师 `.transed` 头里的 `baz` (SA18 = 137.60) **就是** stations.csv 里那个值 -> 他本人用的同一口径 (他头里事件坐标被改写成 "陆上等效源", 所以 SAC `rot to gcp` 算出的是这个角度; 与真实事件大圆 BAZ 差约 171 度)。
- 结论: **NE->RT 用 stations.csv 的 baz 列, 不要换成真实大圆 BAZ**。

## 3. 旧链 Qs 约 320 f^0.80 的归因 (单因子, 波形源固定为导师逐频点成品)
对照 V0 (导师口径) = Qp 203.4 f^1.48 / Qs 682.9 f^0.53 (基准必须回到这里, 差值才可用)。每行只改一个旋钮:

| 旋钮 | Qp Q0/alpha | Qs Q0/alpha |
|---|---|---|
| V1 信号窗改用论文式(1)恒速走时 | 194.3 / 1.42 | 541.6 / 0.51 |
| V3 振幅 peak 替 rms | 202.5 / 1.46 | 668.3 / 0.46 |
| V6 hm 70 -> 80 | 203.4 / 1.48 | 682.9 / 0.53 (只动截距) |
| V7 回归距离改用 stations.csv land_km | 201.6 / 1.48 | 673.3 / 0.53 |
| V8 加 3-sigma 迭代剔除 | 207.3 / 1.44 | 418.9 / 0.73 |
| V2 尾波改 lands | nan | nan (窗落点不过 SNR 闸门) |
| **V10 = V1 + V8** | **194.3 / 1.42** | **311.5 / 0.80** |
| V11 = V10 + peak | 216.1 / 1.32 | 312.7 / 0.79 |

- 旧链当年报 Qp ~211 f^1.27 / Qs ~320 f^0.80 -> **V10 与旧链 Qs 的 Q0 差 3%, alpha 0.80 完全相同**; 再加 peak (V11) 后 Qp 也靠近 (216.1 / 1.32 vs 211 / 1.27)。
- 结论: 旧链 Qs 的 x2 缺口 = **信号窗口径 (论文式(1) 恒速走时替掉 SAC 头 T3/T4) + 3-sigma 迭代剔除离群**, 不是响应格式, 不是旋转, 不是数据。
- 已排除: `H_M` 70 vs 80 只动截距 (Q 不变); `rms` vs `peak` 只让 Qp 的 alpha 1.48 -> 1.46 (Qs 几乎不变)。

## 4. 验证脚本与坑
- `case/zhou_2021_test/5_attenuation_merged/verify_reprocess_from_raw.py` (+ 同名 .log): 吃导师 raw 按配方重处理, 与直接吃他成品做 stage-swap 对照: Qp 逐频点比值 1.000x8 + 1.001, Qs 0.978-1.029。
- 坑 (脚本级, 防重犯):
  - obspy `Trace.slice(t0, t1)` 会少取尾样点 (84001 -> 84000); 对齐要按样点号算 offset (`off = round((his.starttime - raw.starttime)*sr)`), 别拿 slice 后的长度当严格闸门。
  - 缓存键必须带窗口 (起点 + npts): 1999 的 BHZ 逐频点文件是 cut 过的 (54001 点, b=300), BHT 是全长 (84001 点, b=0), 同一台站两个分量窗口不同。
  - SAC 头时间一律相对 reference: `index = (t - b)*sr`; 忽略 b 会让 1999 (b=300) 整批错位。
  - obspy 1.5 的 `read_inventory(..., format="SACPZ")` 不支持; 要比 SACPZ 就按 SAC 的 A0/CONSTANT/SENSITIVITY 约定手算, 或用旧链自己的 `parse_sacpz`。
