# Zhou et al. (2021) Qp/Qs —— 权威配方 + 参数选择分析 + 导师原始代码

**结论（2026-10-12 已复现）**：用导师原始 `.transed_<f>` 文件 + 他四个 MATLAB 函数的原口径重跑：

| 相位 | 复现结果 | 论文 |
|---|---|---|
| Qp（BHZ, V=6.0） | **203.4 f^1.48** | 204 ± 82 · f^(1.48±0.37) |
| Qs（BHT, V=3.6） | **685.4 f^0.53** | 685 ± 359 · f^(0.53±0.48) |

逐频点（2.00→4.00 Hz，步长 0.25，9 个）：
- Qp: 625, 540, 825, 982, 1223, 992, 1615, 1046, 1826
  （SI 图 S3 数字化值 626, 540, 825, 984, 1223, 990, 1615, 1048, 1826 —— **逐点吻合**）
- Qs: 1168, 1125, 1292, 830, 1104, 881, 1397, 1313, 2136

**性质说明**：这是「用作者原始数据 + 原始代码重跑，确认论文数值可复现」，**不是独立复现**
（独立复现要从原始 seed/seed 自己处理、自己算）。它的价值有三：（1）确认 204/685 就是这套口径
产生的、没有不在文档里的隐藏步骤；（2）拿到逐频点真值，可用来校验自建实现的每一步；
（3）定位出我们自己链差在哪（见 §5）。

## 1. 来源（导师工作盘 `/mnt/e/sharedate/T_waves`，Mac 拷贝，**只读勿改**）

| 角色 | 路径 |
|---|---|
| 振幅测量（信号窗/尾波窗/振幅量） | `simulated_T_waves/measure_amplitude_ecnm.m` |
| 调用方（论文三事件） | `simulated_T_waves/compute_Q_ecnm.m` |
| 回归求 Q | `simulated_T_waves/amp2Q_ecnm.m` |
| k→Q 换算 | `simulated_T_waves/k2Q.m` |
| Q(f)→Q0/α 拟合 | `simulated_T_waves/plot_alpha.m` |
| 数据（BHZ/BHT/BHN） | `obs_data/1997-04-28_mw66_Prince_Edward_Islands_Region/*.M.SAC.transed_<f>`；
  `obs_data/1998-05-21_mw63_Prince_Edward_Islands_Region/`；`obs_data/more_eqs/1999-06-20_Mw53_Prince_Edward_Islands_Region/` |

`.transed` = 已旋转（N/E/Z→R/T/Z，故有 BHT）+ 已按各中心频率窄带滤好（`.transed_2.00`…`.transed_4.00`
共 9 个频点各存一份）；SAC 头 `dist` 已被改写成**陆段距离**（SA18: 2250.6 → 629.3 km），
T 标记也为逐台标记（见 §3.8）。

## 2. 导师原始代码（逐字复制）

### 2.1 measure_amplitude_ecnm.m —— 尾波窗定义 / 振幅量
```matlab
% L20-27 默认: switch_plot='off', coda_str 默认 '2S'
switch_plot = 'off';
coda_str = '2S'; % '2S', '2TS' or '3TS'
if nargin >= 5
    switch_plot=varargin{1}; % 'off' or 'fig[ure]'
    if nargin >= 6
        coda_str=varargin{2}; % '2S' or '2TS'
    end
end

% L88-99 信号窗 / 尾波窗
[t_a,sdt] = cut_sac_data(b,ta,twin,dt,sd,'off');
switch coda_str
    case '2S'
        ts2 = sfs.O+2*(sfs.T2-sfs.O); % T2 is the arrival time of S.
    case '2TS'
        ts2 = sfs.T5+2*(sfs.T4-sfs.T5); % T4 is the arrival time of T-S and T5 is the arrival of the T.
    case '3TS'
        ts2 = sfs.T5+3*(sfs.T4-sfs.T5);
end
[t_c,sdc] = cut_sac_data(b,ts2,twin,dt,sd,'off');

% L111-118 振幅量（三种方法择一）
switch method
    case 'rms'
        res = sqrt(sum(sdt.^2))/sqrt(sum(sdc.^2));
    case 'max'
        res = max(abs(sdt))/max(abs(sdc));
    case 'mean'
        res = mean(abs(sdt))/mean(abs(sdc));
end
```
（`twin=[0,20]` 时 `cut_sac_data` 取 `[ta, ta+20]`；`ta` 由 `ref_ap` 指定 '3'→T3、'4'→T4。）

### 2.2 compute_Q_ecnm.m —— 调用（论文三事件 1997/1998/1999）
```matlab
snr = 1; % L8
f_cell={'2.00','2.25','2.50','2.75','3.00','3.25','3.50','3.75','4.00'}; % L9
...
% L21-23 Qp 与 Qs 共用的三事件文件列表
files = {['.../1997-04-28_mw66_Prince_Edward_Islands_Region/*BHZ*',fil];...
    ['.../1998-05-21_mw63_Prince_Edward_Islands_Region/*BHZ*',fil];...
    ['.../more_eqs/1999-06-20_Mw53_Prince_Edward_Islands_Region/*BHZ*',fil]};

amp = measure_amplitude_ecnm(files,'3',[0,20],'rms','nofig'); % L41 Qp: 参考 T3, 窗[0,20], rms
amp_sel=amp(amp.AMPLITUDE>=snr,:);                            % L44 SNR 筛选
amp_sel = sortrows(amp_sel,9);
[Q,Qmin,Qmax,Q_inv_std] = amp2Q_ecnm(amp_sel,str2double(fi),6.0,'fig'); % L72 V_P=6.0

% L91-93 Qs 用 *BHT*（切向）; L108 参考 T4
amp = measure_amplitude_ecnm(files,'4',[0,20],'rms','nofig');
...
[Q,Qmin,Qmax,Q_inv_std] = amp2Q_ecnm(amp_sel,str2double(fi),3.6,'fig'); % L139 V_S=3.6
```
⚠️ **两次调用都只给 5 个参数** → `nargin=5` → `coda_str` 取默认 **`'2S'`**（'2TS'/'3TS' 只是代码里
备着、论文没用）。

### 2.3 amp2Q_ecnm.m —— 几何扩散与回归
```matlab
thk_crust=35;hm=2*thk_crust;                 % L28  hm = 70 km
amp_log=log(amp.*sqrt(dist*hm));             % L29  LHS = ln( (A_dir/A_c) * sqrt(r*hm) )
X=[dist,ones(size(dist))];
[k,bint]=regress(amp_log,X);                 % L41  对 dist 最小二乘
k_std=(bint(1,2)-bint(1,1))/2/1.96;          % L42  斜率 1σ
Q=k2Q(k(1),f,c);                             % L44
Q_std=pi*f*k_std/c/k(1)^2;                   % L46  Q 的 1σ
```
（第 61/69 行画图用 `semilogy(dist, amp.*sqrt(dist*hm))`，y 轴标 `(A/C)*sqrt(r*hm)`；xlim 200–1100 km。）

### 2.4 k2Q.m
```matlab
Q = -pi*f/(c*k);
```

### 2.5 plot_alpha.m —— Q(f) → Q0、α（不加权 log-log 回归）
```matlab
f_log=log(Qp_mat(1,:)');
Q_log=log(Qp_mat(2,:)');
X=[f_log,ones(size(f_log))];
[alpha_Qp,alpha_int]=regress(Q_log,X);       % alpha(1)=α, alpha(2)=ln Q0
Qp=exp(alpha_Qp(2));
Qp_std=exp(alpha_Qp(2))*((alpha_int(2,2)-alpha_int(2,1))/2/1.96);  % 1σ
%% alpha
%alpha_Qp_v=alpha_Qp(1);
alpha_Qp_std=(alpha_int(1,2)-alpha_int(1,1))/2/1.96;
```
⇒ 论文报的 ±82 / ±359 / ±0.37 / ±0.48 = 这个回归系数的 **95% 区间 ÷ 1.96**（即 1σ），不是逐点误差。

## 3. 参数选择逐条分析（取值 / 出处 / 为什么 / 实测后果）

1. **尾波窗 `'2S'` = `O + 2*(T2-O)`，T2 = S_tele（远震 S）**
   - 出处：`compute_Q_ecnm.m` L41/L108 只传 5 参（默认）；`measure_amplitude_ecnm.m` L21、L90-92；
     SI Text S2 "a 20 s time window after two times of travel time of S wave"。
   - 为什么：CNM 的前提是 A_c **与震中距无关**（Aki 1980）——只有充分发育的尾波才满足；
     取 2×S 走时是 Aki & Chouet (1975)/Rautian & Khalturin (1978) 的经典判据。
     `'2TS'` 用 T 到达（T5）当零点，等于把 lapse 定义在**转换点之后**，尾波还没发育 → 不成立。
   - 实测：SA18 的 T2 = 503.5 s（相对发震）→ 窗在发震后 1007–1027 s，**在 T 波到达（T5=1096 s）之前 89 s**、
     T-P_local 之前 175 s。三事件 85 台合并：2S 起点中位 1044 s，比 T 波到达早 81 s（**各台均在 T 到达之前**）。
   - `'2TS'` = `T5+2*(T4−T5)` = **T-S_local 到达后再过 1 个陆段 S 走时**
     （`T4−T5` 即陆段 S 走时：SA18 = 152 s；三事件中位 189 / 140 / 138 s）；`'3TS'` 则为 2 个。
     两者都落在 T 波尾里 → A_dir/A_c < 1 被筛光 → Q 全 nan / 乱值。
   - 这正是 §6 的坑：**不要用「转换点之后 2 倍」或 2×式(1) 恒速走时替代它**。
2. **信号窗 = [T3, T3+20] 与 [T4, T4+20]**（20 s，从到时**起算**，非对称）
   - 出处：`compute_Q_ecnm.m` 双 `[0,20]`；SI Text S2 "The 20 s waveforms after the theoretical arrival times"。
   - 为什么 +20 不对称：T-P_local / T-S_local 是**续至震相**（转换区延展 + 路径弥散），能量集中在
     理论到时之后；对称窗会把 T 到达前的背景框进来。
   - 代价（实测）：包络峰与理论到时差 ±51 s，20 s 窗只抓到包络峰的 **16%–84%**（SA51 仅 16%）→
     窗内测到的是「到时后 20 s 的能量」，不是峰值。
3. **振幅量 = rms**（`sqrt(sum(sdt.^2))/sqrt(sum(sdc.^2))`）
   - 出处：调用传 `'rms'`（`measure_amplitude_ecnm.m` L111-113）。
   - 为什么不用峰值：窄带（±0.1 Hz）波形峰受单点毛刺与相位到达时刻影响大，且 20 s 窗只覆盖
     包络峰的 16–84%；rms 用整窗能量。**SI 图注写的 "peak amplitude" 属用词不严，以代码为准。**
4. **窄带 ±0.1 Hz × 9 个中心频率（2.00–4.00，步长 0.25），且逐频点滤波后另存文件**
   - 出处：SI Text S2 明写 "a narrow band of ±0.1 Hz"；文件命名 `.transed_2.00 … .transed_4.00`。
   - **导师实现（源码级）**：`.transed` 本身 = **去仪器响应**（`bin/z_trans`：
     `transfer from evalresp FNAME $resp to vel freqlim 0.01 0.02 10 20`，第 36 行）；`_<f>` 后缀才是再叠的窄带：
     同目录 2–8 Hz 版由 `bin/sacpreproc` 生成，核心为 SAC **`bp n 4 p 2 c $freq1 $freq2`**
     = **4 极 Butterworth、p 2 双向（零相位）**（第 23 行），随后 `envelope` + `smooth mean halfwidth 128` 出不 `_env`；
     MATLAB 侧同一口径：`iirfilter(sf,'bp','butter','c',[2 8],'p',2,'np',4)`（`fk_map.m` 第 51 行）。
   - 实测频谱（SA18 的 `.transed_2.00/3.00/4.00`，dt=0.05 s、n=120000=6000 s）：能量集中在中心 **±0.10 Hz**
     （±0.15 Hz 已 < −40 dB）→ 与 SI 一致。
   - 为什么：CNM 求 Q(f)，带越窄越接近单频近似。代价 = 振铃长（时间分辨率差），但测的是 20 s 窗
     能量积分，不要求时间分辨率。
   - 注：这是**偏离 CNM 家族惯例**的一步（多数文献用定比例宽带 fH/fL≈2.25，如 2.0–4.5 Hz 中心 3 Hz）；
     复现时别改成宽带。
5. **几何扩散 1/√(r·h_m)**，`thk_crust=35 → hm=70 km`
   - 出处：`amp2Q_ecnm.m` L28-29；SI Text S2 散文（Herrmann & Kijko 1983; Ma'hood et al. 2009）。
   - 为什么：r > h_m（≈2×地壳厚度，远场）→ 体波几何扩散进入柱状段；Lg 经验式给 1/√(r·h_m)。
   - h_m 取 70（源码）还是 80（SI 正文）：**只是常数因子 → 不进斜率 → 不影响 Q**。
6. **V_P=6.0 / V_S=3.6 km/s**
   - 出处：SI Text S2（cite Laske et al. 2013 CRUST1.0）；`compute_Q_ecnm.m` L72/L139。
   - 为什么重要：`Q = -pi*f/(c*k)`，V 线性进 Q —— 取错就线性错。
7. **Qp→BHZ（垂向）、Qs→BHT（切向）**
   - 出处：`compute_Q_ecnm.m` 的 glob `*BHZ*` / `*BHT*`；SI Text S2 "Qp is based on vertical component
     seismograms and Qs from tangential component seismograms"。
   - 为什么：T-P_local 以 P 型（垂向占优）为主，T-S_local 为 S 型（水平横波占优）；切向还顺带压掉
     转换点一侧的 P 泄漏。
   - 落地：导师的 `.transed` 已旋转好（BHT 现成）→ **直接吃他的文件**，避免自写旋转的 180° 反相/噪声坑。
8. **到时用 SAC 头 T 标记**（标记与 `o` 相对，读数必须 `t_i − o`；口径已由 `bin/Ttime` 源码核对）
   - 标记含义：`kt1=P`（P_tele，≈275 s）、**`kt2=S`（S_tele，≈503 s）**、`kt3=T-P`（T-P_local）、
     `kt4=T-S`（T-S_local）、`kt5=T-e`（T 波起始）。SA18 实测（相对发震）：
     275.2 / 503.5 / 1181.5 / 1248.3 / 1096.2 s。
   - T3/T4 由 `bin/Ttime` 算出：`vp=6.0  vs=3.6  vt=1.5`（海段用 **1.5** 而非 1.48），转折点取沿
     事件→台站大圆上**最后一个深于 −1000 m 的点**（陆坡坡折 = SOFAR 声道深度，与 SI Text S1 的
     "both … are converted from regions near the -1000 m isobaths" 一致）：
     `t3 = o + [contdis/vp + (dist−contdis)/vt]`，`t4 = o + [contdis/vs + (dist−contdis)/vt]`。
     ⇒ 形式与正文式(1) 同类（分段恒速），但**几何用转换点（−1000 m）而非海岸线**、海段速度 1.5。
   - T1/T2 **不来自 Ttime**（是 P_tele/S_tele 到时）。若用海岸线 + 1.48 按式(1) 复算，SA18 得
     1224/1293，比标记晚 ~25 s（几何与海段速度双重差异）。
   - 交叉源：`seismic-paper-reproduction/references/zhou-twaves-original-archive.md`
     的《t1/t5 与 T-P/T-S 标记的真实口径》节（含 Ttime 源码与 SA18 逐值校验）。
9. **SNR ≥ 1 = A_TPS_local / A_coda ≥ 1**
   - 出处：`compute_Q_ecnm.m` L8 `snr=1` + L44/L111 `amp(amp.AMPLITUDE>=snr,:)`；SI Text S2
     "SNR >1 (the amplitude ratio of T-P(S)_local to coda)"。
   - 为什么阈值这么低：CNM 要的是「尾波可用、不被自身相位污染」，不是高信噪；阈值抬高会系统性
     剔掉远台 → 距离覆盖变窄 → 斜率不可靠（实测抬高 SNR 会让 α 严重偏离论文）。
   - 无台站黑名单（源码只有这一条筛选）。
10. **三事件合并**（1997-04-28 Mw6.6 / 1998-05-21 Mw6.3 / 1999-06-20 Mw5.3）
    - 出处：`compute_Q_ecnm.m` L21-23；SI "the three earthquakes shown in Figure 1 are used"。
    - 为什么：单事件 30 多台、陆段距离跨 400–1000 km，斜率受个别台（杠杆点）支配；合并增加台-路径
      样本、稳住斜率（实测单事件的 Q(f) 在频间跳 1 个量级）。
11. **回归与误差**：`regress(ln(A·√(r·h_m)), [r, 1])` → 斜率 k → `Q = -pi*f/(c*k)`；
    Q(f) 再对 f 做**不加权 log-log 回归** `ln Q = α ln f + ln Q0`；误差棒 = 回归系数 95% 区间 ÷ 1.96。

## 4. 复现脚本

`/mnt/e/sharedate/case/zhou_2021_test/5_attenuation_merged/mentor_ecnm_repro.py`
—— 直接读导师 `.transed_<f>` 文件 + SAC 头 T 标记，复刻 §2 的两个函数，一次跑
Qp/Qs × '2S'/'2TS'/'3TS' 三套尾波口径对照。

## 5. 走过的弯路（自建链为什么差）

| # | 我们做的 | 导师口径 | 影响 |
|---|---|---|---|
| a | 取窗内**峰值** | **rms** | ❗实测影响很小：单因子只换这一项 → Qp 203→203 (α 1.48→1.46)、Qs 685→675 (0.53→0.45)；逐频点可动 ±20%，拟合值几乎不动。**不是主因** |
| b | 自写 3 阶 Butterworth 双向 ±0.1 Hz | SAC `bp n 4 p 2 c f1 f2`（4 极、双向零相位） | 同类滤波，**不是量级差异** |
| c | 用 `stations.csv` 的 TauP 到时或式(1) 恒速到时 | SAC 头 **T3/T4** | 差 ~1 s（TauP）或 ~25 s（式(1)） |
| d | `hm=80`（SI 正文） | `hm=70`（源码） | 仅常数，不影响 Q |
| e | 按「陆地计时 2×」（=代码里的 `'2TS'`）取尾波 | `'2S'` = 2×S_tele | '2TS' 下 Q 全 nan |

⇒ **要与论文逐位对齐，就直接吃导师的 transed 文件 + 用 SAC 头标记**。

🚫 **反向验证禁令（用户 2026-10-10 明确要求）**：不许扫滤波器阶数/带宽/SNR 去「凑」论文的
204 f^1.48 —— 那是拿答案反推过程，不算复现。参数只能依据 SI 文本 / 导师源码事实设定。

## 6. 注意

- **lapse time 从发震算起**（SI: "the lapse time t_c from origin"）→ 尾波窗 = 2×(S_tele)。
  若改成「转换点之后 2×」（`'2TS'`）→ 窗会落到 T-S_local 上（距 T-S 到达仅 35–55 s，中位 19–22×噪声）
  → A_dir/A_c < 1 → 台被筛光 → Q 全 nan。若改成 2×T-S_local 的**绝对**到时（OLD `ts` 口径）→
  窗落到 T 波到达前的噪声段（SNR 0.5–2.4）。
- **Qp 与 Qs 共用一个尾波窗**（同一 `coda_str`）—— 不是「T-P 一个、T-S 一个」。
- 自建链从 Zenodo 原始数据重处理时，1998 部分台（如 SA03/SA08）BHN 在 T-S 窗附近有 ~1e6 巨值
  （同台 Z 仅 ~1e4）会污染旋转出的 T 分量；导师的 `.transed` 文件无此问题。
- **滤波口径统一是「4 极 Butterworth + 双向」**（SAC `bp n 4 p 2 c f1 f2`；MATLAB `'p',2,'np',4`）。
  早期误记为 "SEIZMO 默认 passes=1" —— 错，他脚本处处显式给 np=4/p=2，默认行为根本不发生。
- **Q0 差一倍（我 320 vs 论文 685）已钉死 = 测量口径，不是数据/响应格式**（2026-10-12 单因子实验）：
  信号窗换成论文式(1) 恒速走时（-21%）+ 3sigma 迭代剔除（-39%）= x0.48 -> 685x0.48=330<=320；
  在他本人逐频点成品波形上能复现出 Qs=311.5 f^0.80（旧链报 ~320 f^0.80，alpha 完全相同）。
  响应格式 SACPZ vs RESP 只差常数 f/2（对 Z 与 T 同等影响，窄带内信号/尾波同乘约掉）-> 已排除；
  旋转用 stations.csv 的 baz 与他 .transed 头同值（signed corr +0.9998），也不是元凶。
  详见同目录 `zhou2021-response-format-and-oldchain-attribution.md`。

## 7. 论文图版面（2026-10-12 从 PDF 提取；出图照此对齐）

**SI Figure S3 (T-P_local) / S4 (T-S_local) = 振幅-距离 9 面板**（= 我们那两张 `fit_Q*_amp_vs_distance`）：
- 横轴 `Continental path distance (km)`，线性 **200-1000**，刻度 200/400/600/800/1000
- 纵轴 `A_TP * sqrt(r*hm) / A_c`（T-S 版写 `A_TS`），**对数 1e1-1e5**，主刻度 10^1..10^5
- 数据点 **空心圆**，拟合线 **黑色实线**；面板标题 = **频带**（如 `1.90-2.10 Hz`，中心 f +- 0.1 Hz）；
  面板内标 `Q=626`；**无图例**；四边全框、刻度朝内
- 论文数字化 Q 值：626 540 825 984 1223 990 1615 1048 1826（= 我们复现值，逐点吻合）

**正文 Figure 5 = Q vs 频率 (a) / 波数 (b)**（= 我们那张 `fit_Q_vs_frequency`）：
- 横轴 (a) `Frequency (Hz)`；(b) 波数 —— **两轴都要对数 (loglog)**：幂律拟合线只在双对数下才是直线，
  半对数会把线画弯。导师就是 `plot_alpha.m:5,31` 的 `loglog` -> 曾被用户指出「拟合线不是直的」。
- 左轴范围 `ylim([3e-4, 1e-2])`（plot_alpha.m:59）；右轴 Q 取倒数（100-3333）。xlim 约 1.7-4.5 Hz。
- 纵轴 **双轴**：左 `Q^-1`、右 `Q`，均对数（左刻度 1e-3/1e-2，右 1e2/1e3）
- 实心点：**红 = Qp、蓝 = Qs**；**无误差棒**；红蓝实线幂律；图例右上角写 `Qp`/`Qs`；面板标签 (a)(b)
- 关键：**数据点与拟合线画在左轴 Q^-1 上**（因此随 f/波数向下衰减），右轴 Q 只是互为倒数的换算刻度。
  画到 Q 轴上会与论文方向相反（Q 随 f 上升）-> 曾被用户一眼看出「和论文相反」。
- (b) 波数 = omega/c = **2*pi*f/V**（V=6.0 P / 3.6 S，plot_alpha.m:73,75,81-82），横轴约 1.7-7.9 km^-1；
  漏掉 2*pi 会差 6.28 倍。图例 `legend('Qp','Qs')` + boxoff，红蓝实心圆点。
- alpha 标注要带 `transform=ax.transAxes`，否则落到数据坐标外、被 bbox_inches='tight' 把画布撑大。

出图脚本 `5_attenuation_merged/plot_fit_attenuation.py`（顶部 `TT_MODE` 切 'sac'/'const'）；
排版自查用 `fk-array-analysis/scripts/check_fig_text_overlap.py` + `fig_ink_margins.py`，目标 相交 0 / 出框 0。
注：面板内 y 刻度与 x 首刻度在左下角易相碰（"200" x "10^1"）-> 给 `tick_params(pad=5.5)` 即可清掉。
