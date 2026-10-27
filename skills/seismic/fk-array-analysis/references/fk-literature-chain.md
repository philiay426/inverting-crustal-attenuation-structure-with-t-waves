# 论文链条逐步推算（Zhou 2021 ← Matsumoto 2016 ← Rost & Thomas 2002）

用途：用户要求"根据论文给的信息 + 它引用的文章的公式原理一步一步推算"时，照本文件给**每一步的
原文出处 + 公式 + 代入本台阵的数值**。不要凭记忆复述文献，也不要只给结论。

## 0. 三篇原文在本机的位置（先找到文件，再引用）

```
论文/T波/用非洲南部T波测定地壳衰减/Zhou et al. - 2021 - Determining Crustal Attenuation With Seismic T.pdf  (+ 附件.pdf = SI)
论文/T波/Geophysical Research Letters - 2016 - Matsumoto - T-phase and tsunami pressure waveforms recorded by near-source IMS.pdf
论文/T波/台阵地震学——方法与应用（Rost & Thomas 2002）/Reviews of Geophysics - 2002 - Rost - ARRAY SEISMOLOGY  METHODS AND APPLICATIONS.pdf  (+ 中文翻译 pdf/docx)
```

⚠️ **在用户的论文库里找一篇文献：按内容 grep，不要只按文件名 find** —— 文件名以**期刊名**开头
（`Geophysical Research Letters - 2016 - Matsumoto - …`），拿作者名去 `-iname '*matsumoto*'` 搜
很容易判成"没下载"而误报"文献缺失"。可靠做法是逐份转文本再搜关键词：

```bash
find "…/Desktop/论文" -iname '*.pdf' -print0 | while IFS= read -r -d '' f; do
  pdftotext -layout "$f" - 2>/dev/null | grep -qi 'Matsumoto' && echo "$f"
done
```

取正文：`cd /tmp && pdftotext -layout "<pdf>" out.txt`，再 `grep -n -i` 定位 Methods 段。

## 1. 步① 无混叠判据（Matsumoto 2016 §3 原文）

> "The criterion for avoiding spatial aliasing in F-K analysis is that **half the wavelength of the
> incoming wave must be larger than the interelement distance**. In the case of the IMS hydrophone
> stations, where the interelement distance is required by the Treaty to be between 1.5 km and 2 km,
> **the frequency must be less than 0.375 Hz**. Above this frequency, spatial aliasing is expected to occur."

```
λ/2 > d  ⟺  f < v / (2·d)          d 取台间距（最坏情况用 d_min）
```

**校验**：他们的 0.375 Hz 就是 `1.5/(2×2) = 0.375` —— 说明判据用法与我们完全一致，也说明
**0.375 Hz 不是经验值，而是按他们自己的台间距算出来的合法上界**。

本台阵代入（`d_min = 100.5 km`）`f_max = v/(2·d_min)`：1.48→0.0074 / 3.6→0.0179 / 4.5→0.0224 /
6.0→0.0299 / 8.12→0.0404 / 11.12→0.0552 Hz —— 与 0.375 Hz 差 7–50 倍，全因台间距差 50 倍。

## 2. 步② 三步预处理（Matsumoto 2016 §3 原文；Zhou 照搬）

> "In the first step, **(1) a sixth-order Butterworth high-pass filter with 1 Hz cutoff** is applied to
> retrieve the T phase and remove the seismic contribution. Successively, **(2) the envelope of the
> high-pass-filtered signal** is calculated to highlight the T-phase contribution. Finally, **(3) a
> sixth-order Butterworth low-pass filter with 0.375 Hz cutoff** is applied to the envelope computed in
> step (2), so that the F-K analysis can be used with the waveform containing the T-phase signal."

动机（原文）：T 波 PSD 只在 **1–100 Hz** 高于地震波（Fig.4，PSD 窗长 5 min），主频 >1 Hz；而三联体
的合法上界只有 0.375 Hz ⇒ 直接用原始波形做 F-K 必然空间混叠 ⇒ 用"高通 → 包络 → 低通"**把 T 波能量
搬进合法频带**。这是回答"为什么用包络、为什么频带这么低"的唯一正确依据。

| 步 | Matsumoto 2016 | Zhou 2021 / 本项目 |
|---|---|---|
| ① 取 T 波 | 6 阶 Butterworth **高通 1 Hz** | **2–8 Hz 带通**（Fig.2c） |
| ② 包络 | \|Hilbert\| | \|Hilbert\|，零平滑 |
| ③ 低通 | 6 阶 Butterworth **低通 0.375 Hz** | **0.25 Hz**（4 阶 filtfilt） |
| ④ 降采样 | — | 20 Hz → 2 Hz（十取一） |

Zhou 原文（§3.2）："applied a **low-pass filter of 0.25 Hz** to the **enveloped vertical component**
seismograms (Figure 2d) highlighting the T phase" —— 同一角色、不同截止。本项目这条链的实测口径与
"包络为什么会有能量"的解释见 `zhou2021-southafrica-fk.md` 的"预处理链原理"节。

## 3. 步③ 目标慢度（Zhou 2021 p.5，别自己预测）

T-P 预测 `13.7 s/°`、T-P 实测 southern 子阵 `10.0 s/°`（11.1 km/s）、T-S `24.7 s/°`（4.5 km/s）。
⚠️ **8.1 / 4.5 / 11.1 都是这两个数的换算（`111.195/u`），不是论文预设的速度**——把它当输入再算回
慢度就是循环；论文两相预测值与式(1) 恒速梯度（6.0/3.6 km/s）也不是同一类量（见
`fk-parameter-selection.md` §6c）。陆地段目标另算：`u_land = 111.195/V`（V = 6.0 / 3.6）⇒ 18.53 / 30.89 s/°。

## 4. 步④ 子阵与平面波假设（Rost & Thomas 2002）

Rost & Thomas §2 [8]："Most array methods assume a plane wave arriving…"；Zhou §3.2 据此把子阵孔径
限制在 **<400 km**（G1 southern 孔径 312 km）。ARF、fk 分析、混叠与孔径-分辨率关系在其 Figure 13/14
与 §[51] 附近 ⇒ 公式化后就是 `fk-parameter-selection.md` §4 的 `Δs = C·111.195/(f·D)`（C = 0.64）、
复制峰间距 `111.195/(f·d_min)` 与几何分辨率地板 `Δs/u ≥ 1.28·d_min/D`。

## 5. 三处必须点明的差异（写方法/汇报时）

1. **截止频率不是普适值**：0.375 Hz 只对 2 km 台间距成立；100.5 km 台间距对应 0.007–0.055 Hz。
2. **传感器与几何**：他们 3 台水听器（三角形、2 km 孔径、压力记录），我们 7 台陆上地震计
   （312 km 孔径、三分量）—— 平面波假设的成立范围与 ARF 形状完全不同，Z/R/T 与三个水听器不可直接类比。
3. **他们的定量结论只有"方位 + 0.6 s/km 量级"**：方位与真实 74.8° 相符；慢度 0.6 s/km ⇔ 1.5 km/s
   （即 11.86 s/°），且只有 3 台设备 ⇒ 论文层面本来就不追求高精度慢度。复现目标应写成
   "峰存在 + 方位 + 落在理论值 ±1 个主瓣内"，而不是"测准到 x%"。

## 6. 本链条的成品

七步（判据 → 预处理 → 目标慢度 → SMAX → 窗 → 频点 → 判读/分辨率）与参数汇总表：
`4_fk_analysis/G1/05_fk_derivation_from_papers.md`。
