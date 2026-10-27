# 波形归一化与纵轴单位（2026-08-24 会话复盘）

用户在多轮追问"单位是速度还是振幅"后定稿。本文记录三套流程的单位规定、差异与踩坑。

## 1. 三个概念必须分开

| 概念 | 内容 | 单位 |
|------|------|------|
| 去响应数据 | `simulate(paz_remove, remove_sensitivity=True)` 输出 | **速度 m/s**（量级 1e-8~1e-7） |
| 绘图归一化 | 除以 P98（全局或分量） | 无量纲，量级 O(1)（±2~±5） |
| 纵轴标签 | 归一化后必须写 | `Norm. Amp.` |

## 2. 判定"这张图是什么单位"的检查顺序

1. 看 SACPZ 头：`INSTGAIN`/`SENSITIVITY` 是否标 `(M/S)`（速度计，如 STS-2、VBB）；标 `(M)` 可能位移。
2. 看去响应 SAC 的实际振幅量级：1e-7 左右 = 速度 m/s；几百~几千 = 原始 counts（未去响应）。
3. 看绘图代码是否除 P98：除过 = 无量纲归一化，标签应 `Norm. Amp.`。

## 3. 归一化因子的三种写法（后果不同）

```python
# A. 分量独立 P98（❌ 三分量图禁用）
p98 = np.percentile(np.abs(cwf[m]), 98)      # 每分量自己的
# → 抹掉三分量相对振幅，无法比较 Z/R/T 强弱

# B. 三分量全局 P98（✅ 波形图标准）
all_abs = np.concatenate([np.abs(z_bp[m]), np.abs(r_bp[m]), np.abs(t_bp[m])])
p98_global = np.percentile(all_abs, 98)
# → 保留三分量相对振幅；波形与包络共用同一因子

# C. 包络 P98 缩放波形（⚠️ 旧模板错误）
ep98 = np.percentile(np.concatenate([env_z, env_r, env_t]), 98)
y / ep98   # 波形也除以包络 P98
# → 波形/包络相对比例失真；2026-08-24 已改为波形全局 P98（B）
```

## 4. 固定 ylim 陷阱

```python
ax.set_ylim(-2.5, 2.5)   # ❌ 固定范围：真实峰值被截断 → "假 ±2"视觉
ax.set_ylim(-y_max, y_max)  # ✅ y_max = 1.05 * max(|波形norm|, |包络norm|)
```

VOI 实测：P 波归一化峰值 ~4.8，固定 ±2.5 时看起来全是 ±2；动态后纵轴 ±5，真实相对振幅可见。

## 5. 三套流程的单位对比

| 流程 | 代表 | 去响应 | 绘图 | 纵轴 | 标签 |
|------|------|--------|------|------|------|
| 新模板 process_image.py | Madagascar 2013-08-17（本次改后） | simulate→速度 | 全局波形 P98 | 动态 | `Norm. Amp.` |
| 旧模板 process_image.py | Southern Africa 2019/2020、Japan 2026 | simulate→速度 | 全局 P98 / 包络 P98 | 固定 ±2.5 | **`Velocity (m/s)`（历史标签错误：已归一化）** |
| 更旧流程 | 2014 Puerto Rico `BPF.png`（H65A） | simulate→速度 | **不归一化** | 真实量级 1e-7 | `Amplitude (2-8 Hz BPF)`（隐含 m/s） |

H65A 实例：TA.H65A BHZ 去响应 max = 4.489e-07，与 BPF.png 图中 ±3.5e-7 吻合 → 未归一化真实速度。

## 6. 处置规则（用户 2026-08-24 明确）

- 归一化统一方向未确认前：**只改当前事件**（"其他事件先不动"）
- 其他事件模板的 `Velocity (m/s)` 标签是错误，统一时批量改为 `Norm. Amp.`
- 用户逐图追问单位是常态，改图后主动说明"数据=速度、图上=归一化 Norm. Amp."
