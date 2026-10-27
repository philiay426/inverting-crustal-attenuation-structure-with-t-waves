# 包络拾取 T-P/T-S 偏振分析标准（2026-07-24 定稿）

## 脚本

- `envelope_polar.py` — 包络拾取版（观测到时±4s）
- `envelope_polar_theory.py` — 理论到时版（理论到时±4s）

## 包络拾取方法

### T-P 拾取
```python
# 理论后0~40s内, 找首个超过50%max的明显波动
m_tp = (t >= ttp) & (t <= ttp + 40)
env_seg = env_z[m_tp]
threshold = 0.5 * np.max(env_seg)
above = np.where(env_seg > threshold)[0]
t_obs_tp = t[m_tp][above[0]] if len(above) > 0 else ttp
```

### T-S 拾取（独立于T-P）
```python
# 从T-P拾取+5s后开始搜索，避免拾取同一个峰
ts_start = max(tts, t_obs_tp + 5)
m_ts = (t >= ts_start) & (t <= ts_start + 30)
env_seg2 = env_z[m_ts]
threshold2 = 0.5 * np.max(env_seg2)
above2 = np.where(env_seg2 > threshold2)[0]
t_obs_ts = t[m_ts][above2[0]] if len(above2) > 0 else tts
```

### 偏振窗口（±4s，带重叠检测）
```python
w_tp = (t_obs_tp - 4, t_obs_tp + 4)
w_ts = (t_obs_ts - 4, t_obs_ts + 4)
# 检查重叠：若T-P窗口与T-S窗口重叠，T-S后延
if w_tp[1] > w_ts[0]:
    shift = w_tp[1] - w_ts[0] + 1
    w_ts = (w_ts[0] + shift, w_ts[1] + shift)
```

## 偏振子窗

| 参数 | 值 |
|------|-----|
| 子窗长度 | **1.0s** |
| 重叠 | 75% |
| 步长 | 0.25s |
| 线性度计算 | Flinn 1965 协方差法 |

**文献依据：**
- Flinn (1965): 子窗覆盖几个完整周期
- Jurkevics (1988): 2-3倍主频周期
- T波 2-8Hz 主频 ~3-5Hz，周期 0.2-0.3s → 1.0s子窗覆盖约3-5个周期

## 2行×3列布局

| T-P Z-R玫瑰 | T-P R-T玫瑰 | T-S R-T玫瑰 |
| 包络+理论/观测到时 | 线性度 | 平面度 |

## 关键参数
- 标题栏：台站名 | Land | SNR | 到时信息 | sub_win=1s
- 包络图：理论到时（虚线）+ 观测到时（实线，50%阈值拾取）
- 线性度/平面度：标注观测到时线
