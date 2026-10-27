# 包络拾取方法 (2026-07-24 定稿)

## 核心方法：50%阈值首次超过

不用 `argmax` 搜绝对最大，改用**首次起跳**检测。

### T-P 拾取

```python
# 理论后0~40s内, 找首个超过50%max的明显波动
m_tp = (t >= ttp) & (t <= ttp + 40)
env_seg = env_z[m_tp]
threshold = 0.5 * np.max(env_seg)
above = np.where(env_seg > threshold)[0]
t_obs_tp = t[m_tp][above[0]] if len(above) > 0 else ttp
```

### T-S 拾取

```python
# 从T-P拾取+5s开始搜，避免拾取同一个峰
ts_start = max(tts, t_obs_tp + 5)
m_ts = (t >= ts_start) & (t <= ts_start + 30)
env_seg2 = env_z[m_ts]
threshold2 = 0.5 * np.max(env_seg2)
above2 = np.where(env_seg2 > threshold2)[0]
t_obs_ts = t[m_ts][above2[0]] if len(above2) > 0 else tts
```

### 偏振窗口（±4s，不重叠）

```python
w_tp = (t_obs_tp - 4, t_obs_tp + 4)
w_ts = (t_obs_ts - 4, t_obs_ts + 4)
# 重叠检测
if w_tp[1] > w_ts[0]:
    shift = w_tp[1] - w_ts[0] + 1
    w_ts = (w_ts[0] + shift, w_ts[1] + shift)
```

## 子窗参数

- **1.0s 子窗，75%重叠**（Flinn 1965 标准）
- T 波主频 3-5Hz，1s 窗口覆盖 3-5 个周期
- 文献参考：Flinn (1965)、Jurkevics (1988)

## 包络计算

```python
from obspy.signal.filter import envelope
from obspy.signal.util import smooth
env = envelope(data)
env_smooth = smooth(env, max(1, int(3 * sr)))  # 3s窗口
```

## 比旧方法改进

| 项目 | 旧 | 新 |
|------|-----|-----|
| 拾取方式 | argmax (绝对峰值) | 50%阈值首次超过 |
| T-S搜索起点 | ttp窗口 | max(tts, T-P+5s) |
| 偏振窗口 | ±5s | ±4s |
| 子窗 | 5.0s | 1.0s |
| 重叠检测 | 无 | 自动后延T-S窗口 |
