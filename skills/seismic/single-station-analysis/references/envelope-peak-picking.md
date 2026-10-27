# 包络拾取标准（2026-07-24 更新）

**核心原则：找首次明显跃升（超过 50% 最大包络值），而非绝对最大值。T-P/T-S 必须独立拾取，不能拾取同一个峰。**

## 拾取方法

```python
# T-P: 理论到时后 0~40s 内, 找首个超过50%最大包络值的点
m_tp = (t >= ttp) & (t <= ttp + 40)
env_seg = env_z[m_tp]
if len(env_seg) > 0:
    threshold = 0.5 * np.max(env_seg)
    above = np.where(env_seg > threshold)[0]
    t_obs_tp = t[m_tp][above[0]] if len(above) > 0 else t[m_tp][np.argmax(env_seg)]
else:
    t_obs_tp = ttp

# T-S: 搜索从 T-P拾取+5s 或 理论到时 开始（避免拾取同一个峰）
ts_search_start = max(tts, t_obs_tp + 5)
m_ts = (t >= ts_search_start) & (t <= ts_search_start + 30)
env_seg2 = env_z[m_ts]
if len(env_seg2) > 0:
    threshold2 = 0.5 * np.max(env_seg2)
    above2 = np.where(env_seg2 > threshold2)[0]
    t_obs_ts = t[m_ts][above2[0]] if len(above2) > 0 else t[m_ts][np.argmax(env_seg2)]
else:
    t_obs_ts = tts
```

## 参数说明

| 参数 | T-P | T-S | 说明 |
|------|-----|-----|------|
| 搜索窗 | 理论后 0~40s | T-P拾取+5s 起 30s | T-S 从 T-P 峰后开始，避免重复拾取 |
| 阈值 | 50% 最大值 | 50% 最大值 | 只拾取明显跃升 |
| 拾取策略 | 首个 > 阈值 | 首个 > 阈值 | 找初动，非最大振幅 |

## 偏振窗口标准（±4s + 防重叠）

```python
w_tp = (t_obs_tp - 4, t_obs_tp + 4)
w_ts = (t_obs_ts - 4, t_obs_ts + 4)
# 检查重叠：若T-P窗口与T-S窗口重叠，T-S后延
if w_tp[1] > w_ts[0]:
    shift = w_tp[1] - w_ts[0] + 1
    w_ts = (w_ts[0] + shift, w_ts[1] + shift)
```

## 设计原则

- **T-S 从 T-P 拾取后开始**：避免 T-P 和 T-S 拾取同一个包络峰（当两者时间接近时常见）
- **±4s 窗口**：±5s 太宽，50%阈值下的初动附近 4s 更精确
- **重叠检测**：如果 ±4s 窗口仍有重叠，T-S 窗口后延至不重叠
