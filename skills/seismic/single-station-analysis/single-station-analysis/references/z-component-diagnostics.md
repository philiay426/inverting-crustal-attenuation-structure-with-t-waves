# Z 分量数据质量诊断

## 快速检查

```python
from obspy import read
tr = read("ZT.VOET..BHZ...SAC")[0]
d = tr.data
print(f"min={d.min():.0f} max={d.max():.0f} mean={d.mean():.0f}")
print(f"负值占比: {(d < 0).sum()/len(d)*100:.1f}%")
```

## 异常判据

| 指标 | 正常 | 异常 |
|:---|:---|:---|
| 负值占比 | ~50% (绕零波动) | **>90% 或 <10%** |
| mean 量级 | ~0 (远小于 std) | 远大于零 |
| |mean| vs |std| | |mean| << std | |mean| ~ std |

## VOET 案例 (南非 Mw6.7)

```
BHZ: min=-734933, max=-44108, mean=-260257 → 100% 负值 → DC 偏移
BHN: min=-92971,  max=116553,  mean=3829     → 正常
BHE: min=-65672,  max=61850,   mean=12252    → 正常
```

## 处理方案

1. **偏振分析**：不受影响（只用 R-T 平面，Z 不参与）
2. **波形剖面**：该站用 `√(R²+T²)` 替代 `√(Z²+R²+T²)`，标注 `(H-only)`
3. **直接排除**：如非关键站，最干净的方案是数据质量剔除
