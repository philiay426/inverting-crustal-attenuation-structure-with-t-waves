# 偏振分析笔记

# T波偏振分析 (polar_analysis.py)

## 脚本位置
- `/home/sy/code/twave/polar_analysis.py` — 综合偏振分析工具
- 调用 process_obs.py 的处理函数

## 功能
对任意台站的三分量数据（Z/NE 或 Z/1/2）做：
1. 2-8 Hz 带通滤波
2. ZNE → ZRT 旋转（BAZ 方向）
3. 动态窗口：T-P/T-S 理论到时的 SNR 截断
   - 噪声窗：P 波之前 (100-250s)
   - 平滑 RMS 包络 (20s 窗)
   - 窗口从 t_arrival - 50s 到 SNR < 1.5 处（最长 120s）
4. 能量比：ER/ET, EZ/EH, log10(ER/ET)
5. 全圆偏振玫瑰图：theta = arctan2(T,R), weight = R²+T²
6. PCA：3×3 协方差矩阵 → linearity, planarity

## 窗口覆盖
台站可通过 `station_overrides` 字典自定义窗口截止时间。

## 自动分量检测
- 优先检测 BHN/BHE（TA 台网）
- 回退 BH1/BH2（IU 台网，如 HRV）
- 文件名 glob 匹配 `{net}.{sta}.{loc}.{comp}*`

## 输出
- 2 行 2 列图：左=T-P，右=T-S
- 第1行：Z/R/T 波形（窗口内归一化，到时虚线标）
- 第2行：振幅加权玫瑰图 + 能量比标注

## 依赖
- `process_obs.py` — SACPZ 解析、到时计算
- `obspy` — 信号处理
- `t_phase_analysis` — tt_t 走时计算

# polar_analysis.py — T-wave Polarization Analysis

Location: `/home/sy/code/twave/polar_analysis.py`

## Purpose

Per-window polarization analysis of T-wave 3-component data:
- ZNE → ZRT rotation
- Energy ratios: ER/ET, EZ/EH, log10(ER/ET)
- Full-circle rose diagrams (R²+T² weighted, 0°/180° = radial, 90°/270° = tangential)
- PCA on 3×3 covariance (linearity = 1−λ₂/λ₁, planarity = 1−2λ₃/(λ₁+λ₂))

## Dynamic SNR-based windowing (original)

- T-P window: [tTP-50, first time SNR<1.5 after tTP], max 120s
- T-S window: [tTS-50, first time SNR<1.5 after tTS], max 120s
- Noise level computed from horizontal RMS in 100-250s window (before P-wave)
- Envelope smoothed with 20s sliding RMS
- If windows overlap (T-P/T-S spacing < 50s), they are still computed independently

## Enhanced Features (June 2026)

See `references/polar-analysis-zr-rose-linplan-20260610.md` for the full update:

- **Z-R rose** for T-P (vertical-radial plane) instead of R-T
- **Fixed window override** for T-S (e.g. 570-600s)
- **T-P window shifted** from tTP-50 → tTP-30
- **Linearity/Planarity bar chart** added as third column
- **English-only labels** (no CJK font dependency)
- **2×3 layout** instead of 2×2

## Usage

```bash
cd /home/sy/code/twave
~/envs/tphase_env/bin/python3 polar_analysis.py STA --network NET
```

Examples:
```bash
# Old style (R-T roses, dynamic windows, 2×2 layout)
python3 polar_analysis.py J59A --network TA

# New style uses hardcoded TS_WIN_START/TS_WIN_END + plane-specific roses
# Edit the constants at top of polar_analysis.py to configure
```

## Output

- Figure: `waveform_contrast/one/polar_analysis_{STA}.png`
- Layout (enhanced):
  - Row 0: T-P Z-R rose / T-S R-T rose / Linearity-Planarity bar chart
  - Row 1: T-P waveform / T-S waveform / Parameter info box

## Key functions

- `load_component(station, network, comp)` — glob-based proc file detection
- `detect_channel_components(station, network)` — auto-detects BHN/BHE vs BH1/BH2
- `rotate_zne_to_zrt(z, n, e, baz)` — standard rotation
- `compute_energy_ratios(zw, rw, tw)` — un-normalized Σ(R²)/Σ(T²) etc.
- `rose_full(theta, weights)` — full-circle histogram, 36 bins
- `pca_3d(z, r, t)` — 3×3 covariance → eigenvalues → linearity/planarity

## Dependencies

- process_obs.py (for `parse_sacpz`, `_compute_arrivals`)
- obspy (bandpass, read)
- t_phase_analysis (tt_t for T-P/T-S arrivals)
- matplotlib (plotting)

# polar_analysis.py — T-wave ZRT Polarization Workflow

Script: `/home/sy/code/twave/polar_analysis.py`

## Usage

```bash
python3 polar_analysis.py J59A --network TA
python3 polar_analysis.py HRV --network IU
```

## Process

All three components MUST be processed through `process_obs.py` first (deconvolve to velocity + 2-8 Hz bandpass). The script auto-detects component naming:
- TA/BH station: BHZ + BHN + BHE
- IU station: BHZ + BH1 + BH2 (via glob, not hardcoded location)

### Steps

1. Load 3 components, identify BAZ from station AZ
2. Apply identical 2-8 Hz bandpass to all three
3. Rotate ZNE → ZRT (R toward source, T = 90° clockwise)
4. Compute noise RMS from pre-P-wave window (100-250s)
5. Compute 20-second sliding RMS envelope of total energy
6. Dynamic window: [T-50, first SNR<1.5] per phase, max 120s
7. Compute energy ratios (un-normalized): ER/ET, EZ/EH, log10(ER/ET)
8. Compute 180°-folded or full-circle rose diagram: θ = arctan2(T, R), weight = R²+T²
9. Compute 3×3 PCA on Z/R/T:
   - Remove mean from each component
   - Covariance = X·Xᵀ/N (3×3)
   - λ1 ≥ λ2 ≥ λ3 (from eigvalsh, ascending order reversed)
   - Linearity = 1 - λ2/λ1
   - Planarity = 1 - 2λ3/(λ1+λ2)

### Station-specific window overrides

In `run()` function, modify `station_overrides` dict:
```python
station_overrides = {
    'HRV': {'tp_end': 1800, 'ts_end': 1800},
}
```

## Output

```
/home/sy/shared/case/waveform_contrast/one/polar_analysis_{STA}.png
```

Two-row figure:
- Row 0, Col 0: Z/R/T waveforms for T-P window (per-window normalization)
- Row 0, Col 1: Z/R/T waveforms for T-S window (per-window normalization)
- Row 1, Col 0: Rose diagram + ER/ET + EZ/EH for T-P window
- Row 1, Col 1: Rose diagram + ER/ET + EZ/EH for T-S window

Arrival times (T-P, T-S) marked as colored dashed vertical lines on waveform panels.

## Interpretation Rules

| ER/ET | Interpretation |
|-------|---------------|
| >1.2  | Radial-dominant (R>T) — P-type or Rayleigh |
| <0.8  | Tangential-dominant (T>R) — S-type or Love |
| 0.8–1.2 | Ambiguous (R≈T) |

| EZ/EH | Interpretation |
|-------|---------------|
| ~0    | Pure Love wave (no vertical) |
| 0.2–0.5 | Mixed Rayleigh+Love |
| >0.8  | Rayleigh-dominant |

| Linearity | Meaning |
|-----------|---------|
| >0.6 | Highly linear particle motion (clean single phase) |
| 0.3–0.6 | Moderate linearity (mixed phases) |
| <0.3 | Near-circular noise (low SNR or multi-phase) |

## Noise Cross-check

For each signal window, a noise window is extracted at `[signal_start-100, signal_start-50]`. If signal and noise have similar ER/ET (Δ < 0.15), the polarization result is unreliable — the signal energy is too low.

## Dependencies

- `process_obs.py` (SACPZ parsing, P/S/T-P/T-S arrival computation via TauP + tt_t)
- ObsPy (signal filtering, SAC I/O)
- `t_phase_analysis` (get_track, tt_t for T-wave travel times)

# Enhanced polar_analysis.py: Z-R Rose, Fixed Windows, Linearity/Planarity Plot

Last updated: 2026-06-10 (updated for P-wave/S-wave mode)
Script: `/home/sy/code/twave/polar_analysis.py`

## What Changed (June 2026 Session)

The original script (June 2024 baseline) had:
- Dynamic SNR-based windows for both T-P and T-S
- R-T rose diagrams for both windows
- 2×2 layout (waveforms × 2, roses × 2)
- Chinese labels in titles and print output
- Linearity/planarity printed to console only

The enhanced version adds:

### 0. Two Analysis Modes

The script now supports **two distinct analysis modes**, controlled by which windows are defined in the `windows` list:

| Mode | Windows | Rose Types | Use Case |
|------|---------|------------|----------|
| **T-wave mode** (original) | T-P window (around T-P), T-S window (around T-S) | T-P: Z-R, T-S: R-T | T-wave polarization discrimination |
| **Body-wave mode** (new) | P-wave window (around P), S-wave window (around S) | P-wave: Z-R, S-wave: R-T | Body-wave polarization at analysis band |

To switch modes, change the `windows` list definition:

**T-wave mode:**
```python
windows = [
    ('T-P', tp_start, tp_end, tTP, 'Z-R'),
    ('T-S', ts_start, ts_end, tTS, 'R-T'),
]
```

**Body-wave mode:**
```python
p_start = max(tP - 30, 0)
p_end   = tP + 60
s_start = TS_WIN_START  # e.g., 570
s_end   = TS_WIN_END    # e.g., 600
windows = [
    ('P-wave', p_start, p_end, tP, 'Z-R'),
    ('S-wave', s_start, s_end, tS, 'R-T'),
]
```

**⚠️ Important:** At 2-8 Hz bandpass, body-wave energy is mostly scattered coda, not coherent P/S. The polarization results at this frequency band are useful for understanding scattered wave properties, not for standard P/S polarization. See "2-8 Hz vs Body Waves" section in the umbrella skill.

When switching to body-wave mode, update ALL references throughout the script:
- Dictionary keys: `'T-P'` → `'P-wave'`, `'T-S'` → `'S-wave'`
- Color scheme: P-wave = `#2ca02c` (green), S-wave = `#d62728` (red, unchanged)
- Phase markers: P-wave window shows P + T-P lines; S-wave window shows S + T-S lines
- Title/judgment labels: all references to 'T-P'/'T-S' must change
- Info text block: window names, ranges, and PCA labels

### 1. Plane-Specific Rose Diagrams

**T-P/P-wave window: Z-R rose (vertical-radial plane)**
- Shows polarization between vertical (Z) and radial (R) components
- `theta = arctan2(R, Z)`, weights = Z² + R²
- Reference lines: Z (blue) at 0°, R (green) at 90°
- Physically: P-wave energy axis in Z-R plane → distinguishes P-type motion from noise

**T-S/S-wave window: R-T rose (horizontal plane)**
- Kept as original: `theta = arctan2(T, R)`, weights = R² + T²
- Reference lines: R (green) at 0°, T (red) at 90°
- Physically: S-wave/Love-wave polarization direction in horizontal plane

Implementation (flexible via `rose_type` per window tuple):
```python
if rose_type == 'Z-R':
    theta = np.arctan2(rw, zw)
    weights = zw**2 + rw**2
    rose_refs = [(0, 'Z', 'blue'), ...]
else:  # 'R-T'
    theta = np.arctan2(tw, rw)
    weights = rw**2 + tw**2
    rose_refs = [(0, 'R', 'green'), ...]
```

### 2. Fixed Window Override

Instead of purely dynamic SNR-based windows, the T-S/S-wave window can be hardcoded:

```python
TS_WIN_START = 570   # fixed S-wave window start (s)
TS_WIN_END   = 600   # fixed S-wave window end (s)
```

When the window is fixed, it no longer depends on arrival time or SNR detection.
Concrete case: S-wave at K62A (S=575s) → fixed window 570-600s captures the S arrival.

### 3. P-wave Window Shift (Backward)

In T-wave mode, T-P window start shifted from `tTP - 50` → `tTP - 30`.
In body-wave mode, P-wave window starts at `tP - 30`, ends at `tP + 60`.
"Backward" in Chinese = later in time (closer to the arrival).

### 4. Linearity & Planarity Bar Chart

New third column in Row 0: grouped bar chart comparing both windows.

```python
# Column 2: Linearity & Planarity
ax_lp.bar(x + offset, values, width, label=metric, color=colors[i])
```

- Green bars = Linearity (1 − λ₂/λ₁, higher = more rectilinear)
- Purple bars = Planarity (1 − 2λ₃/(λ₁+λ₂), higher = more planar/elliptical)
- Eigenvalues λ1/λ2/λ3 shown below bars as text annotation
- Y-range: 0 to 1.05

### 5. English-Only Labels (No CJK Font Dependence)

All print statements, figure titles, axis labels, and annotations are in English.
Required because the user's Linux/WSL system lacks CJK font fallback for matplotlib.

### 6. 2×3 Layout

```
+-----------------------+-----------------------+-----------------------+
| T-P/P Z-R Rose(polar) | T-S/S R-T Rose(polar) | Lin/Plan Bar Chart   |
+-----------------------+-----------------------+-----------------------+
| T-P/P Waveform(Z/R/T) | T-S/S Waveform(Z/R/T) | Parameter Info Box   |
+-----------------------+-----------------------+-----------------------+
```

- figsize: (18, 10)
- height_ratios: [2, 1.5]
- Left column = T-P or P-wave (green theme in body-wave mode), Center = T-S or S-wave (red), Right = metrics + info
- Info box in bottom-right shows: station, event, AZ/BAZ, window ranges, judgments, PCA values

### 7. Waveform Component Labels

Z/R/T traces are all black (`k-`) with text labels on the y-axis margin instead of colored legend:
```python
ax_wf.text(ws - 0.02*(we-ws), offset, 'Z', ...)
ax_wf.text(ws - 0.02*(we-ws), 0, 'R', ...)
ax_wf.text(ws - 0.02*(we-ws), -offset, 'T', ...)
```

### 8. Common P98 Normalization (Matching 08_station_zrt_spec.py)

Since June 2026, the waveform panel uses **global common P98 normalization** instead of per-window local max normalization:

```python
# Collect all data from both windows
all_window_data = []
for wname, wdata in results.items():
    ws, we = wdata['ws'], wdata['we']
    idx = (t >= ws) & (t <= we)
    all_window_data.extend([z_rt[idx], r[idx], t_comp[idx]])

# Single common P98 across ALL 3 components × BOTH windows
common_p98 = max(np.percentile(np.abs(d), 98) for d in all_window_data)
amp = 0.2

# Per-trace scaling:
z_plot = zw_win / common_p98 * amp
r_plot = rw_win / common_p98 * amp
t_plot = tw_win2 / common_p98 * amp
```

This preserves relative amplitudes between P-wave and S-wave windows, so the
user can see that S-wave energy is much stronger than P-wave at 2-8 Hz.

**Why P98 not max:** The 98th percentile is robust against outlier spikes
that might dominate a max-based normalization. P98 captures the typical peak
energy level of each component.

**Why global not per-window:** Previously each window was normalized by its own
`max_win`, making weak P-wave energy appear as large as strong S-wave energy.
Global P98 shows the true amplitude contrast (P ~1e-6 vs S ~1e-5 in velocity units).

## When to Use Each Rose Type

| Window | Rose Plane | Scientific Rationale |
|--------|-----------|---------------------|
| T-P / P-wave | Z-R | P-wave energy is polarized in vertical-radial plane. High Z/R ratio = P-type motion |
| T-S / S-wave | R-T | S-wave/Love-wave energy is horizontal. R/T ratio indicates SH vs SV orientation |

### ⚠️ BHN/BHE `o` Header Missing — Time Axis Misalignment

**Problem:** Proc.sac files for BHZ contain `o=600` (origin time offset), but BHN/BHE proc.sac files downloaded from IRIS may **lack the `o` header entirely**. Without it, `load_component()` defaults to `t_offset = 0` for BHN/BHE, giving them a file-relative time axis while BHZ uses origin-relative time.

**Fix (June 2026):** `load_component()` now accepts a `fallback_o` parameter. When loading 3 components, first read BHZ's `o` value, then pass it to BHN/BHE:

```python
def load_component(station, network, comp, loc='*', fallback_o=None):
    ...
    try:
        o_val = tr.stats.sac.o
    except AttributeError:
        o_val = fallback_o if fallback_o is not None else 0.0
    t_offset = tr.stats.sac.b - o_val
    t = np.arange(len(data)) * delta + t_offset
    return data, delta, t
```

In the `run()` function:
```python
z, delta, t = load_component(station, network, comp_z)
# Get BHZ's o value
bhz_o = None
if matches_z:
    tr_z = read(matches_z[0], headonly=True)[0]
    try: bhz_o = tr_z.stats.sac.o
    except AttributeError: pass
n, _, _ = load_component(station, network, comp_n, fallback_o=bhz_o)
e, _, _ = load_component(station, network, comp_e, fallback_o=bhz_o)
```

Without this fix, BHN/BHE data are shifted by `-o` seconds relative to BHZ (e.g., -600s for `o=600`), causing the rotation to ZRT to mix data from different absolute times.

## Running

```bash
cd /home/sy/code/twave
~/envs/tphase_env/bin/python3 polar_analysis.py J59A --network TA
~/envs/tphase_env/bin/python3 polar_analysis.py K62A --network TA
```

Output: `waveform_contrast/one/polar_analysis_{STA}.png` (200 dpi)

# Polarization Analysis — Flinn (1965) Covariance Matrix Method

## Method

For each sliding window (default 3s, 75% overlap):

1. Form 3×N matrix X = [Z(t), N(t), E(t)]ᵀ (demeaned)
2. Compute covariance matrix C = X·Xᵀ / (N-1) (3×3)
3. Eigenvalue decomposition: C·v = λ·v, λ₁ ≥ λ₂ ≥ λ₃
4. Derive polarization attributes:
   - **Rectilinearity**: R = 1 − (λ₂+λ₃)/(2λ₁) ∈ [0,1]
   - **Planarity**: P = 1 − 2λ₃/(λ₁+λ₂) ∈ [0,1]
   - **Azimuth**: atan2(v_E, v_N), folded to 0-180° (bidirectional axis)
   - **Incidence**: arccos(|v_Z|) ∈ [0,90°]

## Pitfalls Discovered (June 2026)

### Azimuth 0-180° Folding
Without folding, azimuths scatter across 360° because v₁ and -v₁ are equivalent.
Default code folds: `if az >= 180: az -= 180`.

### P-wave Window Width
With 3s sliding window, ±2s around P-wave gives only 5-6 measurement points.
Widening to ±5s gives 13-14 points and dramatically improves azimuth accuracy.
Example: W57A: ±2s → 148° vs BAZ=141° (7° off). ±5s → 142° vs BAZ=141° (1° off).

### Low Linearity in Regional P-waves
At 20-30° distance, P-wave rectilinearity rarely exceeds 0.5-0.6 at TA stations.
Default threshold r_min=0.4 is appropriate; r_min=0.5 filters out too many points.
Only ~4/11 TA stations (W57A, U55A, K60A, J59A) show reliable P-wave azimuths.

### IRIS-downloaded BHN/BHE Headers
FDSN-downloaded SAC files for horizontal components often lack:
- `stla`, `stlo` → read from same-station BHZ file
- `cmpaz`, `cmpinc` → assume standard orientation (BHN=0°, BHE=90°)
- `o`, `nzmsec` → default to 0

### S-wave Polarization Scatter
S-waves have SV+SH components with perpendicular motion, so polarization direction
varies with SV/SH ratio → azimuth naturally scatters. This is expected, not a bug.

### T-wave Polarization
T-wave azimuths are less reliable than P-wave because:
1. T-waves are weaker (lower SNR)
2. Coast-conversion is a distributed source, not a point
3. Multiple conversion points along the coastline produce mixed polarizations

## Script Inventory

| Script | Purpose |
|--------|---------|
| `download_one.py` | Download BHN/BHE + SACPZ from IRIS for stations in raw/ |
| `polarization.py` | Single-station Flinn (generic `--sta`, `--win`, `--r-min`, `--p-win`) |
| `polarization_all.py` | Batch all 11 stations → per-station PNGs + summary map + JSON |
| `plot_polar_map_gmt.py` | GMT terrain map with P/T polarization arrows |
| `run_all_polar.py` | Compute polarization, save JSON (used by GMT script) |

## Key Parameters

- `--win 3.0` — sliding window (s). Smaller = more points, noisier. Larger = smoother, fewer points.
- `--r-min 0.4` — rectilinearity minimum. 0.4 is permissive (catch weak signals). 0.7 is strict (only clean P-waves).
- `--p-win 5.0` — P-wave window half-width (s). Was 2.0s, too narrow.
- `--s-win 3.0` — S-wave window half-width (s).
- `--t-win 3.0` — T-wave window half-width (s).

## GMT Map Guidelines

- Terrain: `@earth_relief_30s`, `cmap='gray'`, `transparency=60`
- Region: `[-86, -64, 14, 48]` (wide enough for all TA stations including N48A at -85°)
- Coast: `0.5p,gray50`, borders: `['1/0.3p,gray60']`
- P arrows: green `#2ca02c` solid, arrow_len=1.5°
- T arrows: blue `#1f77b4` dashed
- Legend: white box at lower-left

# Polarization Analysis Script: polar_analysis.py

Location: `/home/sy/code/twave/polar_analysis.py`

A reusable polarization analysis tool for 3-component T-wave data. 
Designed to be called from command line for any station.

## Usage

```bash
# TA network (default)
python3 /home/sy/code/twave/polar_analysis.py J59A --network TA

# IU network
python3 /home/sy/code/twave/polar_analysis.py HRV --network IU

# Multiple stations
python3 /home/sy/code/twave/polar_analysis.py J59A N63A --network TA
```

## Dependencies

- Imports `parse_sacpz`, `_compute_arrivals`, `apply_bandpass` from `process_obs.py`
- Uses `t_phase_analysis` for tt_t T-P/T-S travel times
- Uses ObsPy for SAC I/O and signal processing

## Output

- PNG figure at `waveform_contrast/one/polar_analysis_{STATION}.png`
- Console output: energy ratios, PCA eigenvalues, judgment

## Figure Layout

**Row 0:** Two-column Z/R/T waveforms
- Left: T-P window (tTP-50 to SNR<1.5 cutoff, max 120s)
- Right: T-S window (tTS-50 to SNR<1.5 cutoff, max 120s)
- T-P/T-S theoretical arrival marked with vertical dashed line + label

**Row 1:** Two-column full-circle rose diagrams
- Left: T-P window amplitude-weighted rose (theta=arctan2(T,R), weight=R2+T2)
- Right: T-S window same
- Reference: R at 0/180 deg, T at 90/270 deg

## Key Methods

- `rotate_zne_to_zrt(z, n, e, baz)` — ZNE to ZRT rotation
- `compute_energy_ratios(zw, rw, tw)` — un-normalized energy ratios
- `rose_full(theta, weights, nbins=36)` — full 360 deg rose histogram
- `pca_3d(z, r, t)` — 3x3 covariance PCA (linearity, planarity)
- `detect_channel_components(station, network)` — auto-detect BHN/BHE vs BH1/BH2

## Dynamic Window Determination

Windows are determined by noise level, not fixed half-widths:
1. Noise RMS computed from pre-P-wave window (100-250s or similar)
2. Total envelope computed as sqrt(Z2 + R2 + T2) with 20s sliding RMS
3. Window start = arrival - 50s
4. Window end = first time after arrival where envelope < 1.5 x noise RMS
5. Capped at max 120s after arrival

## Judgment Logic

- T-P window ER/ET > 1.0 → "偏 radial (R>T)"
- T-S window ER/ET < 1.0 → "偏 tangential (T>R)"
- Otherwise → "ambiguous"

## Known Limitations

- If T-P and T-S are close (<56s, like HRV's 34s), T-S window captures T-P coda → both windows appear radial-dominant
- PCA eigenvalues have large absolute differences between stations (e.g., HRV ~1e-11 vs J59A ~1e-10) — compare relative indices (linearity/planarity), not absolute values
- CJK font warning when Chinese characters are used in titles — cosmetic, matplotlib falls back to a CJK font

# T-P vs T-S Polarization Analysis

Session: 2026-06-04, HRV + M63A + J59A 3-component analysis.

## Workflow

1. **Process all 3 components** through fixed `process_obs.py` (v4+):
   ```bash
   python3 process_obs.py --all  # or specify each file
   ```
2. **Run polarization analysis**:
   ```bash
   python3 /home/sy/code/twave/polar_analysis.py J59A --network TA
   python3 /home/sy/code/twave/polar_analysis.py HRV --network IU
   ```

## Implementation

Script: `/home/sy/code/twave/polar_analysis.py`

### Functions

| Function | Purpose |
|----------|---------|
| `load_component(station, network, comp)` | Load proc file by glob |
| `find_raw_coords(station, network)` | Read STLA/STLO from raw SAC |
| `rotate_zne_to_zrt(z, n, e, baz)` | ZNE → ZRT rotation |
| `compute_energy_ratios(zw, rw, tw)` | ER, ET, EZ, EH from un-normalized data |
| `rose_full(theta, weights)` | Full-circle rose (θ ∈ [-π, π]) |
| `pca_3d(z, r, t)` | PCA on 3×3 covariance (linearity, planarity) |
| `detect_channel_components(station, network)` | Auto-detect BHN/BHE or BH1/BH2 |
| `run(station, network)` | Main analysis |

### Key Design Choices

1. **Per-window normalization**: Each window (T-P, T-S) normalized independently to [-1,1]
2. **Full-circle rose**: No 180° folding — shows directionality (not just axis)
3. **Auto-narrow windows**: If T-P/T-S gap < 50s, ±12s instead of ±25s
4. **Un-normalized energy ratios**: ER/ET, EZ/EH computed on raw (non-normalized) filtered data

## Results

### J59A (TA, AZ=346.9°, gap=56s → ±25s windows)

| Window | s-s | e-e | ER/ET | EZ/EH | Judge |
|--------|-----|-----|-------|-------|-------|
| T-P | 1593 | 1643 | **1.34** | 0.284 | radial (R>T) |
| T-S | 1649 | 1699 | **1.56** | 0.231 | radial (R>T) |

Both windows radial → T-S still contaminated by T-P coda (only 6s gap between windows).

### HRV (IU, AZ=351.4°, gap=34s → ±12s narrow windows)

| Window | s-s | e-e | ER/ET | EZ/EH | Judge |
|--------|-----|-----|-------|-------|-------|
| T-P | 1619 | 1643 | **1.61** | **0.991** | radial, Rayleigh-like |
| T-S | 1653 | 1677 | **1.60** | 0.857 | radial, Rayleigh-like |

EZ/EH≈1 → Rayleigh wave character. Narrow windows still insufficient to separate.

### M63A (TA, gap=22s → ±12s, but SNR too low)

| Window | ER/ET | EZ/EH | Signal E | Noise ER/ET | Reliable? |
|--------|-------|-------|----------|-------------|-----------|
| T-P | 0.72 | 0.380 | 5.3e-15 | 0.68 | ❌ (Δ=0.04) |
| Late peak | 0.55 | 0.472 | 4.5e-15 | 0.67 | ❌ (Δ=0.12) |

M63A signal energy is 6 orders of magnitude below HRV. Noise and signal ER/ET are essentially identical.

## Lessons

1. **T-P and T-S polarization is difficult to separate** with ±25s windows when the gap is <60s.
   Even ±12s narrow windows may not be enough for HRV's 34s gap.
2. **Noise window comparison is essential** for reliability assessment.
3. **Signal energy** is the real limiting factor — M63A is ~10⁶× weaker than HRV.
4. **`simulate_seismometer()` must be used**, not `tr.simulate()`. The latter compresses amplitudes
   by ~200× for IU network data.

# T-wave Polarization Analysis Workflow

## Quick Reference

```bash
# 1. Process 3-component data (if not already done)
cd /path/to/event && python3 process_obs.py --all

# 2. Run polarization analysis
python3 /home/sy/code/twave/polar_analysis.py <STA> --network <NET>
```

## Script Architecture

`/home/sy/code/twave/polar_analysis.py` — reusable single-station T-wave polarization tool.

Key functions (also importable):
- `load_component(station, network, comp)` — glob-based file discovery, handles IU `00/10` locations
- `find_raw_coords(station, network)` — reads SAC header coordinates
- `rotate_zne_to_zrt(z, n, e, baz)` — ZNE → ZRT with back-azimuth
- `compute_energy_ratios(zw, rw, tw)` — EZ, ER, ET, EH
- `rose_full(theta, weights)` — full-circle amplitude-weighted rose histogram
- `pca_3d(z, r, t)` — 3×3 covariance → linearity, planarity

## Window Selection

Dynamic SNR-based windows:
1. Noise RMS from horizontal components (R²+T²) in pre-P window (100s to min(P-50, 250s))
2. 20s sliding-RMS envelope of Z+R+T → SNR = env_smooth / noise_rms
3. T-P window: [tTP-50, first time SNR<1.5 after tTP], capped at 120s
4. T-S window: [tTS-50, first time SNR<1.5 after tTS], capped at 120s

Station-specific overrides via `station_overrides` dict (e.g., `'HRV': {'tp_end': 1800, 'ts_end': 1800}`).

## Output Figure

2×2 grid:
- Row 0: T-P waveform (col 0) | T-S waveform (col 1) — per-window normalization
- Row 1: T-P rose (col 0) | T-S rose (col 1) — full-circle, amplitude-weighted
- Annotations: ER/ET, EZ/EH; T-P/T-S vertical dashed lines on waveforms

## Key Pitfalls

1. **Z-component amplitude mismatch**: If Z is 1000× smaller than horizontals, check all 3 components processed with SAME method. process_obs must use `simulate_seismometer()` (not `tr.simulate()`) consistently.

2. **BH1/BH2 vs BHN/BHE**: IU network uses BH1/BH2; TA uses BHN/BHE. `detect_channel_components()` auto-detects.

3. **File discovery**: Uses glob patterns (`{network}.{station}.*.{comp}.*_proc.sac`), not hard-coded locations.

4. **T-P/T-S proximity**: If interval < 50s, narrow windows (±12s) may still contain T-P coda in T-S window. Report when this happens.

5. **Process before analyze**: Always re-process raw data with the fixed process_obs.py before running polarization analysis, to ensure consistent deconvolution across all components.

# Polarization Analysis for Back-Azimuth Verification

## Purpose

Use P-wave and T-wave polarization directions to verify theoretical back-azimuth (baz) calculated from great-circle path.

## Method

1. Read 3-component waveform (BHZ, BHN, BHE)
2. Bandpass filter (2-8 Hz for T-waves, 0.5-10 Hz for P-waves)
3. Compute covariance matrix in sliding window
4. Eigendecomposition → principal polarization direction
5. Compare observed angle with theoretical baz

## Key Formula

```python
# Covariance matrix of 3-component data
X = np.vstack([z, n, e])
C = np.cov(X)

# Eigendecomposition
eigenvalues, eigenvectors = np.linalg.eigh(C)
v1 = eigenvectors[:, 0]  # principal direction

# Horizontal polarization angle
horizontal_angle = np.degrees(np.arctan2(v1[2], v1[1]))  # atan2(E, N)

# Linearity
linearity = λ1 / (λ2 + λ3)
```

## Data Requirements

- **MUST use raw counts data** (`raw/*.SAC`), NOT response-removed data
- Response-removed data has amplitudes ~10^-7 (too small, unstable polarization)
- Raw counts have amplitudes ~10^4 (good SNR for polarization)
- File extension is **uppercase** `.SAC` in raw directory

## Results: Puerto Rico 2014 Event

### Using Response-Removed Data (OLD - unreliable)

| Station | Phase | Observed | Theory | Δ |
|---------|-------|----------|--------|---|
| Y55A | P-wave | 147.9° | 140.9° | +7° |
| Y55A | T-wave | 117.6° | 140.9° | -23° |

### Using Raw Counts Data (CORRECT - 2026-07-16)

| Station | baz_theory | P_obs | Δbaz | Linearity | Quality |
|---------|-----------|-------|------|-----------|---------|
| Y55A | 140.9° | 113.3° | -27.6° | 1.5 | ✅ Best |
| Y54A | 139.4° | 51.4° | -88.0° | 0.7 | ⚠️ Poor |
| 353A | 131.3° | 61.9° | -69.4° | 1.2 | ⚠️ Poor |
| 256A | 136.5° | -113.2° | -249.7° | 3.0 | ❌ Anomaly |

### Angle Range Note
- `atan2()` returns -180°~180°, theoretical baz is 0°~360°
- Must normalize to same quadrant before comparison:
```python
if p_obs < 0:
    p_obs += 360
```

### Physical Interpretation

**P-wave**: Propagates along great-circle path → polarization direction ≈ theoretical baz ✅

**T-wave**: Converted at ocean-continent boundary → polarization direction points to **conversion point**, not epicenter ⚠️

```
Epicenter ──── Great Circle ────→ Station
     │                              ↑
     └──→ Conversion Point ──→ T-wave ──┘
              ↑
    T-wave polarization points HERE
    (offset from great circle by ~20-30°)
```

## Conclusion

- **P-wave polarization**: Reliable for baz verification (|Δ| < 10°)
- **T-wave polarization**: Points to conversion point, NOT epicenter (|Δ| = 20-30° is normal)
- **For baz verification**: Use P-wave, not T-wave
- **For conversion point**: T-wave polarization direction may help constrain conversion location

## When to Use

| Scenario | Use P-wave | Use T-wave |
|----------|:----------:|:----------:|
| Verify theoretical baz | ✅ | ❌ |
| Find conversion point | ❌ | ✅ |
| Rotate to ZRT | ✅ (use theory baz) | — |
| Station quality check | ✅ | ⚠️ |

## Script Location

`<event_dir>/process_flow/_tmp_polarization/polarization_baz_verify.py`

# Polarization Analysis: PCA + 180° Folded Rose Diagram

Companion to the seismic-twave-figures skill. Covers the polarization analysis
workflow in `polar_analysis.py` at `/home/sy/code/twave/polar_analysis.py`.

## Quick Start

```bash
cd /home/sy/code/twave
python3 polar_analysis.py J59A --network TA
python3 polar_analysis.py HRV --network IU
```

## Input Requirements

- Three-component proc files (BHZ / BHN+BHE or BH1+BH2) processed through
  the fixed `process_obs.py` pipeline (deconvolved to velocity, 2-8 Hz BPF).
- Files in `/home/sy/shared/case/sacprocess/` with naming convention:
  `{NET}.{STA}.{LOC}.{CHN}.M.2014.013.035104_proc.sac`
- The script auto-detects BHN/BHE vs BH1/BH2 and handles both `--` and `00` locations.

## Analysis Steps

### 1. Load & Bandpass (identical for all 3 components)
```python
from obspy.signal.filter import bandpass
z_bp = bandpass(z, 2.0, 8.0, 1/delta, corners=4, zerophase=True)
n_bp = bandpass(n, 2.0, 8.0, 1/delta, corners=4, zerophase=True)
e_bp = bandpass(e, 2.0, 8.0, 1/delta, corners=4, zerophase=True)
```

### 2. ZNE → ZRT Rotation
```python
br = np.radians(baz)  # back-azimuth
r = n * np.cos(br) + e * np.sin(br)
t = -n * np.sin(br) + e * np.cos(br)
# Z is unchanged
```
Where `baz = (180 + 360 - az) % 360`.

### 3. Window Selection
- T-P window: `tTP ± 25s` (or ±12s if overlapping with T-S)
- T-S window: `tTS ± 25s` (or ±12s)
- Auto-detection: if `(tTS - tTP) < 50s`, both windows shrink to ±12s.
- Display range: `tTP - 40s` to `tTS + 40s`

### 4. Energy Ratios (un-normalized)
```python
ez = np.sum(z**2); er = np.sum(r**2); et = np.sum(t**2)
eh = er + et
er_et = er / et
ez_eh = ez / eh
log_er_et = np.log10(er_et)
```

### 5. Rose Diagram (180° Folded)
```python
theta = np.arctan2(t, r)  # arctan2(T, R), range [-π, π]
weights = r**2 + t**2

# 180° folding: 0° and 180° are equivalent, 90° and 270° are equivalent
theta_fold = theta % np.pi  # range [0, π)

# Bin from 0 to π
bins = np.linspace(0, np.pi, 37)  # 36 bins
hist = np.zeros(36)
for i in range(36):
    mask = (theta_fold >= bins[i]) & (theta_fold < bins[i+1])
    hist[i] = np.sum(weights[mask])
hist = hist / np.max(hist)  # normalize
```

**Reference directions on rose:**
- 0°/π → R (radial, green line)
- π/2 → T (tangential, red line)

### 6. PCA on Z/R/T
```python
# Detrend (remove mean)
z0 = z - np.mean(z)
r0 = r - np.mean(r)
t0 = t - np.mean(t)

# 3×N data matrix
X = np.vstack([z0, r0, t0])  # shape (3, N)

# 3×3 covariance matrix
C = X @ X.T / X.shape[1]

# Eigenvalues (ascending → λ3 ≤ λ2 ≤ λ1)
eigvals = np.linalg.eigvalsh(C)
λ3, λ2, λ1 = eigvals

# Polarization indices
linearity = 1 - λ2 / λ1
planarity = 1 - 2 * λ3 / (λ1 + λ2)
```

**Interpretation:**
| linearity | planarity | Interpretation |
|-----------|-----------|---------------|
| ~1        | ~0        | Linearly polarized (one dominant axis) |
| ~0        | ~1        | Planar polarized (energy in a plane) |
| ~0        | ~0        | Spherical (isotropic noise) |
| 0.3-0.6   | 0.4-0.6   | Mixed (typical for T-waves) |

### 7. Noise Reference
Compare each window's ER/ET against a noise window:
```python
noise_start = signal_start - 100
noise_end   = signal_start - 50
```
If signal and noise ER/ET differ by <0.2, the polarization result is unreliable.

## Output Interpretation

| ER/ET | log10(ER/ET) | Interpretation |
|:-----:|:-----------:|---------------|
| >1.2  | >+0.08     | **Radial-dominant** (R>T) — P-wave-like |
| 0.8-1.2| -0.10~+0.08 | **Ambiguous** (R≈T) |
| <0.8  | <-0.10     | **Tangential-dominant** (T>R) — S/Love-wave-like |

| EZ/EH | Interpretation |
|:-----:|---------------|
| ~0    | Pure horizontal (Love wave) |
| 0.3-1.0 | Mixed Rayleigh+Love |
| >1    | Vertical-dominant (P-wave-like) |

## Known Limitations

- T-P and T-S are <60s apart → narrow windows (±12s) reduce sample count,
  making polarization metrics noisier.
- If total energy (λ1) is < 10⁻¹² m²/s², SNR is too low for reliable polarization.
- The script uses `DejaVu Sans` which has no CJK glyphs — title text with
  Chinese characters generates a cosmetic warning. It renders fine via fallback font.

# 包络拾取最新规范 (2026-07-23)

## 拾取策略：分开窗口 + 包络最大值
**结论：分开窗口取 max 是最可靠的方法。** 试过扩窗口、两峰值法、STA/LTA、噪声阈值——都不如原始的"各自窗口取最大包络"稳定。

```python
env_z = smooth(envelope(zf), max(1,int(3*sr)))
# T-P: 理论值 ±60s 内取 max
m_tp = (t>=ttp-60)&(t<=ttp+60)
t_obs_tp = t[m_tp][np.argmax(env_z[m_tp])]
# T-S: 理论值 ±40s 内取 max
m_ts = (t>=tts-40)&(t<=tts+40)
t_obs_ts = t[m_ts][np.argmax(env_z[m_ts])]
```

## 失败的方法记录
| 方法 | 问题 |
|------|------|
| 两最大值全局搜索 | T波包络宽，只有一个主峰 |
| 扩窗口到 +200s | 拾到后面的噪声或二次峰值 |
| 3×背景阈值 + 30%峰值 | 低SNR台站不触发 |
| STA/LTA | 只给初动，用户要的是能量峰值 |

## 偏振结果：澳大利亚 Mw7.0 (2026-07-23)
三站（RK2S/AUNRC/RK1H）包络拾取结果：
| 台站 | T-P实测 | T-S实测 |
|------|---------|---------|
| RK2S | 1004s (+45s) | 1043s (+36s) |
| AUNRC | 965s (+5s) | 1046s (+37s) |
| RK1H | 916s (-45s) | 1049s (+39s) |

T-S 三站一致（+37±2s），T-P差异大但仍在合理范围。
实测到时普遍晚于理论，可能与海水实际路径更长或转换点偏差有关。

# 包络拾取与偏振分析规范

## 包络峰值拾取策略
**结论：分开窗口拾取优于全局峰值搜索。** T 波包络宽且常仅有一个主峰，两个"最大峰值"可能重合。

### 推荐方法
```python
# T-P: 理论值 ±60s 窗口内取 max
m_tp = (t>=ttp-60)&(t<=ttp+60)
t_obs_tp = t[m_tp][np.argmax(env_z[m_tp])]

# T-S: 理论值 ±40s 窗口内取 max  
m_ts = (t>=tts-40)&(t<=tts+40)
t_obs_ts = t[m_ts][np.argmax(env_z[m_ts])]
```

### 失败的尝试
- ❌ 扩窗口到 +200s：拾取到后面的二次峰值或噪声
- ❌ 两个最大值法：T 波包络宽，P/S 初动混在一起
- ❌ 3×背景阈值 + 30%峰值阈值：对低 SNR 台站无效

## 偏振分析双模式
### 文件结构
```
2_single_station_analysis/{STA}/
├── polar_analysis.py          # 统一脚本
├── {STA}_envelope_polar.png   # 实测到时窗口
├── {STA}_theory_polar.png     # 理论到时窗口
├── {STA}_three_filter.png     # 三分量总览（从 waveform_image 复制）
└── ...
```

### 数据源
- 读取 `sacprocess/remove_response/*_resp.sac`（去响应后，非 _proc.sac）
- 带通滤波：2-8 Hz, corners=4
- ZRT 旋转后做偏振

### 两种模式的区别
| 模式 | 偏振窗口 | 标注 |
|------|---------|------|
| envelope_polar | 实测到时 ±5s | 虚线(理论) + 实线(实测) |
| theory_polar | 理论到时 ±5s | 虚线(理论) |

### 网络码注意
- AU 网络：`process_sta('STA', net='AU')`
- S1 网络：`process_sta('STA', net='S1')`
- 务必检查 SAC 文件名中的网络码

# Polarization Window Calibration Workflow

## Issue: Direct P vs PP at Regional Distances

At 15-30° distance, the dominant visible P-wave energy is often **PP** (surface-reflected P),
not the direct P. Always query TauP with `['P','PP','S','SS']` and calibrate windows
against the waveform, not just the first-arrival theoretical time.

### J62A Case Study (Delta=24.55 deg, PR2014)
| Phase | Time | Envelope | Character |
|-------|------|----------|-----------|
| P (direct) | 318s | 1.2e-9 (1.4x noise) | Weak, emergent |
| PP | 350s | 1.8e-7 (180x noise) | **Dominant** |

If polarization windows are centered on tt_P=318s, the PP peak at 350s falls outside.

## Window Calibration Procedure

1. Generate vertical-component waveform + envelope plot (50s major / 10s minor x-ticks)
2. User manually inspects and reports desired window boundaries
3. Update polarization_analysis.py with hardcoded windows
4. Re-run to generate corrected polarization figure

## PR2014 K62A Final Windows
- P: 330-370s, Z-R plane (PP-centered)
- S: 560-620s, R-T plane
- T-P: 1617-1647s, Z-R plane (tt_TP +/-15s)
- T-S: 1660-1700s, R-T plane

## PR2014 J62A (needs recalibration)
- Default windows use tt_P=318s center -> wrong
- Should use PP=350s center or manual calibration

# 单站偏振分析 (2026-07-23)

## 脚本来源
基于 `envelope_polar.py`（BW1H目录），适配任意台站。

## 核心函数

### polarize(t, z, r, tc, ws, we, sub_win=5.0)
- 滑动窗口协方差矩阵特征分解
- 返回：azimuth_ZR, weight, azimuth_RT, weight, linearity_avg, planarity_avg, n_windows
- 子窗口5s，75%重叠

### rose(ax, angles, weights, title, refs)
- 极坐标玫瑰图
- 参考线标注：Z/R/T方向

## 数据要求
- **去响应后的三分量 SAC**（`*_resp.sac`，来自 `sacprocess/remove_response/`）
- 网络+台站必须匹配文件名
- BH1/BH2 自动映射为 BHN/BHE
- 2-8 Hz 带通滤波

## 输出
- `{sta}_envelope_polar.png`: 6-panel图
  - Row1: T-P Z-R玫瑰, T-P R-T玫瑰, T-S R-T玫瑰
  - Row2: 包络+到时, 线性度时间序列, 平面度时间序列

## 偏振解释
- **线性度 ≈ 1**: 纯P波（线性偏振）
- **平面度 ≈ 1**: 纯S波或Love波（平面偏振）
- **T-P窗**: 期望高线性度 + Z-R平面偏振（P波垂直偏振）
- **T-S窗**: 期望高平面度 + R-T平面偏振（S波水平偏振）

## 实测案例
RK2S (澳大利亚Mw7.0):
- T-P: 理论959s, 实测1004s (+45s滞后)
- T-S: 理论1007s, 实测1043s (+36s滞后)
- 实测晚于理论 ~40s，可能海水路径更长或转换点偏差

# T-Wave Polarization: ZRT Rotation + Amplitude-Weighted Rose + Noise Comparison

Alternative to the Flinn (1965) covariance-matrix method. Simpler, works on short windows, and includes a noise-reliability check.

## When to Use

- User has 3-component data and wants to distinguish T-P from T-S/T-Love arrivals
- User wants to check if a late-arriving wave packet has SH/Love character (horizontal > vertical)
- User wants ER/ET and EZ/EH energy ratios
- When the Flinn covariance method is too complex or the window is too short (<5s) for stable eigenvalues

## Workflow

### 1. Load and Rotate ZNE → ZRT

```python
def rotate_zne_to_zrt(z, n, e, baz):
    """R toward source (BAZ direction), T = 90° clockwise"""
    br = np.radians(baz)
    r = n * np.cos(br) + e * np.sin(br)
    t = -n * np.sin(br) + e * np.cos(br)
    return z, r, t
```

BAZ (back azimuth) = |180° - AZ| or (180 + 360 - AZ) % 360.

### 2. Bandpass and Window

Filter all components (typically 2-8 Hz for T-waves). Extract time windows:
- **Signal window**: T-P ±25s (50s window) and/or narrow window (20s)
- **Noise window**: signal_start - 100s to signal_start - 50s (reference)

### 3. Energy Ratios (un-normalized data)

```python
ez = np.sum(zw**2)   # vertical energy
er = np.sum(rw**2)   # radial energy
et = np.sum(tw**2)   # transverse energy
eh = er + et         # horizontal energy

er_et = er / et
ez_eh = ez / eh
log10_er_et = np.log10(er_et)
```

Interpretation:
| ER/ET | EZ/EH | Wave type |
|:-----:|:-----:|-----------|
| >1.5 | ~0 | **Love wave** (T energy, horizontal, radial-dominant) |
| ~1.0 | ~0 | **SH/Love** (transverse-dominant T wave) |
| <0.7 | 0.3-0.5 | **Rayleigh + Love mix** |
| — | >0.5 | Significant vertical → Rayleigh/body wave |

### 4. Amplitude-Weighted Rose Diagram (180° folded)

```python
theta = np.arctan2(t_comp, r)          # arctan2(T, R)
weights = r**2 + t_comp**2             # amplitude-squared weighting
theta_fold = theta % np.pi             # 180° folding (0° = 180°)

# Bin from 0 to π
bins = np.linspace(0, np.pi, nbins + 1)
for i in range(nbins):
    mask = (theta_fold >= bins[i]) & (theta_fold < bins[i + 1])
    hist[i] = np.sum(weights[mask])
```

**180° folding is mandatory** — polarization is axial (0° and 180° are equivalent, 90° and 270° are equivalent). Without folding, the rose shows two opposite lobes and appears split.

### 5. Noise Reliability Check

Compare signal window ER/ET with noise window ER/ET:
- **Signal ≈ Noise** (|Δ| < 0.2) → **Unreliable**. The polarization is dominated by background noise, not the T-wave.
- **|Δ| > 0.5** → Reliable. The polarization direction is signal-specific.
- **Signal energy < noise energy** → Strong warning. Signal window is quieter than the noise window; the "signal" is probably not a coherent wave.

Report all three: signal ER/ET, noise ER/ET, and the difference.

### 6. Plot Layout

Per window (station), 3 panels:

```
+-------------------------------------------+
|  Z/R/T Waveforms (per-component normed)   |
|  Energy ratio text bar (signal + noise)   |
|  Rose diagram (180° folded, polar)        |
+-------------------------------------------+
```

**Waveform normalization**: Per-component max ([−1, 1] each) so the Z trace is visible even when Z << R,T. This is a display choice for the plot; energy ratios are computed from un-normalized data.

## User Preferences (June 2026)

- Rose: amplitude-weighted (R²+T²), NOT count-based
- Theta from R = 0°, T = 90° (Earth-centered reference: green R line at 0°, red T line at 90°)
- Folded 180°: 0° ≡ 180°; only plot [0, π] range
- Reference lines: R (green, solid), T (red, solid)
- Noise window = signal_start - 100s to signal_start - 50s
- Compare signal ER/ET vs noise ER/ET
- If noise comparison shows |Δ| < 0.2, mark result as unreliable
- Plot per-component normalization so Z is visible even when Z ≪ R,T

## Key Pitfalls

1. **ZNE → ZRT requires correct cmpaz** — IU network uses BH1 (cmpaz=0°, N) and BH2 (cmpaz=90°, E). TA network uses BHN/BHE. Always verify with `tr.stats.sac.cmpaz`. Wrong orientation will mix R and T and produce misleading ER/ET ratios.

2. **SACPZ multi-section issue** — IRIS concatenates all location codes into one SACPZ file. Parse only the FIRST response section (stop at second ZEROS or blank line after first CONSTANT). Second-section A0 can be 10¹²× larger, producing wildly wrong amplitudes.

3. **HRV (IU.HRV) horizontal vs vertical amplitude ratio ~1000×** — The horizontal T-wave energy at HRV is ~1000× larger than the vertical. EZ/EH ≈ 0.00 for Love-type T-waves. This is physically correct, not a data error. Plot per-component normalized to see Z waveform shape.

4. **M63A low signal energy** — At ~5×10⁻¹⁵ (velocity m/s), M63A's T-wave energy is ~6 orders of magnitude below HRV. The noise comparison shows signal ER/ET ≈ noise ER/ET → polarization results are unreliable despite SNR > 2.5.

5. **Simulate_seismometer sensitivity** — Always use `simulate_sensitivity=True` in `simulate_seismometer()` to get physically calibrated amplitudes (m/s). Without it, horizontal components may appear 1000× larger than they should be relative to BHZ processed by `tr.simulate()`.

# Large-Scale Polarization Batch Workflow

Full 300+ station 3-C polarization analysis pipeline, developed June 2026 for the 2014 Puerto Rico Mw6.3 event.

## Data Download

```bash
# Download BHN/BHE for ALL stations in raw/ directory
python3 download_all_3c.py --timeout 90
```

Downloads ~640 files from IRIS (350 stations × 2 components). Uses `Client("IRIS")`. Automatically skips existing files. SACPZ instrument responses downloaded alongside waveforms. Total time: ~20-30 minutes.

## Preprocessing

```bash
# Process ALL components (BHZ done already, BHN/BHE new)
python3 process_obs.py --all --skip-existing
```

Output goes to `sacprocess/` with `{NET}.{STA}..{CHN}_proc.sac` naming. Image output to `waveform_image/{STA}_{CHN}/` (component-aware subdirectories to avoid overwriting).

## Polarization Analysis

```bash
# Compute Flinn polarization for all 3-C stations → JSON
python3 run_all_polar_full.py
```

Produces `all_stations_polar.json` with p_az, s_az, tp_az, ts_az for each station. Only processes stations where BHZ+BHN+BHE proc files all exist. Uses:
- TauP iasp91 for P/S arrival times
- Flinn (1965) 3s sliding window, 75% overlap
- Azimuth folded to 0-180° (polarization is bidirectional)
- P-wave window: ±5s around TauP first arrival
- S-wave window: ±3s
- Rectilinearity threshold: 0.4

## T-Wave Polarization (Separate Step)

T-P/T-S arrival times require `tt_t()` which downloads SRTM terrain (~2-3s per station). Run on good stations only:

```bash
python3 compute_twave_polar.py
```

Uses `get_track()` + `tt_t()` for each station. T-P window: ±5s, T-S window: ±5s. For 117 good stations: ~16 minutes.

## GMT Summary Map (Final Style, iterated 5+ versions)

```bash
python3 plot_all_polar_gmt.py
```

**Final map features (June 2026):**
- Light gray SRTM30 terrain (`cmap='gray'`, `transparency=65`), NOT full-color
- **117 stations only** — filtered: `|p_az - (baz%180)| < 35°`. Rejected stations REMOVED (no gray dots).
- **Arrows**: THIN STRAIGHT LINES ONLY (`pen='0.5p'`), NO triangle arrowheads. User said "取消箭头小三角标识，只留一条细直线即可".
- **Arrow lengths**: short — T=`0.8°`, P/S=`1.0-1.2°`
- **Station markers**: black circles `c0.15c` with 4-char labels at `+0.3°/+0.15°`, font `5.5p,Helvetica-Bold`
- **GMT inset** (放大区): `fig.inset(position='jBR+w10c+o0.3c', box='+p0.5p,gray50+gwhite')`. Self-contained map with its own grdimage/coast/plot. BAZ reference dashed lines in gray.
- **Region**: `[-97, -64, 14, 52]` (wide enough for TA station N48A at -85°W)
- **Epicentral distance circles**: `style=f'E{degrees2kilometers(d)}k'`, pen `0.3p,gray60,-`, labeled at ~70% of radius
- **Legend**: white box with colored dots, `0.4-0.55cm` line spacing, font `7-8p,Helvetica`

**GMT arrow rendering CRITICAL PITFALL**: `style='v...'` (GMT vector arrows) does NOT render in PyGMT. Unit `d` not supported for vector length. Use manual line segments:
```python
for r in good:
    rad = np.radians(r['p_az'])
    ex = r['lon'] + alen * np.sin(rad)
    ey = r['lat'] + alen * np.cos(rad)
    fig.plot(x=[r['lon'], ex], y=[r['lat'], ey], pen='0.5p,#2ca02c')
```

**GMT borders PITFALL**: `borders='0.3p,gray60'` fails. Must include feature level prefix: `borders=['1/0.3p,gray60']` (1=national borders).

**Draw order**: T → S → P (T bottom, P top).
- Rejected stations removed entirely
- Black circle markers with station labels
- Thin lines only (no arrowheads)
- Epicentral distance dashed circles every 5°
- Legend at lower-left

## File Inventory

| File | Description |
|------|-------------|
| `all_stations_polar.json` | 312 stations with p_az, s_az, tp_az, ts_az |
| `all_travel_times.json` | 343 stations with TauP P/S times |
| `polarization_results.json` | 11 TA stations only (early subset) |
| `all_polar_map.png/pdf` | Full 312-station map |
| `polar_clean_map.png/pdf` | Filtered 117-station map (final) |
| `polarization_map.png/pdf` | 11-station GMT map |
| `travel_time_map.png/pdf` | T-P arrival time colored map |

# Objective Multi-Station Polarization Summary

## Context

When presenting T-wave polarization results (especially for publication or review), present data **objectively** without force-fitting a specific wave type interpretation. Let the angular distributions and energy ratios speak for themselves.

User preference (2026-06-18): "客观的反应结果，不用刻意贴近love波" — don't bias presentation toward Love wave conclusion.

## Three-Panel Summary Figure

Layout: `figsize=(16, 6)`, GridSpec 3 columns (width_ratios=[1.2, 1, 1])

### Panel 1: Rose Diagram (R-T Plane Direction)

- Data: `RT_dir` from full T-wave window polarization (fold to 0-180° for eigenvector ±ambiguity)
- Bins: 12 bins × 15° (mirror both halves for axial data visual symmetry)
- Axes: 0° labeled **R (Radial)** in green, 90° labeled **T (Transverse)** in red
- No interpretation labels like "Love" or "P/SV" on the figure itself

```python
rt_folded = np.abs(rt_dirs) % 180
theta = np.deg2rad(bin_centers)
theta_full = np.concatenate([theta, theta + np.pi])  # axial symmetry
```

### Panel 2: Bar Chart (Direction Classification)

Three zones based on deviation from R and T:
- **Near Radial**: dev_from_R < 30° (green)
- **Mixed/Intermediate**: 30-60° (gray)
- **Near Transverse**: dev_from_T < 30° (red)

```python
dev_from_R = np.minimum(rt_folded, 180 - rt_folded)  # 0-90 range
dev_from_T = np.abs(rt_folded - 90)                   # 0-90 range
```

Label each bar with count and percentage. No subjective interpretation.

### Panel 3: Energy Partition Boxplot

Boxplot of R%, T%, Z% across all stations:
- R% (green), T% (red), Z% (blue)
- Reference line at 33.3% (equal partition)
- Annotate median values above each box

## PR2014 Reference Values (N=94 stations)

| Metric | Value |
|:-------|:------|
| Near Radial (<30°) | 17 stations (18%) |
| Mixed (30-60°) | 29 stations (31%) |
| Near Transverse (<30°) | 48 stations (51%) |
| Median R% | 33.7% |
| Median T% | 41.3% |
| Median Z% | 23.7% |
| T% > R% | 67/94 (71%) |

## Key Presentation Principles

1. **No wave-type labels on figures** — use "Radial/Transverse" not "P-SV/Love"
2. **Show all stations** — don't cherry-pick high-SNR only for the summary
3. **Include Z component** — non-trivial Z% (~24%) shows the picture is nuanced
4. **Let statistics speak** — 51% near-transverse is a majority but not overwhelming
5. **Avoid loaded language** — "transverse tendency" not "Love wave dominance"

# T-Wave ZRT Polarization with Noise Reliability Check

## When to Use
- User asks for T-wave polarization (ZRT rotation, rose diagram, ER/ET, EZ/EH)
- User wants to distinguish T-P from T-S using polarization
- User provides specific time windows for analysis

## Workflow

### 1. Load & Bandpass (2–8 Hz)
```python
z = bandpass(z_raw, 2, 8, 1/delta, corners=4, zerophase=True)
n = bandpass(n_raw, 2, 8, 1/delta, corners=4, zerophase=True)
e = bandpass(e_raw, 2, 8, 1/delta, corners=4, zerophase=True)
```

### 2. ZNE → ZRT Rotation
```python
baz_rad = np.radians(BAZ)  # BAZ = (180 + 360 - AZ) % 360
R = N * cos(baz_rad) + E * sin(baz_rad)
T = -N * sin(baz_rad) + E * cos(baz_rad)
```

Where `BAZ` = back-azimuth from station to source. For US East Coast from Puerto Rico: BAZ ≈ 188° (adds ~8° due to event at 18.99°N).

### 3. Energy Ratios (un-normalized)
```python
ER = sum(R²), ET = sum(T²), EZ = sum(Z²), EH = ER + ET
# Output: ER/ET, EZ/EH, log10(ER/ET)
```

### 4. Amplitude-Weighted Rose (180° folded)
```python
theta = arctan2(T, R)
theta_fold = theta % np.pi  # fold to [0, π)
weights = R² + T²

# Bin [0, π) into 36 bins, sum weights per bin, normalize by max
for i in range(nbins):
    mask = (theta_fold >= bins[i]) & (theta_fold < bins[i+1])
    hist[i] = sum(weights[mask])
hist /= max(hist)  # normalize
```

Rose plot: bar chart on polar axes, theta=0° = R direction (right), theta=90° = T direction (top).

### 5. Noise Reliability Check
```python
noise_start = signal_start - 100
noise_end = signal_start - 50

diff = abs(ER/ET_signal - ER/ET_noise)
reliable = diff > 0.5      # green
marginal = diff > 0.2      # orange  
unreliable = diff <= 0.2   # red

# Also require signal energy > 2x noise energy
```

## Interpretation

| ER/ET | EZ/EH | Wave type |
|:-----:|:-----:|:-----------|
| >2 | ~0 | Pure Love (T>>R, Z≈0) |
| ~1 | ~1 | Rayleigh (R≈T, Z≈H) |
| >1 | 0.1–0.5 | Mixed radial-dominant |
| <1 | 0.3–0.5 | Mixed transverse-dominant |

## Key Pitfalls

1. **T-S too close to T-P** — At 2000–2800 km range, T-S arrives only 20–60s after T-P. A ±25s window around T-P includes T-S energy. To separate them, need ±5s window and SNR > 5.

2. **Low-SNR = unreliable polarization** — If signal energy is < 2× noise energy, any polarization result is suspect. Report ER/ET for both signal and noise windows and let the user judge.

3. **180° folding is mandatory** — arctan2(T, R) gives [-π, π]. Without folding to [0, π), R+ and R- appear as separate peaks. Fold with `theta % np.pi`.

4. **Use un-normalized data** — Energy ratios require RMS/energy on original velocity data (m/s), NOT normalized traces. Normalization destroys relative amplitude information.

5. **Per-trace normalization in record sections** — For visualizing waveform shape, use per-trace normalization. For polarization amplitude, use global normalization or (better) no normalization.

# T-Wave Polarization Analysis Workflow

## Script: `/home/sy/code/twave/polar_analysis.py`

Single-station 3-component polarization analysis for T-wave discrimination.

### Features
- Auto-detects 3-component data (BHZ/BHN/BHE or BHZ/BH1/BH2)
- Rotates ZNE → ZRT using BAZ
- Computes energy ratios: ER/ET, EZ/EH, log10(ER/ET)
- Full-circle rose diagrams (amplitude-weighted, arctan2(T,R))
- SNR-based dynamic windowing (extends window until SNR < 1.5 with 120s max)
- Station-specific window overrides (e.g., `{'HRV': {'tp_end': 1800, 'ts_end': 1800}}`)
- Theoretical T-P and T-S arrival markers via `tt_t()`
- PCA 3×3 covariance (linearity/planarity kept in code, display optional)

### Usage
```bash
cd /home/sy/code/twave
python3 polar_analysis.py HRV --network IU
python3 polar_analysis.py J59A --network TA
python3 polar_analysis.py X57A --network TA
```

### Output
- PNG to `~/shared/case/waveform_contrast/one/polar_analysis_{STA}.png`
- Console: ER/ET, EZ/EH, log10(ER/ET), PCA eigenvalues, judgment

### Key Findings (2014-01-13 Puerto Rico event)
- **HRV**: ER/ET=1.47-1.61, EZ/EH=0.99 (Rayleigh-type, vertical≈horizontal)
- **J59A**: ER/ET=1.34-1.56, EZ/EH=0.23-0.28 (Love-dominant, low vertical)
- **M63A**: Energy too low for reliable polarization (signal E ~5×10⁻¹⁵ vs HRV ~6×10⁻⁹)
- **X57A**: ER/ET=1.00-1.21, weak linearity (0.02-0.18)

### Pitfalls
- HRV BHZ was initially 1800× too small due to process_obs.py bugs (see `references/obs-deconvolution-pitfalls-june2026.md`)
- EZ/EH=0.000 for HRV was a processing artifact, not physical
- M63A's signal is 6 orders of magnitude weaker than HRV → polarization unreliable
- T-P and T-S windows often overlap (T-S only ~22-56s after T-P), making separation difficult

# Z-T Plane Removal from T-wave Polarization

## Decision (2026-06-25)

Remove Z-T (vertical-transverse) rose diagrams from all T-wave polarization analyses. Only R-T and Z-R planes are used.

## Rationale

- **Z-T is non-standard in literature**: No published T-wave polarization study uses the Z-T plane
- **Physical basis**: Love waves (dominant T-wave mode per De Groot-Hedlin 2020) have near-zero vertical component, making Z-T angles physically meaningless (scattered noise)
- **R-T is the primary diagnostic**: Distinguishes Love (T-dominant, ≈±90°) from Rayleigh (R-dominant, ≈0°/180°)
- **Z-R is supplementary**: Shows Rayleigh wave vertical-radial coupling
- **User confirmed**: "T-coda需要z-t平面吗，我看好像没有人搞这个"

## Affected Functions

| Function | Before | After |
|----------|--------|-------|
| `plot_tcoda_polarization()` | Z-T + R-T roses | R-T + Z-R roses |
| `plot_twave_full_polarization()` | R-T + Z-T + Z-R roses | R-T + Z-R roses + Station Info |

## Layout After Removal

Both functions use the same 2×3 pattern:
```
Row 0:  R-T Rose  |  Z-R Rose  |  Station Info
Row 1:  Z/R/T WF  |  Linearity  |  Planarity
```

The 3rd column naturally fills with station metadata, eliminating wasted blank space.

## Note

`covariance_polarization()` still computes `angles_zt` internally (no performance cost) but these values are no longer used in rose diagrams. They remain available if needed for debugging.
