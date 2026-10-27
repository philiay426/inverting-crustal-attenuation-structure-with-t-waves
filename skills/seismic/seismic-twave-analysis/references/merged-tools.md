# 工具/环境笔记

# Reusable Code Library: /home/sy/code/twave/

## Structure

```
/home/sy/code/twave/
├── polar_analysis.py      # Main T-wave polarization analysis script
├── polar_analysis_自检.txt # Bug/notes log
```

## polar_analysis.py

A self-contained polarization analysis script that calls `process_obs.py`
functions where possible.

### Usage
```bash
cd /home/sy/code/twave
python3 polar_analysis.py J59A --network TA
python3 polar_analysis.py HRV --network IU
python3 polar_analysis.py J59A M63A --network TA   # batch
```

### Key Features
- Auto-detects BHN/BHE vs BH1/BH2 components
- Auto-detects overlapping T-P/T-S windows (shrinks from ±25s to ±12s)
- 180° folded rose diagram with amplitude weighting (R²+T²)
- PCA on 3×3 Z/R/T covariance matrix → linearity, planarity
- Energy ratios: ER/ET, EZ/EH, log10(ER/ET)
- Noise window comparison for reliability check
- Outputs to `waveform_contrast/one/polar_analysis_{STA}.png`

### Key Functions

| Function | Purpose |
|----------|---------|
| `load_component(station, network, comp)` | Load proc file with glob location matching |
| `find_raw_coords(station, network)` | Read SAC header coordinates |
| `rotate_zne_to_zrt(z, n, e, baz)` | ZNE → ZRT rotation |
| `compute_energy_ratios(zw, rw, tw)` | Un-normalized energy: ER/ET, EZ/EH |
| `rose_180(theta, weights)` | 180° folded rose histogram |
| `pca_3d(z, r, t)` | 3×3 covariance → linearity, planarity |
| `detect_channel_components(station, network)` | Auto-detect BHN/BHE vs BH1/BH2 |
| `run(station, network)` | Full analysis entry point |

### Dependencies
- `process_obs.py` (SACPZ parsing, arrivals, bandpass)
- `obspy` (read SAC, filter)
- `t_phase_analysis` (T-P/T-S travel times via `_compute_arrivals`)

# Codex: Reusable T-Wave Analysis Toolkit

Location: `/home/sy/shared/codex/`

A config-driven, event-agnostic Python toolkit for the complete T-wave analysis pipeline.
Designed to be copied and reused for new earthquake events with zero code changes —
only update `config/events.yaml`.

## Directory Structure

```
codex/
├── config/
│   └── events.yaml          # Event parameters (one YAML entry per event)
├── lib/
│   ├── processing.py         # SACPZ parsing, deconvolution, filtering
│   ├── arrivals.py           # TauP + tt_t theoretical arrivals, envelope detection, SNR
│   └── rotation.py           # ZNE→ZRT rotation, energy ratios, rose diagrams
├── scripts/
│   ├── 01_preprocess.py      # Batch preprocessing (BHZ images only)
│   ├── 02_extract_stations.py # Station metadata → Excel
│   ├── 03_map_arrival.py     # GMT arrival-time distribution map
│   └── 04_record_section.py  # Waveform record section (azimuth/time-window params)
├── t_phase_analysis/         # Zhou's T_Phase_Analysis package (moved, not duplicated)
├── process_obs.py            # Legacy single-station manual processing (modified: BHZ-only images)
└── README.md                 # Full 7-step workflow + usage guide
```

## 7-Step Standardized Workflow

| Step | Script | Input | Output |
|------|--------|-------|--------|
| 1 | `01_preprocess.py` | `raw/*.SAC` + `raw/SACPZ/` | `sacprocess/*_proc.sac` + `waveform_image/{STA}/` (BHZ only) |
| 2 | `02_extract_stations.py` | `sacprocess/` + `raw/` | `station_info.xlsx` (NET/STA/LOC/CHA/STLA/STLO/DIST/AZ/BAZ/theoretical arrivals) |
| 3 | `03_map_arrival.py` | `station_info.xlsx` + arrival picks | GMT arrival-time distribution map → `oi/{event_label}/map_*.pdf` |
| 4 | `04_record_section.py` | `sacprocess/` + azimuth/time params | Record section figures (waveform + envelope + arrival markers) |
| 5 | `05_polarization.py` (TODO) | 3-C proc SAC | Multi-window polarization rose diagrams |
| 6 | `06_profile_section.py` (TODO) | `sacprocess/` + azimuth filter | Profile cross-section with elevation + waveform overlay |
| 7 | `07_fk_analysis.py` (TODO) | BHZ proc SAC | F-K spectrum via `t_phase_analysis` |

Steps 5–7 are planned but not yet implemented. Steps 1–4 are production-ready.

## Event Configuration (`config/events.yaml`)

```yaml
puerto_rico:
  label: puerto_rico
  evlo: -66.87
  evla: 18.99
  evdp: 10.0
  date: "2014-01-13"
  raw_dir: /mnt/f/sharedate/2014-01-13-mw63-raw
  case_dir: /mnt/f/sharedate/case
  output_dir: /mnt/f/sharedate/oi/puerto_rico
  ref_el: -1000
  bp_low: 2.0
  bp_high: 8.0

mid_atlantic_2014:
  label: mid_atlantic_2014
  evlo: -45.59
  evla: 23.86
  evdp: 13.7
  date: "2014-07-27"
  raw_dir: /mnt/f/sharedate/2014-07-27-mw61-northern-mid-atlantic-ridge-1
  case_dir: /mnt/f/sharedate/case/2014-07-27-mw61
  output_dir: /mnt/f/sharedate/oi/map_2014_7
  ref_el: -1000
  bp_low: 2.0
  bp_high: 8.0
```

## Key Design Decisions

1. **All paths in YAML** — No hardcoded paths in scripts. `config/events.yaml` is the single source of truth.
2. **Output all to `oi/{event_label}/`** — Figures organized by event, mirroring user's established convention.
3. **BHZ-only images by default** — Preprocessing only generates vertical-component PNGs (raw, remove_response, BPF, spectrogram). Horizontal data stays in proc files for later polarization analysis.
4. **`process_obs.py` preserved as manual entry point** — The original script still works for processing a single station interactively.
5. **`t_phase_analysis/` moved, not copied** — The Zhou package lives at `codex/t_phase_analysis/` and is imported via `sys.path`. No duplication.

## Core Library Modules

### `lib/processing.py`

```python
parse_sacpz(sacpz_path) → (poles, zeros, a0, sensitivity, constant)
  # Stops at second ZEROS section (multi-section SACPZ fix)
  # CONSTANT = A0 × sensitivity — separates correctly

remove_response_obs(data, delta, paz, pre_filt, water_level_db=60.0) → velocity
  # Uses simulate_seismometer(), NOT tr.simulate()
  # water_level in dB (60.0), not ratio
  # np.gradient() for velocity (not tr.differentiate())

bandpass_filter(data, delta, low, high, corners=4) → filtered
  # Zero-phase Butterworth
```

### `lib/arrivals.py`

```python
compute_theoretical_arrivals(evlo, evla, evdp, stlo, stla) → dict
  # Returns: tt_p, tt_s, tt_tp (T-P), tt_ts (T-S), tt_tlove (T-Love)
  # Uses TauP iasp91 for P/S, tt_t() for hybrid-path T-wave
  # T-Love = tt_ocean + land_km / v_love (default 3.6 km/s)

detect_twave_arrival(data, delta, tt_theory, dist_km, window=±120s) → dict
  # Envelope peak in theoretical ±window
  # Celerity filter: rejects if dist_km / t_peak >= 1.8 km/s
  # Returns: t_peak_s, arrival_min, celerity, snr

compute_snr(data, delta, t_signal_start, t_signal_end, t_noise_start, t_noise_end) → float
```

### `lib/rotation.py`

```python
rotate_zne_to_zrt(n, e, z, baz) → (r, t, z)
  # Standard rotation matrix using back-azimuth

compute_energy_ratios(z, r, t) → (er_et, ez_eh, log10_er_et)
  # ER = sum(R²), ET = sum(T²), EZ = sum(Z²), EH = ER+ET

compute_polarization_rose(r, t, window_indices) → (theta_bins, weights)
  # Amplitude-weighted (R²+T²) rose, full 360° circle
  # Returns normalized histogram [0, 2π)

pca_polarization(z, r, t) → (linearity, planarity, azimuth, incidence)
  # Flinn (1965) 3×3 covariance PCA
```

## Usage for a New Event

```bash
# 1. Add event to config/events.yaml
# 2. Run the pipeline
cd /home/sy/shared/codex

python3 scripts/01_preprocess.py --event new_event_label --all
python3 scripts/02_extract_stations.py --event new_event_label
python3 scripts/03_map_arrival.py --event new_event_label
python3 scripts/04_record_section.py --event new_event_label --azimuth 327 350

# Output: oi/{new_event_label}/map_*.pdf, record_section_*.pdf
```

## Processed Events

| Event | Date | Mw | Stations | T-wave SNR | Key Finding |
|-------|------|-----|----------|------------|-------------|
| Puerto Rico | 2014-01-13 | 6.3 | 343 | Good (311 cel-filtered) | Two clear azimuth profiles (327°/350°); strong attenuation signal |
| Mid-Atlantic Ridge | 2014-07-27 | 6.1 | 252 | Weak (13/198 SNR>1.5) | T-wave barely above noise; not suitable for attenuation analysis |

## Known Pitfalls

- **SACPZ multi-section**: IRIS concatenates multiple location codes. `parse_sacpz()` stops at second ZEROS — see `references/process-obs-bugs-20260604.md`.
- **`tr.simulate()` vs `simulate_seismometer()`**: Former produces ~200× smaller amplitudes. Always use `simulate_seismometer()` — see `references/obs-deconvolution-pitfalls-june2026.md`.
- **`tt_t().dist_land` in degrees**: Convert with `degrees2kilometers()` before using as km.
- **STA/LTA on velocity data**: Fails for T-wave detection. Use envelope peak + celerity filter instead.
- **TauP S-wave**: Always query `phase_list=['P','S']` simultaneously — never use `tt_p * 1.7`.

# Cross-Platform Path Handling (2026-07-22)

## Pattern used in process_obs.py and test scripts
```python
import os

# Auto-detect t_phase_analysis path
if os.name == 'nt':  # Windows
    T_PHASE_PATH = r'E:\sharedate\codex\t_phase_analysis'
else:  # WSL / Linux
    T_PHASE_PATH = '/home/sy/codex/t_phase_analysis'
sys.path.insert(0, T_PHASE_PATH)

# Auto-detect BASE_DIR
DEFAULT_BASE_WINDOWS = r'E:\sharedate\case\...'
DEFAULT_BASE_WSL = '/mnt/e/sharedate/case/...'

BASE_DIR = os.environ.get('SAC_BASE_DIR',
                          DEFAULT_BASE_WINDOWS if os.name == 'nt' else DEFAULT_BASE_WSL)
```

User works on both Windows (direct) and WSL (via /mnt/e/). All scripts should support both. Environment variable `SAC_BASE_DIR` allows override for testing.

# GeographicLib ArcDirect Bug: Using Meters Instead of Degrees

## The Bug

When computing azimuth arcs from an epicenter for pyGMT maps, `Geodesic.WGS84.ArcDirect(lat, lon, az, dist)` expects the fourth parameter `dist` to be in **degrees** (arc angle), not meters.

**Wrong code** (produced by thinking of distance in km/m):
```python
g = Geodesic.WGS84.ArcDirect(evla, evlo, az, 4500 * 1000)
```
This interprets 4,500,000 as degrees, wrapping around the Earth ~12,500 times, placing the endpoint at completely wrong coordinates.

**Example** (Barbados 2014 event, 14.67°N, 58.93°W):
- AZ 325° endpoint should be: 45.08°N, 90.65°W (US interior)
- Buggy endpoint was: 14.67°N, 29.88°E (African coast!)

## The Fix

Use `Geodesic.WGS84.Direct(lat, lon, az, dist_m)` where `dist_m` is the distance in **meters**.

```python
g = Geodesic.WGS84.Direct(evla, evlo, az, 4_200_000)  # 4200 km in meters
```

## Correct Function Signature

| Function | Returns | Distance Unit |
|----------|---------|---------------|
| `Geodesic.WGS84.Direct(lat, lon, az, s12)` | Endpoint after s12 **meters** | **Meters** |
| `Geodesic.WGS84.ArcDirect(lat, lon, az, a12)` | Endpoint after a12 **degrees** | **Degrees** |

## Verification

```python
from geographiclib.geodesic import Geodesic
# CORRECT
g = Geodesic.WGS84.Direct(14.6682, -58.9272, 325, 4200 * 1000)
# → 45.08°N, 90.65°W ✓

# WRONG
g = Geodesic.WGS84.ArcDirect(14.6682, -58.9272, 325, 4200 * 1000)
# → Wrong coordinates ✗
```

## Azimuth Arc Helper Function

```python
def azimuth_arc(evlo, evla, az, max_m, npts=200):
    """Compute geodesic azimuth line. max_m = distance in meters."""
    lons, lats = [], []
    for dist_m in np.linspace(0, max_m, npts):
        g = Geodesic.WGS84.Direct(evla, evlo, az, dist_m)
        lons.append(g['lon2'])
        lats.append(g['lat2'])
    return lons, lats
```

## Discovered

2026-06-16, during Barbados 2014 Mw 6.5 arrival map debugging. Cost ~30 min of trial-and-error plus vision-subagent analysis before root cause was found. User spotted the visual bug ("线条出bug了").

# GeographicLib ArcDirect vs Direct: Degrees vs Meters

## The Bug

`Geodesic.WGS84.ArcDirect(lat, lon, az, a12)` — the parameter `a12` is **arc angle in degrees**, NOT distance in meters. This is a critical distinction that cost ~30 minutes of debugging (Barbados 2014 event, June 2026).

## Symptom

Azimuth arcs drawn from an epicenter appear in the **wrong direction**. For example, from 14.67°N, 58.93°W (Barbados area):

| Azimuth | Expected endpoint | Wrong endpoint (ArcDirect with meters) |
|:-------:|:-----------------:|:--------------------------------------:|
| 325° | 45.08°N, 90.65°W (US interior) | 14.67°N, 29.88°E (Africa!) |
| 350° | 50.85°N, 68.64°W (Canada) | 14.67°N, 45.90°W (Atlantic) |

Passing `4500 * 1000 = 4,500,000` as the `a12` parameter to `ArcDirect` is interpreted as **4,500,000 degrees** of arc — roughly 12,500 times around the globe. The endpoint is essentially random.

## Fix

Use `Geodesic.WGS84.Direct(lat, lon, az, s12)` where `s12` is in **meters**:

```python
# CORRECT — distance in meters
g = Geodesic.WGS84.Direct(evla, evlo, az, 4200 * 1000)  # 4200 km

# WRONG — ArcDirect interprets value as degrees
g = Geodesic.WGS84.ArcDirect(evla, evlo, az, 4200 * 1000)  # 4,200,000°
```

## Verification Code

```python
from geographiclib.geodesic import Geodesic

# ArcDirect with 40° arc (correct usage for ArcDirect)
g = Geodesic.WGS84.ArcDirect(14.6682, -58.9272, 325, 40)
print(f"ArcDirect(40° arc): {g['lat2']:.2f}°N, {g['lon2']:.2f}°W")
# → 44.75°N, 90.07°W ✓

# Direct with 4500 km in meters (correct usage for Direct)
g = Geodesic.WGS84.Direct(14.6682, -58.9272, 325, 4500 * 1000)
print(f"Direct(4500km): {g['lat2']:.2f}°N, {g['lon2']:.2f}°W")
# → 45.08°N, 90.65°W ✓
# Both match when used correctly — differences due to ellipsoid vs sphere

# WRONG — ArcDirect with meters
g = Geodesic.WGS84.ArcDirect(14.6682, -58.9272, 325, 4500 * 1000)
print(f"ArcDirect(4,500,000°!): {g['lat2']:.2f}°N, {g['lon2']:.2f}°W")
# → 14.67°N, 29.88°E  ✗ (wrapped globe ~12,500×)
```

## Common Pattern (Azimuth Arc for Maps)

```python
def azimuth_arc(evlo, evla, az, max_m, npts=200):
    """Compute points along a geodesic azimuth line from epicenter.
    max_m = distance in meters."""
    lons, lats = [], []
    for dist_m in np.linspace(0, max_m, npts):
        g = Geodesic.WGS84.Direct(evla, evlo, az, dist_m)
        lons.append(g['lon2']); lats.append(g['lat2'])
    return lons, lats

# 4200 km arc at 325° — use meters!
lons, lats = azimuth_arc(EVLO, EVLA, 325, 4200 * 1000)
fig.plot(x=lons, y=lats, pen='1p,#1f77b4,--')
```

## API Summary

| Function | 4th parameter | Unit | Use case |
|----------|:-------------:|:----:|:--------:|
| `Geodesic.Direct(lat1, lon1, azi1, s12)` | `s12` | **meters** | ✅ Azimuth arcs, distance-based |
| `Geodesic.ArcDirect(lat1, lon1, azi1, a12)` | `a12` | **degrees** | ❌ Rarely needed; arc-angle based |

## Context

This bug was discovered and fixed during the 2014-02-18 Mw 6.5 Barbados event analysis (June 2026). The user noticed "线条出bug了" — the azimuth arcs were going in completely wrong directions.

# GeographicLib `ArcDirect` Pitfall (2026-06-16)

## Symptom

Azimuth arcs drawn from an epicenter at (14.67°N, 58.93°W) using `ArcDirect` with a 4500 km target distance appeared at **completely wrong coordinates**:
- AZ=325° endpoint: **14.67°N, 29.88°E** (eastern hemisphere) — should be **45.08°N, 90.65°W** (US Midwest)
- AZ=350° endpoint: **14.67°N, 45.90°W** (open Atlantic) — should be **50.85°N, 68.64°W** (Canada)

The lines wrapped many times around the globe instead of going toward the US East Coast.

## Root Cause

```python
# WRONG — ArcDirect 4th parameter is ANGULAR DEGREES, not meters
g = Geodesic.WGS84.ArcDirect(evla, evlo, az, 4500 * 1000)
#       4500 * 1000 = 4,500,000 degrees ≈ 12,500× earth circumference!
```

The GeographicLib API:
- `Geodesic.ArcDirect(lat1, lon1, azi1, a12)` — `a12` is **arc length in degrees** (angular distance on the ellipsoid)
- `Geodesic.Direct(lat1, lon1, azi1, s12)` — `s12` is **distance in meters**

## Fix

```python
# CORRECT — Direct with distance in meters
g = Geodesic.WGS84.Direct(evla, evlo, az, 4500 * 1000)  # 4,500,000 meters

# Also correct — ArcDirect with degrees
g = Geodesic.WGS84.ArcDirect(evla, evlo, az, 40)  # 40° arc length (≈ 4440 km)
```

## Validation

```
AZ 325°: (14.67°N, 58.93°W) →            Correct: 45.08°N, 90.65°W
AZ 350°: (14.67°N, 58.93°W) →            Correct: 54.28°N, 70.02°W
AZ 325° ArcDirect(4,500,000 "deg"):  →    Wrong:  14.67°N, 29.88°E
```

## Lesson

When the result of a geodesic computation places a point thousands of km from the expected hemisphere, check whether you're passing **meters to a function that expects degrees**. The `ArcDirect` / `Direct` naming ambiguity is the most common cause.

# GeographicLib Azimuth Arcs in pyGMT

Drawing great-circle azimuth arcs from an epicenter on a pyGMT station map.

## The Bug: ArcDirect Expects Degrees, Not Meters

**Wrong** (causes arcs to wrap around the globe):
```python
from geographiclib.geodesic import Geodesic
g = Geodesic.WGS84.ArcDirect(lat, lon, az, 4500 * 1000)
```
`ArcDirect`'s 4th parameter is **arc angle in degrees**. Passing `4500 * 1000` = 4,500,000° wraps the earth ~12,500 times. A `325°` arc from 14.67°N, 58.93°W erroneously ends at 14.67°N, 29.88°E instead of 45°N, 90°W.

**Correct** (use `Direct()` with meters):
```python
g = Geodesic.WGS84.Direct(lat, lon, az, distance_meters)
```

## Correct Pattern

```python
from geographiclib.geodesic import Geodesic
import numpy as np

def azimuth_arc(evlo, evla, az, max_m, npts=200):
    """Compute points along a geodesic azimuth line from epicenter.
    max_m = distance in meters."""
    lons, lats = [], []
    for dist_m in np.linspace(0, max_m, npts):
        g = Geodesic.WGS84.Direct(evla, evlo, az, dist_m)
        lons.append(g['lon2'])
        lats.append(g['lat2'])
    return lons, lats

# Usage in pyGMT:
fig.plot(x=ax_lons, y=ax_lats, pen='1p,#1f77b4,--')
fig.text(x=ax_lons[-1], y=ax_lats[-1], text=f'{az}°',
         font='8p,Helvetica-Bold,#1f77b4', justify='LM')
```

## Arc Length Constraints

For a map region like `[-90, -55, 10, 50]`, limit arc length to ~4200km to avoid clipping at map edges:
- AZ 325° from 14.67°N, 58.93°W → hits 90°W at ~4500km
- AZ 350° from 14.67°N, 58.93°W → hits 50°N at ~1250km

## Alternative: ArcDirect with Angular Distance

If you want to use `ArcDirect`, pass the angular distance in **degrees**:
```python
arc_deg = 40  # ~4500km in degrees
g = Geodesic.WGS84.ArcDirect(evla, evlo, az, arc_deg)
```

# GeographicLib: Direct() vs ArcDirect() — CRITICAL unit difference

## The Bug (June 2026)

Used `Geodesic.WGS84.ArcDirect(evla, evlo, az, 4500 * 1000)` thinking the 4th parameter was **meters**. It is actually **degrees** of arc angle. Passing 4,500,000 degrees caused the arc to wrap around the earth ~12,500 times, placing the endpoint at completely wrong coordinates.

**Symptom**: AZ 325° arc from (14.67°N, 58.93°W) landed at 29.88°E instead of 90.65°W — `ArcDirect` plus wrong units caused the arc to travel east across Africa instead of NW toward the US.

## Correct API

```python
from geographiclib.geodesic import Geodesic

# CORRECT — Direct() with distance in meters
g = Geodesic.WGS84.Direct(lat, lon, azimuth, distance_meters)

# WRONG — ArcDirect with meters
g = Geodesic.WGS84.ArcDirect(lat, lon, azimuth, distance_meters)  # ❌
```

## ArcDirect Usage

`ArcDirect` is only correct when the 4th parameter is an **angular** arc distance in **degrees**:
```python
# Correct — 40 degrees of arc
g = Geodesic.WGS84.ArcDirect(lat, lon, az, 40)
```

For all seismic workflows (azimuth lines, range rings, epicentral distance arcs), use `Direct()` with meters.

## Verification

```python
# Quick check — endpoint should go NW from Barbados to US
g = Geodesic.WGS84.Direct(14.6682, -58.9272, 325, 4200 * 1000)
print(f"{g['lat2']:.2f}N, {g['lon2']:.2f}W")  # Should be ~45°N, ~91°W
```

## Reference
- GeographicLib Python docs: https://geographiclib.sourceforge.io/Python/doc/

# matplotlib 3.11 中文渲染

## 背景

matplotlib 3.11 (2026-06-12) 重写了字体引擎（libraqm + HarfBuzz + SheenBidi），彻底支持多语言文字。但中文字体仍需系统安装后手动配置。

## 快速启用

在脚本开头加一行：

```python
plt.style.use('/home/sy/tools/mpl_chinese.mplstyle')
```

## mpl_chinese.mplstyle 内容

位置：`/home/sy/tools/mpl_chinese.mplstyle`

```ini
font.family: WenQuanYi Zen Hei
axes.unicode_minus: False
```

## 系统字体

WSL 自带文泉驿正黑：
- 路径：`/usr/share/fonts/truetype/wqy/wqy-zenhei.ttc`
- 字体名：`WenQuanYi Zen Hei`

其他可用中文字体：
- Droid Sans Fallback：`/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf`
- Noto CJK：`sudo apt install fonts-noto-cjk`

## 环境

**统一使用 `/home/sy/miniforge3/envs/obspy/bin/python`**（pygmt 0.19.0 + obspy 1.5 + openpyxl/pillow 齐全）。
备用（勿用于绘图，pygmt 0.18 subplot 有 bug）：
- `~/envs/tphase_env/` — 旧主科研环境
- `~/codex/venv/` — 旧辅助环境

## 注意事项

- 每次清除 matplotlib 缓存后（`rm ~/.cache/matplotlib/fontlist*.json`），需重新运行一次脚本让缓存重建
- `findfont: Failed to find font weight normal` 是已知警告，不影响渲染
- 如果中文仍不显示，检查 `plt.rcParams['font.family']` 是否正确

# PDF Text Extraction on WSL

## Problem
`web_extract` blocks local file paths ("Blocked: URL targets a private or internal network address").

## Solution
Use `pdftotext` from poppler-utils (usually pre-installed on Ubuntu/WSL):

```bash
# Extract text
pdftotext "/mnt/c/Users/username/Desktop/paper.pdf" /tmp/output.txt

# Check if installed
which pdftotext
# If not: sudo apt install poppler-utils
```

## Workflow
1. Windows PDF paths: `/mnt/c/Users/<username>/Desktop/...`
2. Extract to `/tmp/` (cleaned on reboot)
3. Read with `read_file` tool
4. Large PDFs: use offset/limit pagination

## Example
```bash
pdftotext "/mnt/c/Users/75978/Desktop/论文/T波/Zhou2021.pdf" /tmp/zhou2021.txt
wc -l /tmp/zhou2021.txt  # Check size
```

## Why Not Python PDF Libraries?
- `pdfplumber`/`PyPDF2` need `uv pip install --system` (PEP 668)
- `pdftotext` is already installed and fast
- No virtual environment needed

# Spatial Outlier Filter for T-Wave Arrival Picking

## Problem

After STA/LTA-based T-wave arrival detection with a low threshold (≥0.03, cel<1.8 km/s), some stations (~7% in the Caribbean dataset) have their STA/LTA peak at late-time noise (32-42 min) instead of the actual T-wave arrival (20-30 min). These are typically stations with weak T-wave signals where the background noise dominates.

Example from the 2014 Puerto Rico dataset:
- TA.158A: detected at 37.5 min, but 5 nearest neighbors at 22-26 min → deviation +13.9 min
- TA.154A: detected at 37.3 min, neighbors at 22-26 min → deviation +11.6 min

## Solution: Nearest-Neighbor Median Filter

For each detected station, compare its arrival time to the **median of its 5 nearest neighbors** (by geographic distance). Flag as outlier if `|tt_station - median| > 3 min`.

```python
import numpy as np

def spatial_outlier_filter(stations, max_deviation_min=3.0, n_neighbors=5):
    filtered = []
    removed = []
    for r in stations:
        # Compute distances to all other stations
        dists = [(np.sqrt((r['lon']-r2['lon'])**2 + (r['lat']-r2['lat'])**2), r2['tt_min'])
                 for r2 in stations if r2 != r]
        dists.sort()
        # Median of n nearest neighbors
        nn_tts = [t for _, t in dists[:n_neighbors]]
        median = np.median(nn_tts)
        if abs(r['tt_min'] - median) <= max_deviation_min:
            filtered.append(r)
        else:
            removed.append(r)
    return filtered, removed
```

## Results (2014 Puerto Rico, 343 stations)

| Stage | Count | TT range |
|-------|-------|----------|
| Raw detection (STA/LTA≥0.03, cel<1.8) | 311 | 20.6-37.7 min |
| After spatial filter (deviation ≤3 min) | **290** | **22.4-31.8 min** |

The filter removed 21 stations (6.8% of detections), tightening the arrival-time range and eliminating the physically implausible late arrivals.

## When to Use

- Always apply after STA/LTA-based T-wave detection
- After celerity < 1.8 km/s filter but BEFORE figure generation
- The 3-min threshold and 5-neighbor window are calibrated for the USArray TA network (~50 km station spacing). Adjust for sparser networks.

# Station Extraction Pipeline (02_extract_stations.py)

## Correct Parameters
Always use events.yaml config — do NOT hardcode coordinates:
- 2014-01-13 Mw 6.3 Puerto Rico: 18.99°N, -66.87°W, depth=10.0 km

## Pipeline Flow
```
02_extract_stations.py
├── Read config from events.yaml
├── Read SAC files (*.BHZ*_proc.sac)
├── For each station:
│   ├── tt_t() → theoretical T-P, T-S, land_km, ocean_km
│   ├── detect_arrival_envelope() → observed t_obs (envelope peak ±120s)
│   └── compute_snr() → SNR
├── stations.xlsx (BOTH theoretical + observed)
│   Columns: STA, LON, LAT, DIST, LAND, OCEAN, AZ, BAZ, 
│            P_dir, S_dir, T-P_th, T-S_th, T_obs, SNR
└── gmt_stations.txt (observed arrival times for mapping)
    Format: lon lat arrival_min station
```

## Excel Format (header at row 3)
- Row 1: Event info (merged)
- Row 2: Chinese headers
- Row 3: English headers (read with `pd.read_excel(..., header=[2])`)
- Data starts row 4

## Common Mistakes
1. **Wrong header row**: Old Excel was `header=[1]`, new format is `header=[2]`
2. **Wrong column names**: Old had `T-P(s)`, new has `T-P_th` and `T-S_th`
3. **Missing observed column**: Old Excel only had theoretical T-P/T-S; correct has both theoretical AND observed
4. **Wrong parameters**: Old Excel used 19.0428°N/-66.8096°W, depth 20km. Correct is events.yaml: 18.99°N/-66.87°W, depth 10km

## Envelope Detection Issues
- `detect_arrival_envelope()` searches theory_s ± 120s window
- Celerity check uses epicentral distance (not ocean distance) — can fail for stations with short ocean paths
- L63A: observed 1510s vs theoretical 1624s (113s error, SNR=15.9!) — picked noise peak
- **Recommendation**: Use theoretical T-P for analysis, not observed

## Cross-Azimuth T-P Coloring
T-P = ocean_time + land_P_time. Ocean time is 94% of total.
- T-P vs OCEAN: r = 0.981 (nearly perfect correlation)
- T-P vs LAND: r = -0.201
- Within same azimuth group: coastal stations arrive earlier than inland ✓
- Across azimuths: 310° stations (Florida) appear earlier than 350° stations (New England) because ocean path is 1000km shorter
- **This is physically correct** — T-wave reaches Florida coast ~275s before New England
- For "coast→inland gradient" visualization, use land distance coloring instead of total T-P

# Station Info Quick-Extraction (No Preprocessing)

## When to Use

When you receive a new event with raw SAC files and need a fast station inventory
before running the full preprocessing pipeline (which can take hours for 100+ stations).

## Approach

Read SAC **headers only** (`headonly=True`), pick the longest segment per station,
compute distances and TauP arrivals. No deconvolution, no filtering.

## Key Optimizations

### 1. headonly=True
```python
tr = obspy_read(fpath, headonly=True)[0]  # ~100x faster than full read
```

### 2. Longest Segment Selection
For stations with fragmented time segments (multiple SAC files with different start times):
```python
station_candidates = defaultdict(list)
for f in sac_files:
    key = f"{net}.{sta}"
    station_candidates[key].append((f, os.path.getsize(f)))

# Pick largest file (longest continuous segment)
selected[key] = max(station_candidates[key], key=lambda x: x[1])
```

### 3. Skip Warramunga Array Fragments
AU.WB1–WR9 and AU.WC1–WC4 are array elements with ~80 tiny BHZ-only segments each.
These should be excluded from the initial inventory (BHZ-only, no SACPZ needed, not
useful for 3C T-wave analysis at this stage).

## Outputs

- `stations.xlsx`: network, station, lat, lon, dist_km, dist_deg, azimuth, baz, tt_P, tt_S, tp_th_min, ts_th_min, dist_ocean_km, dist_land_km, has_3c, sampling_rate, duration
- `stations.txt`: STATION LON LAT DIST_KM TP_MIN AZIMUTH (GMT-ready, with T-P arrival time for color coding)

## T-P/T-S Computation (Optional)

When `t_phase_analysis` is available, add T-wave travel times:

```python
sys.path.insert(0, '/home/sy/codex/t_phase_analysis')
from t_phase_analysis.traveltime_profile import get_track, tt_t

dist_int = 0.1  # degrees (~11km resolution) — NOT km!
track = get_track(EVLO, EVLA, stlo, stla, dist_int,
                  ref_el=-1000, plot_map=False, plot_elevation=False)
result = tt_t(track, ref_el=-1000, source_depth_km=EVDP, taup_model_name='iasp91')

if result is not None:
    tp_th = result['tt_t']                    # seconds
    ts_th = result['tt_ocean'] + result['tt_land_s']  # seconds
    dist_ocean_km = result['dist_ocean'] * 111.0      # degrees → km
    dist_land_km = result['dist_land'] * 111.0
```

**Critical pitfall**: `get_track()` 5th arg `dist_int` is in **degrees**, not km.
Passing `dist_km` (e.g., 2966) creates a 2-point track → `tt_t()` returns `None`.

See `references/ttt-api-pitfalls-20260707.md` for full API details.

## 3C Detection

Add `has_3c` column by checking for BHE/BH1 files in raw/:

```python
bhe_files = glob.glob(f"{RAW_DIR}/{net}.{sta}..BHE*.SAC") + \
            glob.glob(f"{RAW_DIR}/{net}.{sta}.*.BHE*.SAC") + \
            glob.glob(f"{RAW_DIR}/{net}.{sta}..BH1*.SAC") + \
            glob.glob(f"{RAW_DIR}/{net}.{sta}.*.BH1*.SAC")
has_3c = 'Y' if bhe_files else 'N'
```

## Limitations

- No SNR computation (requires processed waveforms)
- No channel inventory beyond BHZ (only BHZ is scanned for header info)

After this step, the user can decide which stations to preprocess based on
geographic distribution and distance coverage.

# T_Phase_Analysis 官方工具包

## 基本信息

- **来源**：周老师发布版 (Gitee: yzhou89/T_Phase_Analysis)
- **本地路径**：`/home/sy/T_Phase_Analysis/`
- **安装**：`cd /home/sy/T_Phase_Analysis && pip install -e .`

## 模块清单

| 模块 | 功能 |
|:---|:---|
| `traveltime_profile.py` | 走时剖面计算（tt_t, get_track） |
| `snr.py` | 信噪比计算（Welch PSD 频谱能量比） |
| `waveform.py` | 波形处理 |
| `response.py` | 去仪器响应 |
| `pmcc.py` | PMCC 台阵分析 |
| `fdetector.py` | F 统计检测器 |
| `pca.py` | PCA 偏振分析 |
| `denoise.py` | 去噪 |

## 使用注意

这是正式发布版本（带 pyproject.toml），建议 `pip install -e .` 安装后使用。之前手工 copy 的 `traveltime_profile.py` 应统一用此版本的接口。

## 克隆时间

2026-07-09

# T_Phase_Analysis Toolkit Integration (2026-07-16)

## 库路径
- codex 版: `/home/sy/codex/t_phase_analysis`
- gitee 版: `/mnt/f/sharedate/T_Phase_Analysis`

## 模块清单

| 模块 | 用途 | 来源 |
|:---|:---|:---|
| traveltime_profile.py | 走时计算 (tt_t) | codex 版（保持度为单位） |
| waveform.py | 时频谱计算+绘图 | 两版相同 |
| snr.py | 信噪比计算 | 两版相同 |
| pca.py | PCA 偏振时间演化 | gitee 版新增 |
| response.py | MiniSEED 去响应 | gitee 版新增 |
| _validators.py | 输入校验 | gitee 版新增 |
| pmcc.py | 台阵互相关 | 两版相同（暂不接入） |

## API 差异（重要）

codex 版 `tt_t()` 返回：
- `dist_ocean`: 海洋距离（度）
- `dist_land`: 陆地距离（度）

gitee 版 `tt_t()` 返回：
- `dist_ocean_km`: 海洋距离（km）
- `dist_land_km`: 陆地距离（km）

**所有现有脚本依赖 codex 版约定（度），不可替换。**

## 预处理双路径

| 输入格式 | 工具 |
|:---|:---|
| SAC + SACPZ（IRIS TA 台网） | process_obs.py |
| MiniSEED + StationXML | t_phase_analysis.response.remove_response_def() |

## Python 环境

```bash
/home/sy/miniforge3/envs/obspy/bin/python  # 统一环境（pygmt 0.19, obspy 1.5, 全依赖）
```

# Vision Analysis via Subagent (June 2026)

## Use Case
When you cannot view generated figures (headless agent, no GUI), delegate image analysis to a vision subagent.

## Method
```python
delegate_task(
    goal="Analyze figure...",
    context="Figure paths and expected content...",
    toolsets=['vision', 'file', 'terminal']  # MUST include terminal for Python/PIL
)
```

## Key Requirement
The `terminal` toolset is REQUIRED for vision analysis. The subagent uses Python/PIL to read pixel data. Without `terminal`, the subagent cannot run image processing code.

## Limitations
- Not true OCR — works at pixel-level analysis
- Can measure: blank space %, label visibility, color contrast, element overlap
- Cannot assess: aesthetic quality, scientific correctness of labels
- May not work for very large images or complex layouts

# 视觉模型调用规则（2026-07-15 定稿）

## 核心规则

**只在用户明确说要"看看"、"检查"、"确认"图时才调用视觉模型。**

日常画图生成后直接告知文件路径，不主动调用视觉验证。

## 原因

- 视觉模型调用耗时长（20-60秒/张）
- 用户反馈"太慢了"
- 用户可以自己查看图片

## 调用方式（当用户明确要求时）

```bash
bash ~/.hermes/skills/mimo-omni/mimo_api.sh image /path/to/image.png "问题描述" --max-tokens 16384
```

## 不调用的情况

- 生成新图后
- 批量处理后
- 修改参数重新生成后
- 任何用户没有明确要求检查的情况

# Word Report Generation for T-wave Analysis Results

Use `python-docx` to produce `.docx` reports with embedded figures for method documentation and result communication.

## Setup

```bash
# python-docx is usually in the tphase_env
source ~/tphase_env/bin/activate
python3 -c "import docx; print('OK')"
pip install python-docx  # if not installed
```

## Typical Structure

### Method Report (`method_report.docx`)
Covers: data source → processing flow → scripts → parameters → quality control.

### Result Report (`result_report.docx`)
Covers: results → figures → discussion → conclusions.

## Figures in docx

```python
from docx.shared import Inches
img_path = "/path/to/figure.png"
if os.path.exists(img_path):
    doc.add_picture(img_path, width=Inches(5.2))
    last_paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
```

## Chinese Fonts in docx

Set both Latin and East-Asian fonts:

```python
from docx.oxml.ns import qn

style = doc.styles['Normal']
style.font.name = 'Times New Roman'
style.element.rPr.rFonts.set(qn('w:eastAsia'), '宋体')
```

## Styling

```python
from docx.shared import Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH

def add_heading_styled(text, level=1):
    h = doc.add_heading(text, level=level)
    for run in h.runs:
        run.font.color.rgb = RGBColor(0, 51, 102)
    return h
```

## Output Location

- Save to `oi/final_figures/` alongside figures and scripts
- Script names: `write_method_report.py`, `write_result_report.py`

## Pitfalls

1. **Quotes in Chinese strings** — When embedding Chinese text like `这一"Fall Line"界面`, use Unicode escapes (`\u201c`/`\u201d`) or single-quote the outer string to avoid syntax errors from nested double quotes inside doc-strings.
2. **Large files with embedded images** — A docx with three 1+ MB PNG images can be ~10 MB. This is normal for scientific reports.
3. **Table column alignment** — Use `table.style = 'Light Grid Accent 1'` for publication-style tables.

# User Interaction Preferences (HARD RULES)

**These are explicit user corrections — violating them causes frustration.**

## 1. Just Answer Questions, Don't Execute

When user asks a question, answer directly. Don't auto-execute code or start processing unless explicitly asked.

**User's exact words**: "以后我问问题只需要回答就行，不要自己去多余执行"

**Examples**:
- User: "这个放在哪里了？" → Just answer the location, don't run commands
- User: "有没有Linux版？" → Just answer yes/no with alternatives
- User: "端口在哪里看？" → Just explain where to find it

**Violation pattern**: Running terminal commands, checking file sizes, testing scripts when user only asked a question.

## 2. Don't Verify Figures with Vision Model

User will check figures themselves. Only use vision_analyze if user explicitly says "看看图" or "check this image".

**User's exact words**: "以后画图不需要用视觉模型验证分析，我去看就行，除非我指定你看看图"

**Violation pattern**: Running bash ~/.hermes/skills/mimo-omni/mimo_api.sh or vision_analyze after generating figures.

## 3. Progress Reports for Long Tasks

Report progress every 5 minutes during long-running tasks to prevent terminal hibernation.

**User's exact words**: "按时向我汇报进度防止休眠终端"

**Implementation**:
```bash
sleep 300 && echo "=== 5分钟后进度 ===" && tail -20 /tmp/task.log
```

## 4. Clean Up Temp Files

Remove temporary scripts and images after task completion.

## 5. Don't Explain Unless Asked

When user says "只回答就行" or asks a simple question, give a concise answer without executing commands or providing extra context.

**Good**: "代理通了。"
**Bad**: Running additional tests, checking file sizes, providing implementation details.

# Azimuth-Filtered Station Identification from SAC Headers

When you need to find stations within a specific azimuth range from an event,
the fastest approach is reading SAC headers directly — no preprocessing needed.

## Workflow

### 1. Read station coordinates from BHZ SAC headers

```python
from obspy import read

tr = read(f"TA.{sta}..BHZ.*.SAC", headonly=True)[0]
stla = tr.stats.sac.stla
stlo = tr.stats.sac.stlo
```

`headonly=True` is critical — it reads only the header block, not the waveform
data. For 250+ stations this takes <2s instead of minutes.

### 2. Compute azimuth from event to station

⚠️ **Use GeographicLib for BAZ consistency.** The user's polarization and
rotation scripts use GeographicLib for backazimuth. `obspy.geodetics.gps2dist_azimuth`
can differ by ~5.7°. Keep BAZ computation uniform:

```python
from geographiclib.geodesic import Geodesic

evla, evlo = 23.7240, -45.5809  # event coordinates
g = Geodesic.WGS84.Inverse(evla, evlo, stla, stlo)
az = g['azi1']
az = az + 360 if az < 0 else az
baz = (az + 180) % 360
dist_km = g['s12'] / 1000.0
```

`az` = azimuth from event to station (0° = North, clockwise).
`baz` = backazimuth from station to event.

If you prefer a convenience function, use `lib.arrivals.dist_az()` which wraps
the GeographicLib call:
```python
from lib.arrivals import dist_az
dist_km, az, baz = dist_az(evlo, evla, stlo, stla)
```
(Note the argument order: `dist_az(evlo, evla, stlo, stla)` — lon before lat.)

### 3. Filter by target azimuth range

```python
if 315 <= az <= 317:
    # keep this station
```

### 4. Complete example

See `/tmp/find_azimuth_stations.py` (created June 2026 session) or the pattern
below:

```python
import os, glob
from geographiclib.geodesic import Geodesic
from obspy import read

EVENT_DIR = "/mnt/f/sharedate/case/2014-07-27-mw61"
EVLA, EVLO = 23.7240, -45.5809

stations_az = []
for f in glob.glob(os.path.join(EVENT_DIR, "*.BHZ*.SAC")):
    tr = read(f, headonly=True)[0]
    sta, net = tr.stats.station, tr.stats.network
    stla, stlo = tr.stats.sac.stla, tr.stats.sac.stlo
    g = Geodesic.WGS84.Inverse(EVLA, EVLO, stla, stlo)
    az = g['azi1']; az = az + 360 if az < 0 else az
    stations_az.append({'net': net, 'sta': sta, 'az': az,
                        'dist_km': g['s12']/1000, 'lat': stla, 'lon': stlo})

# Filter
target = [s for s in stations_az if 315 <= s['az'] <= 317]
```

## Pitfalls

1. **Multiple networks**: The event directory may contain TA, IU, US, LD, and
   other networks side-by-side. Always glob `*BHZ*` instead of assuming
   `TA.*BHZ*` to discover all stations.

2. **Multiple location codes**: IU.HRV has both `00` and `10` location codes.
   Both will match `*.BHZ*.SAC`. Handle via `tr.stats.location` if needed.

3. **SAC files in root vs `raw/`**: Some events keep SAC in the root dir
   (2014-07-27-mw61), others use `raw/`. Check with `ls` first.

4. **headonly=True errors**: Some SAC files have corrupted headers. Wrap in
   try/except and skip unreadable files.

5. **ArcDirect vs Direct pitfall (azimuth arcs)**: `Geodesic.WGS84.ArcDirect()`
   takes arc angle in **degrees**, NOT meters. When drawing azimuth arcs on maps
   at a given distance, use `Geodesic.WGS84.Direct()` with meters. See
   `references/geographiclib-arcdirect-pitfall.md` for full recipe and example.

## After Identification: Creating Azimuth-Targeted Scripts

Once you have the station list for a specific azimuth range, the next step is
creating **targeted batch scripts** rather than processing all stations.

### Pattern: Two Scripts Per Azimuth Group

For each azimuth range, create two scripts in the event root directory:

| # | Script Name Template | Purpose | Dependency |
|---|---------------------|---------|-----------|
| 1 | `process_{RANGE}.py` | Batch preprocessing (instrument response + BPF) | obspy only — self-contained SACPZ parsing |
| 2 | `02_extract_stations_{RANGE}.py` | Station info → Excel + GMT file | `codex/lib/arrivals` (for T-P/T-S via tt_t) |

### Script 1: Batch Preprocessor (Self-Contained)

```python
# process_315_317.py structure:
TARGET_STATIONS = ["L61A", "K59A", ...]  # hardcoded list
NETWORK = "TA"
CASE_DIR = "/mnt/f/sharedate/case/2014-07-27-mw61"

for sta in TARGET_STATIONS:
    for comp in ['BHZ', 'BHN', 'BHE']:
        tr = read(raw_file)[0]
        tr_resp = remove_response(tr, sacpz_file)    # simulate_seismometer
        tr_bpf = bandpass_filter(tr_resp)             # 2-8 Hz
        SACTrace.from_obspy_trace(tr_bpf).write(proc_path)
```

Key features:
- Reimplements `parse_sacpz()` directly (no external library needed)
- Hardcodes event coordinates + station list
- Can run under any Python env with obspy

### Script 2: Station Info Extraction (Uses Codex Library)

```python
# 02_extract_stations_315.py structure:
import sys; sys.path.insert(0, "/home/sy/shared/codex")
from lib.arrivals import compute_theoretical_arrivals, dist_az

TARGET_STATIONS = [...]  # same list
for fn in bhz_raw_files:
    sta = read(fn).stats.station
    if sta not in TARGET_STATIONS: continue
    arrs = compute_theoretical_arrivals(evlo, evla, evdp, stlo, stla)
    # → P, S, T-P, T-S, land_km, ocean_km
```

Key features:
- Reads only raw SAC headers (no processed data needed for basic station info)
- Uses `compute_theoretical_arrivals()` for T-P/T-S with land/ocean separation
- Writes to Excel (Chinese headers + English abbreviation + SNR≥2 green highlight)
- **Must run from `/home/sy/shared/codex/`** or add to sys.path

### When to Create Azimuth-Targeted Scripts Instead of Codex Pipeline

| Scenario | Chosen Approach |
|----------|----------------|
| Full event analysis needed (all stations, all steps) | Codex pipeline (`01_preprocess.py --all` → `02_extract_stations.py`) |
| First look at a new event, only one azimuth group | Azimuth-targeted scripts (faster, lighter) |
| User says "only process this azimuth range" | Azimuth-targeted scripts |
| Multiple azimuth groups (327° + 350°) | One script pair per group, or codex pipeline + `--az-center` filter |

### Boundary Station Handling

When a station falls just outside the requested azimuth range (e.g. K62A at
317.06° for a 315-317° request), include it in the output with a note so the
user can decide:

```python
# After filtering:
print(f"Boundary stations: K62A (az=317.06°), J60A (az=317.13°)")
```

### Script Naming Convention

Use the azimuth range in the filename for clarity:

```bash
process_315_317.py         # Batch preprocess 315-317° stations
02_extract_stations_315.py # Station info for 315-317° stations
```

The range string should be `{MIN}_{MAX}` (no decimal) for compact filenames.

## Application

This workflow is the first step before any azimuth-filtered analysis. Once
stations are identified and targeted scripts are created, proceed with:
- Run preprocessing script → creates sacprocess/
- Run extraction script → creates stations.xlsx + gmt_stations.txt
- Record section plots with `06_profile_section.py --az-center ...`
- Station arrival maps (pyGMT)
- Polarization analysis (single-station)

# 方位角依赖的海洋路径长度 — 巴哈马效应

## 发现

2014-01-13 Puerto Rico Mw 6.3 事件中，AZ 341–360° 台站的 T-P 到时不是均匀分布——低方位角（341–347°）台站的 T-P 比高方位角（350–354°）台站系统性地早 70–120s。

## 原因

不同方位角的大圆路径穿过的海域不同：

| 方位角 | 路径海域 | 海洋距离 | T-P 范围 |
|--------|----------|----------|----------|
| 341–347° | 巴哈马浅滩（Bahamas Banks） | 2165–2275 km | 1529–1580s |
| 350–354° | 佛罗里达外海深水区 | 2325–2343 km | 1619–1638s |

巴哈马浅滩的大陆架极宽，-1000m 等深线离岸很远，T 波更早"登陆"转换，
海洋段比佛罗里达路径短约 180 km。

## 在剖面图上的表现

按陆地距离排序时，低方位角台站（如 M58A, land=477km, T-P=1529s）
的 T-P 圆点出现在邻近高方位角台站（如 I63A, land=473km, T-P=1638s）
之前 109s。看起来像是"计算错误"，实为真实海底地形效应。

## tt_t() 与 T-Love 包络峰的系统性偏移

`tt_t()` 计算的是 SOFAR 通道传播到 -1000m 等深线的理论走时（T-P 到时）。
而包络自动检测的峰值对应的是陆地段 Love 波到达（T-Love 峰值），
两者之间的差值（~130s 中位数）是正常的陆地段传播时间，不是 bug。

## 标注建议

在波形剖面图上用红色标签标注 az<347° 的台站，并在标题或副标题注明原因。

## 验证数据 (2014-01-13, AZ>341°, SNR≥3)

```
STA      Land(km)  Ocean(km)  T-P(s)   Azimuth   Note
M58A     477       2165       1529     341.7°    巴哈马路径
M59A     415       2212       1554     343.7°    巴哈马路径
L59A     445       2238       1575     344.9°    巴哈马路径
PAL      267       2262       1569     346.0°    巴哈马路径
M61A     283       2275       1580     346.5°    巴哈马路径
I63A     473       2325       1638     353.6°    佛罗里达路径
J62A     384       2341       1637     351.2°    佛罗里达路径
K62A     328       2343       1632     350.1°    佛罗里达路径
L63A     231       2343       1619     350.8°    佛罗里达路径
```

# GMT Map Iterative Refinement Pattern

When producing publication-style GMT station maps for seismic T-wave analysis, expect iterative style corrections. This reference documents the refinement sequence developed during the 2014-01-13 Puerto Rico Mw 6.3 (327° profile) analysis.

## The Iteration Path

### v0 → v1: First draft corrections
| Issue | Correction |
|-------|-----------|
| Colormap too heavy (terrain obscures stations) | Increase transparency `-t40`→`-t60` or switch to contour lines |
| Colorbar inside map area | Move outside: `-DjRM+w8c/0.4c+o1.0c/0c+ml` |
| No azimuth label | Add `327° Profile` text at the profile endpoint |
| No location name | Add "Puerto Rico" next to epicenter |
| Station names cluttering map | Remove unless explicitly requested |

### v1 → v2: User-requested changes
| Issue | Correction |
|-------|-----------|
| Terrain still too light | Switch from custom light CPT to GMT's built-in `-Cglobe` |
| ln(amplitude) is negative on colorbar | Switch to **relative amplitude %**: `rel_pct = amp / max_amp * 100` |
| Station labels still present | Already removed in v1 ✅ |
| Still no 327° label | Add goldenrod `327° Profile` |
| Coastline too bold | Reduce from `1.2p,black` to `0.6p,gray30` |

### v2 → v3: Contour-line experiment (REJECTED — user wanted colorfill back)
| Issue | Correction |
|-------|-----------|
| Terrain resolution too coarse (1 arcmin) | Upgrade to **30 arcsec** (`@earth_relief_30s`) |
| Colorfill terrain obscures stations | Switch to **contour lines with numeric labels** (`grdcontour -C500 -A1000`) |
| Colorbar offset insufficient | Increase offset to `o2.5c` |
| Coastline still slightly heavy | Final: `0.5p,gray30` |

**Result**: User explicitly rejected contour lines and said "还是用地形图吧" (use colorfill instead).

### v3 → v4: Final accepted version
| Issue | Correction |
|-------|-----------|
| Terrain contour lines (rejected) | Back to **color-fill topography** (`gmt grdimage -Cglobe -I+d -t40`) |
| Relative amplitude on colorbar (rejected) | Switch to **actual amplitude ×10⁸** with full range (0~11) |
| Colorbar capped at 5.5 (rejected) | Restore to **full range**, no truncation |
| Colorbar annotation interval too dense | Change to `-Bxa2` (every 2 units) from `-Bxa1` |
| (a) panel label present | Remove panel label for standalone figure |
| No figure-wide title | Add GMT title: `-B+t"..." --FONT_TITLE=16p --MAP_TITLE_OFFSET=1.2c` |
| Anomalous station (V56A) questioned | Verify with multi-band, SNR, and spectral checks; data quality confirmed |

## Reference Template (Final v4 GMT Script Pattern)

See the actual script at `F:\\sharedate\\chubu\\05_summary\\01_gmt_map_327_v4.sh` for the final working version. Key pattern:

```bash
gmt begin figure_map png,eps
  # Title above map
  gmt psbasemap -J${PROJ} -R${R_MAP} \
    -B+t"T-wave Amplitude Decay ..." \
    --FONT_TITLE=16p,Helvetica-Bold,black --MAP_TITLE_OFFSET=1.2c

  # Color-fill topography (30 arcsec)
  gmt grdimage ${TOPO_GRID} -J${PROJ} -R${R_MAP} -Cglobe -I+d -t40 -Q

  # Thin coastline
  gmt pscoast -J${PROJ} -R${R_MAP} -W0.5p,gray30 -Df -A500

  # Distance arcs, great circle path
  gmt psxy arc.txt -J${PROJ} -R${R_MAP} -W1.0p,darkgreen,dashed

  # Stations colored by actual amplitude ×1e8 (full range)
  awk '{print $1, $2, $3}' stations.txt | \
    gmt psxy -J${PROJ} -R${R_MAP} -Sc0.30c -C${CPT} -W0.3p,black

  # Colorbar OUTSIDE map with full-range annotation
  gmt psscale -C${CPT} -J${PROJ} -R${R_MAP} \
    -DjRM+w8c/0.4c+o2.5c/0c+ml -Bxa2+l"T-wave Amplitude (x10@+-8@+)"

  # Epicenter, scale bar, frame (NO panel label)
  echo "${EVLO} ${EVLA}" | gmt psxy -Sa0.6c -Gred -W0.5p,black
  gmt psbasemap -LjBL+c34N+f+w500k+l"500 km"
  gmt psbasemap -Bxa10f5 -Bya5f2.5 -BWSNE
gmt end
```

## Available GMT Earth Relief Resolutions

| Dataset | Resolution | Status | Notes |
|---------|-----------|--------|-------|
| `@earth_relief_01m` | 1 arcmin (~1.85 km) | Always available | Too coarse for publication maps |
| `@earth_relief_30s` | 30 arcsec (~900 m) | Always available | **Default minimum resolution** |
| `@earth_relief_15s` | 15 arcsec (~450 m) | Auto-downloads if cached | Better but larger |
| `@earth_relief_08s` | 8 arcsec (~250 m) | **DOES NOT EXIST** as global grid | Do not attempt to use |
| `@earth_relief_03s` | 3 arcsec (~90 m) | Partial coverage only | Land-only, unpredictable |

## Key Lessons

1. **User rejected both contour lines AND capped colorbar.** Default to colorfill + full-range amplitude scale.
2. **When the user asks about the colorbar**, describe in concrete numbers ("色标显示 T波振幅 × 10⁸，范围 0.5~10.5").
3. **Single-panel standalone map → NO (a)/(b) label → use `-B+t"..."` title instead.**
4. **Anomalous stations should be verified before reporting** — run multi-band/SNR/spectral checks to distinguish geological signal from data error.
