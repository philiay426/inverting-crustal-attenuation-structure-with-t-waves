# 衰减分析笔记

# T-Wave Attenuation Feasibility Figure (4-Panel Composite)

## When to Use

When the user wants a comprehensive amplitude-feasibility figure combining:
- T-wave zoom record section
- Amplitude decay plot (ln envelope peak vs land distance)
- Frequency-dependent amplitude decay (multi-band ln RMS vs land distance)

Example: `figure_350_amplitude_feasibility.py` at `/mnt/f/sharedate/case/`

## Panel Layout

| Panel | Content | Scaling |
|-------|---------|---------|
| (a) | Record section — T-wave zoom window | Global max norm (preserves relative amplitudes) |
| (b) | Record section — same data | Global P98 norm (weak signals become visible) |
| (c) | ln(envelope peak) vs land distance | T-Love ±45s window, color=SNR, station labels |
| (d) | ln(RMS spectral amplitude) vs land distance | 5 sub-bands: 2-3, 3-4, 4-5, 5-6, 6-8 Hz, linear fits |

## Normalization Strategies

Three normalization modes exist, each with a distinct purpose:

| Mode | Formula | When to Use |
|------|---------|-------------|
| **Global max** | `data / max(|all_data|)` | Default for comparing relative amplitudes across stations. Strongest station dominates. |
| **Global P98** | `data / percentile(|all_data|, 98)` | Weak-signal visualization. Suppresses outlier spikes. ~4× larger amplitude than max-norm for typical T-wave data. |
| **Per-trace in window** | `data / max(|data_i|)` per station | Coherence checks. Each station's waveform fills the same vertical space regardless of amplitude. Hides attenuation. |

**User preference**: Global P98 for weak-signal zoom, global max for relative amplitude, per-trace for coherence.

## T-Love Theoretical Arrival

```python
T-Love = tt_ocean + land_distance_km / v_love
```
- `tt_ocean`: ocean path travel time from `tt_t()` (seconds)
- `land_distance_km`: `degrees2kilometers(tr['dist_land'])`
- `v_love`: Love-wave velocity on land, default 3.6 km/s (range 3.4–3.8 km/s)

Distinct from T-P: T-P uses TauP iasp91 P-wave velocity on land (~5.8 km/s), T-Love uses crustal Love-wave velocity (3.4-3.8 km/s). At land distances >400 km, T-Love is 10-30s later than T-P.

## Frequency Bands for Panel (d)

Five sub-bands covering 2-8 Hz with `bandpass(data, fmin, fmax, 1/delta, corners=4, zerophase=True)`:
- 2–3, 3–4, 4–5, 5–6, 6–8 Hz

For each band, compute RMS amplitude in T-Love ±45s window, take ln, fit linear regression vs land distance. Higher frequencies generally show steeper attenuation (more negative b).

## Key Pitfalls

1. **`tt_t().dist_land` is in degrees, not km** — Always convert with `degrees2kilometers()`.
2. **T-Love vs T-P separation**: At land distances <200 km, T-Love and T-P are within 1-2s (indistinguishable on zoom plot). At >500 km, they separate by 15-30s.
3. **Velocity for T-Love**: Use 3.6 km/s as default. Test sensitivity by varying 3.4-3.8 km/s.
4. **Multi-panel spacing**: Use `gridspec.GridSpec(2, 2, hspace=0.30, wspace=0.32)` with sufficient margins.
5. **Station labels in panel (c)**: Use `annotate()` with offset points — don't overlap at dense distances. Draw order: scatter first, then annotate.

# Seismic T-Wave Effective Attenuation Analysis (Zhou-style Method)

Based on Zhou et al. (2021) crustal Qp/Qs inversion methodology, adapted for US East Coast
Love-type seismic T-wave effective attenuation.

## Key Concept

Do NOT assume a single conversion point for the entire US East Coast.
Instead, group stations by azimuth and assign each group an approximate conversion range:

| Profile | Azimuth | Conversion Range | Dist Range | V_land |
|---------|---------|-----------------|------------|--------|
| 350°    | 349-351 | ~2500 km        | 2500-3200 | 3.7 km/s |
| 327°    | 326-328 | ~2000 km        | 2000-2800 | 3.4 km/s |

Reference: De Groot-Hedlin (2020), Figures 3 and 4.

## Methodology (6 Steps)

### Step 0: Directory Structure

```
chubu/
├── 01_az350_selection/        # Azimuth-filtered stations + profile plot
│   ├── stations_az350.txt
│   └── profile_az350.png
├── 02_amplitude_extraction/   # RMS + envelope peak amplitudes
│   └── amplitudes_az350.csv
├── 03_attenuation_fit/        # ln(A) ~ inland_distance regression
│   ├── attenuation_az350.png
│   └── fit_results.txt
└── 04_q_estimation/           # Equivalent Q
    └── q_estimation.md
```

### Step 1: Azimuth Computation and Station Screening

```python
import numpy as np

def calc_azimuth(evla, evlo, stla, stlo):
    """Forward azimuth from source to station (degrees clockwise from north)"""
    dlon = np.radians(stlo - evlo)
    lat1, lat2 = np.radians(evla), np.radians(stla)
    x = np.sin(dlon) * np.cos(lat2)
    y = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)
    az = np.degrees(np.arctan2(x, y))
    return (az + 360) % 360
```

Read each SAC file from `02_filtered_trimmed/`, compute azimuth and distance,
filter by profile parameters, then save to `stations_<profile>.txt`.

### Step 2: Draw Record Section (Profile Plot)

Same as `02_plot_profile.py` in codex — waveforms ordered by distance,
theoretical travel-time lines for 4.0, 3.0, 1.5 km/s using hybrid model:

```python
def travel_time(d, conv_range, v_land, v_ocean=1.5):
    if d <= conv_range:
        return d / v_ocean / 60.0
    else:
        return (conv_range / v_ocean + (d - conv_range) / v_land) / 60.0
```

### Step 3: Compute Inland Distance

```python
inland_distance = station_dist_km - conversion_range
# Only include stations where inland_distance >= 0
```

### Step 4: Extract T-Wave RMS Amplitude

Recommended bandpass: **1-5 Hz** (per De Groot-Hedlin 2020).
Theory arrival: `t_theory_sec = conv_range/1.5 + inland_distance/v_land`
Window: ±120 seconds around theory arrival.

```python
t_theory_sec = conv_range / 1.5 + inland_distance / v_land
idx = np.searchsorted(times, t_theory_sec - 120)
idx2 = np.searchsorted(times, t_theory_sec + 120)
window = data[idx:idx2]
A_rms = np.sqrt(np.mean(window**2))
```

Also compute envelope peak as alternative:
```python
env = envelope(data)
A_peak = np.max(env[idx:idx2])
```

Save to `amplitudes_<profile>.csv`: station, dist_km, inland_km, A_rms, A_peak

### Step 5: Fit ln(A) ~ Inland Distance

```python
from scipy.stats import linregress

mask = amplitudes > 0
x = inland_distances[mask]
y = np.log(amplitudes[mask])
result = linregress(x, y)
b = -result.slope          # effective attenuation coefficient (1/km)
R2 = result.rvalue**2
```

Two versions:
1. **Uncorrected**: `ln(A) = a - b * x`
2. **Geometric spreading corrected**: `A_corr = A * sqrt(x + 1)`, then `ln(A_corr) = a - b * x`

Output two values for comparison.

### Step 6 (Optional): Equivalent Q

```python
Q_eff = np.pi * f_center / (b * v_land)
# f_center = 3.0 Hz for 1-5 Hz band
# v_land = 3.7 km/s for 350° profile
```

⚠️ Label as **"effective Q"** or **"Love-type T-wave equivalent Q"** — do NOT call it Qp or Qs.

### Pitfalls

1. **Do NOT mix azimuths** — each profile (350°, 327°) must be processed independently.
2. **Do NOT call it Qp/Qs** — US East Coast T-waves are Love-type surface waves per De Groot-Hedlin 2020, not body-wave P/S conversions.
3. **Conversion range sensitivity test** — try ±50 km variations to check stability.
4. **Relative amplitude** — if a near-coast reference station exists, use `ln(A_i / A_ref)` to reduce source and conversion efficiency effects.
5. **Amplitude window width** — ±120 s is empirical; adjust based on T-wave duration (shorter for northern group, longer for southern group per paper).
6. **SNR filtering is harmful for de-meaned velocity data** — instrument-response-removed SAC data typically has amplitude ~1e-9 (m/s), with P/S wave energy filling the pre-T-wave noise window (600-1800s), making noise RMS comparable to T-wave RMS. Standard SNR thresholds (≥1.5) reject ALL stations — skip SNR-based rejection and rely on windowed envelope peak extraction around the theoretical arrival. What matters is relative amplitude between stations.
7. **Theoretical arrival window** — the hybrid-path model (1.5 km/s ocean + 3.7 km/s land) only approximates the conversion. The observed envelope peak is typically **~170-300 s earlier** than the hybrid-path prediction for the 327° profile with v_land=3.4 — check timing diagnostics by plotting full waveforms with theoretical lines before extracting amplitudes.
8. **Use envelope peak from a 180s search window, not a fixed RMS window** — The T-wave is a spindle-shaped wave packet (~3-4 min duration per De Groot-Hedlin 2020), not a continuous signal. Search ±180s around theoretical arrival for the envelope maximum, then take RMS in a narrower (±60s) window centered on the peak. This avoids dilution by non-T-wave energy.
9. **Verify celerity** — T-wave criterion: celerity (dist / peak_time) < 1.8 km/s. If celerity > 2.0 km/s, the peak is likely S-wave coda, not T-wave. Check this before trusting the extraction.

## Expected Results

For the 2014-01-13 Puerto Rico event:
- 350° profile: ~40-60 stations (46 found in the wiggle-fill run)
- 2504-3178 km range → inland distance 0-678 km
- Expected b (uncorrected): ~0.001-0.005 1/km (depends on frequency band)
- R² should be > 0.5 for a clear attenuation signal
