# Experiment Results

All numbers below were produced by running the implementation — none are
hand-edited. Regenerate everything with:

```bash
python scripts/generate_demo_results.py
python scripts/run_sd_case_study.py
python scripts/verify_sd_equivalence.py --networks 100
python scripts/benchmark_sd_magnitude.py --runs 50
python scripts/benchmark_sd_direction.py --runs 200
python scripts/benchmark_sd_antenna_height.py --runs 100
```

Environment: Python 3.13, macOS (arm64); package version 0.1.0; 49 unit
tests passing at the commit that produced these results.

Default network for all experiments: 8 stations, 16 baselines, station S00
fixed, noise a = 5 mm, b = 1 ppm, Z-axis scale 1.8, random 3×3 correlations
up to |ρ| ≤ 0.4. Default DIA settings: α_global = 0.05, α_local = 0.001.
Robust defaults: Huber c = 1.5; Hampel (a, b, c) = (1.5, 3.5, 8.0);
a-priori scale.

---

## 1. Clean Gaussian data (seed 42)

Command (equivalent CLI): `python -m gnss_adjust demo --seed 42 --output results/clean`

| method | RMSE [m] | mean pos err [m] | max pos err [m] | flagged | iters | variance factor |
|---|---|---|---|---|---|---|
| WLS | 0.004085 | 0.006653 | 0.011660 | 0 | 1 | 0.757 |
| DIA-reject | 0.004085 | 0.006653 | 0.011660 | 0 | 1 | 0.757 |
| DIA-inflate | 0.004085 | 0.006653 | 0.011660 | 0 | 1 | 0.757 |
| Huber | 0.004085 | 0.006653 | 0.011660 | 0 | 2 | 0.757 |
| Hampel | 0.004085 | 0.006653 | 0.011660 | 0 | 2 | 0.757 |

**Observation — clean-data efficiency.** On clean data every method returns
the identical solution: the global test passes immediately (DIA does
nothing) and no baseline exceeds the robust thresholds (all robust weights
stay at 1). There is no efficiency loss at this noise level.

## 2. Single directional 3D outlier (seed 42)

Bias **(+0.035, −0.025, +0.090) m** injected into baseline 3 (roughly a
10σ, mostly-vertical error).

Command: `python -m gnss_adjust demo --seed 42 --method all --outlier-baseline 3 --outlier-vector 0.035 -0.025 0.090 --output results/directional_outlier`

| method | RMSE [m] | max pos err [m] | precision | recall | F1 | flagged | iters | variance factor |
|---|---|---|---|---|---|---|---|---|
| WLS | 0.008765 | 0.022697 | 0 | 0 | 0 | 0 | 1 | 2.801 |
| DIA-reject | 0.004431 | 0.011212 | 1.0 | 1.0 | 1.0 | 1 | 2 | 0.785 |
| DIA-inflate | 0.004421 | 0.011084 | 1.0 | 1.0 | 1.0 | 1 | 2 | 0.734 |
| Huber | 0.005468 | 0.015519 | 1.0 | 1.0 | 1.0 | 1 | 9 | 1.656 |
| Hampel | 0.004855 | 0.013989 | 1.0 | 1.0 | 1.0 | 1 | 12 | 1.338 |

**Observations.** WLS coordinate error doubles (RMSE 4.1 → 8.8 mm) and its
variance factor jumps to 2.8 — the outlier is smeared over the network. DIA
identifies exactly baseline 3 in one adaptation step and recovers
clean-level accuracy. Huber bounds but does not remove the outlier (final
weight 0.32 at t = 4.6), landing between WLS and DIA. Hampel's descending
segment pushes the weight further down (0.20 at t = 5.0) and lands close to
DIA; the outlier is large enough to be flagged (weight < 0.5) but not large
enough (t < 8) for Hampel's hard rejection.

## 3. Monte Carlo, single 30σ outlier (100 runs, seeds 1000–1099)

Command: `python -m gnss_adjust benchmark --runs 100 --scenario single-outlier --seed 1000 --output results/benchmark_single_outlier`

| method | RMSE mean [m] | RMSE median [m] | F1 | precision | recall | FPR | localization | iters | runtime [s] |
|---|---|---|---|---|---|---|---|---|---|
| WLS | 0.02824 | 0.02435 | 0.00 | 0.00 | 0.00 | 0.000 | 0.00 | 1.0 | 0.0007 |
| DIA-reject | 0.00943 | 0.00657 | 0.96 | 0.96 | 0.96 | 0.003 | 0.95 | 2.0 | 0.0036 |
| DIA-inflate | 0.00944 | 0.00670 | 0.95 | 0.95 | 0.96 | 0.003 | 0.94 | 2.0 | 0.0037 |
| Huber | 0.01035 | 0.00748 | 0.97 | 0.97 | 0.97 | 0.002 | 0.96 | 10.6 | 0.0071 |
| Hampel | 0.00917 | 0.00658 | 0.97 | 0.97 | 0.97 | 0.002 | 0.96 | 4.9 | 0.0038 |

**Observations.** A 30σ outlier triples WLS RMSE. All four countermeasures
recover ~3× accuracy with F1 ≈ 0.95–0.97 and false-positive rates near the
nominal test size. Hampel is both the most accurate and the fastest robust
method here (full rejection converges quickly). DIA-inflate's "converged"
rate (global test eventually passing) was 72%: after inflating the true
outlier, the ×100-inflated block sometimes still leaves the global statistic
marginally above the critical value, ending the loop with status
`no_candidate` — the coordinates are nevertheless fine, which is why its
RMSE matches DIA-reject (the 28 non-"converged" runs in fact average
*better* RMSE, 0.0071 m, than the converged ones, 0.0104 m).

## 4. Outlier-magnitude sweep (30 runs/point, seeds 5000–5029)

Command: `python -m gnss_adjust benchmark --scenario magnitude-sweep --runs 240 --seed 5000 --output results/magnitude_sweep`

Coordinate RMSE (m, mean over 30 runs) vs outlier magnitude (σ units):

| magnitude | WLS | DIA-reject | DIA-inflate | Huber | Hampel |
|---|---|---|---|---|---|
| 0 | 0.0080 | 0.0080 | 0.0080 | 0.0080 | 0.0080 |
| 3 | 0.0077 | 0.0077 | 0.0077 | 0.0076 | 0.0076 |
| 5 | 0.0086 | 0.0078 | 0.0078 | 0.0080 | 0.0080 |
| 8 | 0.0102 | 0.0074 | 0.0074 | 0.0082 | 0.0079 |
| 12 | 0.0129 | 0.0074 | 0.0074 | 0.0082 | 0.0075 |
| 20 | 0.0188 | 0.0074 | 0.0074 | 0.0083 | 0.0074 |
| 30 | 0.0268 | 0.0074 | 0.0074 | 0.0083 | 0.0074 |
| 50 | 0.0431 | 0.0074 | 0.0074 | 0.0083 | 0.0118 |

Detection F1 (mean):

| magnitude | DIA-reject | DIA-inflate | Huber | Hampel |
|---|---|---|---|---|
| 3 | 0.13 | 0.13 | 0.03 | 0.03 |
| 5 | 0.63 | 0.63 | 0.37 | 0.37 |
| 8 | 1.00 | 1.00 | 0.83 | 0.83 |
| 12 | 1.00 | 1.00 | 1.00 | 1.00 |
| 20+ | 1.00 | 0.99–1.00 | 1.00 | 0.92–1.00 |

**Observations.** WLS error grows linearly with outlier size; the protected
methods plateau at the clean-data level once the outlier exceeds ~8σ. Small
(3σ) outliers are hard for everyone — that is a fundamental detectability
limit (they look like noise), not an implementation defect. Huber's RMSE
plateaus slightly above DIA/Hampel because it only bounds, never removes,
the outlier's influence. The Hampel entry at 50σ (0.0118 m, F1 0.92) shows
the known redescending-estimator failure mode: an extreme outlier initially
swamps neighbouring baselines above the hard-rejection threshold in a few
geometries, and the estimator then rejects a clean baseline as well.

## 5. Contamination sweep (30 runs/point, seeds 9000–9029, 20σ outliers)

Command: `python -m gnss_adjust benchmark --scenario contamination-sweep --runs 150 --seed 9000 --output results/contamination_sweep`

Coordinate RMSE (m, mean) vs contamination fraction:

| fraction | WLS | DIA-reject | DIA-inflate | Huber | Hampel |
|---|---|---|---|---|---|
| 0.00 | 0.0075 | 0.0075 | 0.0075 | 0.0075 | 0.0075 |
| 0.05 | 0.0222 | 0.0118 | 0.0117 | 0.0124 | 0.0118 |
| 0.10 | 0.0332 | 0.0178 | 0.0176 | 0.0177 | 0.0160 |
| 0.20 | 0.0390 | 0.0249 | 0.0180 | 0.0244 | 0.0207 |
| 0.30 | 0.0500 | 0.0452 | 0.0304 | 0.0361 | 0.0374 |

**Observations.** All methods degrade with rising contamination, as theory
predicts. Sequential DIA-reject degrades fastest among the protected
methods at 30% contamination — one-at-a-time identification suffers from
masking when many outliers are present, and the rank guard limits how much
can be removed from a 16-baseline network. Covariance inflation and the
robust estimators, which act on all suspicious baselines simultaneously,
hold up better.

## 6. Specific-direction experiments

### Seed-42 directional case study

Command: `python scripts/run_sd_case_study.py`

| quantity | result |
|---|---:|
| planted / detected baseline | 3 / 3 |
| injected vector [m] | (+0.035, -0.025, +0.090) |
| estimated vector [m] | (+0.02762, -0.02431, +0.10582) |
| injected / estimated magnitude [m] | 0.09975 / 0.11204 |
| sign-invariant direction error | 6.736 degrees |
| SD statistic / critical value | 7.535 / 4.033 |
| 3D statistic / critical value | 18.927 / 5.422 |
| `abs(SD^2 - 3T)` | 7.11e-15 |
| SD vs full-3D maximum coordinate difference | 4.55e-12 m |

The SD result identifies the planted baseline, recovers its mostly vertical
direction, and verifies both paper identities. The estimated bias need not
equal the planted vector exactly because it contains the network's realized
Gaussian observation noise as well as the gross error.

### Multi-network equivalence verification

Across 100 networks, 1593 full-rank baseline blocks were tested and seven
zero-redundancy blocks were honestly marked untestable. The maximum
`abs(SD^2 - 3T)` was 2.27e-13; the maximum direction/vector cross-product
norm was 2.40e-17. In 100 coordinate comparisons, directional elimination
and full-3D elimination differed by at most 5.13e-11 m, passing the configured
1e-8 m tolerance.

### Magnitude and orientation recovery

The fixed-direction magnitude sweep used 50 runs per point. Direction
recovery becomes useful at the same scale at which the target becomes
statistically detectable:

| magnitude [sigma] | correct detection | target significant | median angle | 90th-percentile angle |
|---:|---:|---:|---:|---:|
| 1 | 0.00 | 0.00 | 44.99 deg | 74.25 deg |
| 3 | 0.08 | 0.08 | 26.76 deg | 72.35 deg |
| 5 | 0.40 | 0.40 | 14.78 deg | 41.71 deg |
| 8 | 0.88 | 0.90 | 9.26 deg | 19.90 deg |
| 10 | 0.98 | 1.00 | 7.28 deg | 15.06 deg |
| 20 | 0.98 | 1.00 | 3.69 deg | 7.05 deg |
| 50 | 1.00 | 1.00 | 1.37 deg | 2.94 deg |

At fixed 20-sigma magnitude across 200 random sphere orientations, correct
localization was 96.5%, the target was significant in 99.5%, median angular
error was 4.17 degrees, and the 90th percentile was 9.18 degrees. This shows
small but real geometry/noise dependence rather than perfect isotropy.

### Antenna-height-style session proxy

The preliminary experiment applies an 8 cm local-up proxy to up to two
S06 baselines in one synthetic session while other S06 sessions remain clean.
Across 100 runs it obtained 40.5% mean affected-baseline recall, 13.85 degrees
median affected direction error, and 1.13 non-affected significant baselines
per run. This is deliberately **not presented as a reproduction of the
paper's antenna result**: the current generator has block-diagonal baseline
covariance, random local geometry, and no session-level cross-correlation or
raw GNSS solution model. The weak result is recorded as evidence that the
paper's antenna application requires a dedicated correlated session model,
not as evidence against the validated core SD algebra.

## 7. Plots

Each results directory contains the PNGs referenced in the README:
`network.png`, `true_vs_estimated.png`, `coordinate_errors.png`,
`residual_norms.png`, `standardized_residuals.png`, `group_statistics.png`,
`robust_weights.png`, plus `monte_carlo_box.png`,
`rmse_vs_magnitude.png`, `f1_vs_magnitude.png`, and
`rmse_vs_contamination.png` in the benchmark directories. The SD directories
add `direction_recovery_vs_magnitude.png`, `orientation_sweep.png`, and
`incident_direction_error.png` with their raw and aggregate CSV/JSON data.
