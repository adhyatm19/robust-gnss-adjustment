# Real-Data Validation: RINEX → RTKLIB → Baseline Network → Robust Adjustment

This document reports what was done with **actual GNSS observations**, how it
was verified, what was found, and what remains uncertain. Every number below
was produced by the scripts named next to it and can be regenerated with
`scripts/reproduce_real_data.sh` (≈5 min after the 126 MB download).
Companion notes on conventions and statistics: [`REAL_DATA_MATH_NOTES.md`](REAL_DATA_MATH_NOTES.md).

> **Audit status.** This report incorporates the corrections of
> [SCIENTIFIC_AUDIT.md](SCIENTIFIC_AUDIT.md) (approved). Error types are
> reported separately (false alarm, detection, correct identification,
> misidentification); probabilities come from resampled-noise simulations;
> fixed-noise and anomaly-containing results are labelled as conditional.

> **Three kinds of evidence are kept strictly apart in this report.**
> (1) *Unmodified real data* (Experiment A, calibration, loop misclosures,
> EPN comparison): no ground-truth outlier labels exist, so **no detection
> accuracy is claimed**. (2) *Controlled injections into real-data-derived
> baselines* (Experiments B, C): known labels, but the errors are artificial.
> (3) *Physically reprocessed RINEX* with a deliberately wrong antenna height
> (Experiment D): a real processing chain with a known, physically meaningful
> error. Synthetic results from earlier work are in `EXPERIMENTS.md` and are
> not repeated here.

---

## 1. Processing engine

| item | value |
|---|---|
| processor | RTKLIB `rnx2rtkp` 2.4.3 b34 (BSD-2), `https://github.com/tomojitakasu/RTKLIB`, tag `v2.4.3-b34`, commit `180043ee24b6d2b168f98b64be15f69d50046b1a` |
| build | `scripts/build_rtklib.sh` → gitignored `tools/rtklib/` (no system install); `make LDLIBS="-lm -lrt"` because the b34 makefile links an IERS library only needed with `-DIERS_MODEL` |
| variants | `rnx2rtkp-stock` (unmodified) and `rnx2rtkp-hp` (identical source + `third_party/rtklib/rnx2rtkp-output-precision.patch`, `printf` widths only) |
| provenance | `tools/rtklib/build_info.json` (commit, patch SHA-256, compiler, binary SHA-256); every run records its full config file and command |
| invocation | config file only (`-k`); no other command-line processing options, because rnx2rtkp CLI options override the config and would make the recorded configuration incomplete |

### 1.1 What the six uncertainty numbers mean (verified in source)

`src/solution.c` (`outecef`) prints `sdx sdy sdz sdxy sdyz sdzx`.
The first three are standard deviations; the last three are **signed square
roots of covariances**, `sqvar(c) = sign(c)·√|c|` — neither covariances nor
correlation coefficients. Hence `σ_xy = sign(sdxy)·sdxy²`, order (xy, yz, zx).
Internally the covariance is stored as `float32` (`sol_t.qr`), which bounds
any reconstruction at ≈6×10⁻⁸ relative precision.

### 1.2 Output-precision patch and its validation (`scripts/validate_rtklib_patch.py`)

The stock binary prints these fields with `%8.4f` (0.1 mm). Real data, GSI
GEONET 0759→3040 (2005-04-02, GPS, 1 h, 3.34 km, vendored in
`tests/data/gsi_0759_3040/`), four configurations, all 115 epochs each:

| check | result |
|---|---|
| epochs, status, satellites, age, ratio: stock vs hp | identical in all 460 epoch pairs |
| max \|stock − hp\| over all numeric fields | 5.0×10⁻⁵ m = half the stock print resolution (rounding only) |
| my ECEF→ENU rotation of hp covariance vs RTKLIB's *own* ENU output | ≤ 3.2×10⁻⁹ m in the signed-sqrt fields |
| stock-file reconstruction error vs hp reference | σ up to **22 %** wrong (median 0.9 %); correlation up to **0.34** wrong; **4.9 % of epochs flip the sign** of an ENU correlation; never non-PD here |

Verdict: the patch changes formatting only; the full covariance is
reconstructed correctly; **stock `.pos` files are not adequate for 3×3
covariance blocks** at the mm level (`results/rtklib_validation/`). The
pipeline flags any stock-precision covariance whose ENU sigmas or
correlations move by more than 5 % / 0.05 within its rounding box.

### 1.3 Antenna-metadata precedence (verified by experiment, `tests/test_gnss_rtklib_integration.py`)

Source (`postpos.c`): with `antN-anttype=*` RTKLIB takes antenna type **and
height from the RINEX header, overwriting** any `antN-antdel*`; with an
explicit type, the config deltas are used and the header is ignored. Verified
on the GSI pair with a +10 cm height: rover → baseline changes by
−0.10 m·**u**_rover (0.05 mm agreement), base → +0.10 m·**u**_base
(0.06 mm); `*` + config delta → output unchanged; `*` + edited header copy →
perturbation applied; explicit type + edited header → header ignored. The
pipeline therefore always writes an **explicit** type and delta taken from
the parsed header (so the applied values are logged), and warns if RTKLIB
found no ANTEX calibration for the type.

### 1.4 Processing pitfalls found on real data (all now guarded)

1. **`combined` + `out-solstatic=single` can select a weak solution.** RTKLIB
   picks the last epoch with the best *status label*. On Madrid day 001 that
   was an unsmoothed backward-filter epoch from the last minutes of the day
   ("fix", ratio 3.3, σ_U 19 mm, 39 mm off EPN). The pipeline instead runs a
   **forward static filter** and takes its **final epoch** (= estimate from all
   data), records the fix fraction, and never uses the cumulative, mutually
   dependent epoch outputs as separate observations.
2. **Unhealthy satellites used with precise orbits.** On 2024-01-01 G27 had an
   IGS-flagged bad clock (`999999.999999`) *and* broadcast health 63; RTKLIB
   b34 falls back to the broadcast clock and ignores the health flag, so its
   ≈100 m pseudorange error failed the single-point χ² check at **every epoch
   for 6 h, which RTKLIB silently skipped**. `gnss_adjust.gnss.satscreen` now
   excludes, per day and logged, any GPS satellite with a bad SP3 clock or a
   non-zero broadcast health word (G01, G27; G03 on day 005; G08 on day 006).
   A run that solves < 90 % of the common epochs is rejected.
3. **ENU output suppressed in combined mode** (unset smoothed base position
   in `postpos.c`); XYZ output is unaffected and is the only format used.
4. **L1/L2 ambiguity resolution fails beyond ~40 km** with these data and
   RTKLIB (99 km: 1.8 % of epochs fixed with estimated slant ionosphere).

## 2. Real dataset

| item | value |
|---|---|
| source | BKG GNSS data centre `https://igs.bkg.bund.de/root_ftp` (anonymous HTTPS); manifest `configs/real_data/madrid_2024.yaml` |
| stations | MADR (JAVAD TRE_G3TH, AOAD/M_T), CEBR (SEPT POLARX5TR, SEPCHOKE_B3E6), VILL (SEPT POLARX5, SEPCHOKE_B3E6), YEBE (LEICA GR50, LEIAR20 LEIM) — EPN/IGS, Madrid area |
| baseline lengths | MADR–CEBR 10.4, MADR–VILL 25.3, CEBR–VILL 35.3, VILL–YEBE 73.8, MADR–YEBE 99.0, CEBR–YEBE 108.8 km |
| dates | 2024-01-01 … 2024-01-06 (DOY 001–006, GPS week 2295), 24 h daily files |
| sampling / signals | 30 s; RINEX 3.04; processed GPS only, L1C + L2W forced identically at all receivers |
| orbits / antennas | IGS final SP3 (IGS20) for each day ± 1; `igs20_2287.atx` type-mean calibrations |
| reference coordinates | EPN combined weekly solution `EUR0OPSSNX_20233650000_07D_07D_SOL.SNX` (IGS20, epoch 2024-01-03 12:00), an independent multi-analysis-centre product |
| volume | 40 files, 126.3 MB; URL, size, SHA-256 of every file in `results/real_network/madrid_2024_download_log.json` |
| metadata QC | RINEX antenna type and height match SINEX `SITE/ANTENNA` and `SITE/ECCENTRICITY` for every station and day; all pairs share 2880 common epochs (`processing/metadata_qc.json`) |

Within one week the relative station motion is far below 1 mm, so the
weekly reference epoch is adequate.

## 3. Baseline processing (`scripts/process_real_network.py`)

Every pair (6) on every day (6) was processed as one 24 h static session;
36 iono-free solutions plus 18 L1/L2 solutions on the three short pairs, 0
failures, 2859–2867 of 2880 epochs solved per day (after screening).
Configuration `configs/real_data/madrid_2024_processing.yaml`, fixed *a
priori*:

* **main (`iflc_float`)**: ionosphere-free combination, float ambiguities,
  estimated ZTD, precise orbits, solid tides, PCV on, 10° mask — one uniform
  model for all baselines (mixing fixed and float stochastic models in one
  network would confound the calibration);
* **secondary (`l12_fixed_short`)**: L1/L2 with continuous ambiguity
  resolution, broadcast ionosphere, Saastamoinen; short pairs only; used to
  characterise fixed-solution covariances, *not* mixed into the networks.

Base coordinates are the EPN reference values; RTKLIB holds them fixed, so the
reported covariance is that of the baseline **conditional on the base
position** (and, for fixed solutions, on the integer ambiguities). A base
error δ enters the baseline only at second order (≈ L·|δ|/20 000 km).

### 3.1 Formal uncertainty versus repeatability and versus EPN (6 days)

| pair | L [km] | empirical σ E/N/U [mm] | formal σ E/N/U [mm] | mean − EPN E/N/U [mm] |
|---|---|---|---|---|
| MADR→CEBR | 10.4 | 2.2 / 1.6 / 5.6 | 3.5 / 1.6 / 3.8 | +9.7 / −2.4 / +4.1 |
| MADR→VILL | 25.3 | 2.2 / 1.0 / 7.3 | 3.4 / 1.6 / 3.9 | −3.7 / +0.1 / +4.1 |
| CEBR→VILL | 35.3 | 1.1 / 2.1 / 11.5 | 3.4 / 1.6 / 3.7 | −15.7 / +3.9 / +0.1 |
| VILL→YEBE | 73.8 | 4.2 / 2.8 / 15.3 | 3.4 / 1.6 / 3.8 | +7.0 / −2.7 / +0.3 |
| MADR→YEBE | 99.0 | 4.3 / 2.2 / 8.8 | 3.5 / 1.6 / 4.0 | +3.6 / −2.6 / −2.6 |
| CEBR→YEBE | 108.8 | 3.9 / 1.3 / 7.8 | 3.4 / 1.6 / 3.8 | −10.7 / +1.2 / +1.4 |
| *L1/L2 fixed* MADR→CEBR | 10.4 | 1.6 / 1.7 / 4.8 | 0.4 / 0.3 / 0.4 | +5.2 / −5.0 / +11.4 |
| *L1/L2 fixed* MADR→VILL | 25.3 | 5.7 / 2.1 / 10.1 | 0.5 / 0.3 / 0.4 | −10.5 / +0.7 / +35.1 |
| *L1/L2 fixed* CEBR→VILL | 35.3 | 5.7 / 2.5 / 10.9 | 0.5 / 0.3 / 0.4 | −18.6 / +6.3 / +22.5 |

Findings (figure `formal_vs_empirical_sigma.png`, `repeatability_by_day.png`):

* **The RTKLIB covariance is not uniformly optimistic.** Iono-free float:
  horizontal formal σ is realistic or pessimistic (E and N), the vertical is
  optimistic by ≈1.5–4×. Fixed L1/L2: optimistic by ≈10–25× in σ
  (k² ≈ 70–390) — the covariance is conditional on correct integers and
  ignores unmodelled troposphere.
* **Shape mismatch**: a scalar factor cannot fix it (calibration-day variance
  ratios E 0.77, N 0.33, U 1.25).
* **Non-stationarity** (pooled 3D factor k², `calibration.json`):
  days 001–003 **k² = 0.73** (95 % CI 0.48–1.24, 36 dof);
  days 004–006 **k² = 5.17** (3.42–8.72); F = 7.0, p ≈ 5×10⁻⁸, almost all in
  Up (1.2 → 12.5). Days 005/006 show an **elevation-dependent
  processing/model inconsistency** involving VILL (and YEBE on 006): it nearly
  vanishes with a 20° mask (k² 1.47). Not a confirmed station displacement.
  A calibration from quiet days does not transfer.
* **A station-level systematic at CEBR**: every CEBR baseline places CEBR
  ≈ +10 … +16 mm **East** of EPN, on all six days, with 1–2 mm scatter
  (and in both strategies). It is consistent across sessions, so it cancels in
  every loop and is invisible to repeatability and to network testing; only
  the external reference reveals it. **Unresolved systematic discrepancy:**
  the reference is consistent across three EPN analysis centres (±0.7 mm);
  EPN and this pipeline use the same SEPCHOKE_B3E6 calibration; the offset is
  insensitive to orbits, mask, antenna model and satellite selection but
  sensitive to the tracked L2 signal (+4.4 mm with L2L). Reported, not
  corrected, and not identified as a physical error.

## 4. Network construction and cross-baseline correlation

(`gnss_adjust.gnss.network_build`, `scripts/build_real_network.py`)

* **Algebraic independence** — per session, a baseline that closes a loop
  with others of the *same* session is *trivial* (computed from the same
  data). `check_session_independence` rejects such designs. Each evaluation
  session contributes one spanning tree, and the trees differ between days, so
  **every loop closes across sessions**.
* **Stochastic independence is not implied.** Even algebraically independent
  baselines of one session share a receiver (e.g. MADR→CEBR and MADR→VILL on
  day 004) and are correlated. Single-baseline RTKLIB processing provides
  **no cross-baseline covariance and none is invented**: the observation
  covariance is block-diagonal, and that is a stated limitation. Between days
  the independence assumption is also only approximate (sidereally repeating
  multipath, shared antenna models, the station-level effects of §3.1).
* **Empirical probe** (formal covariance, loop statistic w'C⁻¹w ~ χ²(3) if
  the block-diagonal model held; `loop_misclosures.png`):

  | loops | n | values | mean |
  |---|---|---|---|
  | same-session (trivial), day 004 | 3 | 0.26, 1.84, 2.63 | 1.6 |
  | cross-session, quiet days 001–003 | 6 | 1.07, 0.89, 0.63, 2.91, 1.82, 5.08 | 2.1 |
  | cross-session, evaluation days 004–006 | 6 | 5.52, 8.38, 8.94, 7.63, 8.18, 9.01 | 7.9 |

  Same-session loops close better than independence predicts (consistent with
  positive correlation through shared data, but only three loops — weak
  evidence); quiet-day loops agree with χ²(3); evaluation-day loops do not.
  **These loops are not independent of each other**, so the table is a
  diagnostic, not a formal test.

Networks written in the repository format (`results/real_network/networks/`):

| network | sessions | baselines | redundancy r | use |
|---|---|---|---|---|
| `main_formal` / `main_calibrated` | 004–006 | 9 | 18 | Experiment A (formal; formal × 0.73) |
| `quiet_formal` | 001–003 | 9 | 18 | base for Experiments B–D |
| `quiet_2sessions` | 001–002 | 6 | 9 | weaker geometry for B, C |
| `trivial_day004` | 004 only | 6 | 9 | demonstration only (`allow_trivial=True`) |

Each directory has `stations.csv`, `baselines.csv`, `covariance_blocks.csv`
(unchanged format) plus `baseline_metadata.csv` (session, strategy, fix/float,
fix fraction, covariance source/scale, rounding flag),
`reference_coordinates.csv` (EPN, frame, epoch) and `network_design.json`
(sessions, redundancy, loops). Datum: MADR fixed at its EPN coordinate.
Rank(A) = 9 = number of parameters, the network is connected, cond(N) is reported.

## 5. Experiment A — unmodified real data (`results/real_network/experiment_A/`)

`python -m gnss_adjust adjust --network results/real_network/networks/main_formal --reference … --output …`
(α_global 0.05, α_local 0.001, Huber c 1.5, Hampel 1.5/3.5/8, σ₀² = 1).
Baselines 1–3 are day 004, 4–6 day 005, 7–9 day 006. There is **no ground
truth**: nothing here is scored as correct or incorrect detection.

| method | T = v'Pv | critical | σ̂₀² | flagged | iterations |
|---|---|---|---|---|---|
| WLS | 79.29 | 28.87 | 4.41 | – | 1 |
| DIA-reject = SD-reject | 20.87 | 21.03 | 1.74 | 8, 9 | 3 |
| DIA-inflate | 21.85 | 28.87 | 1.21 | 8, 9 | 3 |
| Huber / Hampel | 54.1 / 53.7 | 28.87 | 3.0 | 8 (w ≈ 0.42) | 15 / 20 |

What the data support:

* The day-004–006 network is **inconsistent with the stochastic model**
  (global test rejected; independently recomputed, audit V4).
* Baseline 8 (CEBR→VILL, day 006) is the most inconsistent observation
  (T_b = 39.8; SD direction 79° elevation). This identification is stable
  under a hypothetical within-session correlation model; **the second
  rejection (9) is not** (audit M5).
* **Interpretation.** The day-005/006 inconsistency is an
  **elevation-dependent processing/model inconsistency**: it survives changes
  of orbits, satellite screening, antenna model and troposphere gradients,
  but nearly disappears with a 20° elevation mask (implied VILL height
  deviation +20.4 → +4.2 mm on day 006; `audit_sensitivity/`). It is **not a
  confirmed station displacement** and not shown to be a gross baseline error;
  plausible causes are low-elevation multipath/obstruction or troposphere
  mismodelling, none of which is established.
* Coordinates vs EPN: CEBR ≈ +11.5 mm East for every method (§3.1, §9).

The quiet-day network (001–003) passes the global test (12.37 < 28.87) and
nothing is flagged.

## 6. Experiments B — detection, identification and coordinate accuracy

Four quantities are reported separately throughout:
**false alarm** (a flag when no outlier exists), **detection** (any flag when
one baseline is contaminated), **correct identification** (exactly the
contaminated baseline flagged) and **misidentification** (a different
baseline flagged after a true outlier). Two noise designs are distinguished:
*fixed-noise* (the one observed real realisation plus injected biases) and
*resampled-noise* (formal-covariance noise drawn around the real-network
solution, so truth is known and probabilities are estimable).

### 6.1 False alarms under H₀ (resampled noise, `audit/audit_stats.json`)

2000 H₀ runs per geometry: global test rejects 0.046–0.050 (nominal 0.05);
local 3D tests exceed 0.0000–0.0020 per baseline (nominal 0.001); DIA/SD
flag a baseline in 0.35–0.6 % of runs; Huber/Hampel (weight < 0.5 rule) in
0 %. **Test size is correct.**

### 6.2 Size-matched comparison (resampled noise, `scripts/audit/fair_comparison.py`)

Thresholds calibrated on 3000 H₀ runs (seed A) so that each method flags in
1 % of H₀ runs (DIA/SD: max T_b > 16.4; Huber/Hampel: min weight < 0.71;
r = 9: 15.1 / 0.83); evaluated on **disjoint** simulations (seed B): achieved
false-alarm rate 0.93 % / 0.80 % / 0.80 % (r = 18), 0.90 % / 0.83 % / 0.83 %
(r = 9). 300 single-outlier runs per cell, random direction (binomial SE ≲ 0.03):

| geometry | mm | detection DIA/SD · Huber=Hampel | correct identification | misidentification |
|---|---|---|---|---|
| 3 sessions, r = 18 | 10 | 0.40 · 0.38 | 0.40 · 0.37 | 0.003 · 0.003 |
| | 15 | 0.76 · 0.76 | 0.75 · 0.75 | 0.003 · 0.003 |
| | 20 | 0.93 · 0.92 | 0.92 · 0.92 | 0.010 · 0.003 |
| | ≥ 30 | 1.00 · 1.00 | 1.00 · 1.00 | 0.000 · 0.000 |
| 2 sessions, r = 9 | 10 | 0.42 · 0.35 | 0.35 · 0.29 | 0.07 · 0.06 |
| | 20 | 0.87 · 0.83 | 0.76 · 0.70 | 0.11 · 0.13 |
| | 50 | 1.00 · 1.00 | 0.83 · 0.82 | 0.17 · 0.18 |

* **At matched size, detection and identification are statistically
  indistinguishable between DIA/SD and Huber/Hampel** on these geometries
  (differences ≤ 0.07, mostly within 1–2 SE). The earlier statement that
  robust IRLS is "less sensitive" was largely an artefact of the weight < 0.5
  rule. Huber and Hampel make identical decisions because the calibrated
  thresholds lie where their weight functions coincide. **No unconditional
  ranking is claimed.**
* **Coordinate accuracy** (RMSE of non-datum coordinates vs truth, standard
  procedures, reported separately from detection): at 50 mm, r = 18 —
  WLS 5.66, DIA-reject 1.81, Hampel 1.82, Huber 2.11 mm; at 10 mm all
  ≈ 2.0 mm. Rejection or redescending weights remove the bias; Huber leaves
  part of it.

### 6.3 Statistically indistinguishable baselines (weak geometry)

In the 2-session network station YEBE is connected only by baselines 3
(VILL→YEBE) and 6 (MADR→YEBE). The canonical correlation between their local
test statistics is **1.000**: an outlier on either produces the same residual
pattern, so **no method can tell which baseline is wrong**. With outliers
placed uniformly, the expected misidentification is 2/6 × ½ ≈ 0.17, which is
what all methods show at large magnitudes (0.17–0.18). In the 3-session
networks the largest canonical correlation is 0.52–0.54 (`audit_stats.json`).

### 6.4 Fixed-noise and anomaly-containing conditions (conditional results)

* *Fixed-noise, quiet network* (`experiment_BC/`, 40 trials/cell): conditional
  on the one observed realisation; e.g. DIA correct identification at 10 mm
  0.18 versus 0.40 with resampled noise. Use as illustration only.
* *Fixed-noise, anomaly-containing network* (days 004–006,
  `audit/injection_on_main_network.csv`): the natural day-006 inconsistencies
  are flagged in addition to the injected bias (≈ 1.6 extra flags/run for
  DIA, 0.5–0.8 for Huber/Hampel), so DIA's F1 against injected labels falls
  from 1.00 to 0.57 at 50 mm (Huber 0.73, Hampel 0.80). These extra flags are
  genuine data inconsistencies, not Type I errors, and F1 is not an
  appropriate score in this condition. **Quiet-day results are best-case.**

## 7. Experiment C — specific-direction validation (`experimentC_sd_direction.png`)

* **Identity check on the real network**: SD |w*|², the DIA 3D statistic and
  an independently computed likelihood-ratio statistic agree to
  2×10⁻⁶ relative (absolute-ECEF parameterisation; to 10⁻⁷–10⁻¹⁵ on local
  synthetic networks, `tests/test_sd_dia_independent.py`).
* **Detection follows theory** (resampled noise, 3000 runs,
  `audit/power_resampled.csv`): local SD/3D detection rates agree with the
  non-central χ²(3) power at λ = dᵀP̄ᵢᵢd within binomial SE in every λ bin.
  The fixed-noise figures below are conditional on one realisation.
* **Optimised vs prespecified direction**: a 1D test along the *correct*
  prespecified direction (χ²(1), same λ) detects more often (10 mm: 0.45 vs
  0.35; 15 mm: 0.95 vs 0.85); a test along local *up* applied to random
  directions loses most of its power (15 mm: 0.18). SD's optimised direction
  has exactly the 3D test's power — its added value is the **estimated
  direction**, not detection.
* **Direction recovery** (one outlier, random direction, r = 18): median
  angular error 33.5° (5 mm), 25.3° (10), 11.5° (15), 12.0° (20), 6.3° (30),
  3.3° (50), 2.3° (80); errors fall ≈ ∝ 1/|d|. Below the detection threshold
  the direction is essentially uninformative.
* **Predicted direction uncertainty** from Cov(d̂) = σ₀²P̄ᵢᵢ⁻¹ (largest
  perpendicular σ): median (actual angle / predicted σ) = 0.77 — consistent
  with a conservative bound; a calibrated coverage test is not yet done.
* With two contaminated baselines the second bias contaminates the first
  direction estimate (results kept separate in `sd_direction_summary.csv`).

## 8. Experiment D — reprocessed antenna-height errors (`scripts/real_antenna_experiment.py`)

The RINEX files were not modified; the wrong height was entered in the
RTKLIB configuration and the affected baselines reprocessed.
**Verification**: every reprocessed baseline changed by the predicted
∓Δh·**u** to ≤ 0.0016 mm (station as rover) and ≤ 0.23 mm at Δh = 10 cm
(station as base; `experiment_D/verification.json`).

| case | Δh [cm] | baseline | |w*| (crit 4.03) | vertical 1D |w| (crit 3.29) | SD angle to vertical | DIA flags | Huber/Hampel flags |
|---|---|---|---|---|---|---|---|
| D1 YEBE, day 002 | 1 | MADR→YEBE | 2.98 | 2.76 | 17.7° | – | – |
| | 2 | | 4.89 | 4.76 | 10.3° | 6 | – |
| | 3 | | 6.86 | 6.76 | 7.2° | 6 | 6 |
| | 10 | | 20.79 | 20.76 | 2.3° | 6 | 6 |
| D2 VILL, day 001 | 2 | MADR→VILL / VILL→YEBE | 2.69 / 2.99 | 2.29 / 2.78 | 26.7° / 17.6° | – | – |
| | 3 | | 4.35 / 4.48 | 4.12 / 4.34 | 15.0° / 11.2° | 2, 3 | – |
| | 5 | | 7.90 / 7.55 | 7.78 / 7.46 | 7.9° / 6.5° | 2, 3 | 2, 3 |

* A 2 cm antenna-height error on one baseline is detected by DIA/SD; SD says
  "vertical" (10°). The prespecified vertical test has the larger margin over
  its critical value (4.76/3.29 vs 4.89/4.03).
* The station-level case (D2) is detected as **two separate rejections**;
  none of the methods expresses "one station, one cause". At 1–2 cm the VILL
  height error biases VILL's height by 4–8 mm without any flag.

## 9. Limitations

1. Four stations, six days, one region, GPS only, one processor (RTKLIB b34),
   iono-free float ambiguities for the main network.
2. **Missing cross-baseline covariance is an important statistical
   limitation.** Single-baseline processing provides none; the
   block-diagonal model ignores correlation between baselines sharing a
   receiver in a session. A hypothetical correlation model leaves the global
   result and the identification of baseline 8 intact but changes which
   further baselines would be flagged, and makes the block-diagonal tests
   conservative (audit M5, M8).
3. Formal covariances are conditional (base position, integers), anisotropic
   in error and non-stationary (k² 0.73 vs 5.17, F = 7.0); they are a
   precision model on quiet days only.
4. **VILL** (days 005/006): an elevation-dependent processing/model
   inconsistency, not a confirmed station displacement.
5. **CEBR**: an unresolved ≈ 12 mm East systematic discrepancy relative to a
   reference that is internally consistent across EPN analysis centres; it
   is insensitive to orbits, elevation mask, antenna model and satellite
   selection but **sensitive to the tracked L2 signal** (+4.4 mm with L2L).
   Not a conclusively identified physical error. It cancels in every loop and
   is undetectable by all tests (zero redundancy for constant station errors).
6. Fixed-noise experiments (B, C, D) are conditional on one realisation;
   probabilities are taken from resampled-noise simulations.
7. Absolute-ECEF parameterisation costs ≈10⁻⁶ relative precision in test
   statistics (irrelevant for decisions).
