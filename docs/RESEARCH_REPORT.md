# Outlier Diagnosis in GNSS Baseline Networks: Synthetic Benchmarks and Real-Data Validation

*Research report for Prof. Onkar Dikshit, National Centre for Geodesy, IIT Kanpur.*
*Every quantitative statement cites the experiment or derivation that produces it ([S] = reproducible script/output in this repository). Critical review: `SCIENTIFIC_AUDIT.md`.*

## 1. Research question

How reliably can classical and robust least-squares adjustment **detect**,
**identify** and **mitigate** contaminated baseline observations in GNSS
networks — synthetic and real — and what does specific-direction (SD)
testing add about the direction of a baseline error?

## 2. Mathematical framework

* **Observation model.** $\ell_{ij}=X_j-X_i+e_{ij}$, $\operatorname{Cov}(e_{ij})=Q_{ij}$
  (full 3×3); stacked $\ell+v=Ax$, one fixed datum station,
  $\hat x=(A^\top PA)^{-1}A^\top P\ell$, $v=A\hat x-\ell$,
  $Q_v=Q_\ell-A(A^\top PA)^{-1}A^\top$ (`MATHEMATICS.md` §2–3).
* **Testing.** Global $v^\top Pv/\sigma_0^2\sim\chi^2_r$ with the a-priori
  factor; local 3D test $T_i=(C_i^\top Pv)^\top(C_i^\top PQ_vPC_i)^{-1}(C_i^\top Pv)\sim\chi^2_3$;
  DIA: detect → identify the largest local statistic → reject or inflate →
  re-adjust.
* **Robust estimation.** Huber and Hampel IRLS on 3D blocks, standardised by
  the original weight matrix (`MATHEMATICS.md` §6).
* **Specific direction** (Nie et al. 2019). $\bar P=PQ_vP$,
  $g=-Pv$, $\hat d_i=\bar P_{ii}^{-1}g_i$, $|w^*|^2=g_i^\top\bar P_{ii}^{-1}g_i$,
  direction $\hat d_i/\|\hat d_i\|$ with $\operatorname{Cov}(\hat d_i)=\sigma_0^2\bar P_{ii}^{-1}$.

## 3. What the synthetic prototype established (before this work)

Eight-station synthetic networks with exact truth (`EXPERIMENTS.md`):
under a single 30σ outlier, WLS coordinate RMSE 28.2 mm vs DIA-reject 9.4 mm
and Hampel 9.2 mm, F1 ≈ 0.96–0.97 for all protective methods; DIA degrades
by masking at ≥ 20–30 % contamination; Hampel can swamp at ~50σ [S:
`scripts/generate_demo_results.py`]. These remain synthetic results.

## 4. Real data and processing pipeline

* **Data** [S: `configs/real_data/madrid_2024.yaml`, download log with
  SHA-256]: EPN/IGS stations MADR, CEBR, VILL, YEBE (10–109 km), 1–6 January
  2024, 24 h, 30 s, GPS L1C/L2W, IGS final orbits, IGS20 antenna model;
  reference: EPN combined weekly solution (epoch 2024-01-03). A GSI pair
  (Japan, 2005) verifies the chain.
* **Pipeline** [S: `src/gnss_adjust/gnss/`]: RINEX metadata and overlap
  checks → pinned RTKLIB 2.4.3 b34 (`rnx2rtkp`, config-only, output-precision
  patch) → one static solution per 24 h session (final forward epoch) →
  covariance reconstruction → satellite screening → network assembly that
  forbids same-session loops → adjustment CLI.
* **Verification** [S: `scripts/validate_rtklib_patch.py`, tests]: the patch
  changes formatting only (460 epoch pairs); covariance reconstruction
  matches RTKLIB's internal ENU output to 3×10⁻⁹ m; stock 0.1 mm output
  gives σ errors up to 22 % and correlation sign errors in 4.9 % of epochs;
  antenna-height settings apply as predicted (≤ 0.06 mm on GSI data).
  Three silent RTKLIB pitfalls were found and guarded (status-based static
  selection, unhealthy-satellite use, suppressed output).
* **Networks** [S: `scripts/build_real_network.py`]: one spanning tree per
  day, different trees on different days ⇒ 9 baselines, redundancy 18, six
  cross-session loops; quiet days 001–003, evaluation days 004–006.
* **Reproducibility** [S: audit V7]: a clean clone reproduces every
  numerical result identically; 102 automated tests pass.

## 5. SD versus full-3D testing

* **Identity.** For block-diagonal $Q_\ell$, $|w^*|^2=T_i^{\mathrm{DIA}}=$ the
  likelihood-ratio statistic of a free 3D bias. Verified with independent code
  on synthetic networks (10⁻⁷) and on the real network (all nine statistics
  identical; e.g. 39.785) [S: `tests/test_sd_dia_independent.py`, audit V3].
  Consequently SD-reject and DIA-reject take identical decisions.
* **Correlated case.** Only the general form (the SD formula) is the correct
  3D test; the DIA residual-block statistic deviates by > 5 % [S: same test].
* **Power.** The optimised direction has a χ²(3) null distribution, so SD has
  exactly the 3D test's power; a correctly *prespecified* direction (χ²(1))
  detects more often, a wrong one less [S: Experiment C, `REAL_DATA_MATH_NOTES.md` §6].
* **What SD adds** is the direction estimate: angular error 25° at 10 mm,
  ~12° at 15–20 mm, 3° at 50 mm (fixed-noise, conditional); reprocessed
  antenna-height errors are recovered as vertical within 2–15° once detected
  [S: Experiments C, D]. **No detection advantage is claimed.**

## 6. Validated experiments and findings

**6.1 False alarms under H₀** (resampled noise, real geometry; [S:
`audit_resampled_noise.py`]): global rejection 0.046–0.050 (nominal 0.05),
local 0–0.002 (nominal 0.001), DIA/SD flag in 0.35–0.6 % of runs. The tests
have their nominal size, and power matches the non-central χ²(3) within
sampling error.

**6.2 Detection and identification at equal false-alarm rate** ([S:
`fair_comparison.py`]; thresholds calibrated to 1 % on 3000 H₀ simulations,
evaluated on disjoint simulations, achieved 0.80–0.93 %; 300 runs/cell):

| 3-session network | 10 mm | 15 mm | 20 mm | ≥ 30 mm |
|---|---|---|---|---|
| correct identification DIA/SD | 0.40 | 0.75 | 0.92 | 1.00 |
| correct identification Huber = Hampel | 0.37 | 0.75 | 0.92 | 1.00 |
| misidentification (all) | ≤ 0.003 | ≤ 0.003 | ≤ 0.010 | 0 |

At equal size the methods are indistinguishable in detection and
identification. They differ in **coordinate accuracy** (RMSE at 50 mm: WLS
5.7, Huber 2.1, Hampel 1.8, DIA 1.8 mm). No unconditional ranking follows.

**6.3 Identifiability.** In the 2-session network the tests of baselines 3
and 6 (the only links of YEBE) have canonical correlation **1.000**: an
outlier on either cannot be attributed, giving ≈ 17 % misidentification for
every method [S: `audit_stats.json`, `fair_comparison.csv`].

**6.4 Conditional experiments.** With the one observed noise realisation
the numbers differ (DIA at 10 mm: 0.18 vs 0.40); on the anomaly-containing
days the natural inconsistencies are flagged too, and F1 against injected
labels drops (DIA 1.00 → 0.57 at 50 mm). Quiet-day results are best-case
[S: `experiment_BC/`, `injection_on_main_network.csv`].

**6.5 Unmodified real data** (no ground truth) [S: Experiment A]: the
evaluation-day network fails the global test (79.3 vs 28.9); baseline 8
(day 006) is the most inconsistent observation, a result robust to a
hypothetical cross-covariance model; further rejections are not.

**6.6 Physical errors** [S: Experiment D]: RINEX reprocessed with a wrong
antenna height is detected from 2–3 cm (conditional on one realisation); a
station-level error appears as two separately rejected baselines.

## 7. Covariance miscalibration and model sensitivity

* **Formal vs empirical** [S: `calibration.json`]: iono-free float σ is
  realistic horizontally and optimistic vertically; fixed-ambiguity σ is
  10–25× optimistic.
* **Non-stationarity** [S: `calibration_decomposition.json`]: repeatability
  factor 0.73 (days 001–003) vs 5.17 (004–006), F = 7.0, p ≈ 5×10⁻⁸, mainly
  in height (1.2 → 12.5).
* **Elevation dependence** [S: `audit_sensitivity/`, 252 reprocessed
  sessions]: the day-005/006 inconsistency is unchanged by orbits, satellite
  screening, antenna model and troposphere gradients but nearly vanishes with
  a 20° mask (VILL height deviation +20.4 → +4.2 mm; factor 5.17 → 1.47). It
  is an **elevation-dependent processing/model inconsistency**, not a
  confirmed station displacement.
* **CEBR** [S: `epn_ac_comparison.csv`, `antex_comparison.json`,
  `cebr_l2_signal_test.csv`]: ≈ 12 mm East of EPN on all days; the reference
  agrees across three EPN analysis centres (±0.7 mm); same antenna
  calibration as EPN; insensitive to orbits, mask, antenna model and
  satellite selection, **sensitive to the tracked L2 signal** (+4.4 mm with
  L2L). An **unresolved systematic discrepancy**, not an identified physical
  error; constant over sessions, it is invisible to every test.
* **Cross-baseline covariance** [S: `audit_cross_covariance.py`, a
  *hypothetical* model]: not available from single-baseline processing.
  Plausible within-session correlation would leave the global result and
  baseline 8 unchanged but alter later identifications, and would make the
  block-diagonal tests conservative (size 0.025–0.029 instead of 0.05). This
  remains an important statistical limitation.

## 8. Limitations

Four stations, six days, one region, GPS only, one processor, float
ambiguities; block-diagonal covariance; formal covariances conditional on
base coordinates and integer ambiguities; method comparisons valid for these
geometries and the formal stochastic model only; fixed-noise experiments are
conditional; two real-data effects (VILL, CEBR) unresolved; constant station
errors have no redundancy.

## 9. Next research directions

1. **Stochastic model first:** elevation-aware weighting and component-wise
   variance estimation, judged by repeatability and loop misclosures.
2. **Station-session bias hypotheses:** Baarda's C-matrix test for a bias
   shared by all baselines of one station in one session (equivalent to a
   congruence test), with SD-style direction estimation, evaluated at equal
   false-alarm rate against baseline-wise DIA on a network of ≥ 5 stations,
   with reprocessed antenna perturbations as ground truth. Identifiability
   must be checked first: in the present design only 6 of 12 station-sessions
   per network are distinguishable from baseline outliers
   (`NEXT_RESEARCH_DIRECTION.md`).
3. **Independent processing of CEBR** (another GNSS engine or PPP) to resolve
   the signal-dependent discrepancy.

Neither SD nor the station-session test is a new statistic; the defensible
contribution is a rigorous, reproducible real-data validation of when these
diagnostics identify the cause correctly — and when geometry or the
stochastic model makes that impossible.
