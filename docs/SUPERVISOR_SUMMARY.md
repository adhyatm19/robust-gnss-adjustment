# Robust Adjustment of GNSS Baseline Networks — From Synthetic Benchmarks to Real RINEX Data

*Research summary for discussion with Prof. Onkar Dikshit, National Centre for Geodesy, IIT Kanpur.*
*All figures are reproducible from the repository; details: `REAL_DATA_VALIDATION.md`, critical review: `SCIENTIFIC_AUDIT.md`.*

## 1. Problem

How reliably can classical and robust least-squares adjustment detect,
localise and mitigate contaminated GNSS baseline observations — in synthetic
networks and in a real network processed from RINEX — and what does
specific-direction (SD) detection add about the spatial direction of an
error?

## 2. Methods

* **Model.** Baselines $\ell_{ij}=X_j-X_i+e$ with full 3×3 covariance blocks;
  Gauss–Markov model $\ell+v=Ax$, one fixed datum station, WLS with
  $Q_{\hat x}=(A^\top PA)^{-1}$, $Q_v=Q_\ell-AQ_{\hat x}A^\top$.
* **Testing.** Global χ² test with the a-priori variance factor; component
  w-tests and 3D baseline tests; DIA with rejection or covariance inflation.
* **Robust estimation.** Huber and Hampel IRLS on 3D blocks, standardised by
  the fixed original weight matrix.
* **Specific direction.** $\bar P=PQ_vP$, $\hat d_i=\bar P_{ii}^{-1}g_i$,
  $|w^*|^2=g_i^\top\bar P_{ii}^{-1}g_i$ (Nie et al. 2019).
* **Real-data chain.** RINEX → RTKLIB 2.4.3 b34 (pinned build, output-precision
  patch) → one static solution per 24 h session with formal covariance →
  network assembly that forbids same-session loops → adjustment.
  GNSS processing is delegated to RTKLIB; nothing of it was re-implemented.

## 3. Data

* GSI GEONET 0759/3040 (Japan, 1 h, 3.3 km) — verification of the
  processing chain.
* Four EPN/IGS stations near Madrid (MADR, CEBR, VILL, YEBE; 10–109 km),
  1–6 January 2024, 24 h, 30 s, GPS L1C/L2W, IGS final orbits, IGS20 antenna
  calibrations; reference: EPN combined weekly solution (IGS20, epoch
  2024-01-03). 126 MB, every file checksummed.
* Networks: 9 baselines over three sessions (a different spanning tree per
  day), redundancy 18, six cross-session loops; "quiet" days 001–003 and
  "evaluation" days 004–006.

## 4. Principal results

**Mathematical/statistical (verified independently).**
* SD's statistic equals the 3D outlier test and the likelihood-ratio test of a
  free 3D bias, to numerical precision, on synthetic and on the real network.
  SD therefore has *exactly* the detection power of the 3D test; its
  contribution is the estimated direction with covariance
  $\sigma_0^2\bar P_{ii}^{-1}$.
* With correlated baselines only the general form
  $(C^\top Pv)^\top(C^\top PQ_vPC)^{-1}(C^\top Pv)$ — the SD formula — remains
  correct; the residual-block DIA statistic does not (> 5 % error in a test).
* On the real geometry, resampled-noise simulations confirm the nominal test
  size (global 0.046–0.050 vs 0.05; local ≈ 0.001) and the non-central χ²(3)
  power within sampling error.

**Processing (real data).**
* RTKLIB's default output rounds covariance terms to 0.1 mm; rebuilt 3×3
  blocks had σ errors up to 22 % and correlation sign errors in 4.9 % of
  epochs. The patched output matches RTKLIB's internal values to 3×10⁻⁹ m.
* Three silent RTKLIB pitfalls were found and guarded against: status-based
  selection of a weak static solution (39 mm height error); use of an
  unhealthy satellite with precise orbits, silently discarding 6 h of data;
  suppressed output in one mode.

**Uncertainty (real data).**
* Formal iono-free (float) σ is realistic in the horizontal and optimistic in
  height; fixed-ambiguity σ is 10–25× optimistic.
* The empirical variance factor is not stationary: 0.73 (days 001–003) vs
  5.17 (days 004–006), F = 7.0, p ≈ 5×10⁻⁸, almost entirely vertical. It
  largely disappears with a 20° elevation mask (1.47) — an
  elevation-dependent modelling/environment effect at VILL/YEBE on days
  005–006.
* CEBR is ≈ 12 mm East of the EPN reference on all six days, under every
  processing variant tried; the reference itself is consistent across EPN
  analysis centres (±0.7 mm). The offset is in our processing chain and
  unexplained; it cancels in every loop and is invisible to all tests.

**Outlier diagnosis.**
* *Unmodified real data:* WLS fails the global test (79.3 vs 28.9);
  DIA = SD reject baseline 8 (and, less robustly, 9) from day 006; Huber and
  Hampel down-weight 8. The SD direction of 8 is ~79° elevation (mostly
  vertical). These flags diagnose a model inconsistency; there is no ground
  truth, so no accuracy is claimed. Under a hypothetical within-session
  correlation model, only the identification of 8 is stable.
* *Controlled injections into real-data baselines (resampled noise, 3-session
  network):* correct identification by DIA/SD 0.38 / 0.73 / 0.91 / ~1.0 at
  10 / 15 / 20 / ≥ 30 mm; Huber 0.15 / 0.53 / 0.80 / ≥ 0.98. On a weaker
  2-session network one station is controlled by only two baselines, whose
  tests are perfectly correlated (canonical correlation 1.0) — outliers there
  are not identifiable, causing ~17 % misidentification.
* *Realistic conditions matter:* repeating the injections on the
  anomaly-containing days lowers DIA's F1 at 50 mm from 1.00 to 0.57, while
  Huber/Hampel reach 0.73/0.80 — method rankings depend on whether the
  "clean" data are actually clean.
* *Physical errors:* reprocessing real RINEX with a wrong antenna height
  (applied to ≤ 0.23 mm of the prediction) gives detection from 2–3 cm and an
  SD direction within 2–15° of the vertical; a station-level error is
  rejected as two unrelated baselines.

## 5. Limitations

Four stations, six days, GPS only, one processor, float ambiguities for the
main network; no cross-baseline covariance (single-baseline processing);
formal covariances conditional on the base station and on correct
ambiguities; robust-vs-DIA comparisons not yet at matched size; two
unexplained real-data effects (CEBR offset, VILL elevation-dependent
anomaly). Several of my earlier interpretations were corrected by the audit
(misidentification reported as false positives; single-realisation
statistics; over-interpretation of the VILL flags).

## 6. Most promising direction

What the work supports is a **validated real-data test bed** showing that
(i) on real data the stochastic model — elevation dependence, anisotropy,
non-stationarity — dominates outlier diagnosis more than the choice of
estimator, and (ii) real anomalies act on a *station within a session*
rather than on one baseline. Proposed next steps:

1. **Stochastic model:** elevation-aware weighting and component-wise
   variance estimation, evaluated by repeatability and loop misclosures.
2. **Station-session bias hypotheses** — Baarda's C-matrix test for a bias
   shared by all baselines of one station in one session (equivalent to a
   congruence test), with SD-style direction estimation — evaluated at
   matched size against baseline-wise DIA on a network with ≥ 5 stations,
   using reprocessed antenna perturbations as ground truth
   (`NEXT_RESEARCH_DIRECTION.md`).

Neither SD nor the station-session test is a new statistic; the potential
contribution is a rigorous, real-data validation of when such diagnostics
identify the physical cause correctly.
