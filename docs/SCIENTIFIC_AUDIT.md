# Scientific Audit of the Real-Data Extension

An independent, adversarial re-examination of the RINEX integration and the
real-data experiments (commit `994dae7`), written to decide whether the work
is ready to present. Earlier conclusions were treated as hypotheses to be
tested, not as results to be defended. All new evidence is reproducible:

```bash
python scripts/process_real_network.py configs/real_data/madrid_2024_sensitivity.yaml \
    --output results/real_network/audit_sensitivity      # processing-model sensitivity (252 runs)
python scripts/audit/audit_real_data.py                   # calibration, EPN ACs, ANTEX, L2 test, injections
python scripts/audit/audit_resampled_noise.py             # size / power / identifiability (resampled noise)
python scripts/audit/audit_cross_covariance.py            # HYPOTHETICAL cross-baseline covariance
```

Outputs: `results/real_network/audit/`, `results/real_network/audit_sensitivity/`.

**Severity scale.** *Critical* threatens validity · *Major* changes a research
conclusion · *Minor* worth correcting · *Verified* independently checked.

**Bottom line.** No mathematical error was found in the estimators or test
statistics; their size and power were independently confirmed on the real
geometry. One critical reproducibility defect was found and fixed. Several
interpretations in `REAL_DATA_VALIDATION.md` and `RESEARCH_HANDOFF.md` were
too strong or wrong (§3); they are superseded by this document.

---

## 1. Verified (independently checked)

| # | claim | independent evidence |
|---|---|---|
| V1 | Test **size** is correct on the real geometries | 2000 resampled-noise H₀ realisations per network (`audit_stats.json`): global rejection 0.048 / 0.050 / 0.046 (nominal 0.05); per-baseline local exceedance 0.0000–0.0020 (nominal 0.001); DIA/SD flag anything in 0.35–0.6 % of H₀ runs |
| V2 | **Power** follows the non-central χ²(3) | 3000 resampled single-outlier runs, binned by λ = dᵀP̄d: empirical vs theory within binomial SE in every bin (largest deviation 1.7 SE: 0.123 vs 0.158 at λ∈(5,10]) |
| V3 | SD \|w*\|² = DIA T_b = likelihood-ratio statistic **on the real network** | dense GLS with own design matrix, local origin, LR by augmentation: all nine T_b agree to printed precision (e.g. 39.785, 27.503); previously this was shown independently only on synthetic networks |
| V4 | Experiment A numbers | independent recomputation: v'Pv 79.2903 (= repo), dof 18, χ²₀.₉₅(18) 28.869, coordinates to 7×10⁻⁹ m, DIA second step T(9) = 18.637 (repo 18.64) |
| V5 | Degrees of freedom and critical values | r = 27 − 9 = 18; after one rejection r = 12, critical 21.03; χ²₀.₉₉₉(3) = 16.27, √ = 4.03; two-sided N(0,1) at 0.001: 3.29 |
| V6 | Metric formulas | precision, recall, F1, FPR recomputed from TP/FP/FN for all 20 160 rows: exact (convention: precision = 0 if nothing flagged while outliers exist) |
| V7 | **Clean-checkout reproduction** (after fix C1) | fresh clone + fresh venv + RTKLIB from GitHub + 126 MB re-download: 102 tests pass; all 37 CSV and 13 JSON result files numerically **identical** to the committed ones (only timings, retrieval time and binary hashes differ); 40/40 data files match recorded SHA-256 |
| V8 | RTKLIB pin, patch and conventions | commit `180043e` verified by the build script; patch SHA-256 recorded; validation summary reproduced exactly except binary hashes (m1) |
| V9 | Datum and epoch treatment | EPN weekly epoch 2024-01-03 12:00 within ±3 days of all data; SINEX has no velocities (relative motion of the four stations ≪ 0.1 mm/week); SP3 frame `IGS20`; comparisons are translation-free relative to the datum station; base positions = EPN marker coordinates, antenna heights consistent RINEX↔SINEX |
| V10 | The SD/DIA tautology identified in the Stage 0 audit is resolved | the substantive identity is established by independent code (V3, `tests/test_sd_dia_independent.py`); the 3360-trial "identical decisions" is a consequence of it, not separate evidence (m10) |

---

## 2. Critical

### C1 — Pipeline was not reproducible from a clean checkout (fixed, not yet committed)

* **Evidence.** The first clean-clone run failed:
  `YEBE00ESP_R_20240060000_01D_30S_MO.rnx: first line is not 'RINEX VERSION / TYPE'`.
* **Cause.** Parallel workers processing different baselines of the same
  station-day decompress the same RINEX file into the same work directory;
  `prepare_rinex` wrote in place and accepted any existing file as complete,
  so a worker could read a partially written file. The original runs passed
  only because earlier sequential runs had left complete files behind.
  Whether a run fails is timing-dependent — a silent-corruption risk (a
  truncated-but-parseable file would have produced a wrong baseline without
  an error).
* **Fix.** Atomic publication (`tempfile` + `os.replace`) in
  `gnss_adjust/gnss/rinex.py`; regression test
  `test_prepare_is_safe_under_concurrency` (15 MB file, 32 concurrent workers)
  **fails 3/3 on the old code and passes 3/3 on the new code**. With the fix,
  the clean reproduction succeeds and is numerically identical (V7).
* **Action before release.** Commit the fix (currently uncommitted at the
  user's instruction).

---

## 3. Major

### M1 — The "false-positive rate of 0.22" is a misidentification rate, caused by a non-identifiable pair in the weak geometry

* **Evidence.** In every high-FPR cell the injected outlier *was* detected
  but attributed to the wrong baseline (e.g. outlier on 6 → 3 flagged 21×,
  outlier on 3 → 6 flagged 14×). Under H₀ the procedures' size is correct
  (V1).
* **Reasoning.** α controls Type I errors (alarms under H₀), not Type III
  errors (wrong identification under Hₐ). The separability of two local tests
  is measured by the canonical correlations of
  P̄ᵢᵢ^{-1/2} P̄ᵢⱼ P̄ⱼⱼ^{-1/2}. In the 2-session network the largest is
  **1.000 for baselines 3 and 6**: YEBE is connected only by these two, so an
  outlier on either produces the same residual pattern and the same
  statistic — the hypotheses are not identifiable. Expected error rate:
  P(outlier on 3 or 6) × ½ = 1/6 ≈ 0.17; observed 0.16–0.18 (resampled).
  The 3-session networks have maximum canonical correlation 0.52–0.54.
* **Impact.** The r = 9 results do not show "degraded detection"; they show
  that the experiment included a geometrically impossible identification task.
  The FPR column mixes two different error types.
* **Correction.** Report Type I rate (H₀ only), probability of correct
  identification, and misidentification separately; publish the canonical
  correlations of every geometry; mark non-identifiable pairs instead of
  scoring them; replace "FPR" in Experiment B by these quantities.

### M2 — Experiments B and C used a single realisation of the real noise

* **Evidence.** Fixed-realisation F1 (DIA, r = 18, one outlier): 0.18 at
  10 mm, 0.78 at 15 mm; resampled-noise probability of correct
  identification on the same geometry: 0.38 and 0.73.
* **Reasoning.** With one noise realisation, "detection rate" is conditional
  on that realisation; trials differ only in the injected baseline/direction.
  The earlier statement that empirical detection "tracks theory" was weakly
  supported (it averaged over directions, not noise). V2 now establishes it
  properly.
* **Correction.** Present Experiment B/C as *conditional on the observed
  data*; base statistical claims on resampled-noise Monte Carlo on the real
  geometry (`power_resampled.csv`); keep the fixed-realisation runs as a
  real-data illustration.

### M3 — Quiet-day base network: post-hoc selection makes Experiment B a best case and the method ranking condition-dependent

* **Evidence.** The quiet days (001–003) were chosen as the injection base
  *after* the evaluation days were seen to contain anomalies (its realisation
  is unremarkable: χ² CDF 0.17). The same protocol on the anomaly-containing
  network (`injection_on_main_network.csv`, 30 trials/cell):

  | method | F1, 20 mm | F1, 50 mm | "false" flags / run |
  |---|---|---|---|
  | DIA-reject (quiet → anomalous) | 0.97 → 0.47 | 1.00 → 0.57 | 0 → 1.6–1.7 |
  | Huber | 0.83 → 0.58 | 1.00 → 0.73 | 0 → 0.7–0.8 |
  | Hampel | 0.83 → 0.63 | 1.00 → 0.80 | 0 → 0.5–0.6 |

* **Reasoning.** On the anomalous network most "false positives" are the
  natural day-006 inconsistencies, which F1 against injected labels counts
  as errors even though they are genuine data problems. Both effects are real:
  performance under realistic conditions is much lower than under near-H₀
  conditions, **and** the ranking reverses (robust estimators score higher F1
  because they flag fewer baselines).
* **Correction.** Report both conditions; do not present "DIA best" as a
  general conclusion; when natural anomalies exist, score only injected
  baselines (recall/identification) and report natural flags separately.

### M4 — Robust estimators and DIA are not compared at equal size

* **Evidence.** Under H₀, DIA flags something in 0.35–0.6 % of runs; Huber
  and Hampel in 0 % (V1). Robust "detection" is the heuristic weight < 0.5,
  not a test at level α.
* **Impact.** "Robust IRLS is less sensitive" (Experiment B) mixes a genuine
  effect (fixed-Q⁰ standardisation with local redundancy 1.1–2.2, E[t²] =
  r_b/3) with an arbitrary, more conservative decision threshold.
* **Correction.** Compare at matched size (calibrate the weight threshold to
  the same H₀ flag rate) or with ROC curves; keep the r_b/3 explanation as
  the mechanism, quantified separately.

### M5 — Experiment A: the flagged day-006 baselines reflect an elevation-dependent effect, and the second DIA rejection depends on unknown cross-covariance

* **Evidence (processing sensitivity, 252 reprocessed solutions,
  `sensitivity_summary.csv`).** Implied VILL height deviation on days
  005 / 006 and pooled k² of days 004–006:

  | variant | VILL Up 005 / 006 [mm] | k² (004–006) | CEBR East − EPN [mm] |
  |---|---|---|---|
  | main solution | −11.8 / +20.4 | 5.17 | 12.0 |
  | broadcast orbits | −11.7 / +20.5 | 5.11 | 12.5 |
  | keep G03/G08 | −11.7 / +19.7 | 5.22 | 12.2 |
  | receiver PCV off | −11.9 / +20.2 | 5.13 | 11.9 |
  | ZTD + gradients | −9.1 / +17.5 | 4.53 | 13.4 |
  | 15° mask | −11.9 / +15.9 | 3.34 | 13.5 |
  | **20° mask** | **−1.2 / +4.2** | **1.47** | 12.0 |

  No receiver, antenna, firmware or data-completeness change occurs at VILL
  across the six days.
* **Reasoning.** The anomalies are insensitive to orbits, satellite
  screening, antenna model and troposphere gradients, but almost vanish when
  observations below 20° are removed (the reduction is in the deviations
  themselves, not only in their normalisation). The most plausible causes are
  low-elevation effects at or around VILL on those days (multipath or
  obstruction changes, wet-troposphere mapping, local environment) — not a
  gross baseline error and not an antenna-height-type error. Separately, a
  *hypothetical* shared station-session error model
  (`hypothetical_cross_covariance.json`; correlation ≈ ρ/2 between baselines
  sharing a station within a session, **not measured**) keeps baseline 8 as
  the strongest and the global test rejected, but baseline 9 drops to
  T = 3.0 (ρ = 0.5) / 1.5 (ρ = 1), while baseline 7 rises to 12.1 / 29.9.
* **Impact.** The earlier wording ("session-level vertical anomaly at VILL",
  flags "corroborated") over-interpreted the result. Robust conclusions:
  the day-006 data are inconsistent with the stochastic model, and baseline
  8 is the most inconsistent observation. The identity of the *second*
  rejected baseline is not robust.
* **Correction.** Reword Experiment A accordingly; report the elevation-mask
  dependence as a finding; treat the day-006 flags as diagnosing a modelling
  deficiency (low-elevation stochastic/functional model) rather than
  outliers.

### M6 — The covariance calibration is non-stationary for a physical reason; formal σ is not accuracy

* **Evidence (`calibration_decomposition.json`).** k² = 0.73 (001–003) vs
  5.17 (004–006), variance ratio F = 7.04 with 36/36 dof, p ≈ 5×10⁻⁸
  (under the model's own independence assumptions). By component:
  E 0.77 → 0.97, N 0.33 → 2.76, U 1.23 → 12.46. Leaving out day 006 lowers
  the all-days factor from ≈ 3 to 1.37; without VILL the evaluation factor is
  still 2.59. With a 20° mask the factors become 0.27 / 1.47.
* **Reasoning.** The change is significant, concentrated in the vertical,
  and largely explained by the low-elevation effect of M5 (days 005/006,
  VILL and YEBE). The formal iono-free covariance is a reasonable *precision*
  model on quiet days (k² CI 0.48–1.24) but its shape is wrong (N
  pessimistic, U optimistic) and it ignores elevation-dependent
  mismodelling; fixed-ambiguity covariances are 10–25× optimistic in σ.
  Repeatability is blind to constant errors (M7).
* **Correction.** Do not present a single calibrated scalar as "the"
  stochastic model; present formal σ as conditional precision; recommend
  component-wise (and elevation-aware) variance estimation as research work.

### M7 — CEBR ≈ 12 mm East offset: located in the processing chain, still unexplained

* **Evidence.**
  1. *Reference consistency.* Three EPN analysis centres that processed CEBR
     (BEK, BKG, IGE; `epn_ac_comparison.csv`) agree with the combined
     solution on CEBR−MADR to within ±0.7 mm (E), ±1.2 mm (N). (These
     centres may share processing software, so this shows consistency of
     the reference, not software independence.) The EPN constraint pulls CEBR's
     estimate only 0.5 mm in East from its a priori.
  2. *Antenna model.* CEBR (serial 5661) has no individual calibration in
     `epn_20.atx`; the SEPCHOKE_B3E6 type-mean entry is identical in
     `igs20_2287.atx` and `epn_20.atx` (`antex_comparison.json`). Receiver
     PCV off changes the offset by 0.1 mm.
  3. *Processing sensitivity.* Offset unchanged (11.9–13.5 mm) under broadcast
     orbits, 15°/20° masks, gradients, satellite selection (M5 table).
  4. *Signal dependence.* Using L2L instead of L2W at CEBR moves CEBR East by
     +4.4 mm on average (9.1 → 13.5 mm; `cebr_l2_signal_test.csv`).
* **Reasoning.** Frame, epoch, orbit, elevation and antenna-calibration
  explanations are excluded at the mm level; the offset is constant over six
  days (cancels in loops) and depends on the tracked signal. The most likely
  locus is RTKLIB's single-baseline modelling of this receiver/antenna's
  observations (code/phase biases, signal-specific phase-centre behaviour, or
  multipath) — but this is not demonstrated.
* **Impact.** CEBR coordinates carry a ~1 cm systematic error not reflected
  in any formal σ and undetectable by any of the tests. Accuracy statements
  must exclude or flag CEBR.
* **Correction.** Process CEBR with an independent engine (e.g. PPP with
  GipsyX/PRIDE, or Bernese if available); test other signals and a longer
  span; until then, report as "systematic, processing-related, unexplained".

### M8 — Independent redundancy is smaller than r = 18 suggests

* **Evidence.** Algebraic redundancy is verified (rank 9, 6 cross-session
  loops, no within-session loop). But: (i) a station-level error constant over
  all sessions (CEBR) has **zero** redundancy — it cancels in every loop;
  (ii) the day-006 effect acts on several baselines of one session at once;
  (iii) the lag-1 correlation of daily deviations is −0.31 (expected −0.20
  for independent demeaned data, 90 pairs): no evidence for, but no strong
  evidence against, day-to-day dependence; (iv) if within-session correlation
  exists but is ignored, the block-diagonal tests become **conservative**
  (hypothetical ρ = 0.5 / 1: global size 0.029 / 0.025, local 0.0005 /
  0.0003).
* **Correction.** State redundancy as "18 for session-varying, uncorrelated
  errors"; state explicitly which error classes have no redundancy.

---

## 4. Minor

| # | finding | recommendation |
|---|---|---|
| m1 | RTKLIB binary SHA-256 differs between builds (debug info embeds the build path) | identify the build by source commit + patch hash; do not present binary hashes as reproducible |
| m2 | Python dependencies only lower-bounded | ship the recorded `environment_pip_freeze.txt` as a constraints file |
| m3 | The EPN SINEX does not state its frame explicitly; EPN used ANTEX `IGS20_2290`, the pipeline `igs20_2287` | cite the EPN product documentation for the frame; check the remaining antenna types for calibration differences between releases |
| m4 | CEBR is a constrained (code 1) datum station in the EPN solution | mention; constraint pull is 0.5 mm East, so it does not explain M7 |
| m5 | "Predicted direction uncertainty is realistic (median ratio 0.77)" uses the largest perpendicular σ, i.e. a conservative bound | test angular coverage with resampled noise before calling it calibrated |
| m6 | The real-network identity in Experiment C used deletion with the same solver | superseded by V3 |
| m7 | Absolute-ECEF parameterisation: ~10⁻⁸ m residual and ~10⁻⁶ relative statistic precision | optional: reduce observations to approximate coordinates |
| m8 | Network design (trees, datum, days) hard-coded in `build_real_network.py` | move to the dataset config |
| m14 | Reference coordinates also enter the processing (EPN base positions held fixed in RTKLIB) | effect is second order (≈ L·\|δX\|/20 000 km, < 0.1 mm here); mention for completeness |
| m9 | Robust variance factor counts zero-weight baselines in dof (Stage 0 audit) | still open; affects reported σ̂₀² of Hampel only |
| m10 | "SD and DIA made identical decisions in all 3360 trials" — expected from V3; in the resampled runs one 50 mm r = 9 cell had 0.997 agreement, a numerical tie in the non-identifiable pair | describe as a consequence of the identity, not as empirical support |
| m11 | Experiment D detection thresholds (2–3 cm) are conditional on one noise realisation | resample, or state as conditional |
| m12 | Reproduction needs github.com and igs.bkg.bund.de | document; data integrity is guaranteed by checksums |
| m13 | Experiment B's "FPR ≤ 0.006 … max 0.22" sentence (already corrected once) still mixes error types | replace per M1 |

---

## 5. What the project establishes — classified

1. **Verified mathematical identities:** for block-diagonal Q_ℓ, SD
   \|w*\|² = 3D outlier test = likelihood-ratio statistic (synthetic and real
   geometry); under correlated Q_ℓ only the general form (= SD formula) is the
   correct 3D test; datum invariance of residuals and tests; frame invariance
   of all quadratic forms.
2. **Validated engineering functionality:** RINEX→RTKLIB→baseline chain with
   verified covariance reconstruction, antenna-metadata precedence,
   satellite screening, provenance; clean-checkout reproducibility (after C1);
   truth-free adjustment CLI.
3. **Supported by synthetic / resampled experiments:** correct test size and
   non-central-χ² power on real geometry; identifiability failure of
   near-leaf stations; conservativeness of block-diagonal tests under
   positive within-session correlation (hypothetical model); SD direction
   error ∝ 1/|d|; prespecified-direction power gain (χ²(1) vs χ²(3)).
4. **Supported by real GNSS measurements:** stock RTKLIB output corrupts
   3×3 covariances; RTKLIB processing pitfalls (status-based static
   selection, unhealthy-satellite use, silent epoch loss); formal vs empirical
   precision (float realistic horizontally, optimistic vertically; fixed
   10–25× optimistic); non-stationary, elevation-dependent inconsistency on
   days 005–006; a constant ~12 mm CEBR offset absent from all formal σ and
   invisible to all tests; antenna-height errors of 2–3 cm detectable and
   their direction recovered as vertical (conditional on one realisation).
5. **Tentative interpretations:** low-elevation effects at VILL as the cause
   of the day-005/006 inconsistency; a processing/receiver-signal origin of
   the CEBR offset; the second DIA rejection in Experiment A.
6. **Possible contributions not yet demonstrated:** station-session bias
   hypotheses (see `NEXT_RESEARCH_DIRECTION.md`); elevation-aware or
   component-wise stochastic models; any detection advantage of SD (there is
   none to demonstrate — SD's contribution is directional *interpretation*,
   and its optimised direction cannot beat the equivalent 3D test).

---

## 6. Readiness verdict

**GO WITH CAVEATS.** The mathematics, the statistics and the engineering
pass independent checks, and the real-data findings are genuine and
interesting. The work should not be presented with the earlier
interpretations. Prioritised corrections before presenting:

1. Commit the C1 fix (reproducibility).
2. Replace Experiment B's FPR/F1 summary by Type I size, correct
   identification and misidentification from resampled noise; disclose the
   non-identifiable pair (M1, M2).
3. Present both base conditions (quiet vs anomalous) and drop "DIA best" as a
   general claim; size-match robust thresholds or show ROC (M3, M4).
4. Reword Experiment A: elevation-dependent model inconsistency on day 006;
   only baseline 8 is a robust identification (M5).
5. Present formal σ as conditional precision; state CEBR as an unexplained,
   processing-related systematic and the zero redundancy for constant
   station errors (M6–M8).
