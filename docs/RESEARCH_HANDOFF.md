# Research Handoff — Real-Data Extension

For: the student working on the mathematics of robust GNSS network
adjustment (National Centre for Geodesy, IIT Kanpur). Details and every number:
[`REAL_DATA_VALIDATION.md`](REAL_DATA_VALIDATION.md); conventions and
statistics: [`REAL_DATA_MATH_NOTES.md`](REAL_DATA_MATH_NOTES.md).

## 1. What was completed

* A reproducible **RINEX → RTKLIB → baseline + covariance → network** pipeline
  (`src/gnss_adjust/gnss/`): RINEX 2/3 metadata and overlap checks, pinned
  RTKLIB build (stock + output-precision variant), config-only invocation
  with full logs, `.pos` parsing with verified covariance reconstruction,
  ECEF/ENU transformations, satellite screening, SINEX reference coordinates,
  covariance calibration from repeatability, and network assembly that
  rejects trivial (same-session loop) baselines and reports loop misclosures.
* A truth-free **`python -m gnss_adjust adjust`** command that runs WLS,
  DIA-reject, DIA-inflate, Huber, Hampel and SD-reject on any CSV network and
  reports in local ENU, with optional reference coordinates.
* A **real redundant network**: 4 EPN/IGS stations near Madrid, 10–109 km,
  6 days of 24 h data; evaluation network 9 baselines, r = 18, 6 cross-session
  loops; plus quiet-day, weaker-geometry and deliberately-trivial variants.
* **Experiments A–D** (unmodified data, controlled injections, SD validation,
  reprocessed antenna-height errors), 9 figures, 101 automated tests.
* Repository hygiene (staged separately): generated `.pyc`, `egg-info`,
  `.DS_Store` untracked; local files kept.

## 2. What was verified (and how)

| claim | evidence |
|---|---|
| RTKLIB uncertainty fields are σ and signed √covariance, order xy, yz, zx | source code + reproduction of RTKLIB's own ENU output to 3×10⁻⁹ m |
| precision patch changes formatting only | 460 epoch pairs identical in status/metadata; all fields within ½ print resolution |
| stock 0.1 mm output corrupts 3×3 blocks | 22 % σ errors, 0.34 correlation errors, sign flips in 4.9 % of epochs |
| antenna heights are applied, with documented precedence | ±10 cm tests to ≤ 0.06 mm (GSI); Madrid reprocessing to ≤ 0.23 mm |
| SD statistic = 3D outlier test = likelihood-ratio test | independent dense-GLS code, block-diagonal: 10⁻⁷; real network: 2×10⁻⁶ |
| DIA block statistic ≠ correct test if baselines correlate | > 5 % discrepancy with a correlated Q_ℓ |
| statistics have their nominal null distributions | KS tests: local T ~ χ²(3), global ~ χ²(r), loop ~ χ²(3) |
| residuals and tests do not depend on the fixed station | datum-swap test |
| existing synthetic code unchanged | original 49 tests pass; demo reproduces recorded numbers to 10⁻¹² m |

## 3. What remains incomplete

* No cross-baseline covariance (needs multi-baseline/session processing).
* The CEBR ≈ 11 mm East offset relative to EPN and the day-005/006 VILL height
  anomalies are **unexplained**.
* One processor, GPS only, iono-free float ambiguities, four stations, one week.
* Covariance calibration is a single scalar; the ENU shape mismatch and the
  non-stationarity (k² 0.73 vs 5.17) are reported, not modelled.
* RCS (residual-corrected robust estimation) was not attempted.

## 4. What was actually tested with real data

* Processing of real RINEX (GSI 0759/3040; Madrid 4 stations × 6 days) and
  every verification in §2 involving RTKLIB.
* Formal vs empirical uncertainty, EPN comparison, loop misclosures.
* **Experiment A**: all six methods on the unmodified real network.
  Flags are *interpreted* against independent evidence, never scored.
* **Experiment D**: real RINEX reprocessed with a wrong antenna height
  (a physically meaningful, known error).

## 5. What is synthetic or only partly real

* **Experiments B and C** inject artificial bias vectors into real-data-derived
  baselines (known labels, artificial errors, one fixed noise realisation).
* `tests/test_sd_dia_independent.py`, `test_gnss_network_tools.py` and the
  original benchmark suite (`EXPERIMENTS.md`) are synthetic.

## 6. Commands

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[rinex,dev]"
pytest                                   # 101 tests (RTKLIB tests skip if not built)
scripts/build_rtklib.sh                  # pinned RTKLIB into tools/ (needs git, gcc, make)
scripts/reproduce_real_data.sh python    # everything: validation, download (126 MB),
                                         # processing, networks, Experiments A-D, figures
# single steps
python scripts/validate_rtklib_patch.py
python scripts/fetch_real_data.py configs/real_data/madrid_2024.yaml
python scripts/process_real_network.py configs/real_data/madrid_2024_processing.yaml
python scripts/build_real_network.py
python -m gnss_adjust adjust --network results/real_network/networks/main_formal \
    --reference results/real_network/networks/main_formal/reference_coordinates.csv \
    --output results/real_network/experiment_A/main_formal
python scripts/real_injection_experiments.py --trials 40
python scripts/real_antenna_experiment.py
python scripts/plot_real_results.py
```

## 7. Mathematical issues worth understanding

1. **Formal vs actual covariance.** RTKLIB's covariance is conditional (on
   the base, and on integers when fixed) and assumes white noise per epoch.
   Here it was realistic horizontally but optimistic vertically for float
   solutions, and 10–25× optimistic in σ for fixed ones. Every test, threshold
   and robust weight depends on this matrix — the stochastic model matters
   more than the choice of estimator.
2. **Calibration is itself a statistical problem.** Repeatability captures
   only session-varying errors; constant errors (CEBR) are invisible to it
   *and* to network testing; the factor is non-stationary; a scalar cannot fix
   anisotropy (variance component estimation per ENU component would).
3. **Algebraic vs stochastic independence.** No same-session loop is
   necessary for a well-posed network, not sufficient for a block-diagonal
   covariance. Under correlation only the general form
   $(C^\top Pv)^\top(C^\top PQ_vPC)^{-1}(C^\top Pv)$ is the right test — which
   is the SD formula, not the DIA residual-block formula.
4. **SD = 3D test in detection; the extra is the direction.** The
   optimised direction has a χ²(3) null; a prespecified direction has χ²(1)
   and more power only if it is right. The direction estimate has covariance
   $\sigma_0^2\bar P_{ii}^{-1}$, which was realistic in practice.
5. **Model mismatch.** The anomalies found in real data act on a *station in
   a session* (several baselines at once), not on one baseline. DIA rejected
   two baselines for one physical cause (Experiment D) and resolved a
   day-005/006 conflict by choosing a side (Experiment A).
6. **Robust estimators on sparse networks** are less sensitive than DIA
   because $E[t_b^2]=r_b/3$ with local redundancy $r_b$ ≈ 1–2.
7. **Detectability limits are geometric**: two baselines controlling a
   near-leaf station cannot be told apart (seen at r = 9).

## 8. Next research question for discussion with Prof. Dikshit

**Can a station-session bias hypothesis — a 3D (or prespecified vertical) bias
shared by all baselines of one station in one session, tested with Baarda's
general $C$-matrix statistic and estimated with its direction and
covariance — identify and explain the real anomalies better than
single-baseline DIA/SD, and how does its minimal detectable bias depend on the
session design?**

Why this one: it is directly motivated by what the real data showed (VILL
day 005/006; Experiment D, case D2), it uses only established testing theory
(alternative hypotheses through $C$, MDB/power from the non-centrality), and
the pipeline already provides everything needed (session labels, verified
statistics, reprocessing with known antenna-height errors as ground truth).
A second, complementary question is the stochastic model: component-wise
variance estimation and an explicit session-correlation model, since
calibration — not the estimator — dominated the real-data behaviour.

## 9. The most promising contribution emerging from this work

Not a new estimator. What the results genuinely support is a **validated,
reproducible real-data test bed** — with a verified covariance chain, labelled
physical errors (antenna-height reprocessing) and documented natural
anomalies — showing quantitatively that (a) the SD direction estimate gives
correct and reliably uncertainty-quantified *physical interpretation*
(reprocessed antenna-height errors: SD direction within ≈2–15° of the vertical once detected) without any detection
advantage over the equivalent 3D test, and (b) real anomalies are
station-session effects that single-baseline outlier models represent
poorly. Extending SD-style directional estimation to station-session
hypotheses (§8) is the defensible next step; it should be presented as an
application and validation of established theory, not as a new method.
