# Project Status

Honest classification of every feature. Nothing listed as completed is
planned-but-unwritten; nothing validated is asserted without a test or a
recorded experiment.

## Completed and validated (unit-tested and/or verified in recorded experiments)

- 3D baseline network model with full 3×3 covariance blocks
- Design-matrix construction with fixed-station reduction (sign and
  dimension tests; reversal-consistency test)
- Cholesky-based WLS solver: coordinates, Q_x, residuals, Q_v,
  standardized residuals, variance factors, condition number, dof
  (hand-computed example test; residual orthogonality A'Pv ≈ 0 test;
  noiseless-recovery test; Q_v symmetry/PSD test)
- Covariance construction Q = DRD, validation (symmetry, positive
  definiteness), a/b noise-model unit conversions (tested)
- Synthetic generator: connected redundant graphs, Gaussian and Student-t
  noise (covariance-matched), all outlier injection modes, contamination
  fractions, deterministic seeds (tested via fixtures and scenarios)
- DIA: global χ² detection, component w-tests (with Bonferroni option),
  3D block tests with numerical-rank df, reject and inflate adaptation,
  complete audit trail, rank guard (planted-outlier identification test,
  audit-trail test, bridge-baseline rank-guard test)
- Robust IRLS: Huber and Hampel weight functions (validated against closed
  forms), fixed-reference block standardisation, prior/MAD/a-posteriori
  scale options, non-compounded weights (tested), convergence histories,
  rank guard for hard rejections (tested, including the undecidable
  two-bridge conflict case)
- Metrics: coordinate RMSE / position errors, TP/FP/FN, precision, recall,
  F1, FPR, localization accuracy (used throughout recorded benchmarks)
- CSV input with full validation; CSV round-trip equals in-memory result
  (tested)
- CLI (`demo`, `benchmark` incl. sweeps); full demo pipeline produces all
  artifacts (integration-tested)
- Benchmark framework: 8 scenarios + magnitude/contamination sweeps,
  per-run CSV, aggregate CSV, JSON config/summary, PNG plots, Markdown
  report — all populated from real runs (docs/EXPERIMENTS.md)

## Demonstration only (works, but not exhaustively validated)

- Streamlit app (`app.py`): interactive layer over the tested API; the
  underlying computations are the tested ones, but the UI itself has no
  automated tests
- Plotting functions: exercised by the integration test (files produced,
  correct inputs) but visual correctness is human-reviewed only
- `weak-geometry` and `misscaled-covariance` scenarios: implemented and
  runnable; not part of the recorded headline experiments
- MAD and a-posteriori robust scale options: implemented with safeguards;
  the recorded experiments all use the a-priori scale

## Known behavioural caveats (documented, by design or observed)

- Hampel can reject a clean baseline under extreme (≈50σ) outliers in
  unfavourable geometries (initial swamping) — shown honestly in the
  magnitude sweep
- Sequential DIA degrades at high contamination (≥ 20–30%) due to masking
- DIA-inflate may end with status `no_candidate` while coordinates are
  fine (inflated block keeps the global statistic marginally high)
- Two mutually contradictory observations of a leaf station are
  statistically undecidable; the robust estimator reports both at weight 0
  rather than guessing

## Future work (NOT implemented)

- Real GNSS baseline datasets; RINEX preprocessing; ambiguity resolution
- Minimally constrained / free-network (inner-constraint) adjustment
- Cross-correlated simultaneously observed baselines (session-level Q_l)
- Recursive / sequential adjustment
- Specific-direction (1-D) outlier tests; outlier-pair hypotheses
- Integration with existing / indigenous GNSS processing software

## Explicit assumptions

1. Static network; coordinates do not change during the campaign.
2. Exactly one fixed station defines the datum (translations only).
3. Baseline vectors are already estimated upstream (by a baseline
   processor); this project starts from vectors + covariance blocks.
4. Covariance information is available and, by default, correctly scaled
   (σ₀,prior² = 1); a mis-scaled scenario exists for sensitivity checks.
5. Baselines are mutually uncorrelated (block-diagonal Q_l).
6. Synthetic experiments have exact ground truth; all detection metrics
   rely on it.
7. Results depend materially on tuning parameters (α levels, Huber c,
   Hampel a/b/c, scale option); defaults are conventional, not optimized.
