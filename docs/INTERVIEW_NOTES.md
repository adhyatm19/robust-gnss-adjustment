# Interview Notes

## 30-second explanation

"I built a Python framework that adjusts a 3D GNSS baseline network by
weighted least squares and then compares three ways of surviving bad data:
classical WLS, the geodetic DIA testing procedure (detect–identify–adapt),
and robust Huber/Hampel M-estimation. It runs on synthetic networks with
known ground truth, so I can measure exactly how much each method protects
coordinate accuracy and how reliably it localizes the corrupted baseline.
On a 30σ outlier, WLS coordinate error triples; DIA and Hampel recover
clean-data accuracy with F1 around 0.96."

## 2-minute explanation

Add to the above:

- The observation model: a baseline from i to j observes X_j − X_i plus
  correlated 3D noise with a full 3×3 covariance block. One station is
  fixed for the datum; everything is solved by Cholesky on the normal
  equations, never explicit inverses.
- DIA: global χ² test on v'Pv detects that *something* is wrong;
  per-baseline 3D χ² tests and component w-tests identify *which*
  baseline; adaptation either removes it or inflates its covariance, then
  re-adjusts, with a full audit trail and a rank guard so it can never
  silently disconnect the network.
- Robust: IRLS at the baseline-block level. Huber bounds an outlier's
  influence; Hampel is redescending and rejects extreme outliers entirely.
- Benchmarks: Monte Carlo over seeds, sweeps over outlier magnitude and
  contamination fraction, metrics for both coordinate accuracy and
  detection quality (precision/recall/F1). All numbers come from actually
  running the code; the whole pipeline is reproducible from recorded seeds.
- Honest scope: it consumes baseline vectors + covariances (what a baseline
  processor outputs), not raw RINEX/carrier phase.

## 8-minute technical walkthrough

1. **Model & datum** (1 min). l = Ax with ±I₃ blocks; fixed station moved
   into the observation vector; why exactly 3 datum defects (translations
   only — baseline vectors already fix rotation and scale).
2. **WLS internals** (1.5 min). N = A'PA via Cholesky; Q_x by solving
   against I; v = Ax̂ − l; Q_v = Q_l − AQ_xA'; redundancy r = n − rank(A);
   a-priori vs a-posteriori variance factor and why tests must use the
   a-priori one.
3. **DIA** (2 min). Global test T = v'Pv/σ₀² ~ χ²(r). Identification with
   the 3D block statistic T_b = v_b'Q_v,bb⁺v_b (pseudoinverse + numerical
   rank, since a residual block can be near-singular at low redundancy).
   Adaptation strategies, stop conditions, audit trail, rank guard demo:
   a zero-redundancy bridge baseline is *undetectable* — its residual is
   ~0 no matter how wrong it is — and must never be removed.
4. **Robust IRLS** (2 min). Weight functions and tuning constants; the
   fixed-reference standardisation t_b = √(v_b'P⁰_bv_b/3)/σ and the
   experiment that selected it: I tried three references and two of them
   fail in opposite directions. The effective Q_v of the current weighted
   solve gives a genuine limit cycle (down-weighting inflates the
   observation's own covariance, shrinks its statistic, and the weight
   climbs back — I recorded a period-8 cycle, 0 → 0.08 → 0.48 → … → 1 → 0);
   the original-weight residual covariance Q⁰ − AQ_xA' gives an avalanche
   (the reference under-states a down-weighted residual's dispersion, so
   clean baselines collapse — ten of fifteen dragged below 0.5 in my
   trace). The fixed observation weight P⁰ = (Q⁰)⁻¹ leaves feedback only
   through the solution: monotone separation, no collateral damage.
   Reproducible via `scripts/compare_standardizations.py`. Weights are
   recomputed from scratch each iteration, never compounded.
5. **Results** (1.5 min). Walk the three headline tables in
   EXPERIMENTS.md: clean-data efficiency (identical solutions), the
   directional outlier (WLS 2× worse; DIA exact localization; Huber
   partial, weight 0.32; Hampel stronger, 0.20), Monte Carlo (WLS 28 mm →
   ~9 mm for all protected methods).

## Key concept answers

**What is the design matrix here?** Pure ±I₃ incidence structure — the
model is linear, no linearisation needed. Station i's columns get −I₃,
station j's +I₃, fixed stations contribute to the reduced observation
vector instead. It is the 3D vector generalisation of a levelling network's
incidence matrix.

**Why are baseline components correlated?** All three components are
estimated simultaneously from the same carrier-phase observations through a
shared satellite geometry; the up component is systematically weaker
(geometry is one-sided, tropospheric coupling). Hence a full 3×3 block, and
why the pipeline reweights whole blocks rather than scalar components.

**Why does WLS fail under gross errors?** The quadratic loss gives every
observation influence proportional to its residual — unbounded. The LS
solution smears one gross error over all coordinates (influence function
argument; breakdown point 0).

**DIA vs robust estimation?** DIA is hypothesis testing: hard decisions
with controllable error rates (α), an interpretable audit trail, and the
correct framework when you must *report* which observation was bad. Robust
estimation makes soft, continuous decisions and degrades gracefully — better
with multiple/moderate outliers (no masking cascade), no combinatorial
search. They are complementary: robust estimation is often a good *starting
point* whose weights inform identification.

**Detection vs Identification vs Adaptation?** Detection: "is the whole
model consistent?" (one global χ² test). Identification: "which
observation?" (many local tests, multiplicity handled via Bonferroni
option). Adaptation: "what to do about it?" (remove or down-weight, then
re-test). Keeping the three separate is what makes the procedure auditable.

**Huber vs Hampel?** Huber is monotone: convex problem, unique solution,
never fully rejects — an extreme outlier keeps residual influence bounded
by c. Hampel is redescending: rejects outright beyond c, so extreme
outliers have zero influence, but the objective is non-convex — the result
can depend on the starting point, and an extreme outlier can initially
swamp its neighbours above the rejection threshold (I show this failure
mode honestly in the 50σ sweep entry).

**Meaning of the residual covariance Q_v?** The uncertainty of the
residuals themselves. Its diagonal is what standardises residuals into
w-statistics; Q_vP is the redundancy matrix telling you how visible an
error in each observation is. Low local redundancy ⇒ outliers hide there.

**Meaning of the global χ² test?** Under H₀ (correct model, correct
stochastic model), v'Pv/σ₀² is χ²(r). Exceeding the critical value says
*the model and the data disagree somewhere* — outliers, wrong covariances,
or a wrong functional model; it cannot tell you which, which is why
identification follows.

**Why must one station be fixed?** Baseline vectors are translation
invariant: adding a constant to every station leaves all X_j − X_i
unchanged, so N has a 3-dimensional null space. Fixing one station defines
the datum. (Rotation/scale are observed, unlike in distance-only networks.)

**Limitations of synthetic data?** Real baseline errors are not Gaussian +
occasional clean biases: multipath is time-correlated, half-fixed
ambiguities give discrete error patterns, sessions share satellites so
baselines are cross-correlated, and reported covariances from commercial
processors are notoriously optimistic. Synthetic benchmarks validate the
machinery and rank methods, but the tuning constants would need
re-validation on real data.

**Numerical-stability decisions?** Cholesky for the normal equations
(rejects non-PD early = built-in rank check); covariance only via solves
against I; pseudoinverse + numerical rank for residual blocks;
block-inverse of each 3×3 Q via Cholesky; eigenvalue flooring for random
correlation matrices; explicit condition-number reporting; rank guards in
both DIA and IRLS so no adaptation can silently destroy the network.

**Why doesn't IRLS standardise against Q_v like DIA does?** Because DIA
tests a *fixed* adjustment once, while IRLS feeds its statistic back into
the next solve — and any reference that moves with the weights (or ignores
how down-weighting changes the residual's dispersion) becomes part of the
dynamics. I have the failure traces: the effective-Q_v reference limit-cycles
(a rejected baseline's inflated covariance shrinks its own statistic and it
re-enters at full weight), and the original-weight residual covariance
Q⁰ − AQ_xA' avalanches (it under-states the dispersion of a down-weighted
residual, so t is overestimated and clean baselines collapse one after
another). The fixed observation weight P⁰ is the unique choice of the three
where feedback runs only through the coordinate solution; that is why the
classical equivalent-weight schemes (Danish, IGG) use exactly this form.
`scripts/compare_standardizations.py` reproduces all three behaviours
deterministically — worth running live if the interviewer pushes on it.

## Likely interviewer questions

- *"Your local tests use α = 0.001 — why so small?"* — With 3m component
  tests per adjustment, the family-wise false-alarm rate must stay usable;
  0.001 per test ≈ 5% family-wise for 48 tests, matching the global α. The
  Bonferroni option makes this explicit.
- *"What if two outliers mask each other?"* — Sequential DIA is greedy and
  can fail there (shown in my 30% contamination results). Mitigations:
  robust starting weights, or hypothesis tests for outlier *pairs* — listed
  as future work.
- *"Why block reweighting instead of per-component?"* — Components are
  correlated; a bias in a rotated direction leaks into all components.
  The Mahalanobis block statistic is invariant to that rotation.
- *"How would you handle correlated baselines from one session?"* — Q_l
  stops being block-diagonal; the WLS algebra is unchanged but weight
  assembly and the block tests must operate on session-level groups —
  a planned extension.
- *"Chi-square vs F test?"* — χ² assumes σ₀² known a priori (my default,
  covariances taken as calibrated). If σ₀² were estimated from an
  independent adjustment, ratio tests become F-distributed.

## "Was this integrated into the institute's indigenous GNSS software?"

Honest answer: **No — and it is not represented as such.** This is an
independent research prototype and benchmarking framework. Its purpose is
to validate and compare the adjustment and outlier-handling methods —
demonstrating correctness on data with known ground truth, quantifying the
accuracy/robustness trade-offs, and documenting the algorithms — *before*
any integration into production software. Integration would additionally
require: consuming the production baseline-processor output format,
handling session-correlated covariances, real-data validation of the tuning
constants, and regression tests against the existing adjustment results.
The clean module boundaries (solver / diagnostics / DIA / robust are
independent, typed components) were chosen with that integration in mind.
