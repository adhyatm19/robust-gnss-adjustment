# Next Research Direction: Station-Session Shared-Bias Hypotheses

Evaluation of the proposal in `RESEARCH_HANDOFF.md` §8, after the scientific
audit (`SCIENTIFIC_AUDIT.md`). **Nothing here is implemented yet**; the only
numbers are a feasibility calculation of identifiability on the existing
Madrid design (§3).

## 1. Why consider it

The real data did not behave like isolated baseline blunders:

* the day-005/006 inconsistencies affect **several baselines of one station
  in one session**, and are elevation-dependent (audit M5);
* a deliberate antenna-height error at VILL (Experiment D2) was rejected as
  **two unrelated baselines** by DIA/SD;
* a station error that is **constant over sessions** (CEBR, audit M7)
  cancels in every loop and is undetectable by any hypothesis of this kind.

Baseline-wise DIA, SD and robust estimation all assume that each
observation (baseline) fails on its own. A hypothesis that matches the
physical cause — one station misbehaving during one session — may diagnose
better. Whether it does is an empirical question with clear prerequisites.

## 2. Mathematical formulation

Observation model with a station-session bias ∇ ∈ ℝ³ for station k in
session s:

$$
\ell = A x + C_{ks}\,\nabla + e,\qquad e\sim N(0,\sigma_0^2 Q_\ell),
$$

where the 3m×3 matrix $C_{ks}$ has, for every baseline i of session s,
block $+I_3$ if k is the *to*-station of i, $-I_3$ if k is the *from*-station,
and 0 otherwise (all other sessions: 0). Physically ∇ is a session-specific
displacement of station k's antenna phase centre (wrong height, local
multipath/obstruction, setup error).

* **Test statistic** (Baarda; Teunissen 2000), valid for any positive-definite
  $Q_\ell$:
  $$T_{ks}=(C^\top P v)^\top (C^\top P Q_v P C)^{-1}(C^\top P v)\sim\chi^2_q,$$
  with $q=\operatorname{rank}(C^\top P Q_v P C)\le 3$; the estimate is
  $\hat\nabla=(C^\top P Q_v P C)^{-1}C^\top P(-v)$ with covariance
  $\sigma_0^2(C^\top PQ_vPC)^{-1}$ — exactly the SD algebra with $C_{ks}$
  replacing a single-baseline selector, so the SD direction
  $\hat\nabla/\|\hat\nabla\|$ carries over unchanged.
* **Prespecified direction** (e.g. local up u for an antenna-height
  hypothesis): $C_{ks}u$, a 1-dof test with more power *if* the direction is
  right (audit V2/Experiment C).
* **Equivalent view.** Adding ∇ is the same as giving station k a separate
  coordinate for session s. The test is therefore a **congruence (stability)
  test** of station k between session s and the rest, as in deformation
  analysis. This is established theory; the research contribution would be
  its *application and validation* for diagnosing GNSS baseline networks, not
  a new statistic.
* **Random-effects counterpart.** If station-session errors are random rather
  than fixed, the same $C_{ks}$ define covariance components
  $Q_\ell=Q_0+\sum \sigma_{ks}^2 C_{ks}C_{ks}^\top$ — the hypothetical
  cross-covariance model used in the audit (M5, M8). Bias testing and
  variance-component estimation are the two ends of the same model and
  should be studied together.

## 3. Identifiability requirements

1. **Rank.** $C_{ks}$ must not lie in range(A); otherwise $Q_vPC=0$. This
   requires station k to be observed in **at least one other session** (a
   station observed in a single session can absorb any bias in its
   coordinates). Holds for all stations of the Madrid design (rank 3 for all
   12 station-sessions of each network, including the datum station).
2. **Separability from baseline hypotheses.** If station k has only **one**
   baseline in session s, $C_{ks}$ equals that baseline's outlier selector:
   canonical correlation 1.000 — the two explanations are indistinguishable.
   With two baselines, the Madrid design gives canonical correlations
   0.64–0.69 with the nearest single-baseline hypothesis: separable, but only
   for larger biases. In a 4-station spanning-tree session at most two
   stations have two baselines, so **in each 3-session network only 6 of the
   12 station-sessions can be discriminated from baseline outliers, and only
   weakly.**
3. **Separability between stations.** Two stations whose only connection in
   session s is one baseline produce indistinguishable signatures there; the
   canonical-correlation matrix across all candidate hypotheses (baseline and
   station-session) is the design criterion, as for the near-leaf failure in
   the audit (M1).
4. **Constant station errors** (same in all sessions) remain unidentifiable
   by construction (they are coordinate changes); only external information
   (reference coordinates, other software) can reveal them.
5. **Stochastic model.** The test is only as good as $Q_\ell$. The audit
   shows elevation-dependent inconsistency and non-stationary variance
   factors; a station-session test run on a mis-specified $Q_\ell$ will
   confound "bias" with "variance misfit". Elevation-aware weighting or a
   20° mask should be part of any evaluation.

## 4. How to test whether it improves diagnosis

**Design.** More stations per session so that most stations have ≥ 2
baselines per session (≥ 5–6 stations; non-tree session designs only if the
trivial-baseline problem is handled by session-level processing). Report the
full canonical-correlation matrix of all hypotheses *before* running
experiments.

**Competing procedures** (all at matched Type I size, calibrated by H₀
simulation on the real geometry):
(a) baseline-wise DIA = SD;
(b) station-session tests alone;
(c) a combined procedure choosing among both families (e.g. maximum of
size-normalised statistics, or an information criterion), with an explicit
"undecidable" outcome when candidates are not separable.

**Ground truth, in increasing realism.**
1. Resampled-noise Monte Carlo on the real geometry with injected
   station-session biases and single-baseline biases (known labels, both
   kinds mixed).
2. Real RINEX reprocessed with known station-level perturbations
   (antenna height, horizontal antenna offset) — the Experiment D mechanism,
   which already exists and is verified.
3. Unmodified data (VILL days 005/006), interpreted only — no scores.

**Metrics.** Type I size; probability of correct identification of the
*cause* (station-session vs baseline, and which one); misidentification
matrix; minimal detectable / identifiable bias per hypothesis; bias in
coordinates after adaptation; direction error of $\hat\nabla$.

**Success criterion.** At matched size, (c) identifies the correct cause
substantially more often than (a) for station-session errors, without losing
accuracy for single-baseline errors, *and* this holds on reprocessed real
data — not only in simulation.

## 5. Is it the most promising next step?

**Promising, but not first.** The audit shows that the dominant problem on
real data is the **stochastic/functional model** (elevation-dependent
inconsistency, non-stationary and anisotropic variance factor, ~1 cm
undetectable station systematics), not the choice of outlier test. A
station-session test on a mis-specified covariance would mostly re-diagnose
the same model deficiency. Recommended order:

1. **Stochastic model first:** elevation-dependent weighting / mask study and
   component-wise variance estimation, evaluated by repeatability and loop
   misclosures (uses existing data and code).
2. **Then station-session hypotheses**, together with their random-effects
   (variance-component) counterpart, on a redesigned network with ≥ 5
   stations and verified identifiability, using reprocessed antenna
   perturbations as ground truth.

This sequence keeps the work within established estimation and testing
theory, gives the station-session idea a fair test, and addresses the
largest error source the real data revealed.
