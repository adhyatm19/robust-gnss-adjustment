# Mathematical and Statistical Notes for the Real-Data Pipeline

Conventions, transformations and testing assumptions introduced or made
explicit by the RINEX integration. The adjustment model itself is unchanged
and documented in [`MATHEMATICS.md`](MATHEMATICS.md); the results are in
[`REAL_DATA_VALIDATION.md`](REAL_DATA_VALIDATION.md).

## 1. From RTKLIB output to an observation block

RTKLIB (relative, static) estimates the rover position $X_r$ with the base
position $X_b^{0}$ held fixed. The repository observation is

$$
\ell_{br} = X_r - X_b^{0}, \qquad \ell_{br} + v = X_r - X_b \;\;(\text{model } l_{ij}=X_j-X_i),
$$

so **from = base, to = rover**. Its covariance block is RTKLIB's rover
covariance, i.e. $\operatorname{Cov}(\hat X_r \mid X_b = X_b^0)$ — and, for a
fixed solution, conditional on the integer ambiguities. It excludes the base
coordinate uncertainty and ambiguity uncertainty, and it is computed under
RTKLIB's own stochastic model (white measurement noise per epoch,
elevation-dependent weights), which ignores time correlation over the
session. It is therefore a **formal** covariance and is calibrated
empirically (§4), not assumed correct.

**Field conversion** (verified in `solution.c`): with printed
$(s_x,s_y,s_z,s_{xy},s_{yz},s_{zx})$,

$$
Q=\begin{bmatrix} s_x^2 & \operatorname{sgn}(s_{xy})s_{xy}^2 & \operatorname{sgn}(s_{zx})s_{zx}^2\\
\cdot & s_y^2 & \operatorname{sgn}(s_{yz})s_{yz}^2\\ \cdot&\cdot& s_z^2\end{bmatrix}.
$$

**Validity is checked, not assumed**: symmetry (by construction), positive
definiteness (Cholesky/eigenvalues), condition number, and *rounding
sensitivity*: every printed value is only known to ±r/2 (r = print
resolution); sampling that box shows how far the ENU standard deviations and
correlations can move. With r = 0.1 mm and |ρ| ≈ 0.75 the set includes
matrices whose ENU correlations differ by > 0.3 and change sign
(`results/rtklib_validation/`). Hence the patched high-precision output.

**One solution per session.** RTKLIB's epoch outputs are a *cumulative*
filter: epoch $k$ uses all data up to $k$, so successive outputs are almost
perfectly correlated. Only the final forward epoch (the estimate from all data
of the static session) is used; treating epochs as independent observations
would overstate the information by orders of magnitude.

## 2. Frames

ECEF→local ENU at a point with geodetic latitude φ and longitude λ (WGS84):

$$
R=\begin{bmatrix}-\sin\lambda & \cos\lambda & 0\\ -\sin\varphi\cos\lambda & -\sin\varphi\sin\lambda & \cos\varphi\\
\cos\varphi\cos\lambda & \cos\varphi\sin\lambda & \sin\varphi\end{bmatrix},\qquad
d_{\mathrm{ENU}} = R\,d,\quad Q_{\mathrm{ENU}} = R\,Q\,R^\top .
$$

$R$ is orthonormal with det = +1, so every quadratic form is invariant:
$d_{\mathrm{ENU}}^\top Q_{\mathrm{ENU}}^{-1} d_{\mathrm{ENU}} = d^\top Q^{-1} d$.
**All tests are frame-invariant**; ENU is used only for interpretation
(residuals, directions, coordinate differences). The adjustment is done in
ECEF; residuals/directions of baseline $i\to j$ are rotated at station $i$.
Verified against RTKLIB's own ENU output to ≤ 3×10⁻⁹ m (signed-sqrt fields).

## 3. Datum and the fixed reference station

One station (MADR) is held at its EPN coordinate. Consequences:

* **Residuals, $Q_v$, the global test and every local test are invariant**
  to which station is fixed (minimal constraint) — verified in
  `tests/test_adjust_cli.py::test_datum_invariance_of_residuals_and_tests`.
* $Q_{\hat x}$ describes precision **relative to the datum station**, not
  absolute accuracy; the datum station's own uncertainty
  (EPN σ ≈ 0.05–0.17 mm) is not propagated.
* A bias *at the datum station* or a **consistent station-level bias**
  (same in all sessions) is not detectable by any of the tests: it is either
  absorbed by the datum or cancels around every loop. Comparisons with the
  reference therefore use translation-free differences
  $(\hat X_i-\hat X_{\text{datum}})-(X_i^{\text{ref}}-X_{\text{datum}}^{\text{ref}})$,
  which is how the CEBR East offset becomes visible.
* The fixed base coordinates in RTKLIB enter the baseline only at second
  order, $\delta b \approx L\,|\delta X_b| / \rho$ with $\rho\approx 2\times10^4$ km.

**Numerics.** Parameters are absolute ECEF coordinates (~6.4×10⁶ m). Double
precision then leaves ≈10⁻⁸ m in residuals and ≈10⁻⁶ relative in test
statistics (verified: coordinates agree with a local-origin solution to
10⁻⁷ m). Negligible for decisions; reducing observations to approximate
coordinates would remove it.

## 4. Covariance calibration (`gnss_adjust.gnss.calibration`)

For a pair observed in $D$ sessions with estimates $b_d$, formal $Q_d$ and
true covariance $k^2 Q_d$:

$$
\bar b = \Big(\sum P_d\Big)^{-1}\sum P_d b_d,\qquad
S=\sum_d (b_d-\bar b)^\top P_d (b_d-\bar b)\sim k^2\chi^2_{3(D-1)},
$$

pooled $\hat k^2=\sum S/\sum \mathrm{dof}$ with the exact χ² interval. The
estimator is unbiased (tested by simulation). It measures only **session-to-
session** variability: constant errors are invisible, and positive
correlation between sessions biases $\hat k^2$ low. It is computed on
calibration days (001–003) and applied to separate evaluation days to avoid
circularity — and on these data it does **not** transfer (0.73 vs 5.17).

## 5. Cross-baseline correlation

Three distinct things:

1. **Within-baseline (3×3) covariance** — provided by RTKLIB, used in full.
2. **Cross-covariance between baselines** — e.g. two baselines of one session
   sharing a receiver: both double differences contain that receiver's
   undifferenced errors, so with equal per-receiver noise the correlation of
   their observation errors is of order 0.5. Single-baseline processing does
   not compute it; recovering it would require multi-baseline (session)
   processing. **Not available, not invented**: $Q_\ell$ is block-diagonal.
3. **Correlation induced by common processing choices** (same base
   coordinates, orbits, antenna models, troposphere model) — systematic and
   not described by any covariance here.

*Algebraic* independence (no within-session loop: enforced by
`check_session_independence`) is necessary but not sufficient for
*stochastic* independence. If the true $Q_\ell$ had off-diagonal blocks, the
correct local test for baseline $i$ is the general form

$$
T_i=(C_i^\top P v)^\top\,(C_i^\top P Q_v P C_i)^{-1}\,(C_i^\top P v),
$$

which **equals the SD quantity** $g_i^\top\bar P_{ii}^{-1}g_i$ for any
positive-definite $Q_\ell$, whereas the DIA block statistic
$v_b^\top Q_{v,bb}^{-1}v_b$ equals it **only** for block-diagonal $Q_\ell$
(tested: > 5 % discrepancy under correlation,
`tests/test_sd_dia_independent.py`). Supporting a full $Q_\ell$ in the
solver is straightforward once a defensible model exists; the DIA module
would then need the general form.

## 6. Test statistics, distributions and what "SD adds"

With a-priori $\sigma_0^2$ (not the a-posteriori factor):

* global: $v^\top P v/\sigma_0^2 \sim \chi^2_r$, $r=n-\operatorname{rank}A$;
* local 3D (DIA): $T_b\sim\chi^2_3$; SD: $|w^*|^2 = g_i^\top\bar P_{ii}^{-1}g_i = T_b$
  (block-diagonal case) — verified against an independent likelihood-ratio
  computation $\Omega_0-\Omega_a$ with a free 3D bias, and the H₀
  distribution against χ²(3) by Kolmogorov–Smirnov.
* The SD *direction* $u^*=\hat d/\|\hat d\|$ is chosen by maximising over all
  directions using the same data, so $|w^*|^2$ follows χ²(3), **not** χ²(1).
  A direction specified *before* seeing the data (e.g. local up for an
  antenna-height hypothesis) gives $w_u = u^\top g/\sqrt{u^\top\bar P u}\sim N(0,1)$
  under H₀, with non-centrality $\lambda_u=(u^\top\bar P d)^2/(u^\top\bar P u)\le\lambda=d^\top\bar P d$,
  equality iff $u\parallel d$. Correctly prespecified → more power (χ²(1) at
  the same λ); wrongly prespecified → power lost. SD never has more power than
  the 3D test; its contribution is the direction estimate with covariance
  $\sigma_0^2\bar P_{ii}^{-1}$ (empirically realistic: median actual/predicted
  angular error 0.77).
* Sequential procedures (DIA, SD snooping) and maximum-over-baselines
  selection do not have nominal overall error rates (Teunissen 2018); α_local
  is a per-test level. Bonferroni is available but off by default.

## 7. Robust estimation on real geometry

IRLS standardises with the fixed $P^0_b$ (see `MATHEMATICS.md` §6). Under H₀,
$E[t_b^2]=\operatorname{tr}(Q_{v,bb}P^0_b)/3 = r_b/3$ (local redundancy
number $r_b$ ∈ [0,3]). In the Madrid networks $r_b$ = 1.1–2.2, so a Huber
constant $c$ acts like $c\sqrt{3/r_b}$ ≈ 1.2–1.6 $c$ in properly standardised
units: the robust estimators are systematically **less sensitive** than DIA
on sparse real networks, which Experiment B confirms.
