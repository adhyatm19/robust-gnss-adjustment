# Mathematical Documentation

This document collects the full mathematical model implemented in
`gnss_adjust`, with the exact conventions used in the code.

## 1. Observation model

Station $i$ has Cartesian coordinates $X_i = [X_i, Y_i, Z_i]^\top$ (metres).
A GNSS baseline observed **from** station $i$ **to** station $j$ is

$$
\ell_{ij} = X_j - X_i + e_{ij}, \qquad \operatorname{Cov}(e_{ij}) = Q_{ij}
\in \mathbb{R}^{3\times3},
$$

with a **full** covariance block $Q_{ij}$: the X/Y/Z components of a
processed baseline are correlated because they are estimated jointly from
the same carrier-phase data through a shared satellite geometry.

### Terminology

| symbol | meaning |
|---|---|
| $Q$ | covariance matrix (m²) |
| $Q/\sigma_0^2$ | cofactor matrix |
| $P = Q^{-1}$ | weight matrix |
| $\sigma_{0,\text{prior}}^2$ | a-priori variance factor (default 1) |
| $\hat\sigma_0^2$ | a-posteriori variance factor |

## 2. Linearised network model

With $m$ baselines and $u$ unknown stations, stack all observations
($n = 3m$ scalars):

$$
\ell + v = A x, \qquad Q_\ell = \operatorname{blockdiag}(Q_1,\dots,Q_m),
\qquad P = Q_\ell^{-1}.
$$

**Design matrix.** Each baseline contributes three rows. For a baseline
$i \to j$: $-I_3$ in the columns of station $i$ if unknown, $+I_3$ in the
columns of station $j$ if unknown. A fixed endpoint's coordinates are moved
into the observation vector:

- $i$ fixed: $\ell' = \ell + X_i^{\text{fixed}}$,
- $j$ fixed: $\ell' = \ell - X_j^{\text{fixed}}$.

**Datum.** The model is invariant under a common translation of all
stations, so without a fixed station $\operatorname{rank}(A) = 3(u{+}1) - 3$
and $N$ is singular. Fixing exactly one station removes the 3 translation
degrees of freedom. (Rotation and scale are already determined because
baseline *vectors*, not just distances, are observed.)

## 3. Weighted least squares

$$
N = A^\top P A, \quad u = A^\top P \ell, \quad \hat x = N^{-1} u
$$

solved via Cholesky factorisation of $N$ (`scipy.linalg.cho_factor` /
`cho_solve`), never through an explicit inverse. $Q_{\hat x} = N^{-1}$ is
produced by solving $N Q_{\hat x} = I$ with the same factorisation.

**Residuals** (sign convention used throughout code, tests and docs):

$$
v = A\hat x - \ell.
$$

$$
Q_v = Q_\ell - A Q_{\hat x} A^\top, \qquad
r = n - \operatorname{rank}(A), \qquad
\hat\sigma_0^2 = \frac{v^\top P v}{r}.
$$

$Q_v$ is symmetric positive semidefinite; its rank equals the redundancy
$r$. The **redundancy matrix** $R = Q_v P$ distributes total redundancy over
observations; a low local redundancy means an outlier there is hard to see.

**Standardised residuals**: $w_i = v_i / (\sigma_{0,\text{prior}}
\sqrt{(Q_v)_{ii}})$, which are $\mathcal N(0,1)$ under the null hypothesis.

## 4. DIA

### 4.1 Detection (global model test)

$$
T = \frac{v^\top P v}{\sigma_{0,\text{prior}}^2} \sim \chi^2(r)
\quad\text{under } H_0;
\qquad \text{reject if } T > \chi^2_{1-\alpha_g}(r).
$$

The **a-priori** factor is used deliberately: an outlier inflates the
a-posteriori $\hat\sigma_0^2$ and would mask itself if the test were
normalised by it.

### 4.2 Identification

- **Component w-test**: $|w_i| > \Phi^{-1}(1-\alpha_\ell/2)$, with optional
  Bonferroni correction $\alpha_\ell / (3m)$.
- **3D baseline group test**: for the residual block $v_b$ with covariance
  block $Q_{v,bb}$,

$$
T_b = \frac{v_b^\top Q_{v,bb}^{+} v_b}{\sigma_{0,\text{prior}}^2}
\sim \chi^2(\text{df}), \qquad
\text{df} = \operatorname{rank}(Q_{v,bb}) \ (\text{normally } 3),
$$

using the Moore–Penrose pseudoinverse and the numerical rank so the test
stays valid when a residual block is (nearly) constrained.

The identified baseline is the one with the largest ratio
$T_b / \chi^2_{1-\alpha_b}(\text{df})$ among those exceeding it.

### 4.3 Adaptation

Either the identified baseline is **removed**, or its covariance block is
**inflated** by a configurable factor (soft rejection). The loop re-adjusts
and repeats until: the global test passes; no candidate exceeds the local
threshold; the maximum number of adaptations is reached; or a removal would
disconnect the network or exhaust redundancy (**rank guard**). All decisions
are recorded in an audit trail.

## 5. Robust M-estimation (IRLS)

Instead of minimising $\sum_b v_b^\top P_b v_b$, a robust M-estimator
minimises $\sum_b \rho(t_b)$ for a slower-growing loss $\rho$. The IRLS
implementation reweights whole 3×3 blocks:

$$
t_b = \frac{1}{\sigma_{\text{scale}}}
\sqrt{\frac{v_b^\top P^0_b v_b}{3}}, \qquad
P^{\text{eff}}_b = w(t_b)\, P^0_b ,
$$

where $P^0_b$ is the **original** weight block. Multipliers are recomputed
from scratch each iteration (never compounded).

**Why the fixed reference $P^0_b$?** Standardising against the residual
covariance of the current *weighted* solve couples the statistic to the
weights: down-weighting inflates the observation's effective covariance,
shrinks its own statistic, and the weight oscillates (limit cycles observed
experimentally in this project). With a fixed reference, feedback runs only
through the coordinate solution — rejecting an outlier moves the solution
toward the clean data, growing the outlier's residual and shrinking the
others', a monotone separation. This matches classical equivalent-weight
practice (Danish method, IGG schemes). The trade-off — mild conservatism for
high-redundancy baselines since $\operatorname{Var}(v_b) \preceq Q^0_b$ — is
absorbed by the robust scale options.

### Weight functions

Huber ($c > 0$, default 1.5):

$$
w(t) = \begin{cases} 1 & t \le c \\ c/t & t > c \end{cases}
$$

Hampel ($0 < a < b < c$, defaults 1.5 / 3.5 / 8.0):

$$
w(t) = \begin{cases}
1 & 0 \le t \le a\\
a/t & a < t \le b\\
\dfrac{a\,(c-t)}{(c-b)\,t} & b < t \le c\\
0 & t > c
\end{cases}
$$

Huber bounds the influence of an outlier but never removes it; Hampel is
*redescending* — beyond $c$ the observation is rejected outright, so an
extreme outlier has exactly zero influence.

### Scale options

- `prior`: $\sigma_{\text{scale}} = \sigma_{0,\text{prior}}$;
- `mad`: median of the raw $t_b$ divided by the median of
  $\sqrt{\chi^2_3/3}$ (≈ 0.888), a MAD-style consistent estimate under
  Gaussian noise;
- `aposteriori`: $\sqrt{\hat\sigma_0^2}$ capped from above (an
  outlier-inflated factor must not mask the outlier) and floored from below.

### Safeguards

Zero residuals give $t = 0 \Rightarrow w = 1$; scales are floored at
`min_scale`; and a **rank guard** floors hard rejections at a tiny positive
weight whenever removing them would disconnect the network, so the normal
matrix stays positive definite while the baselines remain reported at
weight 0. Two mutually contradictory observations of a leaf station are
statistically undecidable — the estimator then reports both at weight 0 with
`converged` reflecting the actual outcome, rather than silently choosing one.

## 6. Synthetic covariance construction

Baseline sigma from the standard GNSS specification, with careful unit
handling ($a$ in mm, $b$ in ppm, $L$ and the result in metres):

$$
\sigma = \sqrt{ (a\cdot10^{-3})^2 + (b\cdot10^{-6} L)^2 } .
$$

Per-component sigmas $D = \operatorname{diag}(\sigma\,s_X, \sigma\,s_Y,
\sigma\,s_Z)$ (Z typically 1.5–2× worse), a random valid correlation matrix
$R$ (eigenvalue-floored and re-normalised), and

$$
Q_b = D\,R\,D \succ 0 .
$$

Student-t noise with $\nu$ dof is scaled by $\sqrt{(\nu-2)/\nu}$ so its
covariance still equals $Q_b$, isolating the effect of heavy tails from a
scale change.
