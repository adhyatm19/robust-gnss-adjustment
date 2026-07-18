# Robust GNSS Baseline Network Adjustment

**Weighted Least Squares, DIA Outlier Detection, and Robust M-Estimation**

A Python research prototype that adjusts a redundant 3D GNSS baseline network and
benchmarks how classical and robust estimation methods behave when observations
are contaminated by gross errors and non-Gaussian noise.

> **Scope statement (read this first).** This is a *synthetic, baseline-level*
> research and benchmarking prototype. Inputs are baseline vectors
> (`ΔX, ΔY, ΔZ`) with full 3×3 covariance blocks — the kind of quantities a
> GNSS baseline processor produces *after* carrier-phase processing. It is
> **not** a raw RINEX processor, **not** a carrier-phase ambiguity-resolution
> engine, **not** a replacement for commercial GNSS software, and it is **not
> yet integrated** into any production or indigenous GNSS software. It is an
> independent framework for validating adjustment and outlier-handling methods
> before such integration.

---

## 1. Motivation

GNSS baseline networks are adjusted by weighted least squares (WLS), which is
optimal under Gaussian noise but has **zero breakdown point**: a single gross
error (wrong antenna height, half-fixed ambiguities, multipath-corrupted
session) contaminates *every* estimated coordinate, because least squares
smears the error over the whole network. Two families of countermeasures
exist:

* **Statistical testing (DIA)** — detect that something is wrong (global
  chi-square test), identify which observation is wrong (local w-tests /
  3D block tests), and adapt (remove or down-weight it), iteratively.
* **Robust M-estimation** — replace the quadratic loss with one that grows
  more slowly for large residuals, implemented as iteratively reweighted
  least squares (IRLS) with Huber or Hampel weight functions.

This project implements both on the same network model and measures — with
known synthetic ground truth — coordinate accuracy, detection quality, and
clean-data efficiency.

## 2. Mathematical model

Station `i` has unknown coordinates `X_i = [X, Y, Z]ᵀ`. A baseline observed
from `i` to `j` is

```
l_ij = X_j − X_i + e_ij,     Cov(e_ij) = Q_ij  (full 3×3)
```

Stacking all `m` baselines (`n = 3m` scalar observations) gives the linear
Gauss–Markov model

```
l + v = A x,      Q_l = blockdiag(Q_1, …, Q_m),      P = Q_l⁻¹
```

One station is held **fixed** to define the datum; its contribution moves into
the reduced observation vector. Each baseline contributes three rows to `A`:
`−I₃` in the columns of an unknown from-station and `+I₃` in the columns of an
unknown to-station.

### Weighted least squares

```
N = AᵀPA,   u = AᵀPl,   x̂ = N⁻¹u        (solved by Cholesky, never an explicit inverse)

v   = A x̂ − l                            (residual sign convention used everywhere)
Q_x = N⁻¹                                (obtained by solving N Q_x = I)
Q_v = Q_l − A Q_x Aᵀ
r   = n − rank(A)                        (degrees of freedom)
σ̂₀² = vᵀPv / r                           (a-posteriori variance factor)
```

The solver reports network connectivity, `rank(A)`, `cond(N)`, symmetry and
positive (semi)definiteness of the covariance matrices, degrees of freedom,
and the datum station.

### DIA — Detection, Identification, Adaptation

* **Detection** (global test): `T = vᵀPv / σ₀,prior²  ~  χ²(r)` under H₀;
  reject when `T > χ²₁₋α(r)`. The **a-priori** variance factor is used —
  the a-posteriori one is inflated by the very outlier being tested.
* **Identification**: component w-tests
  `w_i = v_i / (σ₀,prior √(Q_v)_ii)` compared with `N(0,1)` quantiles
  (optionally Bonferroni-corrected), and per-baseline 3D group tests
  `T_b = v_bᵀ Q_v,bb⁺ v_b / σ₀,prior²  ~  χ²(df)` with `df` = numerical rank
  of the residual block (normally 3).
* **Adaptation**: either **reject** the identified baseline or **inflate**
  its covariance block, then re-adjust; loop until the global test passes,
  no candidate exceeds the local threshold, a maximum number of adaptations
  is reached, or removal would disconnect the network / exhaust redundancy
  (rank guard). A complete audit trail is kept; nothing is deleted silently.

### Robust M-estimation (IRLS on 3D blocks)

Each iteration computes a block-standardised Mahalanobis measure against the
*fixed original* weight block `P0_b = Q0_b⁻¹`,

```
t_b = √( v_bᵀ P0_b v_b / 3 ) / σ_scale
```

and multiplies the baseline's original weight block by `w(t_b)` (multipliers
are recomputed from scratch each iteration — never compounded). The fixed
reference is deliberate: standardising against the residual covariance of the
current *weighted* solve couples the statistic to the weights and produces
IRLS limit cycles, whereas a fixed reference feeds back only through the
coordinate solution and separates outliers monotonically (classic
equivalent-weight practice; the rigorous `Q_v`-based statistics remain in the
DIA pipeline):

**Huber** (monotone, bounded influence):

```
w(t) = 1        t ≤ c
w(t) = c/t      t > c            (default c = 1.5)
```

**Hampel** (redescending, three-part, 0 < a < b < c):

```
w(t) = 1                          0 ≤ t ≤ a
w(t) = a/t                        a < t ≤ b
w(t) = a(c−t) / ((c−b) t)         b < t ≤ c
w(t) = 0                          t > c          (defaults a=1.5, b=3.5, c=8)
```

Robust scale options: known a-priori σ₀, MAD-based scale (with the χ(3)
consistency constant), or a safeguarded a-posteriori scale.

## 3. Synthetic data generation

* Random true coordinates in a configurable box; random spanning tree plus
  extra edges → connected, redundant graph; one fixed station.
* Per-baseline sigma from the classic model `σ = √(a_mm² + (b_ppm·L)²)`,
  with per-axis scale factors (Z noisier) and a random valid correlation
  matrix: `Q_b = D R D`.
* Noise: multivariate Gaussian or Student-t (heavy-tailed, covariance-matched).
* Gross errors: single component, explicit 3D vector, random/fixed direction
  with magnitude in metres or sigma units, multiple baselines, or a
  contamination fraction. Reported covariances can be deliberately mis-scaled.
* Every experiment records its seed and full configuration.

## 4. Benchmark metrics

Coordinate: RMSE, mean/max 3D station position error, per-station error.
Detection: TP/FP/FN, precision, recall, F1, false-positive rate, exact
localization. Numerical: iteration counts, convergence rate, runtime,
`cond(N)`, final variance factor, number of flagged baselines. Clean-data
efficiency: robust-vs-WLS RMSE with no outliers. **All reported numbers are
produced by running the code** (see `scripts/generate_demo_results.py`).

## 5. Installation

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[app,dev]"          # package + streamlit + pytest
pytest                                # 35 deterministic tests
```

Requires Python ≥ 3.11. Core dependencies: NumPy, SciPy, Pandas, Matplotlib,
NetworkX, PyYAML; Streamlit only for the demo app.

## 6. Command-line interface

```bash
# one experiment, all methods, plots + report
python -m gnss_adjust demo --seed 42 --output results/demo

# directional 3D outlier on baseline 3
python -m gnss_adjust demo --seed 42 --method all \
    --outlier-baseline 3 --outlier-vector 0.035 -0.025 0.090 \
    --output results/directional_outlier

# Monte Carlo benchmark
python -m gnss_adjust benchmark --runs 200 --scenario single-outlier \
    --output results/benchmark

# sweeps
python -m gnss_adjust benchmark --scenario magnitude-sweep --runs 80 --output results/mag
python -m gnss_adjust benchmark --scenario contamination-sweep --runs 50 --output results/cont
```

Scenarios: `clean`, `single-small-outlier`, `single-outlier`,
`directional-outlier`, `multiple-outliers`, `student-t`, `weak-geometry`,
`misscaled-covariance`, plus `magnitude-sweep` and `contamination-sweep`.

## 7. Streamlit app

```bash
streamlit run app.py
```

Interactive control of seed, network size/density, noise level, outlier
placement (vector or sigma-magnitude), significance levels, and robust tuning
constants; displays the network diagram, method comparison, residual
diagnostics, robust weights, and the full DIA audit trail.

## 8. Example output

`python -m gnss_adjust demo --seed 42 --outlier-baseline 3 --outlier-vector
0.035 -0.025 0.090 --output results/directional_outlier` produces
`method_comparison.csv`, `summary.json`, `report.md`, the network CSVs, and
eight PNG figures (network map, true-vs-estimated, coordinate errors,
residual norms, w-tests, group statistics, robust weights). Actual numbers
from the committed runs are in [docs/EXPERIMENTS.md](docs/EXPERIMENTS.md).

## 9. Repository structure

```
robust-gnss-adjustment/
├── app.py                     # Streamlit demo
├── configs/                   # YAML experiment configs
├── src/gnss_adjust/
│   ├── types.py               # dataclasses & conventions
│   ├── network.py             # design matrix, datum, connectivity
│   ├── covariance.py          # Q = DRD construction & validation
│   ├── simulation.py          # synthetic generator & outlier injection
│   ├── wls.py                 # Cholesky-based WLS solver
│   ├── diagnostics.py         # global test, w-tests, 3D block tests
│   ├── dia.py                 # DIA loop with audit trail & rank guard
│   ├── robust.py              # Huber/Hampel IRLS
│   ├── metrics.py             # coordinate & detection metrics
│   ├── benchmark.py           # scenarios, Monte Carlo, sweeps
│   ├── plotting.py            # matplotlib figures
│   ├── io.py                  # CSV/JSON/YAML I/O & validation
│   └── cli.py                 # argparse CLI
├── tests/                     # 35 deterministic pytest tests
├── examples/                  # CSV input format examples
├── scripts/generate_demo_results.py
└── docs/                      # MATHEMATICS, EXPERIMENTS, INTERVIEW_NOTES, PROJECT_STATUS
```

## 10. Input CSV format

`stations.csv`: `station_id,x,y,z,is_fixed` — `baselines.csv`:
`baseline_id,from_station,to_station,dx,dy,dz` — `covariance_blocks.csv`:
`baseline_id,q_xx,q_xy,q_xz,q_yy,q_yz,q_zz` (metres / m²). Loading validates
station references, duplicates, covariance symmetry/positive-definiteness,
and graph connectivity.

## 11. Limitations

* Baseline vectors are treated as *given* — no raw observation processing,
  no ambiguity resolution, no troposphere/ionosphere modelling.
* Baselines are assumed mutually uncorrelated (block-diagonal `Q_l`);
  simultaneously observed sessions violate this.
* Static network, single fixed-station datum (no free-network /
  minimally-constrained variants yet).
* Synthetic ground truth only; tuning constants and significance levels
  materially influence results and would need re-validation on real data.

## 12. Future work

Real GNSS baseline datasets; RINEX preprocessing; ambiguity resolution;
minimally constrained and free-network (inner-constraint) adjustment;
correlated simultaneously-observed baselines; recursive/sequential
adjustment; specific-direction (one-dimensional) outlier tests; integration
with existing GNSS processing software.

## 13. References

1. Baarda, W. (1968). *A Testing Procedure for Use in Geodetic Networks.*
   Netherlands Geodetic Commission, Publ. on Geodesy 2(5).
2. Teunissen, P. J. G. (2000). *Testing Theory: an Introduction.* Delft
   University Press.
3. Teunissen, P. J. G. (2018). Distributional theory for the DIA method.
   *Journal of Geodesy* 92, 59–80.
4. Huber, P. J. (1964). Robust estimation of a location parameter.
   *Annals of Mathematical Statistics* 35, 73–101.
5. Hampel, F. R. et al. (1986). *Robust Statistics: The Approach Based on
   Influence Functions.* Wiley.
6. Koch, K.-R. (1999). *Parameter Estimation and Hypothesis Testing in
   Linear Models.* 2nd ed., Springer.
7. Strang, G. & Borre, K. (1997). *Linear Algebra, Geodesy, and GPS.*
   Wellesley-Cambridge Press.
