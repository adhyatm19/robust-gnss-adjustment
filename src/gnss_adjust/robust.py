"""Robust M-estimation via iteratively reweighted least squares (IRLS).

Reweighting acts on whole 3D baseline blocks. At each iteration:

1. run WLS with effective weights ``w_b * P_b`` (``P_b`` from the *original*
   covariance block — multipliers are recomputed each iteration from scratch,
   never compounded);
2. for each baseline compute the block-standardised Mahalanobis measure

       t_b = sqrt( v_b' P0_b v_b / 3 ) / sigma_scale

   against the *fixed original* weight block ``P0_b = Q0_b^{-1}``, so that
   t_b ~ 1 for clean data and the classic tuning constants (Huber c=1.5
   etc.) apply — see :func:`_block_t_values` for why the reference must not
   move with the current robust weights;
3. map ``t_b`` through the Huber or Hampel weight function;
4. repeat until both coordinates and weights converge.

A rank guard floors hard rejections at a tiny positive weight whenever they
would disconnect the network, so the solver never fails on a singular normal
matrix; the affected baselines remain reported at weight zero.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
from scipy import stats

from . import network as net
from .covariance import weight_block
from .types import AdjustmentResult, BaselineObservation, RobustResult, Station
from .wls import AdjustmentError, solve_wls

logger = logging.getLogger(__name__)

# Median of sqrt(chi2_3 / 3): consistency constant for the MAD-style scale
# of 3-dof block-standardised residuals under Gaussian noise.
_CHI3_MEDIAN = float(np.sqrt(stats.chi2.ppf(0.5, 3) / 3.0))


def huber_weight(t: float, c: float = 1.5) -> float:
    """Huber weight function: 1 for |t| <= c, c/|t| beyond."""
    if c <= 0:
        raise ValueError("Huber tuning constant c must be positive")
    t = abs(float(t))
    return 1.0 if t <= c else c / t


def hampel_weight(t: float, a: float = 1.5, b: float = 3.5, c: float = 8.0) -> float:
    """Hampel three-part redescending weight function (0 <= a < b < c)."""
    if not (0 < a < b < c):
        raise ValueError(f"Hampel thresholds must satisfy 0 < a < b < c, got {a}, {b}, {c}")
    t = abs(float(t))
    if t <= a:
        return 1.0
    if t <= b:
        return a / t
    if t <= c:
        return (a * (c - t)) / ((c - b) * t)
    return 0.0


@dataclass
class RobustConfig:
    method: str = "huber"                    # "huber" | "hampel"
    huber_c: float = 1.5
    hampel_a: float = 1.5
    hampel_b: float = 3.5
    hampel_c: float = 8.0
    scale: str = "prior"                     # "prior" | "mad" | "aposteriori"
    sigma0_prior_sq: float = 1.0
    max_iterations: int = 50
    coord_tol_m: float = 1e-6                # convergence: max |delta x|
    weight_tol: float = 1e-4                 # convergence: max |delta w|
    min_scale: float = 1e-3                  # safeguard against collapse
    downweight_threshold: float = 0.5        # report weight below this as "strong"
    zero_weight_floor: float = 1e-6          # rank guard: floor for hard rejections
                                             # when they would disconnect the network


def _weight_fn(cfg: RobustConfig):
    if cfg.method == "huber":
        return lambda t: huber_weight(t, cfg.huber_c)
    if cfg.method == "hampel":
        return lambda t: hampel_weight(t, cfg.hampel_a, cfg.hampel_b, cfg.hampel_c)
    raise ValueError(f"unknown robust method {cfg.method!r}")


def _block_t_values(
    result: AdjustmentResult,
    observations: list[BaselineObservation],
    weight_blocks: dict[int, np.ndarray],
    sigma_scale: float,
) -> dict[int, float]:
    """Per-baseline block-standardised residual (before robust scaling):

        t_b = sqrt( v_b' P0_b v_b / 3 ) / sigma_scale

    where ``P0_b`` is the weight matrix of the *original* covariance block —
    a fixed reference that never changes with the current robust weights.

    Rationale: standardising against the residual covariance ``Qv_bb`` of the
    current weighted solve couples the statistic to the weights themselves
    (down-weighting an observation inflates its effective covariance and
    shrinks its own statistic), which produces limit cycles in the IRLS
    iteration. The other fixed candidate, the original-weight residual
    covariance ``Q0 - A Qx A'``, fails in the opposite direction: it
    subtracts the redundancy term and therefore under-states the dispersion
    of a down-weighted observation's residual, overestimating t and
    avalanching clean baselines toward zero weight. With the fixed
    observation weight ``P0`` the only feedback runs through the coordinate
    solution: down-weighting an outlier moves the solution toward the clean
    observations, which *grows* the outlier's residual and *shrinks* the
    clean residuals — a monotone separation that converges. This is the
    classic equivalent-weight practice (Danish method, IGG schemes); both
    failure modes are reproduced by ``scripts/compare_standardizations.py``.
    The price is slight conservatism for high-redundancy baselines
    (Var(v_b) <= Q0_b), which the robust scale options can absorb. The
    rigorous ``Qv``-standardised statistics remain available in the DIA
    pipeline (:mod:`gnss_adjust.diagnostics`), which tests a *fixed*
    adjustment rather than iterating on its own output.
    """
    t: dict[int, float] = {}
    for k, obs in enumerate(observations):
        v_b = result.residual_block(k)
        if not np.all(np.isfinite(v_b)):
            t[obs.baseline_id] = 0.0
            continue
        stat = float(v_b @ weight_blocks[obs.baseline_id] @ v_b)
        t[obs.baseline_id] = float(np.sqrt(max(stat, 0.0) / 3.0)) / sigma_scale
    return t


def _robust_scale(cfg: RobustConfig, result: AdjustmentResult, t_raw: dict[int, float]) -> float:
    """Scale used to standardise t_b, per the configured strategy, safeguarded."""
    if cfg.scale == "prior":
        s = float(np.sqrt(cfg.sigma0_prior_sq))
    elif cfg.scale == "mad":
        vals = np.array([x for x in t_raw.values() if np.isfinite(x)])
        s = float(np.median(vals) / _CHI3_MEDIAN) if vals.size else 1.0
    elif cfg.scale == "aposteriori":
        vf = result.variance_factor
        s = float(np.sqrt(vf)) if np.isfinite(vf) and vf > 0 else 1.0
        s = min(s, 10.0)  # an outlier-inflated factor must not mask itself
    else:
        raise ValueError(f"unknown scale option {cfg.scale!r}")
    return max(s, cfg.min_scale)


def robust_adjust(
    stations: list[Station],
    observations: list[BaselineObservation],
    config: RobustConfig | None = None,
) -> RobustResult:
    """Huber or Hampel IRLS adjustment at the 3D baseline-block level."""
    cfg = config or RobustConfig()
    wfun = _weight_fn(cfg)

    weights: dict[int, float] = {o.baseline_id: 1.0 for o in observations}
    prev_coords: dict[str, np.ndarray] | None = None
    obj_hist: list[float] = []
    dx_hist: list[float] = []
    dw_hist: list[float] = []
    scale_hist: list[float] = []
    converged = False
    result: AdjustmentResult | None = None
    t_scaled: dict[int, float] = {}

    def _guarded_solve(w: dict[int, float], validate: bool) -> AdjustmentResult:
        """Solve WLS; if hard rejections would break the network rank, floor them."""
        solve_w = w
        positive = [o for o in observations if w[o.baseline_id] > 0]
        if len(positive) < len(observations) and not net.is_connected(stations, positive):
            solve_w = {b: max(x, cfg.zero_weight_floor) for b, x in w.items()}
            logger.warning("robust rank guard: zero weights floored at %g to keep "
                           "the network connected", cfg.zero_weight_floor)
        try:
            return solve_wls(stations, observations, sigma0_prior_sq=cfg.sigma0_prior_sq,
                             weight_multipliers=solve_w, validate_covariances=validate)
        except AdjustmentError:
            if solve_w is not w:
                raise
            solve_w = {b: max(x, cfg.zero_weight_floor) for b, x in w.items()}
            logger.warning("robust rank guard: normal matrix singular, retrying "
                           "with zero weights floored at %g", cfg.zero_weight_floor)
            return solve_wls(stations, observations, sigma0_prior_sq=cfg.sigma0_prior_sq,
                             weight_multipliers=solve_w, validate_covariances=validate)

    # Fixed standardisation reference: original per-baseline weight blocks.
    P0 = {o.baseline_id: weight_block(o.covariance) for o in observations}

    for it in range(1, cfg.max_iterations + 1):
        result = _guarded_solve(weights, validate=(it == 1))
        obj_hist.append(result.vTPv)

        # raw t (scale 1), then robust scale, then scaled t
        t_raw = _block_t_values(result, observations, P0, sigma_scale=1.0)
        sigma_scale = _robust_scale(cfg, result, t_raw)
        scale_hist.append(sigma_scale)
        t_scaled = {bid: t / sigma_scale for bid, t in t_raw.items()}

        new_weights = {bid: wfun(t) for bid, t in t_scaled.items()}
        max_dw = max(abs(new_weights[b] - weights[b]) for b in weights)
        dw_hist.append(max_dw)

        max_dx = float("inf")
        if prev_coords is not None:
            max_dx = max(
                float(np.linalg.norm(result.coordinates[sid] - prev_coords[sid]))
                for sid in result.coordinates
            )
            dx_hist.append(max_dx)
        prev_coords = {k: v.copy() for k, v in result.coordinates.items()}
        weights = new_weights

        if max_dx <= cfg.coord_tol_m and max_dw <= cfg.weight_tol:
            converged = True
            break

    assert result is not None
    # Final solve with the converged weights so reported result matches them.
    final = _guarded_solve(weights, validate=False)
    downweighted = sorted(b for b, w in weights.items() if w < cfg.downweight_threshold)
    if downweighted:
        logger.info("%s IRLS down-weighted baselines: %s", cfg.method, downweighted)

    return RobustResult(
        final=final,
        method=cfg.method,
        weights=dict(weights),
        t_values=dict(t_scaled),
        iterations=len(obj_hist),
        converged=converged,
        objective_history=obj_hist,
        max_coord_change_history=dx_hist,
        max_weight_change_history=dw_hist,
        scale_history=scale_hist,
        downweighted=downweighted,
    )
