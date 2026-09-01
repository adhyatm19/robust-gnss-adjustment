"""Paper-faithful specific-direction calculations for one GNSS baseline.

Nie, Yang, and Shen (2019), Equations (8), (14), (16), and (17):

``P_bar_ii d_hat_i = g_i``, ``u_i* = d_hat_i / ||d_hat_i||``, and
``|w_i*| = sqrt(d_hat_i.T P_bar_ii d_hat_i) / sigma0``.

No explicit matrix inverse is formed.  The full-rank 3x3 reliability block is
solved with a symmetric positive-definite solver.  The implementation uses
the project's residual sign convention as documented in ``reliability.py``.
"""

from __future__ import annotations

import warnings

import numpy as np
import scipy.linalg
from scipy import stats

from ..types import AdjustmentResult
from .reliability import baseline_reliability
from .results import SpecificDirectionResult


def angular_error_degrees(
    estimated_direction: np.ndarray,
    reference_direction: np.ndarray,
    *,
    sign_invariant: bool = True,
) -> float:
    """Return the angle between two directions in degrees.

    The paper's maximizing direction is defined up to sign, so the default
    uses ``abs(dot)`` and reports an angle in the range 0 to 90 degrees.
    """
    estimated = np.asarray(estimated_direction, dtype=float).reshape(3)
    reference = np.asarray(reference_direction, dtype=float).reshape(3)
    estimated_norm = float(np.linalg.norm(estimated))
    reference_norm = float(np.linalg.norm(reference))
    if (
        estimated_norm == 0.0
        or reference_norm == 0.0
        or not np.all(np.isfinite(estimated))
        or not np.all(np.isfinite(reference))
    ):
        return float("nan")
    cosine = float(estimated @ reference / (estimated_norm * reference_norm))
    if sign_invariant:
        cosine = abs(cosine)
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def _critical_values(alpha: float) -> tuple[float, float]:
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must lie strictly between 0 and 1")
    chi2_critical = float(stats.chi2.ppf(1.0 - alpha, df=3))
    return float(np.sqrt(chi2_critical)), chi2_critical / 3.0


def test_specific_direction(
    result: AdjustmentResult,
    baseline_index: int,
    *,
    alpha: float = 0.001,
    rank_tolerance: float = 1e-10,
    zero_tolerance: float = 1e-12,
    warn_equivalence_tolerance: float = 1e-10,
) -> SpecificDirectionResult:
    """Evaluate the SD and equivalent 3D statistics for one baseline.

    The known a-priori variance-factor form from the paper is used, matching
    the rest of this project's diagnostics.  Blocks with numerical rank below
    three are not covered by the paper's positive-definite derivation and are
    returned as untestable.
    """
    p_ii, g_i = baseline_reliability(result, baseline_index)
    p_ii = 0.5 * (p_ii + p_ii.T)
    eigenvalues = np.linalg.eigvalsh(p_ii)
    scale = max(float(eigenvalues[-1]), 0.0)
    tolerance = max(rank_tolerance * scale, np.finfo(float).eps)
    rank = int(np.sum(eigenvalues > tolerance))
    positive = eigenvalues[eigenvalues > tolerance]
    condition = (
        float(positive.max() / positive.min()) if len(positive) else float("inf")
    )
    critical_sd, critical_3d = _critical_values(alpha)

    if rank < 3:
        nan3 = np.full(3, np.nan)
        return SpecificDirectionResult(
            baseline_id=result.baseline_ids[baseline_index],
            baseline_index=baseline_index,
            reliability_block=p_ii,
            g=g_i,
            outlier_vector=nan3.copy(),
            outlier_magnitude=float("nan"),
            direction=nan3.copy(),
            sd_statistic=float("nan"),
            statistic_3d=float("nan"),
            critical_sd=critical_sd,
            critical_3d=critical_3d,
            significant=False,
            testable=False,
            reliability_rank=rank,
            reliability_condition=condition,
            equivalence_error=float("nan"),
        )

    d_hat = scipy.linalg.solve(p_ii, g_i, assume_a="pos", check_finite=True)
    magnitude = float(np.linalg.norm(d_hat))
    # With zero evidence every spatial direction gives the same statistic;
    # return undefined rather than fabricating a preferred unit vector.
    direction = (
        np.full(3, np.nan) if magnitude <= zero_tolerance else d_hat / magnitude
    )
    quadratic = float(d_hat @ p_ii @ d_hat)
    # Round-off may produce a tiny negative value for a theoretically PSD form.
    quadratic = max(quadratic, 0.0)
    sigma0_sq = float(result.sigma0_prior_sq)
    if not np.isfinite(sigma0_sq) or sigma0_sq <= 0.0:
        raise ValueError("sigma0_prior_sq must be finite and positive")
    sd_statistic = float(np.sqrt(quadratic / sigma0_sq))
    statistic_3d = quadratic / (3.0 * sigma0_sq)
    equivalence_error = abs(sd_statistic ** 2 - 3.0 * statistic_3d)
    equivalence_scale = max(sd_statistic ** 2, 1.0)
    if equivalence_error > warn_equivalence_tolerance * equivalence_scale:
        warnings.warn(
            "specific-direction identity |w*|^2 = 3*T failed beyond tolerance",
            RuntimeWarning,
            stacklevel=2,
        )

    return SpecificDirectionResult(
        baseline_id=result.baseline_ids[baseline_index],
        baseline_index=baseline_index,
        reliability_block=p_ii,
        g=g_i,
        outlier_vector=d_hat,
        outlier_magnitude=magnitude,
        direction=direction,
        sd_statistic=sd_statistic,
        statistic_3d=statistic_3d,
        critical_sd=critical_sd,
        critical_3d=critical_3d,
        significant=bool(sd_statistic > critical_sd),
        testable=True,
        reliability_rank=rank,
        reliability_condition=condition,
        equivalence_error=equivalence_error,
    )
