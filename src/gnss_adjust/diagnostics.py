"""Statistical diagnostics: global model test, component w-tests, 3D baseline tests.

All tests use the *a-priori* variance factor ``sigma0_prior_sq`` (default 1),
never the a-posteriori one — mixing them is a classic DIA implementation bug,
because an outlier inflates the a-posteriori factor and masks itself.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
from scipy import stats

from .types import AdjustmentResult, BaselineTest

logger = logging.getLogger(__name__)


@dataclass
class GlobalTest:
    statistic: float           # v'Pv / sigma0_prior^2
    critical: float            # chi2.ppf(1 - alpha, r)
    dof: int
    alpha: float
    passed: bool


def global_model_test(result: AdjustmentResult, alpha: float = 0.05) -> GlobalTest:
    """Chi-square overall model test on v'Pv."""
    r = result.degrees_of_freedom
    if r <= 0:
        return GlobalTest(statistic=float("nan"), critical=float("nan"),
                          dof=r, alpha=alpha, passed=True)
    T = result.vTPv / result.sigma0_prior_sq
    crit = float(stats.chi2.ppf(1.0 - alpha, r))
    return GlobalTest(statistic=float(T), critical=crit, dof=r, alpha=alpha,
                      passed=bool(T <= crit))


def component_critical(alpha_local: float, n_tests: int = 1, bonferroni: bool = False) -> float:
    """Two-sided normal critical value for the w-test, optionally Bonferroni-adjusted."""
    a = alpha_local / n_tests if bonferroni else alpha_local
    return float(stats.norm.ppf(1.0 - a / 2.0))


def baseline_tests(
    result: AdjustmentResult,
    alpha_local: float = 0.001,
    alpha_baseline: float = 0.001,
    bonferroni: bool = False,
    true_outliers: set[int] | None = None,
) -> list[BaselineTest]:
    """Per-baseline identification statistics.

    Component statistic:  w_i = v_i / (sigma0_prior * sqrt(Qv_ii))
    Group statistic:      T_b = v_b' pinv(Qv_bb) v_b / sigma0_prior^2,
    compared against chi2(df_eff) where df_eff is the numerical rank of Qv_bb
    (normally 3, less if the residual block is constrained).
    """
    m = len(result.baseline_ids)
    n_comp_tests = 3 * m
    w_crit = component_critical(alpha_local, n_comp_tests, bonferroni)
    sigma0_sq = result.sigma0_prior_sq

    out: list[BaselineTest] = []
    for k, bid in enumerate(result.baseline_ids):
        v_b = result.residual_block(k)
        Qv_bb = 0.5 * (result.residual_cov_block(k) + result.residual_cov_block(k).T)
        w_b = result.standardized_residuals[3 * k: 3 * k + 3]
        # numerical rank of the residual block
        eigvals = np.linalg.eigvalsh(Qv_bb)
        tol = max(eigvals.max(), 0.0) * 1e-10 + 1e-30
        df_eff = int(np.sum(eigvals > tol))
        if df_eff > 0:
            T_b = float(v_b @ np.linalg.pinv(Qv_bb, rcond=1e-10) @ v_b) / sigma0_sq
            crit_b = float(stats.chi2.ppf(1.0 - alpha_baseline, df_eff))
        else:
            T_b, crit_b = float("nan"), float("nan")
        max_w = float(np.nanmax(np.abs(w_b))) if np.any(np.isfinite(w_b)) else float("nan")
        out.append(BaselineTest(
            baseline_id=bid,
            residual=v_b.copy(),
            residual_norm=float(np.linalg.norm(v_b)),
            standardized_components=w_b.copy(),
            max_abs_component_stat=max_w,
            component_critical=w_crit,
            group_statistic=T_b,
            group_critical=crit_b,
            group_df=df_eff,
            flagged_component=bool(np.isfinite(max_w) and max_w > w_crit),
            flagged_group=bool(np.isfinite(T_b) and T_b > crit_b),
            is_true_outlier=(bid in true_outliers) if true_outliers is not None else None,
        ))
    return out
