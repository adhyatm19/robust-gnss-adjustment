"""Detection - Identification - Adaptation (DIA) outlier handling loop.

The loop:

1. run WLS;
2. **Detection** — global chi-square test on v'Pv / sigma0_prior^2;
3. if the global test passes: stop ("clean" or "adapted");
4. **Identification** — the baseline with the largest ratio of its 3D group
   statistic to its critical value, among baselines exceeding the local
   threshold and not adapted before;
5. **Adaptation** — either remove the baseline or inflate its covariance
   block, then re-adjust;
6. stop on: global test passing, no candidate above the local threshold,
   maximum adaptations reached, or a removal that would disconnect the
   network / destroy redundancy (rank guard).

Nothing is deleted silently: every iteration is recorded in the audit trail.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, replace

import numpy as np

from . import diagnostics as dg
from . import network as net
from .types import BaselineObservation, DIAIteration, DIAResult, Station
from .wls import solve_wls

logger = logging.getLogger(__name__)


@dataclass
class DIAConfig:
    alpha_global: float = 0.05
    alpha_local: float = 0.001          # per component w-test
    alpha_baseline: float = 0.001       # per 3D baseline group test
    bonferroni: bool = False
    sigma0_prior_sq: float = 1.0
    adaptation: str = "reject"          # "reject" | "inflate"
    inflation_factor: float = 100.0     # covariance multiplier for "inflate"
    max_adaptations: int = 10
    min_redundancy: int = 1             # never adapt below this many dof


def run_dia(
    stations: list[Station],
    observations: list[BaselineObservation],
    config: DIAConfig | None = None,
    true_outliers: set[int] | None = None,
) -> DIAResult:
    """Run the full DIA loop and return final adjustment plus audit trail."""
    cfg = config or DIAConfig()
    if cfg.adaptation not in ("reject", "inflate"):
        raise ValueError(f"unknown adaptation strategy {cfg.adaptation!r}")

    active = list(observations)
    flagged: list[int] = []
    iterations: list[DIAIteration] = []
    prev_coords: dict[str, np.ndarray] | None = None
    status = "clean"

    for it in range(cfg.max_adaptations + 1):
        result = solve_wls(stations, active, sigma0_prior_sq=cfg.sigma0_prior_sq)
        gtest = dg.global_model_test(result, alpha=cfg.alpha_global)
        tests = dg.baseline_tests(
            result,
            alpha_local=cfg.alpha_local,
            alpha_baseline=cfg.alpha_baseline,
            bonferroni=cfg.bonferroni,
            true_outliers=true_outliers,
        )
        max_dx = None
        if prev_coords is not None:
            max_dx = max(
                float(np.linalg.norm(result.coordinates[sid] - prev_coords[sid]))
                for sid in result.coordinates
            )
        prev_coords = {k: v.copy() for k, v in result.coordinates.items()}

        record = DIAIteration(
            iteration=it,
            global_statistic=gtest.statistic,
            global_critical=gtest.critical,
            global_passed=gtest.passed,
            identified_baseline=None,
            identified_statistic=None,
            local_critical=None,
            action="none",
            active_baselines=[o.baseline_id for o in active],
            max_coordinate_change=max_dx,
            baseline_tests=tests,
        )

        if gtest.passed:
            record.action = "stopped:global_test_passed"
            iterations.append(record)
            status = "adapted" if flagged else "clean"
            break

        # ---- Identification: strongest group exceedance among un-adapted baselines
        candidates = [
            t for t in tests
            if t.flagged_group and t.baseline_id not in flagged
        ]
        if not candidates:
            record.action = "stopped:no_candidate_above_local_threshold"
            iterations.append(record)
            status = "no_candidate"
            break
        best = max(candidates, key=lambda t: t.group_statistic / t.group_critical)
        record.identified_baseline = best.baseline_id
        record.identified_statistic = best.group_statistic
        record.local_critical = best.group_critical

        if it >= cfg.max_adaptations:
            record.action = "stopped:max_adaptations"
            iterations.append(record)
            status = "max_iterations"
            break

        # ---- Adaptation
        if cfg.adaptation == "reject":
            remaining = [o for o in active if o.baseline_id != best.baseline_id]
            dof_after = 3 * len(remaining) - result.diagnostics.n_parameters
            if (not net.is_connected(stations, remaining)) or dof_after < cfg.min_redundancy:
                record.action = "stopped:rank_guard"
                iterations.append(record)
                status = "rank_guard"
                logger.warning(
                    "DIA: refusing to remove baseline %s (would disconnect the "
                    "network or exhaust redundancy)", best.baseline_id)
                break
            active = remaining
            record.action = f"reject:baseline_{best.baseline_id}"
        else:
            active = [
                replace(o, covariance=o.covariance * cfg.inflation_factor)
                if o.baseline_id == best.baseline_id else o
                for o in active
            ]
            record.action = f"inflate:baseline_{best.baseline_id}:x{cfg.inflation_factor:g}"
        flagged.append(best.baseline_id)
        iterations.append(record)
        logger.info("DIA iteration %d: %s (T=%.2f > %.2f)", it, record.action,
                    best.group_statistic, best.group_critical)
    else:  # loop exhausted without break (should not happen; safety)
        status = "max_iterations"

    final = solve_wls(stations, active, sigma0_prior_sq=cfg.sigma0_prior_sq)
    return DIAResult(final=final, iterations=iterations, flagged_baselines=flagged,
                     status=status, adaptation=cfg.adaptation)
