"""Specific-direction detect-identify-adapt-readjust workflows."""

from __future__ import annotations

import logging
from dataclasses import dataclass, replace

import numpy as np

from .. import diagnostics as dg
from .. import network as net
from ..types import AdjustmentResult, BaselineObservation, Station
from ..wls import solve_wls
from .detection import detect_specific_direction
from .results import SDIteration, SDSnoopingResult, SpecificDirectionResult

logger = logging.getLogger(__name__)


@dataclass
class SDConfig:
    """Settings for sequential SD data snooping.

    ``reject`` follows the paper's instruction to eliminate the whole vector
    once the SD/3D-equivalent statistic is significant.  ``correct`` subtracts
    the estimated directional bias and is useful for demonstrating Equation
    (19); it should not be interpreted as independent new information.
    """

    alpha_global: float = 0.05
    alpha_local: float = 0.001
    bonferroni: bool = False
    sigma0_prior_sq: float = 1.0
    adaptation: str = "reject"  # "reject" | "correct"
    max_adaptations: int = 10
    min_redundancy: int = 1


def correct_observation(
    observations: list[BaselineObservation],
    test: SpecificDirectionResult,
) -> list[BaselineObservation]:
    """Subtract the estimated additive SD bias from one observation vector."""
    if not test.testable or not np.all(np.isfinite(test.outlier_vector)):
        raise ValueError("cannot correct an untestable specific-direction result")
    found = False
    corrected: list[BaselineObservation] = []
    for observation in observations:
        if observation.baseline_id == test.baseline_id:
            corrected.append(replace(
                observation,
                vector=observation.vector - test.outlier_vector,
            ))
            found = True
        else:
            corrected.append(observation)
    if not found:
        raise ValueError(f"baseline {test.baseline_id} is not present")
    return corrected


def eliminate_along_specific_direction(
    stations: list[Station],
    observations: list[BaselineObservation],
    test: SpecificDirectionResult,
    *,
    sigma0_prior_sq: float = 1.0,
) -> AdjustmentResult:
    """Readjust after estimating the one SD nuisance parameter.

    At the maximizing direction, ``u*d`` equals the complete 3D outlier
    estimate.  Subtracting it before readjustment is algebraically equivalent
    to augmenting the observation equation with that nuisance parameter.
    """
    return solve_wls(
        stations,
        correct_observation(observations, test),
        sigma0_prior_sq=sigma0_prior_sq,
    )


def eliminate_full_3d(
    stations: list[Station],
    observations: list[BaselineObservation],
    baseline_id: int,
    *,
    sigma0_prior_sq: float = 1.0,
) -> AdjustmentResult:
    """Readjust after eliminating all three components of one baseline."""
    remaining = [obs for obs in observations if obs.baseline_id != baseline_id]
    if len(remaining) == len(observations):
        raise ValueError(f"baseline {baseline_id} is not present")
    return solve_wls(stations, remaining, sigma0_prior_sq=sigma0_prior_sq)


def run_sd_snooping(
    stations: list[Station],
    observations: list[BaselineObservation],
    config: SDConfig | None = None,
) -> SDSnoopingResult:
    """Run global detection, SD identification, adaptation, and readjustment."""
    cfg = config or SDConfig()
    if cfg.adaptation not in ("reject", "correct"):
        raise ValueError(f"unknown adaptation strategy {cfg.adaptation!r}")

    active = list(observations)
    flagged: list[int] = []
    iterations: list[SDIteration] = []
    status = "clean"

    for iteration in range(cfg.max_adaptations + 1):
        adjustment = solve_wls(
            stations, active, sigma0_prior_sq=cfg.sigma0_prior_sq)
        global_test = dg.global_model_test(adjustment, alpha=cfg.alpha_global)
        detection = detect_specific_direction(
            adjustment, alpha=cfg.alpha_local, bonferroni=cfg.bonferroni)
        candidate = detection.identified
        record = SDIteration(
            iteration=iteration,
            global_statistic=global_test.statistic,
            global_critical=global_test.critical,
            global_passed=global_test.passed,
            identified_baseline=candidate.baseline_id if candidate else None,
            identified_statistic=candidate.sd_statistic if candidate else None,
            local_critical=candidate.critical_sd if candidate else None,
            estimated_outlier_vector=(candidate.outlier_vector.copy() if candidate else None),
            action="none",
            active_baselines=[obs.baseline_id for obs in active],
            baseline_tests=detection.tests,
        )

        if global_test.passed:
            record.action = "stopped:global_test_passed"
            iterations.append(record)
            status = "adapted" if flagged else "clean"
            break
        if candidate is None or candidate.baseline_id in flagged:
            record.action = "stopped:no_candidate_above_local_threshold"
            iterations.append(record)
            status = "no_candidate"
            break
        if iteration >= cfg.max_adaptations:
            record.action = "stopped:max_adaptations"
            iterations.append(record)
            status = "max_iterations"
            break

        if cfg.adaptation == "reject":
            remaining = [
                obs for obs in active if obs.baseline_id != candidate.baseline_id
            ]
            dof_after = 3 * len(remaining) - adjustment.diagnostics.n_parameters
            if (not net.is_connected(stations, remaining)) or dof_after < cfg.min_redundancy:
                record.action = "stopped:rank_guard"
                iterations.append(record)
                status = "rank_guard"
                break
            active = remaining
            record.action = f"reject:baseline_{candidate.baseline_id}"
        else:
            active = correct_observation(active, candidate)
            record.action = f"correct:baseline_{candidate.baseline_id}"

        flagged.append(candidate.baseline_id)
        iterations.append(record)
        logger.info(
            "SD iteration %d: %s (|w*|=%.3f > %.3f)",
            iteration,
            record.action,
            candidate.sd_statistic,
            candidate.critical_sd,
        )
    else:
        status = "max_iterations"

    final = solve_wls(stations, active, sigma0_prior_sq=cfg.sigma0_prior_sq)
    return SDSnoopingResult(
        final=final,
        iterations=iterations,
        flagged_baselines=flagged,
        status=status,
        adaptation=cfg.adaptation,
    )

