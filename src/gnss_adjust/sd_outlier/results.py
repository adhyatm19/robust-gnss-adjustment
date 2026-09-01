"""Structured outputs for specific-direction outlier detection."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..types import AdjustmentResult


@dataclass(frozen=True)
class SpecificDirectionResult:
    """Paper-equation results for one tested 3D baseline.

    ``outlier_vector`` estimates the additive bias in the observation, so it
    has the same sign as a bias planted by :class:`~gnss_adjust.OutlierSpec`.
    A rank-deficient reliability block is reported as ``testable=False``;
    vector, direction, and statistics are then NaN rather than misleading.
    """

    baseline_id: int
    baseline_index: int
    reliability_block: np.ndarray
    g: np.ndarray
    outlier_vector: np.ndarray
    outlier_magnitude: float
    direction: np.ndarray
    sd_statistic: float
    statistic_3d: float
    critical_sd: float
    critical_3d: float
    significant: bool
    testable: bool
    reliability_rank: int
    reliability_condition: float
    equivalence_error: float


@dataclass(frozen=True)
class SpecificDirectionDetection:
    """All local SD tests from one WLS adjustment and the best candidate."""

    adjustment: AdjustmentResult
    tests: list[SpecificDirectionResult]
    identified_baseline: int | None
    identified: SpecificDirectionResult | None
    alpha: float
    bonferroni: bool


@dataclass
class SDIteration:
    """Audit record for one detect-identify-adapt iteration."""

    iteration: int
    global_statistic: float
    global_critical: float
    global_passed: bool
    identified_baseline: int | None
    identified_statistic: float | None
    local_critical: float | None
    estimated_outlier_vector: np.ndarray | None
    action: str
    active_baselines: list[int]
    baseline_tests: list[SpecificDirectionResult] = field(default_factory=list)


@dataclass
class SDSnoopingResult:
    """Outcome of sequential SD data snooping."""

    final: AdjustmentResult
    iterations: list[SDIteration]
    flagged_baselines: list[int]
    status: str
    adaptation: str
