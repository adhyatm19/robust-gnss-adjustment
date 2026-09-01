"""Test every active baseline and identify the strongest SD candidate."""

from __future__ import annotations

from ..types import AdjustmentResult
from .results import SpecificDirectionDetection
from .specific_direction import test_specific_direction


def detect_specific_direction(
    result: AdjustmentResult,
    *,
    alpha: float = 0.001,
    bonferroni: bool = False,
) -> SpecificDirectionDetection:
    """Return all SD tests and the largest significant local statistic.

    When ``bonferroni=True``, ``alpha`` is divided by the number of baseline
    vectors tested.  Untestable zero-redundancy blocks remain in the audit
    output but cannot become candidates.
    """
    n_tests = len(result.baseline_ids)
    local_alpha = alpha / n_tests if bonferroni and n_tests else alpha
    tests = [
        test_specific_direction(result, k, alpha=local_alpha)
        for k in range(n_tests)
    ]
    candidates = [test for test in tests if test.testable and test.significant]
    identified = max(candidates, key=lambda test: test.sd_statistic, default=None)
    return SpecificDirectionDetection(
        adjustment=result,
        tests=tests,
        identified_baseline=identified.baseline_id if identified else None,
        identified=identified,
        alpha=alpha,
        bonferroni=bonferroni,
    )

