"""Core vector, direction, and statistic identities of the SD method."""

from dataclasses import replace

import numpy as np
import pytest

from gnss_adjust.sd_outlier import (
    detect_specific_direction,
    test_specific_direction as evaluate_specific_direction,
)
from gnss_adjust.types import BaselineObservation
from gnss_adjust.wls import solve_wls


def _inject(observations, baseline_id, vector):
    return [
        replace(obs, vector=obs.vector + np.asarray(vector, dtype=float))
        if obs.baseline_id == baseline_id else obs
        for obs in observations
    ]


def test_noiseless_network_has_zero_sd_evidence(square_noiseless):
    stations, observations = square_noiseless
    detection = detect_specific_direction(solve_wls(stations, observations))

    assert detection.identified_baseline is None
    for test in detection.tests:
        assert test.testable
        assert test.sd_statistic == pytest.approx(0.0, abs=1e-7)
        assert not test.significant
        assert np.isnan(test.direction).all()


def test_planted_outlier_has_unit_parallel_direction(simulated_outlier):
    stations, observations, _ = simulated_outlier
    detection = detect_specific_direction(solve_wls(stations, observations))
    test = next(item for item in detection.tests if item.baseline_id == 3)

    assert detection.identified_baseline == 3
    assert test.significant
    assert np.linalg.norm(test.direction) == pytest.approx(1.0, abs=1e-12)
    assert np.linalg.norm(np.cross(test.direction, test.outlier_vector)) < 1e-12
    assert test.sd_statistic ** 2 == pytest.approx(3.0 * test.statistic_3d)
    assert test.equivalence_error < 1e-10


def test_reversing_bias_reverses_sd_direction(square_noiseless):
    stations, observations = square_noiseless
    bias = np.array([0.25, -0.10, 0.40])
    positive = solve_wls(stations, _inject(observations, 2, bias))
    negative = solve_wls(stations, _inject(observations, 2, -bias))
    t_pos = evaluate_specific_direction(positive, positive.baseline_ids.index(2))
    t_neg = evaluate_specific_direction(negative, negative.baseline_ids.index(2))

    np.testing.assert_allclose(t_pos.direction, -t_neg.direction, atol=1e-10)
    assert t_pos.sd_statistic == pytest.approx(t_neg.sd_statistic, rel=1e-10)


def test_unredundant_baseline_is_not_fabricated():
    from gnss_adjust.types import Station

    stations = [
        Station("A", np.zeros(3), is_fixed=True),
        Station("B", np.array([1.0, 0.0, 0.0])),
    ]
    observations = [
        BaselineObservation(7, "A", "B", np.array([1.2, 0.0, 0.0]), np.eye(3))
    ]
    result = solve_wls(stations, observations)
    test = evaluate_specific_direction(result, 0)

    assert not test.testable
    assert test.reliability_rank == 0
    assert not test.significant
    assert np.isnan(test.direction).all()
