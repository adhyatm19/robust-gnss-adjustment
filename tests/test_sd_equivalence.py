"""Numerical checks of the paper's SD-to-3D equivalence."""

import numpy as np

from gnss_adjust.sd_outlier import (
    detect_specific_direction,
    eliminate_along_specific_direction,
    eliminate_full_3d,
)
from gnss_adjust.wls import solve_wls


def test_sd_statistic_equals_three_dimensional_statistic(simulated_clean):
    stations, observations, _ = simulated_clean
    detection = detect_specific_direction(solve_wls(stations, observations))

    for test in detection.tests:
        if test.testable:
            np.testing.assert_allclose(
                test.sd_statistic ** 2,
                3.0 * test.statistic_3d,
                rtol=1e-12,
                atol=1e-12,
            )


def test_sd_elimination_matches_full_3d_elimination(simulated_outlier):
    stations, observations, _ = simulated_outlier
    detection = detect_specific_direction(solve_wls(stations, observations))
    assert detection.identified is not None

    sd_adjusted = eliminate_along_specific_direction(
        stations, observations, detection.identified)
    full_3d_adjusted = eliminate_full_3d(
        stations, observations, detection.identified.baseline_id)

    for station_id in sd_adjusted.coordinates:
        np.testing.assert_allclose(
            sd_adjusted.coordinates[station_id],
            full_3d_adjusted.coordinates[station_id],
            atol=2e-9,
        )

