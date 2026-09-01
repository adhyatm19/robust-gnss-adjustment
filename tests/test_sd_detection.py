"""End-to-end SD detection and sequential snooping."""

import numpy as np

from gnss_adjust.metrics import coordinate_metrics
from gnss_adjust.sd_outlier import SDConfig, detect_specific_direction, run_sd_snooping
from gnss_adjust.wls import solve_wls


def test_sd_localizes_large_planted_outlier(simulated_outlier):
    stations, observations, _ = simulated_outlier
    detection = detect_specific_direction(solve_wls(stations, observations))

    assert detection.identified_baseline == 3
    assert detection.identified is not None
    assert detection.identified.outlier_vector.shape == (3,)


def test_sd_snooping_rejects_outlier_and_improves_coordinates(simulated_outlier):
    stations, observations, truth = simulated_outlier
    fixed = {station.station_id for station in stations if station.is_fixed}
    wls = solve_wls(stations, observations)
    sd = run_sd_snooping(stations, observations, SDConfig(adaptation="reject"))
    wls_rmse = coordinate_metrics(
        truth.true_coordinates, wls.coordinates, exclude=fixed).rmse_m
    sd_rmse = coordinate_metrics(
        truth.true_coordinates, sd.final.coordinates, exclude=fixed).rmse_m

    assert sd.flagged_baselines == [3]
    assert sd.status == "adapted"
    assert sd_rmse < wls_rmse
    assert sd.iterations[0].action == "reject:baseline_3"
    assert sd.iterations[-1].global_passed


def test_sd_correction_records_estimated_vector(simulated_outlier):
    stations, observations, _ = simulated_outlier
    sd = run_sd_snooping(
        stations,
        observations,
        SDConfig(adaptation="correct", max_adaptations=1),
    )

    assert sd.flagged_baselines == [3]
    assert sd.iterations[0].estimated_outlier_vector is not None
    assert np.isfinite(sd.iterations[0].estimated_outlier_vector).all()
