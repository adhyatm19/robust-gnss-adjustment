"""Weighted least-squares correctness."""

import numpy as np
import pytest

from gnss_adjust import network as net
from gnss_adjust.types import BaselineObservation, Station
from gnss_adjust.wls import solve_wls


def test_noiseless_network_recovers_truth(square_noiseless):
    stations, obs = square_noiseless
    result = solve_wls(stations, obs)
    for s in stations:
        np.testing.assert_allclose(result.coordinates[s.station_id], s.coords, atol=1e-9)
    assert result.vTPv == pytest.approx(0.0, abs=1e-12)
    assert result.degrees_of_freedom == 3 * 6 - 3 * 3  # 18 obs, 9 params


def test_hand_computed_example():
    """Two stations, one unknown, two observations of the same baseline.

    A fixed at origin. Observations of (B - A): (1, 0, 0) with sigma^2 = 1
    and (1.3, 0, 0) with sigma^2 = 0.5 (all components independent).
    Weighted mean for X: (1*1 + 2*1.3) / (1 + 2) = 1.2.
    """
    stations = [Station("A", np.zeros(3), is_fixed=True), Station("B", np.zeros(3))]
    obs = [
        BaselineObservation(0, "A", "B", np.array([1.0, 0.0, 0.0]), np.eye(3)),
        BaselineObservation(1, "A", "B", np.array([1.3, 0.0, 0.0]), 0.5 * np.eye(3)),
    ]
    result = solve_wls(stations, obs)
    np.testing.assert_allclose(result.coordinates["B"], [1.2, 0.0, 0.0], atol=1e-12)
    # residuals: v = A x - l -> obs0: 0.2, obs1: -0.1 (X component)
    assert result.residuals[0] == pytest.approx(0.2)
    assert result.residuals[3] == pytest.approx(-0.1)
    # v'Pv = 0.2^2 * 1 + 0.1^2 * 2 = 0.06 ; dof = 6 - 3 = 3
    assert result.vTPv == pytest.approx(0.06)
    assert result.degrees_of_freedom == 3
    # Q_x for B: (P1 + P2)^{-1} = diag(1/3) on X etc.
    np.testing.assert_allclose(np.diag(result.covariance), [1 / 3] * 3, atol=1e-12)


def test_residual_orthogonality(simulated_clean):
    stations, obs, _ = simulated_clean
    result = solve_wls(stations, obs)
    system = net.build_design_system(stations, obs)
    grad = system.A.T @ result.P @ result.residuals
    scale = max(abs(result.residuals).max(), 1.0)
    assert np.abs(grad).max() < 1e-8 * scale / 1e-3  # tight absolute check
    assert np.abs(grad).max() < 1e-6


def test_residual_covariance_symmetric_psd(simulated_clean):
    stations, obs, _ = simulated_clean
    result = solve_wls(stations, obs)
    Qv = result.residual_covariance
    np.testing.assert_allclose(Qv, Qv.T, atol=1e-12)
    eig = np.linalg.eigvalsh(0.5 * (Qv + Qv.T))
    assert eig.min() > -1e-10 * max(eig.max(), 1e-30)


def test_variance_factor_near_one_on_clean_data(simulated_clean):
    stations, obs, _ = simulated_clean
    result = solve_wls(stations, obs)
    # correctly scaled covariances -> a-posteriori factor ~ 1 (chi2 spread)
    assert 0.3 < result.variance_factor < 3.0


def test_zero_weight_excludes_observation(square_noiseless):
    stations, obs = square_noiseless
    # corrupt baseline 0 badly, then exclude it via zero weight
    bad = list(obs)
    bad[0] = BaselineObservation(0, bad[0].from_station, bad[0].to_station,
                                 bad[0].vector + np.array([9.9, 0, 0]),
                                 bad[0].covariance)
    result = solve_wls(stations, bad, weight_multipliers={0: 0.0})
    for s in stations:
        np.testing.assert_allclose(result.coordinates[s.station_id], s.coords, atol=1e-9)
