"""Huber / Hampel weight functions and IRLS behaviour."""

import numpy as np
import pytest

from gnss_adjust.robust import RobustConfig, hampel_weight, huber_weight, robust_adjust
from gnss_adjust.wls import solve_wls


def test_huber_weight_function():
    assert huber_weight(0.0) == 1.0
    assert huber_weight(1.5, c=1.5) == 1.0
    assert huber_weight(3.0, c=1.5) == pytest.approx(0.5)
    assert huber_weight(-3.0, c=1.5) == pytest.approx(0.5)
    with pytest.raises(ValueError):
        huber_weight(1.0, c=0.0)


def test_hampel_weight_function():
    a, b, c = 1.5, 3.5, 8.0
    assert hampel_weight(1.0, a, b, c) == 1.0
    assert hampel_weight(2.0, a, b, c) == pytest.approx(0.75)
    # descending segment continuous at b and reaching 0 at c
    assert hampel_weight(b, a, b, c) == pytest.approx(a / b)
    assert hampel_weight(c, a, b, c) == pytest.approx(0.0)
    assert hampel_weight(100.0, a, b, c) == 0.0
    with pytest.raises(ValueError):
        hampel_weight(1.0, a=3.5, b=1.5, c=8.0)


def test_huber_converges_and_downweights(simulated_outlier):
    stations, obs, truth = simulated_outlier
    result = robust_adjust(stations, obs, RobustConfig(method="huber"))
    assert result.converged
    assert result.weights[3] < 0.2               # 0.5 m outlier -> strong down-weight
    assert all(w > 0 for w in result.weights.values())  # Huber never hits zero
    wls = solve_wls(stations, obs)
    def rmse(coords):
        d = [coords[s] - truth.true_coordinates[s] for s in coords if s != "S00"]
        return float(np.sqrt(np.mean(np.square(d))))
    assert rmse(result.final.coordinates) < rmse(wls.coordinates)


def test_hampel_suppresses_extreme_outlier(simulated_outlier):
    stations, obs, truth = simulated_outlier
    result = robust_adjust(stations, obs, RobustConfig(method="hampel"))
    assert result.converged
    assert result.weights[3] == 0.0              # redescending: full rejection
    assert 3 in result.downweighted
    # non-outliers keep essentially full weight
    others = [w for b, w in result.weights.items() if b != 3]
    assert np.median(others) > 0.8


def test_robust_close_to_wls_on_clean_data(simulated_clean):
    stations, obs, truth = simulated_clean
    wls = solve_wls(stations, obs)
    for method in ("huber", "hampel"):
        rob = robust_adjust(stations, obs, RobustConfig(method=method))
        assert rob.converged
        for sid in wls.coordinates:
            # clean-data efficiency: within 2 mm of WLS
            assert np.linalg.norm(rob.final.coordinates[sid] - wls.coordinates[sid]) < 2e-3


def test_weights_not_compounded(simulated_outlier):
    """Multipliers are recomputed from scratch: final weight equals w(t_final)."""
    stations, obs, _ = simulated_outlier
    cfg = RobustConfig(method="huber")
    result = robust_adjust(stations, obs, cfg)
    from gnss_adjust.robust import huber_weight as hw
    for bid, t in result.t_values.items():
        assert result.weights[bid] == pytest.approx(hw(t, cfg.huber_c), abs=1e-9)


def _chain_stations():
    from gnss_adjust.types import Station
    return [
        Station("A", np.array([0.0, 0.0, 0.0]), is_fixed=True),
        Station("B", np.array([1000.0, 0.0, 0.0])),
        Station("C", np.array([2000.0, 0.0, 0.0])),
    ]


def test_hampel_rejects_bad_bridge_with_majority():
    """With 2 clean + 1 gross bridge observation, only the bad one is zeroed."""
    from gnss_adjust.types import BaselineObservation
    q = np.eye(3) * 1e-4
    obs = [
        BaselineObservation(0, "A", "B", np.array([1000.0, 0.0, 0.0]), q),
        BaselineObservation(1, "A", "B", np.array([1000.005, 0.0, 0.0]), q),
        BaselineObservation(2, "B", "C", np.array([1000.0, 0.0, 0.0]), q),
        BaselineObservation(3, "B", "C", np.array([1000.003, 0.0, 0.0]), q),
        BaselineObservation(4, "B", "C", np.array([1000.0, 0.0, 0.2]), q),  # gross
    ]
    result = robust_adjust(_chain_stations(), obs, RobustConfig(method="hampel"))
    assert result.weights[4] == 0.0
    assert all(result.weights[b] > 0.5 for b in (0, 1, 2, 3))
    assert np.linalg.norm(result.final.coordinates["C"] - np.array([2000.0, 0.0, 0.0])) < 0.02


def test_hampel_rank_guard_survives_ambiguous_bridges():
    """Two contradictory bridge observations both get rejected; the rank guard
    must keep the solver alive (no crash, finite coordinates) even though the
    conflict is statistically undecidable (no third observation to arbitrate)."""
    from gnss_adjust.types import BaselineObservation
    q = np.eye(3) * 1e-4
    obs = [
        BaselineObservation(0, "A", "B", np.array([1000.0, 0.0, 0.0]), q),
        BaselineObservation(1, "A", "B", np.array([1000.005, 0.0, 0.0]), q),
        BaselineObservation(2, "B", "C", np.array([1000.0, 0.0, 0.0]), q),
        BaselineObservation(3, "B", "C", np.array([1000.0, 0.0, 3.0]), q),
    ]
    result = robust_adjust(_chain_stations(), obs, RobustConfig(method="hampel"))
    assert np.all(np.isfinite(result.final.coordinates["C"]))
    # both conflicting bridges flagged as fully down-weighted -> visible to user
    assert result.weights[2] == 0.0 and result.weights[3] == 0.0
    assert {2, 3}.issubset(set(result.downweighted))


def test_convergence_history_recorded(simulated_outlier):
    stations, obs, _ = simulated_outlier
    result = robust_adjust(stations, obs, RobustConfig(method="huber"))
    assert result.iterations == len(result.objective_history)
    assert len(result.max_weight_change_history) == result.iterations
    assert result.max_weight_change_history[-1] <= 1e-4
    assert len(result.scale_history) == result.iterations
