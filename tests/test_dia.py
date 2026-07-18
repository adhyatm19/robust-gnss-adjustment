"""DIA detection, identification, adaptation, audit trail, rank guard."""

import numpy as np

from gnss_adjust.dia import DIAConfig, run_dia
from gnss_adjust.diagnostics import global_model_test
from gnss_adjust.types import BaselineObservation, Station
from gnss_adjust.wls import solve_wls


def test_clean_data_passes_global_test(simulated_clean):
    stations, obs, _ = simulated_clean
    result = run_dia(stations, obs, DIAConfig(alpha_global=0.001))
    assert result.status == "clean"
    assert result.flagged_baselines == []
    assert len(result.iterations) == 1
    assert result.iterations[0].global_passed


def test_dia_identifies_planted_outlier(simulated_outlier):
    stations, obs, truth = simulated_outlier
    result = run_dia(stations, obs, DIAConfig(adaptation="reject"))
    assert result.flagged_baselines == [3]
    assert result.status == "adapted"
    # final adjustment must pass the global test
    assert result.iterations[-1].global_passed


def test_dia_inflate_identifies_planted_outlier(simulated_outlier):
    stations, obs, truth = simulated_outlier
    result = run_dia(stations, obs, DIAConfig(adaptation="inflate", inflation_factor=1000))
    assert 3 in result.flagged_baselines
    assert result.status == "adapted"


def test_dia_improves_coordinates(simulated_outlier):
    stations, obs, truth = simulated_outlier
    wls = solve_wls(stations, obs)
    dia = run_dia(stations, obs, DIAConfig(adaptation="reject"))
    def rmse(coords):
        d = [coords[s] - truth.true_coordinates[s]
             for s in coords if not any(st.station_id == s and st.is_fixed for st in stations)]
        return float(np.sqrt(np.mean(np.square(d))))
    assert rmse(dia.final.coordinates) < rmse(wls.coordinates)


def test_audit_trail_complete(simulated_outlier):
    stations, obs, _ = simulated_outlier
    result = run_dia(stations, obs, DIAConfig(adaptation="reject"),
                     true_outliers={3})
    assert len(result.iterations) >= 2
    first = result.iterations[0]
    assert not first.global_passed
    assert first.identified_baseline == 3
    assert first.action == "reject:baseline_3"
    assert 3 in first.active_baselines
    # true-outlier bookkeeping recorded on the per-baseline tests
    t3 = next(t for t in first.baseline_tests if t.baseline_id == 3)
    assert t3.is_true_outlier is True
    last = result.iterations[-1]
    assert 3 not in last.active_baselines
    assert last.global_passed
    # every iteration recorded stats and critical values
    for it in result.iterations:
        assert np.isfinite(it.global_statistic)
        assert np.isfinite(it.global_critical)


def test_rank_guard_prevents_disconnection():
    """A bridge baseline with a huge outlier must not be silently removed."""
    stations = [
        Station("A", np.array([0.0, 0.0, 0.0]), is_fixed=True),
        Station("B", np.array([1000.0, 0.0, 0.0])),
        Station("C", np.array([2000.0, 0.0, 0.0])),
    ]
    q = np.eye(3) * 1e-4
    # A-B measured twice (redundant), B-C measured once (a bridge) with a gross error
    obs = [
        BaselineObservation(0, "A", "B", np.array([1000.0, 0.0, 0.0]), q),
        BaselineObservation(1, "A", "B", np.array([1000.01, 0.0, 0.0]), q),
        BaselineObservation(2, "B", "C", np.array([1000.0, 0.0, 5.0]), q),  # 5 m error, undetectable
    ]
    result = run_dia(stations, obs, DIAConfig(adaptation="reject"))
    # baseline 2 has zero redundancy: its residual is ~0, cannot be identified;
    # whatever happens, it must never be removed (network would disconnect)
    active_final = result.iterations[-1].active_baselines
    assert 2 in active_final
    assert all("reject:baseline_2" != it.action for it in result.iterations)
