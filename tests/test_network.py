"""Design matrix structure, datum handling, connectivity."""

import numpy as np
import pytest

from gnss_adjust import network as net
from gnss_adjust.types import BaselineObservation, Station


def test_design_matrix_shape_and_signs(square_noiseless):
    stations, obs = square_noiseless
    system = net.build_design_system(stations, obs)
    m, u = len(obs), 3  # 3 unknown stations
    assert system.A.shape == (3 * m, 3 * u)
    assert system.unknown_ids == ["S01", "S02", "S03"]

    # baseline S01 -> S02 must have -I3 at S01 columns, +I3 at S02 columns
    k = next(i for i, o in enumerate(obs)
             if o.from_station == "S01" and o.to_station == "S02")
    rows = slice(3 * k, 3 * k + 3)
    col = {sid: 3 * j for j, sid in enumerate(system.unknown_ids)}
    np.testing.assert_array_equal(system.A[rows, col["S01"]:col["S01"] + 3], -np.eye(3))
    np.testing.assert_array_equal(system.A[rows, col["S02"]:col["S02"] + 3], np.eye(3))
    # all other columns zero
    np.testing.assert_array_equal(system.A[rows, col["S03"]:col["S03"] + 3], np.zeros((3, 3)))


def test_fixed_station_moved_to_observations(square_noiseless):
    stations, obs = square_noiseless
    system = net.build_design_system(stations, obs)
    # baseline S00(fixed) -> S01: l_reduced = l + X_fixed = X_S01 exactly (noiseless)
    k = next(i for i, o in enumerate(obs)
             if o.from_station == "S00" and o.to_station == "S01")
    np.testing.assert_allclose(system.l_reduced[3 * k:3 * k + 3],
                               np.array([1000.0, 0.0, 10.0]), atol=1e-12)


def test_reversed_baseline_consistency(square_noiseless):
    """Reversing a baseline flips its observed vector and its A rows."""
    stations, obs = square_noiseless
    k = next(i for i, o in enumerate(obs)
             if o.from_station == "S01" and o.to_station == "S02")
    reversed_obs = list(obs)
    reversed_obs[k] = obs[k].reversed()
    s1 = net.build_design_system(stations, obs)
    s2 = net.build_design_system(stations, reversed_obs)
    rows = slice(3 * k, 3 * k + 3)
    np.testing.assert_allclose(s2.A[rows], -s1.A[rows])
    np.testing.assert_allclose(s2.l_reduced[rows], -s1.l_reduced[rows])
    # unchanged rows unaffected
    other = slice(0, 3 * k)
    np.testing.assert_allclose(s2.A[other], s1.A[other])


def test_no_fixed_station_raises():
    stations = [Station("A", np.zeros(3)), Station("B", np.ones(3))]
    obs = [BaselineObservation(0, "A", "B", np.ones(3), np.eye(3) * 1e-4)]
    with pytest.raises(net.NetworkError, match="datum"):
        net.build_design_system(stations, obs)


def test_unknown_endpoint_raises(square_stations):
    obs = [BaselineObservation(0, "S00", "SXX", np.ones(3), np.eye(3) * 1e-4)]
    with pytest.raises(net.NetworkError, match="unknown station"):
        net.validate_observations(square_stations, obs)


def test_duplicate_baseline_id_raises(square_stations):
    q = np.eye(3) * 1e-4
    obs = [BaselineObservation(0, "S00", "S01", np.ones(3), q),
           BaselineObservation(0, "S01", "S02", np.ones(3), q)]
    with pytest.raises(net.NetworkError, match="duplicate baseline"):
        net.validate_observations(square_stations, obs)


def test_connectivity_detection(square_stations):
    q = np.eye(3) * 1e-4
    connected = [BaselineObservation(0, "S00", "S01", np.ones(3), q),
                 BaselineObservation(1, "S01", "S02", np.ones(3), q),
                 BaselineObservation(2, "S02", "S03", np.ones(3), q)]
    assert net.is_connected(square_stations, connected)
    assert not net.is_connected(square_stations, connected[:2])
    assert not net.is_connected_without(square_stations, connected, removed_baseline=1)
