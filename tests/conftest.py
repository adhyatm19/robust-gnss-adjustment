"""Shared fixtures: small deterministic networks."""

from __future__ import annotations

import numpy as np
import pytest

from gnss_adjust.simulation import NetworkConfig, NoiseModel, OutlierSpec, simulate_network
from gnss_adjust.types import BaselineObservation, Station


@pytest.fixture
def square_stations() -> list[Station]:
    """Four stations on a 1 km square; S00 fixed."""
    return [
        Station("S00", np.array([0.0, 0.0, 0.0]), is_fixed=True),
        Station("S01", np.array([1000.0, 0.0, 10.0])),
        Station("S02", np.array([1000.0, 1000.0, -5.0])),
        Station("S03", np.array([0.0, 1000.0, 20.0])),
    ]


def make_noiseless_observations(stations: list[Station]) -> list[BaselineObservation]:
    """All 6 pairwise baselines, exact vectors, simple diagonal covariance."""
    coords = {s.station_id: s.coords for s in stations}
    Q = np.diag([1e-4, 1e-4, 4e-4])  # (10 mm, 10 mm, 20 mm) sigmas
    obs = []
    bid = 0
    ids = [s.station_id for s in stations]
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            obs.append(BaselineObservation(
                baseline_id=bid, from_station=ids[i], to_station=ids[j],
                vector=coords[ids[j]] - coords[ids[i]], covariance=Q.copy()))
            bid += 1
    return obs


@pytest.fixture
def square_noiseless(square_stations):
    return square_stations, make_noiseless_observations(square_stations)


@pytest.fixture
def simulated_clean():
    cfg = NetworkConfig(n_stations=8, n_baselines=16, seed=7)
    return simulate_network(cfg)


@pytest.fixture
def simulated_outlier():
    """One deterministic huge outlier (0.5 m in Z) on baseline 3."""
    cfg = NetworkConfig(
        n_stations=8, n_baselines=16, seed=7,
        outliers=[OutlierSpec(baseline_id=3, vector=(0.0, 0.0, 0.5))])
    return simulate_network(cfg)
