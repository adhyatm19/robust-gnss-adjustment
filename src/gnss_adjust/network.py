"""Network construction: design matrix, datum handling, connectivity checks."""

from __future__ import annotations

import logging

import networkx as nx
import numpy as np

from .types import BaselineObservation, DesignSystem, Station

logger = logging.getLogger(__name__)

I3 = np.eye(3)


class NetworkError(ValueError):
    """Raised for structural problems (unknown stations, disconnection, rank loss)."""


def station_map(stations: list[Station]) -> dict[str, Station]:
    m: dict[str, Station] = {}
    for s in stations:
        if s.station_id in m:
            raise NetworkError(f"duplicate station id {s.station_id!r}")
        m[s.station_id] = s
    return m


def validate_observations(stations: list[Station], observations: list[BaselineObservation]) -> None:
    """Check endpoint existence, duplicate baseline ids, degenerate baselines."""
    smap = station_map(stations)
    seen_ids: set[int] = set()
    for obs in observations:
        if obs.baseline_id in seen_ids:
            raise NetworkError(f"duplicate baseline id {obs.baseline_id}")
        seen_ids.add(obs.baseline_id)
        for sid in (obs.from_station, obs.to_station):
            if sid not in smap:
                raise NetworkError(f"baseline {obs.baseline_id} references unknown station {sid!r}")
        if obs.from_station == obs.to_station:
            raise NetworkError(f"baseline {obs.baseline_id} connects a station to itself")


def build_graph(stations: list[Station], observations: list[BaselineObservation]) -> nx.MultiGraph:
    g: nx.MultiGraph = nx.MultiGraph()
    g.add_nodes_from(s.station_id for s in stations)
    for obs in observations:
        g.add_edge(obs.from_station, obs.to_station, key=obs.baseline_id)
    return g


def is_connected(stations: list[Station], observations: list[BaselineObservation]) -> bool:
    g = build_graph(stations, observations)
    return g.number_of_nodes() > 0 and nx.is_connected(g)


def is_connected_without(
    stations: list[Station],
    observations: list[BaselineObservation],
    removed_baseline: int,
) -> bool:
    """Would the network stay connected if ``removed_baseline`` were dropped?"""
    remaining = [o for o in observations if o.baseline_id != removed_baseline]
    return is_connected(stations, remaining)


def build_design_system(
    stations: list[Station],
    observations: list[BaselineObservation],
) -> DesignSystem:
    """Assemble ``A`` and the reduced observation vector.

    For a baseline from ``i`` to ``j`` with expectation ``X_j - X_i``:

    * station ``i`` unknown  ->  ``-I3`` in its columns;
    * station ``j`` unknown  ->  ``+I3`` in its columns;
    * a fixed endpoint's contribution is moved into the observation vector,
      producing ``l_reduced`` with  E[l_reduced] = A x.
    """
    validate_observations(stations, observations)
    smap = station_map(stations)
    fixed = {s.station_id: s.coords.copy() for s in stations if s.is_fixed}
    if not fixed:
        raise NetworkError("no fixed station: the datum is undefined")
    unknown_ids = sorted(s.station_id for s in stations if not s.is_fixed)
    col = {sid: 3 * k for k, sid in enumerate(unknown_ids)}

    m = len(observations)
    A = np.zeros((3 * m, 3 * len(unknown_ids)))
    l_red = np.zeros(3 * m)
    baseline_ids: list[int] = []
    for k, obs in enumerate(observations):
        r = 3 * k
        l_red[r:r + 3] = obs.vector
        i, j = obs.from_station, obs.to_station
        if smap[i].is_fixed and smap[j].is_fixed:
            logger.warning("baseline %s connects two fixed stations; it constrains nothing", obs.baseline_id)
        if smap[i].is_fixed:
            l_red[r:r + 3] += smap[i].coords     # E[l] = X_j - X_i_fixed
        else:
            A[r:r + 3, col[i]:col[i] + 3] = -I3
        if smap[j].is_fixed:
            l_red[r:r + 3] -= smap[j].coords     # E[l] = X_j_fixed - X_i
        else:
            A[r:r + 3, col[j]:col[j] + 3] = I3
        baseline_ids.append(obs.baseline_id)

    return DesignSystem(A=A, l_reduced=l_red, unknown_ids=unknown_ids,
                        baseline_ids=baseline_ids, fixed=fixed)


def datum_station(stations: list[Station]) -> str | None:
    fixed = [s.station_id for s in stations if s.is_fixed]
    return fixed[0] if fixed else None
