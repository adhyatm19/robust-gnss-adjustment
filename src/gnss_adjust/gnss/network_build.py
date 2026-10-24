"""Assemble processed baselines into the repository's network format, with design checks.

Algebraic versus stochastic independence
----------------------------------------
* **Algebraic (functional) independence**: within one observing session, a
  set of baselines is algebraically independent iff it contains no cycle
  (it is a forest). A baseline that closes a cycle inside one session is a
  *trivial* baseline: up to processing differences it is a linear
  combination of the others, computed from the same observations, so it adds
  no new information while the block-diagonal model would count it as
  independent. :func:`check_session_independence` rejects such designs
  unless explicitly allowed (e.g. to *demonstrate* the effect).
* **Stochastic independence** is a different, stronger property. Even
  algebraically independent baselines of one session share receiver data
  (e.g. a common base) and are therefore correlated; single-baseline
  processing (RTKLIB) does not provide that cross-covariance, and none is
  invented here: the observation covariance stays block-diagonal and the
  omission is reported as a limitation. Redundancy must come from genuinely
  separate sessions (different days), whose independence is itself an
  assumption that the loop-misclosure and repeatability tests can probe.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

from .. import io as gio
from ..types import BaselineObservation, Station


class NetworkDesignError(ValueError):
    """The requested baseline set violates the network design rules."""


@dataclass
class Loop:
    baselines: list[int]               # baseline ids
    signs: list[int]                   # +1 if traversed from->to
    sessions: list[str]
    trivial: bool                      # all baselines from one session
    misclosure: np.ndarray             # sum of signed vectors (m)
    covariance: np.ndarray             # block-diagonal prediction (m^2)
    chi2: float                        # w' C^-1 w, ~chi2(3) if model holds
    length_m: float


@dataclass
class NetworkDesign:
    stations: list[Station]
    observations: list[BaselineObservation]
    metadata: pd.DataFrame
    sessions: dict[str, list[int]]
    redundancy: int
    loops: list[Loop] = field(default_factory=list)


def check_session_independence(edges: list[tuple[str, str, str, int]],
                               allow_trivial: bool = False) -> dict[str, list[int]]:
    """Each session's baselines must form a forest (no within-session loop).

    ``edges`` are ``(session, from, to, baseline_id)``. Returns the baseline ids
    per session; raises :class:`NetworkDesignError` on a trivial baseline.
    """
    sessions: dict[str, list[int]] = {}
    graphs: dict[str, nx.Graph] = {}
    for sess, a, b, bid in edges:
        g = graphs.setdefault(sess, nx.Graph())
        if g.has_node(a) and g.has_node(b) and nx.has_path(g, a, b) and not allow_trivial:
            raise NetworkDesignError(
                f"baseline {bid} ({a}->{b}) closes a loop within session {sess}: it is a trivial "
                "(algebraically dependent) baseline computed from the same data")
        g.add_edge(a, b)
        sessions.setdefault(sess, []).append(bid)
    return sessions


def fundamental_loops(observations: list[BaselineObservation], session_of: dict[int, str]) -> list[Loop]:
    """A cycle basis of the baseline multigraph with misclosures (block-diagonal model).

    Repeated baselines between the same pair give 2-loops (repeatability).
    """
    g = nx.MultiGraph()
    for o in observations:
        g.add_edge(o.from_station, o.to_station, key=o.baseline_id)
    tree = nx.minimum_spanning_tree(nx.Graph(g))
    by_id = {o.baseline_id: o for o in observations}
    tree_edges = set()
    for a, b in tree.edges():
        bid = min(k for k in g[a][b])  # first baseline of the pair spans the tree
        tree_edges.add(bid)
    loops = []
    for o in observations:
        if o.baseline_id in tree_edges:
            continue
        # path in tree from o.to back to o.from closes the loop o.from -> o.to -> ... -> o.from
        path = nx.shortest_path(tree, o.to_station, o.from_station)
        ids, signs = [o.baseline_id], [1]
        for a, b in zip(path[:-1], path[1:]):
            bid = min(k for k in g[a][b] if k in tree_edges)
            e = by_id[bid]
            ids.append(bid)
            signs.append(1 if (e.from_station, e.to_station) == (a, b) else -1)
        w = sum(s * by_id[i].vector for i, s in zip(ids, signs))
        c = sum(by_id[i].covariance for i in ids)
        sess = [session_of[i] for i in ids]
        loops.append(Loop(baselines=ids, signs=signs, sessions=sess,
                          trivial=len(set(sess)) == 1, misclosure=w, covariance=c,
                          chi2=float(w @ np.linalg.solve(c, w)),
                          length_m=float(sum(np.linalg.norm(by_id[i].vector) for i in ids))))
    return loops


def build_network(
    records: list[dict],
    selection: list[tuple[str, str, str]],
    datum_station: str,
    approx_xyz: dict[str, np.ndarray],
    datum_xyz: np.ndarray,
    covariance_scale: float = 1.0,
    covariance_label: str = "formal",
    allow_trivial: bool = False,
) -> NetworkDesign:
    """Select ``(session, from, to)`` solutions and build a validated network.

    ``approx_xyz`` gives non-datum station coordinates (only used for
    reporting/rotations -- the linear adjustment does not need them);
    ``datum_xyz`` the coordinates held fixed. ``covariance_scale`` multiplies
    every formal block (e.g. an independently calibrated ``k^2``).
    """
    index = {(r["session_id"], r["from_station"], r["to_station"]): r for r in records}
    obs: list[BaselineObservation] = []
    meta = []
    edges = []
    for bid, key in enumerate(selection, start=1):
        if key not in index:
            raise NetworkDesignError(f"no processed solution for {key}")
        r = index[key]
        q = covariance_scale * np.asarray(r["covariance_ecef"], dtype=float)
        obs.append(BaselineObservation(baseline_id=bid, from_station=r["from_station"],
                                       to_station=r["to_station"],
                                       vector=np.asarray(r["vector_ecef"], dtype=float),
                                       covariance=q))
        edges.append((r["session_id"], r["from_station"], r["to_station"], bid))
        meta.append({"baseline_id": bid, "session": r["session_id"], "from": r["from_station"],
                     "to": r["to_station"], "strategy": r.get("strategy", ""),
                     "status": r["status"], "fix_fraction": r["fix_fraction"],
                     "length_m": r["length_m"], "covariance_source": r["covariance_source"],
                     "covariance_scale": covariance_scale, "covariance_label": covariance_label,
                     "rounding_sensitive": r["covariance_check"]["rounding_sensitive"]})
    sessions = check_session_independence(edges, allow_trivial=allow_trivial)
    names = sorted({o.from_station for o in obs} | {o.to_station for o in obs})
    if datum_station not in names:
        raise NetworkDesignError(f"datum station {datum_station} not in the network")
    stations = [Station(n, datum_xyz if n == datum_station else approx_xyz[n],
                        is_fixed=(n == datum_station)) for n in names]
    g = nx.MultiGraph()
    g.add_nodes_from(names)
    g.add_edges_from((o.from_station, o.to_station) for o in obs)
    if not nx.is_connected(g):
        raise NetworkDesignError("network graph is not connected")
    redundancy = 3 * len(obs) - 3 * (len(names) - 1)
    session_of = {bid: s for s, ids in sessions.items() for bid in ids}
    return NetworkDesign(stations=stations, observations=obs, metadata=pd.DataFrame(meta),
                         sessions=sessions, redundancy=redundancy,
                         loops=fundamental_loops(obs, session_of))


def write_network(design: NetworkDesign, out_dir: Path, reference: dict[str, np.ndarray],
                  reference_info: dict) -> None:
    """Repository CSVs + sidecars: baseline_metadata.csv, reference_coordinates.csv,
    network_design.json (sessions, redundancy, loops)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    gio.write_network_csv(design.stations, design.observations, out_dir)
    design.metadata.to_csv(out_dir / "baseline_metadata.csv", index=False)
    pd.DataFrame([{"station_id": k, "x": v[0], "y": v[1], "z": v[2], **reference_info}
                  for k, v in sorted(reference.items())]).to_csv(
        out_dir / "reference_coordinates.csv", index=False)
    (out_dir / "network_design.json").write_text(json.dumps({
        "sessions": design.sessions,
        "redundancy": design.redundancy,
        "n_stations": len(design.stations),
        "n_baselines": len(design.observations),
        "datum_station": next(s.station_id for s in design.stations if s.is_fixed),
        "cross_baseline_covariance": "not available from single-baseline processing; "
                                     "block-diagonal model used (limitation)",
        "loops": [{"baselines": lp.baselines, "signs": lp.signs, "sessions": lp.sessions,
                   "trivial": lp.trivial, "misclosure_mm": (1e3 * lp.misclosure).tolist(),
                   "chi2_3": lp.chi2, "length_km": lp.length_m / 1e3} for lp in design.loops],
    }, indent=1))
