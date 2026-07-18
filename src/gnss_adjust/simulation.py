"""Reproducible synthetic GNSS baseline-network generator.

Generates true station coordinates, a connected redundant baseline graph,
realistic full 3x3 covariance blocks, Gaussian or Student-t noise, and
optional gross-error injection. Every function takes an explicit seed or
``numpy.random.Generator`` — there is no global random state.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field, asdict
from typing import Any

import numpy as np

from . import covariance as cov
from .types import BaselineObservation, SimulationTruth, Station

logger = logging.getLogger(__name__)


@dataclass
class NoiseModel:
    """GNSS baseline noise model: sigma = sqrt(a_mm^2 + (b_ppm * L)^2)."""

    a_mm: float = 5.0
    b_ppm: float = 1.0
    axis_scale: tuple[float, float, float] = (1.0, 1.0, 1.8)  # Z usually worse
    max_abs_correlation: float = 0.4
    student_t_dof: float | None = None    # None = Gaussian; else heavy-tailed
    covariance_scale_error: float = 1.0   # simulate mis-scaled reported Q (1 = honest)


@dataclass
class OutlierSpec:
    """One gross error injected into a single baseline observation.

    Exactly one of ``vector``, (``component``, magnitude), or (``direction``,
    magnitude) must be given. Magnitude may be in metres (``magnitude_m``) or
    in units of the baseline's mean component sigma (``magnitude_sigma``).
    """

    baseline_id: int
    vector: tuple[float, float, float] | None = None       # explicit bias (m)
    component: int | None = None                           # 0=X, 1=Y, 2=Z
    direction: tuple[float, float, float] | None = None    # unit-normalised
    magnitude_m: float | None = None
    magnitude_sigma: float | None = None

    def resolve(self, Q: np.ndarray, rng: np.random.Generator) -> np.ndarray:
        """Return the 3-vector bias in metres for a baseline with covariance Q."""
        mean_sigma = float(np.sqrt(np.trace(Q) / 3.0))
        if self.magnitude_m is not None:
            mag = float(self.magnitude_m)
        elif self.magnitude_sigma is not None:
            mag = float(self.magnitude_sigma) * mean_sigma
        else:
            mag = None
        if self.vector is not None:
            return np.asarray(self.vector, dtype=float).reshape(3)
        if self.component is not None:
            if mag is None:
                raise ValueError("component outlier needs magnitude_m or magnitude_sigma")
            e = np.zeros(3)
            e[self.component] = mag
            return e
        if mag is None:
            raise ValueError("outlier needs a vector, or a magnitude with component/direction")
        if self.direction is not None:
            d = np.asarray(self.direction, dtype=float).reshape(3)
        else:
            d = rng.standard_normal(3)
        n = np.linalg.norm(d)
        if n == 0:
            raise ValueError("outlier direction must be non-zero")
        return mag * d / n


@dataclass
class NetworkConfig:
    """Everything needed to reproduce a synthetic experiment."""

    n_stations: int = 8
    n_baselines: int | None = None       # default: ~2x spanning tree
    extent_m: float = 10_000.0           # horizontal box size
    height_range_m: float = 300.0
    fixed_station: str | None = None     # default: first station
    noise: NoiseModel = field(default_factory=NoiseModel)
    outliers: list[OutlierSpec] = field(default_factory=list)
    contamination_fraction: float | None = None   # alternative to explicit outliers
    contamination_magnitude_sigma: float = 20.0
    seed: int = 0

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        return d


def generate_true_stations(cfg: NetworkConfig, rng: np.random.Generator) -> list[Station]:
    """Random true coordinates in a box; station ids S00, S01, ..."""
    if cfg.n_stations < 3:
        raise ValueError("need at least 3 stations for a meaningful network")
    ids = [f"S{k:02d}" for k in range(cfg.n_stations)]
    fixed_id = cfg.fixed_station or ids[0]
    if fixed_id not in ids:
        raise ValueError(f"fixed station {fixed_id!r} not among generated ids {ids}")
    xy = rng.uniform(0.0, cfg.extent_m, size=(cfg.n_stations, 2))
    z = rng.uniform(-cfg.height_range_m / 2, cfg.height_range_m / 2, size=cfg.n_stations)
    return [
        Station(station_id=sid, coords=np.array([x, y, h]), is_fixed=(sid == fixed_id))
        for sid, (x, y), h in zip(ids, xy, z)
    ]


def generate_baseline_pairs(
    n_stations: int,
    n_baselines: int | None,
    rng: np.random.Generator,
) -> list[tuple[int, int]]:
    """A connected redundant graph: random spanning tree + extra random edges."""
    n_tree = n_stations - 1
    max_edges = n_stations * (n_stations - 1) // 2
    target = n_baselines if n_baselines is not None else min(max_edges, 2 * n_tree)
    if target < n_tree:
        raise ValueError(f"{target} baselines cannot connect {n_stations} stations")
    if target > max_edges:
        raise ValueError(f"at most {max_edges} distinct pairs exist for {n_stations} stations")

    order = rng.permutation(n_stations)
    edges: set[tuple[int, int]] = set()
    for k in range(1, n_stations):
        a = int(order[k])
        b = int(order[rng.integers(0, k)])       # attach to a random earlier node
        edges.add((min(a, b), max(a, b)))
    all_pairs = [(i, j) for i in range(n_stations) for j in range(i + 1, n_stations)
                 if (i, j) not in edges]
    rng.shuffle(all_pairs)
    for pair in all_pairs:
        if len(edges) >= target:
            break
        edges.add(pair)
    out = sorted(edges)
    rng.shuffle(out)
    return out


def simulate_network(cfg: NetworkConfig) -> tuple[list[Station], list[BaselineObservation], SimulationTruth]:
    """Generate stations, noisy baseline observations, and ground truth.

    Returns ``(stations, observations, truth)``. ``stations`` carry the *true*
    coordinates (the adjustment only reads the fixed station's coordinates
    and, for iterative use, nothing else — the model is linear).
    """
    rng = np.random.default_rng(cfg.seed)
    stations = generate_true_stations(cfg, rng)
    ids = [s.station_id for s in stations]
    pairs = generate_baseline_pairs(cfg.n_stations, cfg.n_baselines, rng)

    # --- resolve contamination fraction into concrete outlier specs
    outliers = list(cfg.outliers)
    if cfg.contamination_fraction is not None:
        n_bad = int(round(cfg.contamination_fraction * len(pairs)))
        bad_ids = rng.choice(len(pairs), size=n_bad, replace=False)
        for bid in bad_ids:
            outliers.append(OutlierSpec(
                baseline_id=int(bid),
                magnitude_sigma=cfg.contamination_magnitude_sigma,
            ))
    outlier_by_id: dict[int, list[OutlierSpec]] = {}
    for spec in outliers:
        outlier_by_id.setdefault(spec.baseline_id, []).append(spec)

    coords = {s.station_id: s.coords for s in stations}
    observations: list[BaselineObservation] = []
    injected: dict[int, np.ndarray] = {}
    for bid, (i, j) in enumerate(pairs):
        sid_i, sid_j = ids[i], ids[j]
        true_vec = coords[sid_j] - coords[sid_i]
        length = float(np.linalg.norm(true_vec))
        sigma = cov.baseline_sigma_m(cfg.noise.a_mm, cfg.noise.b_ppm, length)
        sigmas = cov.component_sigmas(sigma, cfg.noise.axis_scale)
        R = cov.random_correlation(rng, cfg.noise.max_abs_correlation)
        Q_true = cov.build_covariance_block(sigmas, R)

        L = np.linalg.cholesky(Q_true)
        if cfg.noise.student_t_dof is not None:
            dof = float(cfg.noise.student_t_dof)
            if dof <= 2:
                raise ValueError("student_t_dof must exceed 2 for finite covariance")
            # multivariate t with the same covariance: scale so Cov = Q_true
            g = rng.chisquare(dof) / dof
            e = (L @ rng.standard_normal(3)) * np.sqrt((dof - 2) / dof) / np.sqrt(g)
        else:
            e = L @ rng.standard_normal(3)

        obs_vec = true_vec + e
        if bid in outlier_by_id:
            bias = np.zeros(3)
            for spec in outlier_by_id[bid]:
                bias = bias + spec.resolve(Q_true, rng)
            obs_vec = obs_vec + bias
            injected[bid] = bias

        Q_reported = Q_true * cfg.noise.covariance_scale_error
        observations.append(BaselineObservation(
            baseline_id=bid, from_station=sid_i, to_station=sid_j,
            vector=obs_vec, covariance=Q_reported,
        ))

    unresolved = set(outlier_by_id) - {o.baseline_id for o in observations}
    if unresolved:
        raise ValueError(f"outlier specs reference non-existent baselines: {sorted(unresolved)}")

    truth = SimulationTruth(
        true_coordinates={s.station_id: s.coords.copy() for s in stations},
        outlier_baselines=sorted(injected),
        outlier_vectors=injected,
        seed=cfg.seed,
        config=cfg.to_dict(),
    )
    logger.info("simulated network: %d stations, %d baselines, %d outliers (seed=%d)",
                cfg.n_stations, len(observations), len(injected), cfg.seed)
    return stations, observations, truth
