"""CSV / JSON / YAML input-output and validation.

CSV formats (all headers required, coordinates and observations in metres,
covariances in m^2):

stations.csv:            station_id, x, y, z, is_fixed
baselines.csv:           baseline_id, from_station, to_station, dx, dy, dz
covariance_blocks.csv:   baseline_id, q_xx, q_xy, q_xz, q_yy, q_yz, q_zz
"""

from __future__ import annotations

import dataclasses
import json
import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from . import covariance as cov
from . import network as net
from .types import BaselineObservation, Station

logger = logging.getLogger(__name__)

STATION_COLUMNS = ["station_id", "x", "y", "z", "is_fixed"]
BASELINE_COLUMNS = ["baseline_id", "from_station", "to_station", "dx", "dy", "dz"]
COV_COLUMNS = ["baseline_id", "q_xx", "q_xy", "q_xz", "q_yy", "q_yz", "q_zz"]


class InputError(ValueError):
    """Raised for malformed or inconsistent input files."""


def _require_columns(df: pd.DataFrame, columns: list[str], name: str) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise InputError(f"{name}: missing columns {missing}")


def read_stations_csv(path: str | Path) -> list[Station]:
    df = pd.read_csv(path)
    _require_columns(df, STATION_COLUMNS, "stations.csv")
    stations = []
    for _, row in df.iterrows():
        stations.append(Station(
            station_id=str(row["station_id"]),
            coords=np.array([row["x"], row["y"], row["z"]], dtype=float),
            is_fixed=bool(row["is_fixed"]),
        ))
    if sum(s.is_fixed for s in stations) == 0:
        raise InputError("stations.csv: no station has is_fixed=1 (datum undefined)")
    return stations


def read_baselines_csv(baselines_path: str | Path, covariance_path: str | Path) -> list[BaselineObservation]:
    bdf = pd.read_csv(baselines_path)
    cdf = pd.read_csv(covariance_path)
    _require_columns(bdf, BASELINE_COLUMNS, "baselines.csv")
    _require_columns(cdf, COV_COLUMNS, "covariance_blocks.csv")
    cov_by_id: dict[int, np.ndarray] = {}
    for _, r in cdf.iterrows():
        bid = int(r["baseline_id"])
        if bid in cov_by_id:
            raise InputError(f"covariance_blocks.csv: duplicate baseline_id {bid}")
        Q = np.array([
            [r["q_xx"], r["q_xy"], r["q_xz"]],
            [r["q_xy"], r["q_yy"], r["q_yz"]],
            [r["q_xz"], r["q_yz"], r["q_zz"]],
        ], dtype=float)
        cov_by_id[bid] = Q

    observations: list[BaselineObservation] = []
    seen: set[int] = set()
    for _, r in bdf.iterrows():
        bid = int(r["baseline_id"])
        if bid in seen:
            raise InputError(f"baselines.csv: duplicate baseline_id {bid}")
        seen.add(bid)
        if bid not in cov_by_id:
            raise InputError(f"baseline {bid}: no covariance block provided")
        Q = cov_by_id[bid]
        try:
            cov.validate_covariance_block(Q, name=f"Q(baseline {bid})")
        except cov.CovarianceError as exc:
            raise InputError(str(exc)) from exc
        observations.append(BaselineObservation(
            baseline_id=bid,
            from_station=str(r["from_station"]),
            to_station=str(r["to_station"]),
            vector=np.array([r["dx"], r["dy"], r["dz"]], dtype=float),
            covariance=Q,
        ))
    extra = set(cov_by_id) - seen
    if extra:
        logger.warning("covariance_blocks.csv has blocks for unknown baselines: %s", sorted(extra))
    return observations


def load_network_from_csv(
    stations_path: str | Path,
    baselines_path: str | Path,
    covariance_path: str | Path,
) -> tuple[list[Station], list[BaselineObservation]]:
    """Load and fully validate a network from the three CSV files."""
    stations = read_stations_csv(stations_path)
    observations = read_baselines_csv(baselines_path, covariance_path)
    net.validate_observations(stations, observations)
    if not net.is_connected(stations, observations):
        raise InputError("network graph is not connected")
    return stations, observations


def write_network_csv(
    stations: list[Station],
    observations: list[BaselineObservation],
    out_dir: str | Path,
) -> None:
    """Write stations/baselines/covariance CSVs in the documented format."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([
        {"station_id": s.station_id, "x": s.coords[0], "y": s.coords[1],
         "z": s.coords[2], "is_fixed": int(s.is_fixed)}
        for s in stations
    ]).to_csv(out / "stations.csv", index=False)
    pd.DataFrame([
        {"baseline_id": o.baseline_id, "from_station": o.from_station,
         "to_station": o.to_station, "dx": o.vector[0], "dy": o.vector[1],
         "dz": o.vector[2]}
        for o in observations
    ]).to_csv(out / "baselines.csv", index=False)
    pd.DataFrame([
        {"baseline_id": o.baseline_id,
         "q_xx": o.covariance[0, 0], "q_xy": o.covariance[0, 1],
         "q_xz": o.covariance[0, 2], "q_yy": o.covariance[1, 1],
         "q_yz": o.covariance[1, 2], "q_zz": o.covariance[2, 2]}
        for o in observations
    ]).to_csv(out / "covariance_blocks.csv", index=False)


class _NumpyEncoder(json.JSONEncoder):
    def default(self, o: Any) -> Any:
        if isinstance(o, Path):
            return str(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, (np.bool_,)):
            return bool(o)
        if dataclasses.is_dataclass(o) and not isinstance(o, type):
            return dataclasses.asdict(o)
        return super().default(o)


def save_json(obj: Any, path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, cls=_NumpyEncoder)


def load_yaml(path: str | Path) -> dict[str, Any]:
    with open(path) as fh:
        return yaml.safe_load(fh)
