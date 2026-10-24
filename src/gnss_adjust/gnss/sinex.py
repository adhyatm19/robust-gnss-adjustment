"""Minimal SINEX 2.x reader: station positions, their sigmas, antennas and eccentricities.

Used to obtain independent reference coordinates (e.g. the EPN combined
weekly solution) and to cross-check RINEX antenna metadata. Positions refer
to the station marker; ``SITE/ECCENTRICITY`` gives the antenna reference
point offset (up, north, east) from the marker.
"""

from __future__ import annotations

import datetime as dt
import gzip
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class SinexStation:
    code: str
    xyz: np.ndarray
    sigma_xyz: np.ndarray
    ref_epoch: dt.datetime
    antenna_type: str | None = None          # 20-char IGS code incl. radome, stripped
    eccentricity_une: np.ndarray | None = None


@dataclass
class SinexSolution:
    path: str
    header: str
    stations: dict[str, SinexStation] = field(default_factory=dict)


def sinex_epoch(s: str) -> dt.datetime:
    """``YY:DOY:SSSSS`` -> datetime (YY < 50 -> 20YY)."""
    yy, doy, sod = (int(x) for x in s.strip().split(":"))
    year = 2000 + yy if yy < 50 else 1900 + yy
    return dt.datetime(year, 1, 1) + dt.timedelta(days=doy - 1, seconds=sod)


def read_sinex(path: str | Path, stations: list[str] | None = None) -> SinexSolution:
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    wanted = {s.upper() for s in stations} if stations else None
    est: dict[str, dict[str, tuple[float, float, dt.datetime]]] = {}
    ant: dict[str, str] = {}
    ecc: dict[str, np.ndarray] = {}
    block = None
    header = ""
    with opener(path, "rt", encoding="ascii", errors="replace") as fh:
        for line in fh:
            if line.startswith("%=SNX"):
                header = line.strip()
                continue
            if line.startswith("+"):
                block = line[1:].strip()
                continue
            if line.startswith("-"):
                block = None
                continue
            if line.startswith("*") or block is None:
                continue
            # site code: columns 2-5 in SITE/* blocks, 15-18 in SOLUTION/ESTIMATE
            code = (line[14:18] if block == "SOLUTION/ESTIMATE" else line[1:5]).upper()
            if wanted is not None and code not in wanted:
                continue
            if block == "SOLUTION/ESTIMATE" and line[7:11] in ("STAX", "STAY", "STAZ"):
                est.setdefault(code, {})[line[10]] = (
                    float(line[47:68]), float(line[69:80]), sinex_epoch(line[27:39]))
            elif block == "SITE/ANTENNA":
                ant[code] = line[42:62].strip()
            elif block == "SITE/ECCENTRICITY":
                ecc[code] = np.array([float(x) for x in line[46:].split()[:3]])
    sol = SinexSolution(path=str(path), header=header)
    for code, comp in est.items():
        if set(comp) != {"X", "Y", "Z"}:
            continue
        sol.stations[code] = SinexStation(
            code=code,
            xyz=np.array([comp[c][0] for c in "XYZ"]),
            sigma_xyz=np.array([comp[c][1] for c in "XYZ"]),
            ref_epoch=comp["X"][2],
            antenna_type=ant.get(code),
            eccentricity_une=ecc.get(code),
        )
    return sol
