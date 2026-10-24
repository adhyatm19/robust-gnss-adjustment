"""RINEX observation-file metadata, decompression and session-overlap checks.

Only header metadata and the epoch structure are read -- observations are
processed by RTKLIB, never here. Supported: RINEX 2.x and 3.x/4.x observation
files, optionally Hatanaka-compressed (``.crx``/``.??d``) and/or gzip/Unix
compressed (``.gz``/``.Z``). Source files are never modified: decompressed
copies are written to a caller-chosen work directory.
"""

from __future__ import annotations

import datetime as dt
import gzip
import hashlib
import logging
import math
import os
import re
import tempfile
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)


class RinexError(ValueError):
    """Malformed or unsupported RINEX file."""


@dataclass
class RinexObsMeta:
    path: str
    version: float
    marker_name: str
    marker_number: str
    receiver_type: str
    antenna_type: str             # 20-char IGS antenna+radome code, stripped
    approx_xyz: np.ndarray | None
    antenna_delta_hen: np.ndarray  # (H, E, N) metres, as in the header
    header_interval_s: float | None
    time_first: dt.datetime | None
    time_last: dt.datetime | None
    obs_types: dict[str, list[str]]  # system letter (v3) or "*" (v2) -> codes
    n_epochs: int
    epoch_first: dt.datetime | None
    epoch_last: dt.datetime | None
    median_interval_s: float | None
    systems_in_data: list[str]
    epochs: list[dt.datetime] = field(default_factory=list, repr=False)
    warnings: list[str] = field(default_factory=list)

    @property
    def station_id(self) -> str:
        return self.marker_name[:4].upper()

    def summary(self) -> dict:
        d = asdict(self)
        d.pop("epochs")
        for k in ("approx_xyz", "antenna_delta_hen"):
            d[k] = None if d[k] is None else [float(x) for x in d[k]]
        for k in ("time_first", "time_last", "epoch_first", "epoch_last"):
            d[k] = None if d[k] is None else d[k].isoformat()
        return d


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _decompressed_name(name: str) -> str:
    for suffix in (".gz", ".Z"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    if name.endswith(".crx"):
        return name[:-4] + ".rnx"
    m = re.match(r"^(.*\.\d\d)d$", name, flags=re.IGNORECASE)
    if m:
        return m.group(1) + "o"
    return name


def prepare_rinex(src: str | Path, workdir: str | Path) -> Path:
    """Return a plain-text copy of ``src`` in ``workdir`` (decompressing if needed).

    Plain files are copied byte-for-byte so RTKLIB never reads from the
    archive directory and the original is never touched.
    """
    src = Path(src)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    out = workdir / _decompressed_name(src.name)
    if out.exists() and out.stat().st_mtime >= src.stat().st_mtime:
        return out
    data = src.read_bytes()
    name = src.name
    needs_hatanaka = name.endswith((".crx", ".crx.gz", ".crx.Z")) or bool(
        re.search(r"\.\d\dd(\.gz|\.Z)?$", name, flags=re.IGNORECASE))
    if name.endswith(".Z") or needs_hatanaka:
        try:
            import hatanaka  # optional dependency: pip install gnss-adjust[rinex]
        except ImportError as exc:  # pragma: no cover - depends on environment
            raise RinexError(
                f"{src}: Hatanaka/Unix-compress decompression needs the 'hatanaka' "
                "package (pip install -e '.[rinex]')") from exc
        data = hatanaka.decompress(data)
    elif name.endswith(".gz"):
        data = gzip.decompress(data)
    # Atomic publish: concurrent workers decompressing the same file must never
    # see a partially written output (a race observed in a clean checkout).
    fd, tmp = tempfile.mkstemp(dir=workdir, prefix=f".{out.name}.", suffix=".part")
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
        os.replace(tmp, out)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise
    return out


def _parse_time_v3(s: str) -> dt.datetime:
    parts = s.split()
    sec = float(parts[5])
    whole = int(math.floor(sec))
    return dt.datetime(int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3]),
                       int(parts[4]), whole, int(round((sec - whole) * 1e6)) % 1_000_000)


def _epoch_v2(line: str) -> dt.datetime:
    yy = int(line[1:3])
    year = yy + (2000 if yy < 80 else 1900)
    sec = float(line[15:26])
    whole = int(math.floor(sec + 1e-7))
    micro = max(0, int(round((sec - whole) * 1e6)))
    return dt.datetime(year, int(line[4:6]), int(line[7:9]), int(line[10:12]),
                       int(line[13:15]), whole) + dt.timedelta(microseconds=micro)


def read_obs_meta(path: str | Path) -> RinexObsMeta:
    """Read header metadata and scan the epoch structure of a plain RINEX obs file."""
    path = Path(path)
    lines = path.read_text(encoding="ascii", errors="replace").splitlines()
    if not lines or "RINEX VERSION / TYPE" not in lines[0]:
        raise RinexError(f"{path}: first line is not 'RINEX VERSION / TYPE'")
    version = float(lines[0][:9])
    if lines[0][20] not in ("O", "o"):
        raise RinexError(f"{path}: not an observation file (type {lines[0][20]!r})")

    hdr: dict[str, str] = {}
    obs_types: dict[str, list[str]] = {}
    last_sys = None
    warnings: list[str] = []
    i = 0
    for i, line in enumerate(lines):
        label = line[60:].strip()
        if label == "END OF HEADER":
            break
        if label in ("MARKER NAME", "MARKER NUMBER", "REC # / TYPE / VERS", "ANT # / TYPE",
                     "APPROX POSITION XYZ", "ANTENNA: DELTA H/E/N", "INTERVAL",
                     "TIME OF FIRST OBS", "TIME OF LAST OBS"):
            hdr.setdefault(label, line[:60])
        elif label == "ANTENNA: DELTA X/Y/Z":
            warnings.append("antenna eccentricity given as DELTA X/Y/Z (not H/E/N)")
        elif label == "# / TYPES OF OBSERV":       # RINEX 2
            obs_types.setdefault("*", []).extend(line[6:60].split())
        elif label == "SYS / # / OBS TYPES":        # RINEX 3/4
            if line[0] != " ":
                last_sys = line[0]
                obs_types[last_sys] = line[7:60].split()
            elif last_sys:
                obs_types[last_sys].extend(line[7:60].split())
    else:
        raise RinexError(f"{path}: END OF HEADER not found")
    body = lines[i + 1:]

    def fnum(label: str, n: int) -> np.ndarray | None:
        if label not in hdr:
            return None
        try:
            vals = [float(x) for x in hdr[label].split()[:n]]
        except ValueError:
            return None
        return np.array(vals) if len(vals) == n else None

    approx = fnum("APPROX POSITION XYZ", 3)
    if approx is not None and np.linalg.norm(approx) < 1.0:
        warnings.append("APPROX POSITION XYZ is zero")
        approx = None
    delta = fnum("ANTENNA: DELTA H/E/N", 3)
    if delta is None:
        warnings.append("ANTENNA: DELTA H/E/N missing; assumed zero")
        delta = np.zeros(3)
    interval = fnum("INTERVAL", 1)

    def ftime(label: str) -> dt.datetime | None:
        if label not in hdr:
            return None
        try:
            return _parse_time_v3(hdr[label][:43])
        except (ValueError, IndexError):
            return None

    epochs: list[dt.datetime] = []
    systems: set[str] = set()
    if version >= 3.0:
        for line in body:
            if line.startswith(">"):
                flag = int(line[31:32] or 0)
                if flag <= 1:
                    epochs.append(_parse_time_v3(line[2:29]))
            elif line[:1].isalpha():
                systems.add(line[0])
    else:
        ntypes = len(obs_types.get("*", []))
        if ntypes == 0:
            raise RinexError(f"{path}: no '# / TYPES OF OBSERV' in RINEX 2 header")
        lines_per_sat = math.ceil(ntypes / 5)
        k = 0
        while k < len(body):
            line = body[k]
            if len(line) < 32 or not line[28:29].strip().isdigit():
                k += 1
                continue
            flag = int(line[28])
            nsat = int(line[29:32])
            if flag > 1:                 # event: nsat header records follow
                k += 1 + nsat
                continue
            sats = line[32:68]
            extra = math.ceil(nsat / 12) - 1
            for e in range(extra):
                sats += body[k + 1 + e][32:68]
            for j in range(nsat):
                c = sats[3 * j: 3 * j + 1]
                systems.add("G" if c in (" ", "") else c)
            epochs.append(_epoch_v2(line))
            k += 1 + extra + nsat * lines_per_sat

    med = None
    if len(epochs) > 1:
        diffs = np.diff([e.timestamp() for e in epochs])
        med = float(np.median(diffs))
        if np.any(diffs <= 0):
            warnings.append("non-increasing epoch times in file")

    antenna = hdr.get("ANT # / TYPE", " " * 40)[20:40].strip()
    if not antenna:
        warnings.append("antenna type missing from header")
    return RinexObsMeta(
        path=str(path),
        version=version,
        marker_name=hdr.get("MARKER NAME", "").strip(),
        marker_number=hdr.get("MARKER NUMBER", "").strip(),
        receiver_type=hdr.get("REC # / TYPE / VERS", " " * 40)[20:40].strip(),
        antenna_type=antenna,
        approx_xyz=approx,
        antenna_delta_hen=delta,
        header_interval_s=float(interval[0]) if interval is not None else None,
        time_first=ftime("TIME OF FIRST OBS"),
        time_last=ftime("TIME OF LAST OBS"),
        obs_types=obs_types,
        n_epochs=len(epochs),
        epoch_first=epochs[0] if epochs else None,
        epoch_last=epochs[-1] if epochs else None,
        median_interval_s=med,
        systems_in_data=sorted(systems),
        epochs=epochs,
        warnings=warnings,
    )


@dataclass
class OverlapReport:
    start: dt.datetime | None
    end: dt.datetime | None
    overlap_s: float
    common_epochs: int
    fraction_of_shorter: float
    ok: bool
    reason: str


def check_overlap(a: RinexObsMeta, b: RinexObsMeta, min_overlap_s: float = 1800.0,
                  min_fraction: float = 0.5) -> OverlapReport:
    """Check that two stations observed simultaneously for long enough.

    Common epochs are matched to the nearest 0.1 s (receivers may timestamp
    with sub-millisecond offsets, as in the GSI test data).
    """
    if not a.epochs or not b.epochs:
        return OverlapReport(None, None, 0.0, 0, 0.0, False, "one file has no epochs")
    start = max(a.epochs[0], b.epochs[0])
    end = min(a.epochs[-1], b.epochs[-1])
    overlap = max(0.0, (end - start).total_seconds())
    ka = {round(e.timestamp(), 1) for e in a.epochs}
    kb = {round(e.timestamp(), 1) for e in b.epochs}
    common = len(ka & kb)
    shorter = min(len(ka), len(kb))
    frac = common / shorter if shorter else 0.0
    if overlap < min_overlap_s:
        return OverlapReport(start, end, overlap, common, frac, False,
                             f"overlap {overlap:.0f} s < required {min_overlap_s:.0f} s")
    if frac < min_fraction:
        return OverlapReport(start, end, overlap, common, frac, False,
                             f"only {frac:.0%} of epochs are common")
    return OverlapReport(start, end, overlap, common, frac, True, "ok")
