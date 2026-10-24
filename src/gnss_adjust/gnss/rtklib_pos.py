"""Parser for RTKLIB ``.pos`` solution files (ECEF/XYZ format) and covariance reconstruction.

Verified against RTKLIB 2.4.3 b34 source (``src/solution.c``, ``outecef``):

* columns after the time: ``x y z Q ns sdx sdy sdz sdxy sdyz sdzx age ratio``;
* ``sdx, sdy, sdz`` are standard deviations ``sqrt(c_xx)`` etc.;
* ``sdxy, sdyz, sdzx`` are **signed square roots of covariances**,
  ``sqvar(c) = sign(c) * sqrt(|c|)`` -- not covariances, not correlation
  coefficients. Hence ``c_xy = sign(sdxy) * sdxy**2``;
* the internal covariance (``sol_t.qr``) is stored in single precision
  (float32), which bounds any reconstruction at ~6e-8 relative precision;
* the stock binary prints these with ``%8.4f`` (0.1 mm resolution); the
  repository's output-precision patch prints ``%16.10f``;
* the ENU format (``out-solformat=enu``, ``outenu``) prints the baseline
  rotated to the base station's local frame with fields
  ``e n u Q ns sde sdn sdu sden sdnu sdue age ratio`` -- the same
  signed-square-root convention, ordered (en, nu, ue).

``Q`` (solution status): 1 fix, 2 float, 3 sbas, 4 dgps, 5 single, 6 ppp.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

Q_STATUS = {1: "fix", 2: "float", 3: "sbas", 4: "dgps", 5: "single", 6: "ppp"}
ECEF_COLUMNS = ["x", "y", "z", "Q", "ns", "sdx", "sdy", "sdz", "sdxy", "sdyz", "sdzx", "age", "ratio"]


class PosFormatError(ValueError):
    """Raised when a .pos file is not an RTKLIB ECEF solution file."""


@dataclass
class PosRecord:
    time: str                       # as printed (GPST unless the header says otherwise)
    xyz: np.ndarray                 # rover ECEF position (m); ENU baseline for frame "enu"
    q_status: int
    n_sat: int
    sd: np.ndarray                  # (sdx, sdy, sdz) m
    sd_off: np.ndarray              # signed sqrt (sdxy, sdyz, sdzx) m
    age: float
    ratio: float
    resolution_m: float             # print resolution of the sd fields

    @property
    def covariance(self) -> np.ndarray:
        return covariance_from_rtklib(self.sd, self.sd_off)

    @property
    def status_name(self) -> str:
        return Q_STATUS.get(self.q_status, f"unknown({self.q_status})")


@dataclass
class PosFile:
    path: Path
    program: str | None
    input_files: list[str]
    obs_start: str | None
    obs_end: str | None
    ref_pos: np.ndarray | None      # base position used by RTKLIB (m, ECEF; xyz format only)
    frame: str                      # "xyz" (rover ECEF position) or "enu" (baseline in base ENU)
    options: dict[str, str]         # echoed processing options (out-outopt=on)
    records: list[PosRecord] = field(default_factory=list)


def signed_square(s: np.ndarray | float) -> np.ndarray:
    """Inverse of RTKLIB ``sqvar``: ``sign(s) * s**2``."""
    s = np.asarray(s, dtype=float)
    return np.sign(s) * s * s


def covariance_from_rtklib(sd: np.ndarray, sd_off: np.ndarray) -> np.ndarray:
    """3x3 ECEF covariance from RTKLIB ``(sdx, sdy, sdz)`` and ``(sdxy, sdyz, sdzx)``."""
    sx, sy, sz = np.asarray(sd, dtype=float).reshape(3)
    cxy, cyz, czx = signed_square(np.asarray(sd_off, dtype=float).reshape(3))
    return np.array([
        [sx * sx, cxy, czx],
        [cxy, sy * sy, cyz],
        [czx, cyz, sz * sz],
    ])


def _decimals(token: str) -> int:
    return len(token.split(".")[1]) if "." in token else 0


def parse_pos_file(path: str | Path, frame: str = "xyz") -> PosFile:
    """Parse an RTKLIB ``.pos`` file written with ``out-solformat=xyz`` (or ``enu``)."""
    if frame not in ("xyz", "enu"):
        raise ValueError(f"unsupported frame {frame!r}")
    marker = ("x-ecef(m)", "sdxy(m)") if frame == "xyz" else ("e-baseline(m)", "sden(m)")
    path = Path(path)
    program = obs_start = obs_end = None
    ref_pos = None
    inputs: list[str] = []
    options: dict[str, str] = {}
    columns_ok = False
    records: list[PosRecord] = []
    with open(path, encoding="ascii", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            if line.startswith("%"):
                body = line[1:].strip()
                if ":" in body and not body.startswith("("):
                    key, _, value = body.partition(":")
                    key, value = key.strip(), value.strip()
                    if key == "program":
                        program = value
                    elif key == "inp file":
                        inputs.append(value)
                    elif key == "obs start":
                        obs_start = value
                    elif key == "obs end":
                        obs_end = value
                    elif key == "ref pos" and frame == "xyz":
                        ref_pos = np.array([float(t) for t in value.split()[:3]])
                    elif key:
                        options[key] = value
                if marker[0] in body and marker[1] in body:
                    columns_ok = True
                continue
            if not columns_ok:
                raise PosFormatError(
                    f"{path}: solution lines before a {frame} column header; "
                    f"RTKLIB must be run with out-solformat={frame} and out-outhead=on")
            tokens = line.split()
            if len(tokens) != 2 + len(ECEF_COLUMNS):
                raise PosFormatError(f"{path}: unexpected column count in line: {line!r}")
            t = " ".join(tokens[:2])
            v = tokens[2:]
            records.append(PosRecord(
                time=t,
                xyz=np.array([float(v[0]), float(v[1]), float(v[2])]),
                q_status=int(v[3]),
                n_sat=int(v[4]),
                sd=np.array([float(x) for x in v[5:8]]),
                sd_off=np.array([float(x) for x in v[8:11]]),
                age=float(v[11]),
                ratio=float(v[12]),
                resolution_m=10.0 ** (-min(_decimals(x) for x in v[5:11])),
            ))
    if not columns_ok:
        raise PosFormatError(f"{path}: no RTKLIB {frame} column header found")
    return PosFile(path=path, program=program, input_files=inputs, obs_start=obs_start,
                   obs_end=obs_end, ref_pos=ref_pos, frame=frame, options=options,
                   records=records)


@dataclass
class CovarianceCheck:
    """Validity of a reconstructed covariance and its sensitivity to print rounding."""

    symmetric: bool
    positive_definite: bool
    min_eigenvalue: float
    condition_number: float
    max_abs_correlation: float
    resolution_m: float
    # Envelope over all matrices consistent with the printed (rounded) fields:
    rounding_nonpd_fraction: float
    rounding_max_rel_sigma_change_enu: float
    rounding_max_abs_corr_change_enu: float
    rounding_sensitive: bool


def check_covariance(record: PosRecord, ref_xyz: np.ndarray, n_samples: int = 2000,
                     seed: int = 0, sensitive_rel: float = 0.05,
                     sensitive_corr: float = 0.05) -> CovarianceCheck:
    """Validate ``record.covariance`` and quantify rounding sensitivity.

    Every printed field is only known to within +-resolution/2. Sampling that
    box (deterministically seeded) gives the set of covariance matrices the
    printed numbers cannot distinguish; the matrix is flagged
    ``rounding_sensitive`` if, within that set, any local ENU standard
    deviation changes by more than ``sensitive_rel`` (relative) or any ENU
    correlation by more than ``sensitive_corr`` (absolute), or if the set
    contains non-positive-definite matrices.
    """
    from .frames import covariance_ecef_to_enu  # local import keeps module light

    q = record.covariance
    eig = np.linalg.eigvalsh(q)
    pd = bool(eig[0] > 0.0)
    sig = np.sqrt(np.diag(q))
    corr = q / np.outer(sig, sig)
    max_corr = float(np.max(np.abs(corr[np.triu_indices(3, 1)])))

    def enu_stats(qm: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        qe = covariance_ecef_to_enu(qm, ref_xyz)
        s = np.sqrt(np.clip(np.diag(qe), 0.0, None))
        c = qe / np.outer(s, s)
        return s, c[np.triu_indices(3, 1)]

    s0, c0 = enu_stats(q)
    rng = np.random.default_rng(seed)
    half = 0.5 * record.resolution_m
    nonpd = 0
    max_rel = 0.0
    max_dc = 0.0
    for _ in range(n_samples):
        sd = record.sd + rng.uniform(-half, half, 3)
        so = record.sd_off + rng.uniform(-half, half, 3)
        qs = covariance_from_rtklib(np.abs(sd), so)
        if np.linalg.eigvalsh(qs)[0] <= 0.0:
            nonpd += 1
            continue
        s, c = enu_stats(qs)
        max_rel = max(max_rel, float(np.max(np.abs(s / s0 - 1.0))))
        max_dc = max(max_dc, float(np.max(np.abs(c - c0))))
    frac = nonpd / n_samples
    return CovarianceCheck(
        symmetric=True,  # symmetric by construction
        positive_definite=pd,
        min_eigenvalue=float(eig[0]),
        condition_number=float(eig[-1] / eig[0]) if pd else float("inf"),
        max_abs_correlation=max_corr,
        resolution_m=record.resolution_m,
        rounding_nonpd_fraction=frac,
        rounding_max_rel_sigma_change_enu=max_rel,
        rounding_max_abs_corr_change_enu=max_dc,
        rounding_sensitive=bool(frac > 0 or max_rel > sensitive_rel or max_dc > sensitive_corr),
    )
