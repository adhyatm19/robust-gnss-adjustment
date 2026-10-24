"""RINEX pair -> one validated static baseline vector with covariance and provenance.

Orientation convention: RTKLIB's base receiver is the **from**-station and
its rover the **to**-station, so the baseline is

    b = X_rover - X_base     (ECEF metres)

which is exactly the repository's observation model ``l_ij = X_j - X_i``.
The covariance is RTKLIB's formal rover-position covariance with the base
position held fixed, i.e. the formal covariance of ``b`` *conditional on*
the base coordinates and (for fixed solutions) on the fixed integer
ambiguities. It is not an externally validated accuracy.
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from . import frames
from .rinex import RinexObsMeta, check_overlap, prepare_rinex, read_obs_meta
from .rtklib import AntennaSetting, ProcessingConfig, RunRecord, run_rnx2rtkp
from .rtklib_pos import CovarianceCheck, check_covariance, parse_pos_file

logger = logging.getLogger(__name__)

USABLE_STATUS = {1: "fix", 2: "float"}


class BaselineError(RuntimeError):
    """A baseline could not be produced or failed validation."""


@dataclass
class BaselineSolution:
    session_id: str
    from_station: str                  # base
    to_station: str                    # rover
    vector_ecef: np.ndarray            # X_to - X_from (m)
    covariance_ecef: np.ndarray        # formal (m^2)
    status: str                        # fix | float
    ratio: float
    n_sat: int
    base_xyz: np.ndarray
    rover_xyz: np.ndarray
    covariance_source: str
    covariance_check: CovarianceCheck
    overlap_s: float
    common_epochs: int
    selection: str                     # "single" (RTKLIB solstatic) or "last" (final forward epoch)
    epochs_solved: int                 # solution epochs RTKLIB wrote
    fix_fraction: float                # share of written epochs with status fix
    rover_antenna: AntennaSetting
    base_antenna: AntennaSetting
    config: dict
    run: RunRecord
    rover_meta: dict
    base_meta: dict
    warnings: list[str] = field(default_factory=list)

    @property
    def length_m(self) -> float:
        return float(np.linalg.norm(self.vector_ecef))

    def vector_enu(self) -> np.ndarray:
        return frames.vector_ecef_to_enu(self.vector_ecef, self.base_xyz)

    def covariance_enu(self) -> np.ndarray:
        return frames.covariance_ecef_to_enu(self.covariance_ecef, self.base_xyz)

    def to_record(self) -> dict:
        """Flat JSON-serialisable record (for the processing log)."""
        d = asdict(self)
        for k in ("vector_ecef", "covariance_ecef", "base_xyz", "rover_xyz"):
            d[k] = np.asarray(d[k]).tolist()
        d["length_m"] = self.length_m
        d["vector_enu"] = self.vector_enu().tolist()
        d["covariance_enu"] = self.covariance_enu().tolist()
        return d


def antenna_from_header(meta: RinexObsMeta, delta_up_offset_m: float = 0.0) -> AntennaSetting:
    """Explicit RTKLIB antenna setting reproducing the RINEX header.

    RINEX gives ``ANTENNA: DELTA H/E/N``; RTKLIB's ``antdel`` is ordered
    (east, north, up) with up = H. ``delta_up_offset_m`` deliberately
    perturbs the height (for controlled antenna-height experiments) without
    touching the RINEX file.
    """
    h, e, n = (float(x) for x in meta.antenna_delta_hen)
    return AntennaSetting(mode="explicit", antenna_type=meta.antenna_type,
                          delta_enu=(e, n, h + delta_up_offset_m))


def process_baseline(
    *,
    session_id: str,
    rover_obs: Path,
    base_obs: Path,
    aux_files: list[Path],
    base_xyz: np.ndarray,
    workdir: Path,
    cfg: ProcessingConfig,
    rover_station: str | None = None,
    base_station: str | None = None,
    rover_antenna: AntennaSetting | None = None,
    base_antenna: AntennaSetting | None = None,
    variant: str = "hp",
    min_overlap_s: float = 1800.0,
    allow_float: bool = True,
    selection: str = "single",
    min_epoch_fraction: float = 0.9,
) -> BaselineSolution:
    """Decompress inputs, check metadata/overlap, run RTKLIB, parse and validate.

    ``selection="single"`` expects RTKLIB ``out-solstatic=single`` (one line).
    ``selection="last"`` expects all epochs (``out-solstatic=all``) of a
    *forward* static filter and takes the final epoch: for a static position
    that is the estimate from all data. (With ``combined`` +
    ``solstatic=single`` RTKLIB picks the last epoch with the best *status*,
    which near the session end can be an unsmoothed backward solution built
    from minutes of data -- observed on the Madrid data.) Every written epoch
    is cumulative, so epochs are never used as separate observations.
    ``min_epoch_fraction`` guards against RTKLIB silently skipping epochs.
    """
    if selection not in ("single", "last"):
        raise ValueError(f"unknown selection {selection!r}")
    if selection == "last" and (cfg.solstatic != "all" or cfg.solution_type != "forward"):
        raise ValueError("selection='last' requires solstatic='all' and solution_type='forward'")
    workdir = Path(workdir)
    rnx_dir = workdir / "rinex"
    rov = prepare_rinex(rover_obs, rnx_dir)
    bas = prepare_rinex(base_obs, rnx_dir)
    aux = [prepare_rinex(a, rnx_dir) if str(a).endswith((".gz", ".Z")) else Path(a)
           for a in aux_files]
    mr, mb = read_obs_meta(rov), read_obs_meta(bas)
    rover_station = rover_station or mr.station_id
    base_station = base_station or mb.station_id
    warnings = [f"{rover_station}: {w}" for w in mr.warnings]
    warnings += [f"{base_station}: {w}" for w in mb.warnings]

    ov = check_overlap(mr, mb, min_overlap_s=min_overlap_s)
    if not ov.ok:
        raise BaselineError(f"{base_station}->{rover_station} ({session_id}): {ov.reason}")

    rover_antenna = rover_antenna or antenna_from_header(mr)
    base_antenna = base_antenna or antenna_from_header(mb)
    out_pos = workdir / "solutions" / f"{session_id}_{base_station}_{rover_station}_{variant}.pos"
    run = run_rnx2rtkp(rov, bas, aux, out_pos, cfg, base_xyz, rover_antenna, base_antenna,
                       variant=variant)
    pos = parse_pos_file(out_pos)
    if not pos.records:
        raise BaselineError(f"{out_pos}: RTKLIB wrote no solution ({run.stderr_tail})")
    if selection == "single" and len(pos.records) != 1:
        raise BaselineError(f"{out_pos}: expected one static solution line "
                            f"(out-solstatic=single), got {len(pos.records)}")
    rec = pos.records[-1]
    n_written = len(pos.records)
    fix_fraction = float(np.mean([r.q_status == 1 for r in pos.records]))
    if selection == "last" and n_written < min_epoch_fraction * ov.common_epochs:
        raise BaselineError(f"{out_pos}: RTKLIB solved only {n_written} of {ov.common_epochs} "
                            "common epochs (check excluded/unhealthy satellites)")
    if rec.q_status not in USABLE_STATUS:
        raise BaselineError(f"{out_pos}: solution status {rec.status_name} is not a "
                            "carrier-phase relative solution")
    if rec.q_status == 2 and not allow_float:
        raise BaselineError(f"{out_pos}: ambiguities not fixed (float) and allow_float=False")
    if pos.ref_pos is None or np.max(np.abs(pos.ref_pos - base_xyz)) > max(rec.resolution_m, 1e-6):
        raise BaselineError(f"{out_pos}: RTKLIB reference position {pos.ref_pos} differs from "
                            f"the supplied base position {base_xyz}")
    if cfg.antex_file and cfg.rcv_pcv:
        # RTKLIB echoes the matched receiver antenna ("antenna1/2" header lines);
        # an empty type means no calibration was found and none was applied.
        for key, ant in (("antenna1", rover_antenna), ("antenna2", base_antenna)):
            echoed = pos.options.get(key, "").split("(")[0].strip()
            if ant.mode == "explicit" and not echoed:
                warnings.append(f"{key}: no ANTEX calibration found for {ant.antenna_type!r}; "
                                "receiver PCO/PCV not applied")
    check = check_covariance(rec, base_xyz)
    if not check.positive_definite:
        raise BaselineError(f"{out_pos}: reconstructed covariance is not positive definite")
    source = ("rtklib_formal_float32" if rec.resolution_m <= 1e-9
              else f"rtklib_formal_rounded_{rec.resolution_m:g}m")
    if check.rounding_sensitive:
        warnings.append("covariance is sensitive to print rounding "
                        f"(max ENU sigma change {check.rounding_max_rel_sigma_change_enu:.1%}, "
                        f"max ENU corr change {check.rounding_max_abs_corr_change_enu:.2f})")
    return BaselineSolution(
        session_id=session_id,
        from_station=base_station,
        to_station=rover_station,
        vector_ecef=rec.xyz - np.asarray(base_xyz, dtype=float),
        covariance_ecef=rec.covariance,
        status=USABLE_STATUS[rec.q_status],
        ratio=rec.ratio,
        n_sat=rec.n_sat,
        base_xyz=np.asarray(base_xyz, dtype=float),
        rover_xyz=rec.xyz,
        covariance_source=source,
        covariance_check=check,
        overlap_s=ov.overlap_s,
        common_epochs=ov.common_epochs,
        selection=selection,
        epochs_solved=n_written,
        fix_fraction=fix_fraction,
        rover_antenna=rover_antenna,
        base_antenna=base_antenna,
        config=cfg.to_dict(),
        run=run,
        rover_meta=mr.summary(),
        base_meta=mb.summary(),
        warnings=warnings,
    )
