"""Invocation of RTKLIB ``rnx2rtkp`` for static relative baseline processing.

Design decisions (verified against RTKLIB 2.4.3 b34 source):

* **Configuration only.** All processing options are written to a config
  file passed with ``-k``; no other command-line processing options are
  used, because rnx2rtkp command-line options take precedence over the
  config file and would make the recorded configuration incomplete.
* **Antenna metadata precedence** (``postpos.c``, ``setpcv``): when
  ``antN-anttype = *`` the antenna type *and* the ENU antenna delta are
  taken from the RINEX header, **overwriting** any ``antN-antdel*`` in the
  config. With an explicit antenna type the config ``antN-antdel*`` values
  are used and the RINEX header delta is ignored. This module therefore
  writes an *explicit* type and delta by default (taken from the parsed
  RINEX header, so the applied values are known and logged); a perturbation
  experiment changes only the delta it writes.
* **Receiver numbering:** ``ant1`` is the rover (first obs file), ``ant2``
  the base (second obs file). The base position is always given explicitly
  (``ant2-postype = xyz``); RTKLIB treats it as error-free, so the reported
  covariance is that of the rover position *given* the base, i.e. of the
  baseline vector.
* **Static output:** ``out-solstatic = single`` makes RTKLIB print one line:
  the last epoch with the best solution status (fix > float > ...). The
  printed time is the *first* epoch of the session (postpos.c), so it must
  not be read as the solution epoch.
"""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

NAVSYS_BITS = {"G": 1, "S": 2, "R": 4, "E": 8, "J": 16, "C": 32, "I": 64}


class RTKLIBError(RuntimeError):
    """Raised when rnx2rtkp is unavailable or fails."""


def _repo_tools_dir() -> Path | None:
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "pyproject.toml").exists() and (parent / "src" / "gnss_adjust").exists():
            return parent / "tools" / "rtklib" / "bin"
    return None


def find_rnx2rtkp(variant: str = "hp") -> Path:
    """Locate the rnx2rtkp binary.

    ``variant`` is ``"hp"`` (output-precision patch) or ``"stock"``. Search
    order: ``$GNSS_ADJUST_RTKLIB_BIN/rnx2rtkp-<variant>``, then
    ``<repo>/tools/rtklib/bin/rnx2rtkp-<variant>`` (built by
    ``scripts/build_rtklib.sh``). A system ``rnx2rtkp`` is accepted only for
    ``variant="system"``, so an unknown build is never used silently.
    """
    if variant == "system":
        found = shutil.which("rnx2rtkp")
        if not found:
            raise RTKLIBError("no rnx2rtkp on PATH")
        return Path(found)
    name = f"rnx2rtkp-{variant}"
    candidates = []
    if os.environ.get("GNSS_ADJUST_RTKLIB_BIN"):
        candidates.append(Path(os.environ["GNSS_ADJUST_RTKLIB_BIN"]) / name)
    tools = _repo_tools_dir()
    if tools is not None:
        candidates.append(tools / name)
    for c in candidates:
        if c.is_file() and os.access(c, os.X_OK):
            return c
    raise RTKLIBError(
        f"{name} not found (looked in {[str(c) for c in candidates]}); "
        "build it with scripts/build_rtklib.sh or set GNSS_ADJUST_RTKLIB_BIN")


def binary_info(binary: Path) -> dict:
    """Version string reported by the binary plus its SHA-256."""
    h = hashlib.sha256(binary.read_bytes()).hexdigest()
    path = str(binary)
    tools = _repo_tools_dir()
    if tools is not None:
        try:  # record repository-relative paths so logs are machine-independent
            path = str(binary.resolve().relative_to(tools.parents[2]))
        except ValueError:
            pass
    info = {"path": path, "sha256": h, "variant": binary.name}
    build_info = binary.parent.parent / "build_info.json"
    if build_info.exists():
        import json
        info["build_info"] = json.loads(build_info.read_text())
    return info


@dataclass
class AntennaSetting:
    """Antenna model applied by RTKLIB for one receiver.

    ``mode="explicit"`` writes ``antenna_type`` and ``delta_enu`` (metres,
    east/north/up from marker to antenna reference point) to the config.
    ``mode="rinex"`` writes ``*`` and lets RTKLIB read both from the header.
    """

    mode: str = "explicit"
    antenna_type: str = ""
    delta_enu: tuple[float, float, float] = (0.0, 0.0, 0.0)


@dataclass
class ProcessingConfig:
    pos_mode: str = "static"
    frequency: str = "l1+2"
    solution_type: str = "combined"      # forward | backward | combined
    elevation_mask_deg: float = 10.0
    nav_systems: tuple[str, ...] = ("G",)
    iono: str = "brdc"                   # brdc (short baselines) | est-stec | dual-freq
    tropo: str = "saas"                  # saas | est-ztd
    ephemeris: str = "brdc"              # brdc | precise (needs SP3 among inputs)
    tide_correction: str = "off"         # off | on
    ar_mode: str = "continuous"          # continuous | fix-and-hold | instantaneous | off
    ar_threshold: float = 3.0
    sat_pcv: bool = False                # requires satellite ANTEX
    rcv_pcv: bool = False                # requires receiver ANTEX
    phase_windup: bool = True
    antex_file: str | None = None
    solstatic: str = "single"            # single | all
    extra: dict[str, str] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)


def build_config_text(cfg: ProcessingConfig, base_xyz: np.ndarray,
                      rover_antenna: AntennaSetting, base_antenna: AntennaSetting) -> str:
    navsys = sum(NAVSYS_BITS[s] for s in cfg.nav_systems)
    lines = {
        "pos1-posmode": cfg.pos_mode,
        "pos1-frequency": cfg.frequency,
        "pos1-soltype": cfg.solution_type,
        "pos1-elmask": f"{cfg.elevation_mask_deg:g}",
        "pos1-dynamics": "off",
        "pos1-tidecorr": cfg.tide_correction,
        "pos1-ionoopt": cfg.iono,
        "pos1-tropopt": cfg.tropo,
        "pos1-sateph": cfg.ephemeris,
        "pos1-posopt1": "on" if cfg.sat_pcv else "off",
        "pos1-posopt2": "on" if cfg.rcv_pcv else "off",
        "pos1-posopt3": "on" if cfg.phase_windup else "off",
        "pos1-navsys": str(navsys),
        "pos2-armode": cfg.ar_mode,
        "pos2-arthres": f"{cfg.ar_threshold:g}",
        "out-solformat": "xyz",
        "out-outhead": "on",
        "out-outopt": "on",
        "out-timesys": "gpst",
        "out-timeform": "hms",
        "out-timendec": "3",
        "out-solstatic": cfg.solstatic,
        "ant2-postype": "xyz",
        "ant2-pos1": f"{base_xyz[0]:.6f}",
        "ant2-pos2": f"{base_xyz[1]:.6f}",
        "ant2-pos3": f"{base_xyz[2]:.6f}",
    }
    for n, ant in (("1", rover_antenna), ("2", base_antenna)):
        if ant.mode == "rinex":
            lines[f"ant{n}-anttype"] = "*"
        elif ant.mode == "explicit":
            lines[f"ant{n}-anttype"] = ant.antenna_type
            lines[f"ant{n}-antdele"] = f"{ant.delta_enu[0]:.6f}"
            lines[f"ant{n}-antdeln"] = f"{ant.delta_enu[1]:.6f}"
            lines[f"ant{n}-antdelu"] = f"{ant.delta_enu[2]:.6f}"
        else:
            raise ValueError(f"unknown antenna mode {ant.mode!r}")
    if cfg.antex_file:
        lines["file-satantfile"] = cfg.antex_file
        lines["file-rcvantfile"] = cfg.antex_file
    lines.update(cfg.extra)
    return "".join(f"{k:<18}={v}\n" for k, v in lines.items())


@dataclass
class RunRecord:
    command: list[str]
    config_file: str
    config_text: str
    output_file: str
    returncode: int
    elapsed_s: float
    stderr_tail: str
    binary: dict


def run_rnx2rtkp(rover_obs: Path, base_obs: Path, aux_files: list[Path], out_pos: Path,
                 cfg: ProcessingConfig, base_xyz: np.ndarray,
                 rover_antenna: AntennaSetting, base_antenna: AntennaSetting,
                 variant: str = "hp", timeout_s: float = 3600.0) -> RunRecord:
    """Run one static relative baseline solution; raise on failure."""
    binary = find_rnx2rtkp(variant)
    out_pos = Path(out_pos)
    out_pos.parent.mkdir(parents=True, exist_ok=True)
    conf = out_pos.with_suffix(".conf")
    text = build_config_text(cfg, np.asarray(base_xyz, dtype=float),
                             rover_antenna, base_antenna)
    conf.write_text(text)
    cmd = [str(binary), "-k", str(conf), "-o", str(out_pos),
           str(rover_obs), str(base_obs), *map(str, aux_files)]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s)
    elapsed = time.perf_counter() - t0
    # rnx2rtkp writes a carriage-return progress line to stderr; keep the tail.
    tail = proc.stderr.replace("\r", "\n").strip().splitlines()[-3:]
    record = RunRecord(command=cmd, config_file=str(conf), config_text=text,
                       output_file=str(out_pos), returncode=proc.returncode,
                       elapsed_s=elapsed, stderr_tail="\n".join(tail),
                       binary=binary_info(binary))
    if proc.returncode != 0:
        raise RTKLIBError(f"rnx2rtkp failed ({proc.returncode}): {record.stderr_tail}")
    if not out_pos.exists() or out_pos.stat().st_size == 0:
        raise RTKLIBError(f"rnx2rtkp produced no output: {record.stderr_tail}")
    return record
