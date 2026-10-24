"""Calibrate covariances and assemble the real Madrid networks from processed baselines.

Input:  results/real_network/processing/baselines_all.json (scripts/process_real_network.py)
Output: results/real_network/calibration.json and networks/<name>/ in the
        repository CSV format plus sidecars (see gnss_adjust.gnss.network_build).

Design (fixed before looking at adjustment results):
* calibration sessions: days 001-003 (all six pairs, repeatability);
* evaluation sessions: days 004-006, one *different* spanning tree per day,
  so that every loop closes across sessions (no trivial baselines);
* covariance variants: formal RTKLIB covariance, and formal x k^2 with k^2
  estimated on the calibration days only (hold-out, no circularity);
* a deliberately trivial network (all six pairs of day 004) is built only to
  measure what same-session loops look like.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gnss_adjust.gnss.calibration import repeatability_calibration  # noqa: E402
from gnss_adjust.gnss.network_build import build_network, write_network  # noqa: E402
from gnss_adjust.gnss.sinex import read_sinex  # noqa: E402

CAL_DAYS = ["2024001", "2024002", "2024003"]
EVAL_TREES = {                       # (from, to) per evaluation session; each is a spanning tree
    "2024004": [("MADR", "CEBR"), ("MADR", "VILL"), ("VILL", "YEBE")],
    "2024005": [("CEBR", "VILL"), ("MADR", "CEBR"), ("MADR", "YEBE")],
    "2024006": [("CEBR", "YEBE"), ("CEBR", "VILL"), ("MADR", "VILL")],
}
# Injection base network (Experiments B-D): the same tree design on the quiet
# calibration days, whose repeatability is consistent with the formal covariance.
QUIET_TREES = {
    "2024001": [("MADR", "CEBR"), ("MADR", "VILL"), ("VILL", "YEBE")],
    "2024002": [("CEBR", "VILL"), ("MADR", "CEBR"), ("MADR", "YEBE")],
    "2024003": [("CEBR", "YEBE"), ("CEBR", "VILL"), ("MADR", "VILL")],
}
DATUM = "MADR"
ALL_PAIRS = [("MADR", "CEBR"), ("MADR", "VILL"), ("CEBR", "VILL"), ("VILL", "YEBE"),
             ("MADR", "YEBE"), ("CEBR", "YEBE")]


def main() -> None:
    proc = json.loads((ROOT / "results/real_network/processing/baselines_all.json").read_text())
    records = proc["records"]
    m = yaml.safe_load((ROOT / "configs/real_data/madrid_2024.yaml").read_text())
    snx_path = ROOT / proc["reference_sinex"]
    snx = read_sinex(snx_path, m["stations"])
    reference = {k: v.xyz for k, v in snx.stations.items()}
    ref_info = {"source": snx_path.name, "frame": "IGS20 (EPN combined weekly)",
                "epoch": snx.stations[DATUM].ref_epoch.isoformat()}
    # approximate (non-reference) coordinates for non-datum stations: RINEX headers
    approx = {}
    for r in records:
        approx.setdefault(r["to_station"], np.array(r["rover_meta"]["approx_xyz"]))
        approx.setdefault(r["from_station"], np.array(r["base_meta"]["approx_xyz"]))

    out = ROOT / "results" / "real_network"
    cal: dict = {"calibration_days": CAL_DAYS, "evaluation_days": list(EVAL_TREES)}
    by = lambda strat, days: [r for r in records if r["strategy"] == strat and r["session_id"] in days]
    for strat in ("iflc_float", "l12_fixed_short"):
        res = {}
        for label, days in (("calibration_days", CAL_DAYS), ("evaluation_days", list(EVAL_TREES)),
                            ("all_days", CAL_DAYS + list(EVAL_TREES))):
            r = repeatability_calibration(by(strat, days))
            res[label] = {"k2": r.k2, "k2_ci95": r.k2_ci95, "dof": r.dof, "n_pairs": r.n_pairs,
                          "n_sessions": r.n_sessions, "enu_variance_ratio": r.enu_variance_ratio,
                          "per_pair": r.per_pair}
        cal[strat] = res
    k2 = cal["iflc_float"]["calibration_days"]["k2"]
    (out / "calibration.json").write_text(json.dumps(cal, indent=1))

    iflc = [r for r in records if r["strategy"] == "iflc_float"]
    selection = [(s, a, b) for s, tree in EVAL_TREES.items() for a, b in tree]
    variants = {
        "main_formal": dict(selection=selection, covariance_scale=1.0, covariance_label="formal"),
        "main_calibrated": dict(selection=selection, covariance_scale=k2,
                                covariance_label="formal x k2 (days 001-003 repeatability)"),
        "quiet_formal": dict(selection=[(s, a, b) for s, tree in QUIET_TREES.items()
                                        for a, b in tree],
                             covariance_scale=1.0, covariance_label="formal"),
        "quiet_2sessions": dict(selection=[(s, a, b) for s, tree in QUIET_TREES.items()
                                           if s != "2024003" for a, b in tree],
                                covariance_scale=1.0, covariance_label="formal"),
        "trivial_day004": dict(selection=[("2024004", a, b) for a, b in ALL_PAIRS],
                               covariance_scale=1.0, covariance_label="formal",
                               allow_trivial=True),
    }
    summary = {}
    for name, kw in variants.items():
        design = build_network(iflc, datum_station=DATUM, approx_xyz=approx,
                               datum_xyz=reference[DATUM], **kw)
        write_network(design, out / "networks" / name, reference, ref_info)
        summary[name] = {
            "baselines": len(design.observations), "redundancy": design.redundancy,
            "loops": len(design.loops), "trivial_loops": sum(lp.trivial for lp in design.loops),
            "loop_chi2": [round(lp.chi2, 3) for lp in design.loops],
        }
    print(json.dumps({"k2_iflc_calibration": cal["iflc_float"]["calibration_days"]["k2"],
                      "k2_ci95": cal["iflc_float"]["calibration_days"]["k2_ci95"],
                      "k2_iflc_evaluation": cal["iflc_float"]["evaluation_days"]["k2"],
                      "enu_ratio_cal": cal["iflc_float"]["calibration_days"]["enu_variance_ratio"],
                      "k2_fixed_calibration": cal["l12_fixed_short"]["calibration_days"]["k2"],
                      "networks": summary}, indent=1))


if __name__ == "__main__":
    main()
