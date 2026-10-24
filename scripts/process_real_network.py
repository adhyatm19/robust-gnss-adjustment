"""Process every station pair of a real dataset, day by day, with RTKLIB.

Reads a processing YAML (configs/real_data/*_processing.yaml) and its dataset
manifest; expects the data fetched by scripts/fetch_real_data.py. Each
(day, pair) is one independent static session solution (one line from
RTKLIB with out-solstatic=single); per-epoch filter outputs are never used
as separate observations because they are cumulative and strongly dependent.

Outputs (results/real_network/processing/):
  baselines_all.json   full records (vector, covariance, status, provenance)
  baselines_all.csv    one row per solution (ENU sigmas/correlations, status)
  metadata_qc.json     RINEX header/SINEX cross-checks per station and day
Intermediate files go to data/work/<dataset>/ (gitignored).
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gnss_adjust.gnss import frames  # noqa: E402
from gnss_adjust.gnss.baseline import BaselineError, process_baseline  # noqa: E402
from gnss_adjust.gnss.rinex import prepare_rinex, read_obs_meta  # noqa: E402
from gnss_adjust.gnss.rtklib import ProcessingConfig  # noqa: E402
from gnss_adjust.gnss.satscreen import rinex2_gps_unhealthy, sp3_bad_clock_satellites  # noqa: E402
from gnss_adjust.gnss.sinex import read_sinex  # noqa: E402


def day_inputs(m: dict, raw: Path, work: Path, doy: int) -> tuple[dict[str, Path], list[Path], list[str]]:
    year = m["year"]
    obs = {}
    for sta in m["stations"]:
        name = m["obs_template"].format(long=m["station_long_names"][sta], year=year, doy=doy)
        obs[sta] = prepare_rinex(raw / "obs" / name, work / "rinex")
    nav = prepare_rinex(raw / "nav" / m["nav_template"].format(doy=doy, yy=year % 100), work / "rinex")
    orbits = []
    excluded = set(rinex2_gps_unhealthy(nav))
    date = dt.date(year, 1, 1) + dt.timedelta(days=doy - 1)
    for y, d in m["orbit_days"]:
        if abs((y - year) * 365 + d - doy) <= 1:  # day and its neighbours
            f = prepare_rinex(raw / "orbit" / m["orbit_template"].format(year=y, doy=d), work / "orbit")
            orbits.append(f)
            excluded |= {s for s in sp3_bad_clock_satellites(f, (date.year, date.month, date.day))
                         if s.startswith("G")}
    return obs, [nav, *orbits], sorted(excluded)


def _one(job: dict) -> dict:
    try:
        sol = process_baseline(**job["kwargs"])
        return {"ok": True, "record": sol.to_record() | {"strategy": job["strategy"]}}
    except (BaselineError, RuntimeError) as exc:
        return {"ok": False, "job": job["label"], "error": str(exc)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("config", type=Path)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--days", type=int, nargs="*", default=None)
    ap.add_argument("--pairs", nargs="*", default=None, help="e.g. MADR-CEBR")
    ap.add_argument("--output", type=Path, default=ROOT / "results" / "real_network" / "processing")
    args = ap.parse_args()
    pc = yaml.safe_load(args.config.read_text())
    m = yaml.safe_load((ROOT / pc["dataset"]).read_text())
    raw = ROOT / "data" / "raw" / m["name"]
    work = ROOT / "data" / "work" / m["name"]
    snx_path = raw / "products" / Path(m["reference_sinex"]["url"]).name
    snx = read_sinex(snx_path, m["stations"])
    antex = raw / "products" / Path(m["antex"]["url"]).name
    strategies = {name: {**pc["processing"], **over} for name, over in pc["strategies"].items()}
    days = args.days or m["days"]
    pairs = [tuple(p) for p in pc["pairs"]]
    if args.pairs:
        pairs = [p for p in pairs if f"{p[0]}-{p[1]}" in args.pairs]

    qc: list[dict] = []
    screening: dict[int, list[str]] = {}
    jobs: list[dict] = []
    for doy in days:
        obs, aux, excluded = day_inputs(m, raw, work, doy)
        screening[doy] = excluded
        for sta, path in obs.items():
            meta = read_obs_meta(path)
            ref = snx.stations[sta]
            hen = meta.antenna_delta_hen
            qc.append({
                "doy": doy, "station": sta, "rinex_antenna": meta.antenna_type,
                "sinex_antenna": ref.antenna_type,
                "antenna_match": meta.antenna_type == ref.antenna_type,
                "rinex_delta_h_m": float(hen[0]),
                "sinex_ecc_up_m": float(ref.eccentricity_une[0]),
                "eccentricity_match": bool(abs(hen[0] - ref.eccentricity_une[0]) < 1e-4
                                           and np.allclose(hen[1:], ref.eccentricity_une[[2, 1]],
                                                           atol=1e-4)),
                "receiver": meta.receiver_type, "epochs": meta.n_epochs,
                "interval_s": meta.median_interval_s, "systems": meta.systems_in_data,
                "has_L1C_L2W": {"L1C", "L2W"} <= set(meta.obs_types.get("G", [])),
                "warnings": meta.warnings,
            })
        for strat, params in strategies.items():
            p = dict(params)
            only = p.pop("pairs", None)
            # automatic screening (default) or a fixed list for sensitivity analysis
            fixed_excl = p.pop("fixed_exclusions", None)
            excl = excluded if fixed_excl is None else sorted(fixed_excl)
            p["extra"] = {**p.get("extra", {}), "pos1-exclsats": " ".join(excl)}
            cfg = ProcessingConfig(**{**p, "nav_systems": tuple(p["nav_systems"]),
                                      "antex_file": str(antex)})
            for base, rover in pairs:
                if only is not None and [base, rover] not in only:
                    continue
                label = f"{strat} {m['year']}-{doy:03d} {base}->{rover}"
                jobs.append({"label": label, "strategy": strat, "kwargs": dict(
                    session_id=f"{m['year']}{doy:03d}", rover_obs=obs[rover], base_obs=obs[base],
                    aux_files=aux, base_xyz=snx.stations[base].xyz,
                    workdir=work / strat / f"d{doy:03d}", cfg=cfg, rover_station=rover,
                    base_station=base, variant=pc["rtklib_variant"], selection="last")})
    bad_qc = [q for q in qc if not (q["antenna_match"] and q["eccentricity_match"] and q["has_L1C_L2W"])]
    if bad_qc:
        raise SystemExit(f"metadata QC failed: {bad_qc}")

    results = []
    with cf.ProcessPoolExecutor(max_workers=args.workers) as ex:
        for res in ex.map(_one, jobs):
            results.append(res)
            tag = (f"{res['record']['strategy']} {res['record']['session_id']} "
                   f"{res['record']['from_station']}->{res['record']['to_station']}"
                   if res["ok"] else res["job"])
            print(("ok   " if res["ok"] else "FAIL ") + tag +
                  (f"  {res['record']['status']} fix%={100 * res['record']['fix_fraction']:.0f}"
                   f" ratio={res['record']['ratio']:.1f}"
                   f" t={res['record']['run']['elapsed_s']:.0f}s" if res["ok"] else f"  {res['error']}"),
                  flush=True)

    args.output.mkdir(parents=True, exist_ok=True)
    records = [r["record"] for r in results if r["ok"]]
    failures = [r for r in results if not r["ok"]]
    # repository-relative paths so the log is machine-independent
    text = json.dumps({"config": pc, "reference_sinex": str(snx_path.relative_to(ROOT)),
                       "satellite_screening": screening, "records": records,
                       "failures": failures}, indent=1, default=str)
    (args.output / "baselines_all.json").write_text(text.replace(str(ROOT) + "/", ""))
    rows = []
    for r in records:
        qe = np.array(r["covariance_enu"])
        s = np.sqrt(np.diag(qe))
        ref_b = snx.stations[r["to_station"]].xyz - snx.stations[r["from_station"]].xyz
        d_enu = np.array(r["vector_enu"]) - frames.vector_ecef_to_enu(
            ref_b, snx.stations[r["from_station"]].xyz)
        rows.append({
            "strategy": r["strategy"],
            "session": r["session_id"], "from": r["from_station"], "to": r["to_station"],
            "length_km": r["length_m"] / 1e3, "status": r["status"],
            "fix_fraction": r["fix_fraction"], "epochs_solved": r["epochs_solved"],
            "ratio": r["ratio"],
            "diff_ref_e_mm": 1e3 * d_enu[0], "diff_ref_n_mm": 1e3 * d_enu[1],
            "diff_ref_u_mm": 1e3 * d_enu[2],
            "n_sat": r["n_sat"], "sigma_e_mm": 1e3 * s[0], "sigma_n_mm": 1e3 * s[1],
            "sigma_u_mm": 1e3 * s[2], "corr_en": qe[0, 1] / s[0] / s[1],
            "corr_nu": qe[1, 2] / s[1] / s[2], "corr_ue": qe[2, 0] / s[2] / s[0],
            "covariance_source": r["covariance_source"], "elapsed_s": r["run"]["elapsed_s"],
            "warnings": "; ".join(r["warnings"]),
        })
    pd.DataFrame(rows).sort_values(["strategy", "session", "from", "to"]).to_csv(
        args.output / "baselines_all.csv", index=False)
    (args.output / "metadata_qc.json").write_text(
        json.dumps(qc, indent=1, default=str).replace(str(ROOT) + "/", ""))
    print(f"{len(records)} solutions, {len(failures)} failures -> {args.output}")


if __name__ == "__main__":
    main()
