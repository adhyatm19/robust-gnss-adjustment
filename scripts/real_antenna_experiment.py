"""Experiment D: physically meaningful antenna-height errors, by *reprocessing* real RINEX.

A wrong antenna height is entered for one station in the RTKLIB
configuration (explicit antenna type + perturbed ``antdelu``; the RINEX
files are never modified), the affected baselines are reprocessed, and the
perturbed solutions replace the originals in the quiet-days network.

Cases (quiet network: days 001-003, one spanning tree per day):
  D1  YEBE, day 002: YEBE occurs in one baseline that day (MADR->YEBE)
      -> a single-baseline physical outlier along YEBE's local vertical;
  D2  VILL, day 001: VILL occurs in two baselines (MADR->VILL as rover,
      VILL->YEBE as base) -> a station-level error violating the
      single-baseline outlier model.

Checks: (i) RTKLIB applied the perturbation (baseline change vs the
predicted -dh*u_rover / +dh*u_base), (ii) what each method flags, (iii) SD
direction vs local vertical, (iv) prespecified vertical 1D test vs the
optimised 3D/SD test, (v) the induced coordinate error.
Outputs: results/real_network/experiment_D/
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import build_real_network as brn  # noqa: E402
import process_real_network as prn  # noqa: E402
from gnss_adjust import benchmark as bm  # noqa: E402
from gnss_adjust.gnss import frames  # noqa: E402
from gnss_adjust.gnss.baseline import antenna_from_header, process_baseline  # noqa: E402
from gnss_adjust.gnss.network_build import build_network  # noqa: E402
from gnss_adjust.gnss.rinex import read_obs_meta  # noqa: E402
from gnss_adjust.gnss.rtklib import ProcessingConfig  # noqa: E402
from gnss_adjust.gnss.sinex import read_sinex  # noqa: E402
from gnss_adjust.sd_outlier import SDConfig, detect_specific_direction, run_sd_snooping  # noqa: E402
from gnss_adjust.wls import solve_wls  # noqa: E402

CASES = {"D1": ("YEBE", 2), "D2": ("VILL", 1)}
DH_M = [0.01, 0.02, 0.03, 0.05, 0.10]
ALPHA = 0.001


def main() -> None:
    pc = yaml.safe_load((ROOT / "configs/real_data/madrid_2024_processing.yaml").read_text())
    m = yaml.safe_load((ROOT / pc["dataset"]).read_text())
    raw, work = ROOT / "data/raw" / m["name"], ROOT / "data/work" / m["name"]
    proc = json.loads((ROOT / "results/real_network/processing/baselines_all.json").read_text())
    records = [r for r in proc["records"] if r["strategy"] == "iflc_float"]
    snx = read_sinex(ROOT / proc["reference_sinex"], m["stations"])
    ref = {k: v.xyz for k, v in snx.stations.items()}
    antex = raw / "products" / Path(m["antex"]["url"]).name
    params = {**pc["processing"], **pc["strategies"]["iflc_float"]}
    selection = [(s, a, b) for s, tree in brn.QUIET_TREES.items() for a, b in tree]
    approx = {r["to_station"]: np.array(r["rover_meta"]["approx_xyz"]) for r in records}
    approx |= {r["from_station"]: np.array(r["base_meta"]["approx_xyz"]) for r in records}
    settings = bm.MethodSettings()
    clean_net = build_network(records, selection, brn.DATUM, approx, ref[brn.DATUM])
    clean = solve_wls(clean_net.stations, clean_net.observations)
    z1 = stats.norm.ppf(1 - ALPHA / 2)

    rows, verify = [], []
    for case, (station, doy) in CASES.items():
        obs, aux, excluded = prn.day_inputs(m, raw, work, doy)
        session = f"{m['year']}{doy:03d}"
        cfg = ProcessingConfig(**{**params, "nav_systems": tuple(params["nav_systems"]),
                                  "antex_file": str(antex),
                                  "extra": {**params["extra"], "pos1-exclsats": " ".join(excluded)}})
        affected = [(s, a, b) for (s, a, b) in selection if s == session and station in (a, b)]
        meta = read_obs_meta(obs[station])
        up = frames.enu_rotation_at(ref[station])[2]
        for dh in DH_M:
            new = {}
            for (_, a, b) in affected:
                kw = {"rover_antenna": antenna_from_header(meta, dh)} if b == station else \
                     {"base_antenna": antenna_from_header(meta, dh)}
                sol = process_baseline(
                    session_id=session, rover_obs=obs[b], base_obs=obs[a], aux_files=aux,
                    base_xyz=ref[a], workdir=work / "antenna" / f"{case}_{dh:.2f}", cfg=cfg,
                    rover_station=b, base_station=a, variant=pc["rtklib_variant"],
                    selection="last", **kw)
                rec = sol.to_record() | {"strategy": "iflc_float"}
                orig = next(r for r in records if (r["session_id"], r["from_station"],
                                                    r["to_station"]) == (session, a, b))
                delta = np.array(rec["vector_ecef"]) - np.array(orig["vector_ecef"])
                predicted = (-dh if b == station else dh) * up
                verify.append({"case": case, "dh_m": dh, "baseline": f"{a}->{b}",
                               "station_role": "rover" if b == station else "base",
                               "applied_antdelu_m": (rec["rover_antenna"] if b == station
                                                     else rec["base_antenna"])["delta_enu"][2],
                               "delta_mm": (1e3 * delta).tolist(),
                               "predicted_mm": (1e3 * predicted).tolist(),
                               "max_abs_deviation_mm": float(1e3 * np.max(np.abs(delta - predicted)))})
                new[(session, a, b)] = rec
            recs = [new.get((r["session_id"], r["from_station"], r["to_station"]), r) for r in records]
            net = build_network(recs, selection, brn.DATUM, approx, ref[brn.DATUM])
            ids_aff = {i + 1 for i, key in enumerate(selection) if key in new}
            outc = bm.run_all_methods(net.stations, net.observations, settings)
            sdr = run_sd_snooping(net.stations, net.observations, SDConfig(alpha_local=ALPHA))
            wls = outc["WLS"].detail
            sd = detect_specific_direction(wls, alpha=ALPHA)
            g = -(wls.P @ wls.residuals)
            flags = {k: sorted(v.flagged) for k, v in outc.items()} | {"SD-reject": sdr.flagged_baselines}
            for bid in sorted(ids_aff):
                k = bid - 1
                t = sd.tests[k]
                P = t.reliability_block
                gk = g[3 * k:3 * k + 3]
                w_up = float(up @ gk / np.sqrt(up @ P @ up))
                elev = float(np.degrees(np.arcsin(min(1.0, abs(t.direction @ up)))))
                rows.append({
                    "case": case, "station": station, "session": session, "dh_m": dh,
                    "baseline_id": bid, "baseline": f"{selection[k][1]}->{selection[k][2]}",
                    "affected_ids": " ".join(map(str, sorted(ids_aff))),
                    "global_T": wls.vTPv, "global_crit": float(stats.chi2.ppf(0.95, wls.degrees_of_freedom)),
                    "sd_statistic": t.sd_statistic, "sd_critical": t.critical_sd,
                    "sd_significant": t.significant, "sd_identified_max": sd.identified_baseline == bid,
                    "sd_angle_to_vertical_deg": 90.0 - elev,
                    "sd_bias_norm_mm": 1e3 * t.outlier_magnitude,
                    "w_up_1d": w_up, "w_up_critical": z1, "up_test_significant": abs(w_up) > z1,
                    **{f"flags_{name}": " ".join(map(str, f)) for name, f in flags.items()},
                    "station_coord_error_wls_mm": 1e3 * float(np.linalg.norm(
                        wls.coordinates[station] - clean.coordinates[station])),
                    "station_coord_error_dia_mm": 1e3 * float(np.linalg.norm(
                        outc["DIA-reject"].coordinates[station] - clean.coordinates[station])),
                    "station_up_error_wls_mm": 1e3 * float(up @ (wls.coordinates[station]
                                                                 - clean.coordinates[station])),
                })
    out = ROOT / "results/real_network/experiment_D"
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / "summary.csv", index=False)
    (out / "verification.json").write_text(json.dumps(verify, indent=1))
    print(pd.DataFrame(verify)[["case", "dh_m", "baseline", "station_role", "applied_antdelu_m",
                                "max_abs_deviation_mm"]].to_string(index=False))
    cols = ["case", "dh_m", "baseline", "sd_statistic", "sd_significant", "sd_angle_to_vertical_deg",
            "w_up_1d", "up_test_significant", "flags_DIA-reject", "flags_Huber", "flags_Hampel",
            "station_up_error_wls_mm"]
    print(pd.DataFrame(rows)[cols].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
