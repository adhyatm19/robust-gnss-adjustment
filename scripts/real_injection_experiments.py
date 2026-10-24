"""Experiments B and C: controlled outlier injection into real-data-derived baselines.

IMPORTANT: these are *controlled injections* of known bias vectors into
baselines estimated by RTKLIB from real RINEX observations. They are not
naturally occurring errors; detection metrics are computed only against the
injected labels. The base network (quiet days 001-003, formal covariance)
passes the global test unmodified, and contains one fixed realisation of
real noise, so variability across trials comes from the injected baseline,
direction and magnitude -- not from re-sampled noise.

Experiment B: detection/localisation/coordinate impact for WLS, DIA-reject,
DIA-inflate, Huber, Hampel and SD-reject, for 1 or 2 contaminated baselines,
random / vertical / horizontal directions, 5-80 mm magnitudes, on two
geometries (3 sessions, r=18; 2 sessions, r=9).

Experiment C: specific-direction validation -- angular error of the SD
bias direction versus magnitude and non-centrality; the predicted direction
uncertainty; the numerical SD = 3D = likelihood-ratio identity on the real
network; and the detection rate of the optimised-direction (SD/3D) test
versus a *prespecified*-direction 1D test (correct up-direction, wrong
direction), against the theoretical non-central chi-square power.

Outputs: results/real_network/experiment_BC/
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gnss_adjust import benchmark as bm  # noqa: E402
from gnss_adjust import io as gio  # noqa: E402
from gnss_adjust import metrics as met  # noqa: E402
from gnss_adjust.gnss import frames  # noqa: E402
from gnss_adjust.network_report import read_reference_csv  # noqa: E402
from gnss_adjust.sd_outlier import (SDConfig, angular_error_degrees,  # noqa: E402
                                    detect_specific_direction, run_sd_snooping)
from gnss_adjust.wls import solve_wls  # noqa: E402

METHODS = [*bm.ALL_METHODS, "SD-reject"]
ALPHA_LOCAL = 0.001


def load(net: str):
    d = ROOT / "results" / "real_network" / "networks" / net
    st, obs = gio.load_network_from_csv(d / "stations.csv", d / "baselines.csv",
                                        d / "covariance_blocks.csv")
    return st, obs, read_reference_csv(d / "reference_coordinates.csv")


def unit_direction(mode: str, rng: np.random.Generator, ref_xyz: np.ndarray) -> np.ndarray:
    r = frames.enu_rotation_at(ref_xyz)
    if mode == "random":
        v = rng.normal(size=3)
        return v / np.linalg.norm(v)
    if mode == "up":
        return r[2] * rng.choice([-1.0, 1.0])
    if mode == "horizontal":
        az = rng.uniform(0, 2 * np.pi)
        return np.cos(az) * r[0] + np.sin(az) * r[1]
    raise ValueError(mode)


def run_methods(st, obs, settings, sd_cfg):
    out = {}
    t0 = time.perf_counter()
    oc = bm.run_all_methods(st, obs, settings)
    for k, v in oc.items():
        out[k] = v
    t1 = time.perf_counter()
    sd = run_sd_snooping(st, obs, sd_cfg)
    out["SD-reject"] = bm.MethodOutcome(
        method="SD-reject", coordinates=sd.final.coordinates, flagged=set(sd.flagged_baselines),
        weights=None, iterations=len(sd.iterations), converged=sd.status in ("clean", "adapted"),
        runtime_s=time.perf_counter() - t1, condition_number=sd.final.diagnostics.condition_number_N,
        variance_factor=sd.final.variance_factor, detail=sd)
    return out


def experiment(geometry: str, net: str, n_trials: int, seed: int, magnitudes_mm, modes, n_out_list):
    st, obs, ref = load(net)
    settings = bm.MethodSettings()
    sd_cfg = SDConfig(alpha_local=ALPHA_LOCAL)
    smap = {s.station_id: s for s in st}
    fixed = {s.station_id for s in st if s.is_fixed}
    datum = next(iter(fixed))
    clean = solve_wls(st, obs)
    clean_coords = clean.coordinates
    clean_sd = detect_specific_direction(clean, alpha=ALPHA_LOCAL)
    pbar = {t.baseline_id: t.reliability_block for t in clean_sd.tests}   # geometry only
    ids = [o.baseline_id for o in obs]
    rng = np.random.default_rng(seed)
    crit3 = stats.chi2.ppf(1 - ALPHA_LOCAL, 3)
    z1 = stats.norm.ppf(1 - ALPHA_LOCAL / 2)
    runs, dirs = [], []
    run_id = 0
    for n_out in n_out_list:
        for mode in modes:
            for mag in magnitudes_mm:
                for _ in range(n_trials):
                    run_id += 1
                    chosen = sorted(rng.choice(ids, size=n_out, replace=False).tolist())
                    inj = {}
                    for b in chosen:
                        o = next(x for x in obs if x.baseline_id == b)
                        inj[b] = 1e-3 * mag * unit_direction(mode, rng, smap[o.from_station].coords)
                    cont = [replace(o, vector=o.vector + inj[o.baseline_id]) if o.baseline_id in inj
                            else o for o in obs]
                    outc = run_methods(st, cont, settings, sd_cfg)
                    for name, oc in outc.items():
                        dm = met.detection_metrics(set(chosen), oc.flagged, set(ids))
                        d_clean = np.array([oc.coordinates[s] - clean_coords[s]
                                            for s in oc.coordinates if s not in fixed])
                        d_ref = np.array([(oc.coordinates[s] - oc.coordinates[datum])
                                          - (ref[s] - ref[datum]) for s in oc.coordinates
                                          if s not in fixed])
                        runs.append({
                            "geometry": geometry, "run_id": run_id, "n_outliers": n_out,
                            "direction_mode": mode, "magnitude_mm": mag, "method": name,
                            "true_outliers": " ".join(map(str, chosen)),
                            "flagged": " ".join(map(str, sorted(oc.flagged))),
                            "tp": dm.true_positives, "fp": dm.false_positives,
                            "fn": dm.false_negatives, "precision": dm.precision,
                            "recall": dm.recall, "f1": dm.f1, "fpr": dm.false_positive_rate,
                            "localization_correct": dm.localization_correct,
                            "rmse_vs_clean_mm": 1e3 * float(np.sqrt(np.mean(d_clean ** 2))),
                            "rmse_vs_epn_mm": 1e3 * float(np.sqrt(np.mean(d_ref ** 2))),
                            "iterations": oc.iterations, "runtime_s": oc.runtime_s,
                            "converged": oc.converged,
                            "min_weight_true": (min(oc.weights[b] for b in chosen)
                                                if oc.weights else np.nan),
                            "min_weight_clean": (min(w for b, w in oc.weights.items()
                                                     if b not in chosen) if oc.weights else np.nan),
                        })
                    # ---- Experiment C: local SD test on the contaminated WLS adjustment
                    wls = outc["WLS"].detail
                    sd = detect_specific_direction(wls, alpha=ALPHA_LOCAL)
                    g_all = -(wls.P @ wls.residuals)
                    for b in chosen:
                        k = ids.index(b)
                        t = sd.tests[k]
                        d = inj[b]
                        u_true = d / np.linalg.norm(d)
                        o = obs[k]
                        up = frames.enu_rotation_at(smap[o.from_station].coords)[2]
                        P = pbar[b]
                        lam = float(d @ P @ d)
                        g = g_all[3 * k:3 * k + 3]
                        w_true = float(u_true @ g / np.sqrt(u_true @ P @ u_true))
                        w_up = float(up @ g / np.sqrt(up @ P @ up))
                        cov_d = np.linalg.inv(t.reliability_block)
                        dh = t.outlier_vector
                        proj = np.eye(3) - np.outer(dh, dh) / (dh @ dh)
                        sig_ang = float(np.degrees(np.arctan2(
                            np.sqrt(np.linalg.eigvalsh(proj @ cov_d @ proj)[-1]), np.linalg.norm(dh))))
                        dirs.append({
                            "geometry": geometry, "run_id": run_id, "n_outliers": n_out,
                            "direction_mode": mode, "magnitude_mm": mag, "baseline_id": b,
                            "noncentrality": lam,
                            "power_3d_theory": float(stats.ncx2.sf(crit3, 3, lam)),
                            "power_1d_theory": float(stats.ncx2.sf(z1 ** 2, 1, lam)),
                            "sd_statistic": t.sd_statistic, "sd_significant": t.significant,
                            "identified_as_max": sd.identified_baseline == b,
                            "angle_deg": angular_error_degrees(t.direction, u_true),
                            "angle_sigma_pred_deg": sig_ang,
                            "estimated_norm_mm": 1e3 * t.outlier_magnitude,
                            "w_true_direction": w_true, "detect_1d_true": abs(w_true) > z1,
                            "w_up_direction": w_up, "detect_1d_up": abs(w_up) > z1,
                            "true_elevation_deg": float(np.degrees(np.arcsin(abs(up @ u_true)))),
                        })
    return pd.DataFrame(runs), pd.DataFrame(dirs)


def identity_check(net: str) -> dict:
    """SD |w*|^2 vs DIA T_b vs likelihood-ratio (deletion) on the real network."""
    from gnss_adjust.diagnostics import baseline_tests
    st, obs, _ = load(net)
    full = solve_wls(st, obs)
    sd = detect_specific_direction(full)
    dia = baseline_tests(full)
    rows = []
    for k, o in enumerate(obs):
        rest = [x for x in obs if x.baseline_id != o.baseline_id]
        lr = full.vTPv - solve_wls(st, rest).vTPv   # block-diagonal: augmentation == deletion
        rows.append({"baseline_id": o.baseline_id, "sd_sq": sd.tests[k].sd_statistic ** 2,
                     "dia_T": dia[k].group_statistic, "lr": lr})
    df = pd.DataFrame(rows)
    return {"network": net, "rows": rows,
            "max_rel_sd_vs_lr": float((abs(df.sd_sq - df.lr) / df.lr).max()),
            "max_rel_dia_vs_lr": float((abs(df.dia_T - df.lr) / df.lr).max())}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--trials", type=int, default=30)
    ap.add_argument("--seed", type=int, default=20240101)
    ap.add_argument("--output", type=Path, default=ROOT / "results/real_network/experiment_BC")
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    mags = [5, 10, 15, 20, 30, 50, 80]
    modes = ["random", "up", "horizontal"]
    runs, dirs = [], []
    for geometry, net, seed_off in (("3 sessions (r=18)", "quiet_formal", 0),
                                    ("2 sessions (r=9)", "quiet_2sessions", 1)):
        r, d = experiment(geometry, net, args.trials, args.seed + seed_off, mags, modes, [1, 2])
        runs.append(r)
        dirs.append(d)
    runs = pd.concat(runs)
    dirs = pd.concat(dirs)
    runs.to_csv(args.output / "runs.csv", index=False)
    dirs.to_csv(args.output / "sd_directions.csv", index=False)
    agg = runs.groupby(["geometry", "n_outliers", "direction_mode", "magnitude_mm", "method"]).agg(
        f1=("f1", "mean"), precision=("precision", "mean"), recall=("recall", "mean"),
        fpr=("fpr", "mean"), localization=("localization_correct", "mean"),
        rmse_vs_clean_mm=("rmse_vs_clean_mm", "mean"), rmse_vs_epn_mm=("rmse_vs_epn_mm", "mean"),
        iterations=("iterations", "mean"), runtime_ms=("runtime_s", lambda x: 1e3 * x.mean()),
        runs=("run_id", "nunique")).reset_index()
    agg.to_csv(args.output / "aggregate.csv", index=False)
    keys = ["geometry", "n_outliers", "direction_mode", "magnitude_mm"]
    dagg = dirs.groupby(keys).agg(
        n=("angle_deg", "size"), noncentrality_median=("noncentrality", "median"),
        median_angle_deg=("angle_deg", "median"),
        p90_angle_deg=("angle_deg", lambda x: np.nanpercentile(x, 90)),
        median_angle_over_pred=("angle_deg", "median"),
        sd_detect_rate=("sd_significant", "mean"), power_3d_theory=("power_3d_theory", "mean"),
        detect_1d_true_rate=("detect_1d_true", "mean"), power_1d_theory=("power_1d_theory", "mean"),
        detect_1d_up_rate=("detect_1d_up", "mean")).reset_index()
    dagg["median_angle_over_pred"] = dirs.assign(r=dirs.angle_deg / dirs.angle_sigma_pred_deg) \
        .groupby(keys).r.median().values
    dagg.to_csv(args.output / "sd_direction_summary.csv", index=False)
    ident = [identity_check("quiet_formal"), identity_check("main_formal")]
    (args.output / "identity_check.json").write_text(json.dumps(ident, indent=1))
    (args.output / "config.json").write_text(json.dumps(
        {"trials_per_cell": args.trials, "seed": args.seed, "magnitudes_mm": mags,
         "direction_modes": modes, "n_outliers": [1, 2], "alpha_local": ALPHA_LOCAL,
         "methods": METHODS, "note": "controlled injections into real-data-derived baselines"},
        indent=1))
    print(agg[(agg.n_outliers == 1) & (agg.direction_mode == "random")]
          .pivot_table(index=["geometry", "magnitude_mm"], columns="method", values="f1").round(2))
    print(dagg[(dagg.direction_mode == "random") & (dagg.n_outliers == 1)].round(3).to_string())
    print(json.dumps({k: v for k, v in ident[0].items() if k != "rows"}),
          json.dumps({k: v for k, v in ident[1].items() if k != "rows"}))


if __name__ == "__main__":
    main()
