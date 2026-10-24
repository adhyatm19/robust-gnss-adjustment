"""Size-matched comparison of DIA/SD, Huber and Hampel on the real Madrid geometries.

Noise is resampled from the formal block-diagonal covariance around the WLS
solution of the real network (truth known by construction).

1. CALIBRATION (seed A, H0 only): each method gets one decision statistic
   and a threshold chosen so that it flags a baseline in exactly ALPHA of
   H0 runs:
     DIA/SD   : max_b T_b (first adjustment; T_b = |w*|^2) > c_DIA
     Huber    : min_b w_b (converged IRLS weights)        < w_Huber
     Hampel   : min_b w_b                                  < w_Hampel
   The flagged baseline is the argmax/argmin.
2. EVALUATION (seed B, disjoint): fresh H0 runs (achieved false-alarm rate)
   and single-outlier runs (random direction, 10-50 mm), scored as
   detection (any flag), correct identification, misidentification.
3. COORDINATE ACCURACY of the *standard* procedures (WLS, DIA-reject with
   its default alpha, Huber, Hampel) vs the simulated truth, reported
   separately from detection.

Output: results/real_network/audit/fair_comparison.{json,csv}
"""
from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from gnss_adjust import io as gio  # noqa: E402
from gnss_adjust.diagnostics import baseline_tests  # noqa: E402
from gnss_adjust.dia import DIAConfig, run_dia  # noqa: E402
from gnss_adjust.robust import RobustConfig, robust_adjust  # noqa: E402
from gnss_adjust.wls import solve_wls  # noqa: E402

ALPHA = 0.01
N_CAL, N_H0_TEST, N_TEST = 3000, 3000, 300
MAGS = [10, 15, 20, 30, 50]
OUT = ROOT / "results/real_network/audit"


def load(net):
    d = ROOT / "results/real_network/networks" / net
    return gio.load_network_from_csv(d / "stations.csv", d / "baselines.csv", d / "covariance_blocks.csv")


def simulate(rng, obs, x0, bias=None):
    out = []
    for o in obs:
        v = x0[o.to_station] - x0[o.from_station] + rng.multivariate_normal(np.zeros(3), o.covariance)
        if bias and o.baseline_id in bias:
            v = v + bias[o.baseline_id]
        out.append(replace(o, vector=v))
    return out


def evaluate(st, obs):
    w = solve_wls(st, obs)
    T = {t.baseline_id: t.group_statistic for t in baseline_tests(w)}
    hub = robust_adjust(st, obs, RobustConfig(method="huber"))
    ham = robust_adjust(st, obs, RobustConfig(method="hampel"))
    dia = run_dia(st, obs, DIAConfig())
    return {"T": T, "w_huber": hub.weights, "w_hampel": ham.weights,
            "coords": {"WLS": w.coordinates, "DIA-reject": dia.final.coordinates,
                       "Huber": hub.final.coordinates, "Hampel": ham.final.coordinates}}


def decisions(e, thr):
    b_t = max(e["T"], key=e["T"].get)
    out = {"DIA/SD": b_t if e["T"][b_t] > thr["DIA/SD"] else None}
    for name, key in (("Huber", "w_huber"), ("Hampel", "w_hampel")):
        b = min(e[key], key=e[key].get)
        out[name] = b if e[key][b] < thr[name] else None
    return out


def main():
    res, rows = {}, []
    for net in ("quiet_formal", "quiet_2sessions"):
        st, obs = load(net)
        x0 = solve_wls(st, obs).coordinates
        fixed = {s.station_id for s in st if s.is_fixed}
        ids = [o.baseline_id for o in obs]
        # ---- calibration (seed A)
        rng = np.random.default_rng(1001)
        cal = [evaluate(st, simulate(rng, obs, x0)) for _ in range(N_CAL)]
        thr = {"DIA/SD": float(np.quantile([max(e["T"].values()) for e in cal], 1 - ALPHA)),
               "Huber": float(np.quantile([min(e["w_huber"].values()) for e in cal], ALPHA)),
               "Hampel": float(np.quantile([min(e["w_hampel"].values()) for e in cal], ALPHA))}
        # ---- evaluation (seed B, disjoint)
        rng = np.random.default_rng(2002)
        h0 = [decisions(evaluate(st, simulate(rng, obs, x0)), thr) for _ in range(N_H0_TEST)]
        fa = {m: float(np.mean([d[m] is not None for d in h0])) for m in thr}
        res[net] = {"alpha_target": ALPHA, "thresholds": thr, "achieved_false_alarm_test": fa,
                    "n_calibration": N_CAL, "n_h0_test": N_H0_TEST}
        print(net, json.dumps(res[net]))
        for mag in MAGS:
            for _ in range(N_TEST):
                b = int(rng.choice(ids))
                u = rng.normal(size=3)
                u /= np.linalg.norm(u)
                e = evaluate(st, simulate(rng, obs, x0, {b: 1e-3 * mag * u}))
                dec = decisions(e, thr)
                row = {"net": net, "mag_mm": mag, "true": b}
                for m, f in dec.items():
                    row[f"{m}_detect"] = f is not None
                    row[f"{m}_correct"] = f == b
                    row[f"{m}_misid"] = f is not None and f != b
                for m, c in e["coords"].items():
                    d = np.array([c[s] - x0[s] for s in c if s not in fixed])
                    row[f"{m}_rmse_mm"] = 1e3 * float(np.sqrt(np.mean(d ** 2)))
                rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "fair_comparison_runs.csv", index=False)
    agg = df.groupby(["net", "mag_mm"]).mean(numeric_only=True).drop(columns="true").round(3)
    agg.to_csv(OUT / "fair_comparison.csv")
    res["summary"] = agg.reset_index().to_dict(orient="records")
    (OUT / "fair_comparison.json").write_text(json.dumps(res, indent=1))
    print(agg.to_string())


if __name__ == "__main__":
    main()
