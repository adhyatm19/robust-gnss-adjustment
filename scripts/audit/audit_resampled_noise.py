"""Audit simulations on the REAL Madrid geometries with RESAMPLED noise.

Truth = WLS solution of the unmodified network; observations regenerated as
truth + e, e ~ N(0, Q_formal) (block-diagonal, exactly the assumed model).
This isolates the statistical procedures from the single real noise
realisation used in Experiment B.
"""
import json, sys
from dataclasses import replace
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from gnss_adjust import benchmark as bm, io as gio
from gnss_adjust.wls import solve_wls
from gnss_adjust.diagnostics import baseline_tests, global_model_test
from gnss_adjust.sd_outlier import detect_specific_direction, run_sd_snooping, SDConfig
from gnss_adjust.gnss import frames

OUT = ROOT / "results/real_network/audit"; OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(777)
ALPHA = 0.001

def load(net):
    d = ROOT / "results/real_network/networks" / net
    return gio.load_network_from_csv(d / "stations.csv", d / "baselines.csv", d / "covariance_blocks.csv")

def regenerate(st, obs, x0, extra=None):
    out = []
    for o in obs:
        e = rng.multivariate_normal(np.zeros(3), o.covariance)
        v = x0[o.to_station] - x0[o.from_station] + e
        if extra and o.baseline_id in extra:
            v = v + extra[o.baseline_id]
        out.append(replace(o, vector=v))
    return out

def methods(st, obs):
    oc = bm.run_all_methods(st, obs, bm.MethodSettings())
    sd = run_sd_snooping(st, obs, SDConfig(alpha_local=ALPHA))
    fl = {k: set(v.flagged) for k, v in oc.items()}
    fl["SD-reject"] = set(sd.flagged_baselines)
    return oc["WLS"].detail, fl

res = {}
# ---------------- S1: size under H0
for net, n_mc in (("quiet_formal", 2000), ("quiet_2sessions", 2000), ("main_formal", 2000)):
    st, obs = load(net)
    x0 = solve_wls(st, obs).coordinates
    glob = 0; anyflag = {m: 0 for m in [*bm.ALL_METHODS, "SD-reject"]}; local_exc = np.zeros(len(obs))
    for _ in range(n_mc):
        o2 = regenerate(st, obs, x0)
        w, fl = methods(st, o2)
        glob += not global_model_test(w).passed
        for m in anyflag: anyflag[m] += bool(fl[m])
        local_exc += np.array([t.flagged_group for t in baseline_tests(w)])
    res[f"H0_{net}"] = {"n_mc": n_mc, "global_reject_rate": glob / n_mc,
                        "any_flag_rate": {m: v / n_mc for m, v in anyflag.items()},
                        "local_exceedance_rate_per_baseline": (local_exc / n_mc).round(4).tolist(),
                        "expected_local": ALPHA, "expected_global": 0.05}
    print(net, json.dumps(res[f"H0_{net}"]))

# ---------------- S2: power & misidentification with resampled noise (single outlier, random dir)
rows = []
crit3 = stats.chi2.ppf(1 - ALPHA, 3)
for net in ("quiet_formal", "quiet_2sessions"):
    st, obs = load(net)
    clean = solve_wls(st, obs); x0 = clean.coordinates
    pbar = {t.baseline_id: t.reliability_block for t in detect_specific_direction(clean).tests}
    smap = {s.station_id: s for s in st}
    ids = [o.baseline_id for o in obs]
    for mag in (10, 15, 20, 30, 50):
        for _ in range(300):
            b = int(rng.choice(ids)); u = rng.normal(size=3); u /= np.linalg.norm(u)
            d = 1e-3 * mag * u
            o2 = regenerate(st, obs, x0, {b: d})
            w, fl = methods(st, o2)
            sd = detect_specific_direction(w, alpha=ALPHA)
            k = ids.index(b)
            lam = float(d @ pbar[b] @ d)
            rows.append({"net": net, "mag": mag, "b": b, "lam": lam,
                         "power_theory": float(stats.ncx2.sf(crit3, 3, lam)),
                         "local_detect": sd.tests[k].significant,
                         "dia_correct": fl["DIA-reject"] == {b},
                         "dia_wrong": bool(fl["DIA-reject"] - {b}),
                         "dia_any": bool(fl["DIA-reject"]),
                         "sd_eq_dia": fl["SD-reject"] == fl["DIA-reject"],
                         "huber_correct": fl["Huber"] == {b}, "huber_wrong": bool(fl["Huber"] - {b}),
                         "angle": float(np.degrees(np.arccos(min(1, abs(sd.tests[k].direction @ u))))),
                         })
df = pd.DataFrame(rows); df.to_csv(OUT / "power_resampled.csv", index=False)
g = df.groupby(["net", "mag"]).agg(n=("lam", "size"), power_theory=("power_theory", "mean"),
    local_detect=("local_detect", "mean"), dia_correct=("dia_correct", "mean"),
    dia_wrong=("dia_wrong", "mean"), huber_correct=("huber_correct", "mean"),
    sd_eq_dia=("sd_eq_dia", "mean"), median_angle=("angle", "median")).round(3)
print(g.to_string())
res["power_resampled"] = g.reset_index().to_dict(orient="records")
# binned by lambda (theory vs empirical, both nets)
bins = [0, 5, 10, 15, 20, 30, 50, 100, 1e9]
df["lbin"] = pd.cut(df.lam, bins)
gb = df.groupby("lbin", observed=True).agg(n=("lam", "size"), theory=("power_theory", "mean"),
                                           empirical=("local_detect", "mean")).round(3)
gb["binom_se"] = np.sqrt(gb.theory * (1 - gb.theory) / gb.n).round(3)
print(gb.to_string()); res["power_by_lambda"] = gb.reset_index().astype(str).to_dict(orient="records")

# ---------------- S3: confusability (canonical correlations of local test vectors)
for net in ("quiet_formal", "quiet_2sessions", "main_formal"):
    st, obs = load(net)
    w = solve_wls(st, obs)
    Pb = w.P @ w.residual_covariance @ w.P
    m = len(obs); worst = []
    for i in range(m):
        for j in range(i + 1, m):
            Pii = Pb[3*i:3*i+3, 3*i:3*i+3]; Pjj = Pb[3*j:3*j+3, 3*j:3*j+3]; Pij = Pb[3*i:3*i+3, 3*j:3*j+3]
            Li = np.linalg.cholesky(np.linalg.inv(Pii)); Lj = np.linalg.cholesky(np.linalg.inv(Pjj))
            s = np.linalg.svd(Li.T @ Pij @ Lj, compute_uv=False)
            worst.append((round(float(s.max()), 4), obs[i].baseline_id, obs[j].baseline_id))
    worst.sort(reverse=True)
    res[f"canonical_corr_{net}"] = worst[:5]
    print(net, "max canonical correlations between local tests:", worst[:5])

# ---------------- S4: where does the quiet realisation sit?
for net in ("quiet_formal", "main_formal"):
    st, obs = load(net); w = solve_wls(st, obs)
    p = float(stats.chi2.cdf(w.vTPv, w.degrees_of_freedom))
    res[f"realisation_{net}"] = {"vTPv": w.vTPv, "dof": w.degrees_of_freedom, "chi2_cdf": p}
    print(net, "v'Pv", round(w.vTPv, 2), "cdf", round(p, 4))
(OUT / "audit_stats.json").write_text(json.dumps(res, indent=1, default=str))
