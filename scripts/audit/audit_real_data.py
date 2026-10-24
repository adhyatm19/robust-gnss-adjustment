"""Scientific-audit checks on the real-data results (docs/SCIENTIFIC_AUDIT.md).

Prerequisites: scripts/reproduce_real_data.sh has run, and the processing
sensitivity analysis:
  python scripts/process_real_network.py configs/real_data/madrid_2024_sensitivity.yaml \
      --output results/real_network/audit_sensitivity
Downloads (logged with SHA-256): ten EPN analysis-centre weekly SINEX files
for GPS week 2295 and epn_20.atx (~62 MB) into data/raw/audit_epn_ac/.

Writes results/real_network/audit/:
  calibration_decomposition.json   k^2 by period/component/station/day, F-test, lag-1 corr
  sensitivity_summary.csv          processing-model sensitivity (VILL, CEBR, k^2)
  epn_ac_comparison.csv            CEBR/VILL/YEBE - MADR per EPN analysis centre
  antex_comparison.json            SEPCHOKE_B3E6 entries, IGS20 vs EPN ATX
  cebr_l2_signal_test.csv          CEBR position with L2W vs L2L
  injection_on_main_network.csv    Experiment-B protocol on the anomaly-containing network
  audit_download_log.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import fetch_real_data as frd  # noqa: E402
import process_real_network as prn  # noqa: E402
import real_injection_experiments as rie  # noqa: E402
from gnss_adjust.gnss import frames  # noqa: E402
from gnss_adjust.gnss.baseline import process_baseline  # noqa: E402
from gnss_adjust.gnss.calibration import repeatability_calibration  # noqa: E402
from gnss_adjust.gnss.rtklib import ProcessingConfig  # noqa: E402
from gnss_adjust.gnss.sinex import read_sinex  # noqa: E402

OUT = ROOT / "results/real_network/audit"
CAL = ["2024001", "2024002", "2024003"]
EVAL = ["2024004", "2024005", "2024006"]
STATIONS = ["MADR", "CEBR", "VILL", "YEBE"]
ACS = ["BEK0", "BEV0", "BKG0", "COD0", "GFZ0", "IGE0", "LPT0", "NKG0", "RGA0", "WUT0"]
BKG = "https://igs.bkg.bund.de/root_ftp"


def k2_component(recs, comp=None):
    """Pooled repeatability factor, 3D (comp=None) or one ENU component (0,1,2)."""
    if comp is None:
        return repeatability_calibration(recs).k2
    g: dict = {}
    for r in recs:
        g.setdefault((r["from_station"], r["to_station"]), []).append(r)
    s = dof = 0.0
    for rs in g.values():
        if len(rs) < 2:
            continue
        x = np.array([r["vector_enu"][comp] for r in rs])
        w = np.array([1 / np.array(r["covariance_enu"])[comp, comp] for r in rs])
        xb = (w * x).sum() / w.sum()
        s += (w * (x - xb) ** 2).sum()
        dof += len(rs) - 1
    return s / dof


def calibration_decomposition(records):
    R = [r for r in records if r["strategy"] == "iflc_float"]
    sel = lambda days, ex=(): [r for r in R if r["session_id"] in days
                              and not set(ex) & {r["from_station"], r["to_station"]}]
    a, b = k2_component(sel(CAL)), k2_component(sel(EVAL))
    dof = repeatability_calibration(sel(CAL)).dof
    F = b / a
    out = {"k2_cal": a, "k2_eval": b, "dof_each": dof, "F": F,
           "F_p_two_sided": float(2 * min(stats.f.sf(F, dof, dof), stats.f.cdf(F, dof, dof))),
           "by_component": {c: {"cal": k2_component(sel(CAL), i), "eval": k2_component(sel(EVAL), i)}
                            for i, c in enumerate("ENU")},
           "eval_without_VILL": k2_component(sel(EVAL, ["VILL"])),
           "all_days_leave_one_out": {d: k2_component(sel([x for x in CAL + EVAL if x != d]))
                                      for d in CAL + EVAL}}
    dev: dict = {}
    for r in R:
        dev.setdefault((r["from_station"], r["to_station"]), {})[r["session_id"]] = np.array(r["vector_enu"])
    pairs = []
    for v in dev.values():
        mtx = np.array([v[d] for d in CAL + EVAL])
        mtx -= mtx.mean(axis=0)
        for c in range(3):
            pairs += list(zip(mtx[:-1, c], mtx[1:, c]))
    x = np.array(pairs)
    out["lag1_correlation_daily_deviations"] = float(np.corrcoef(x[:, 0], x[:, 1])[0, 1])
    out["lag1_expected_iid_demeaned_6days"] = -1 / 5
    return out


def sensitivity_summary():
    path = ROOT / "results/real_network/audit_sensitivity/baselines_all.json"
    recs = json.loads(path.read_text())["records"]
    d = pd.read_csv(path.with_suffix(".csv"))
    d["session"] = d.session.astype(str)
    rows = []
    for st, g in d.groupby("strategy"):
        rec = [r for r in recs if r["strategy"] == st]
        g = g.copy()
        for c in "enu":
            g[f"dev_{c}"] = g[f"diff_ref_{c}_mm"] - g.groupby(["from", "to"])[f"diff_ref_{c}_mm"].transform("mean")
        row = {"strategy": st,
               "k2_cal": repeatability_calibration([r for r in rec if r["session_id"] in CAL]).k2,
               "k2_eval": repeatability_calibration([r for r in rec if r["session_id"] in EVAL]).k2}
        for day in ("2024005", "2024006"):
            v = g[(g.session == day) & ((g["from"] == "VILL") | (g["to"] == "VILL"))]
            row[f"VILL_up_{day[-3:]}_mm"] = float(np.mean(np.where(v["to"] == "VILL", 1, -1) * v.dev_u))
        c = g[(g["from"] == "CEBR") | (g["to"] == "CEBR")]
        e = np.where(c["to"] == "CEBR", 1, -1) * c.diff_ref_e_mm
        row.update(CEBR_east_vs_EPN_mm=float(e.mean()), CEBR_east_sd_mm=float(e.std()),
                   mean_formal_sigma_u_mm=float(g.sigma_u_mm.mean()))
        rows.append(row)
    return pd.DataFrame(rows)


def epn_ac_comparison(log):
    dest = ROOT / "data/raw/audit_epn_ac"
    files = []
    for ac in ACS:
        name = f"{ac}EPNFIN_20233650000_07D_07D_SOL.SNX.gz"
        try:
            log.append(frd.fetch(f"{BKG}/EUREF/products/2295/{name}", dest / name))
            files.append(dest / name)
        except RuntimeError as exc:
            log.append({"url": name, "error": str(exc)})
    eur = read_sinex(ROOT / "data/raw/madrid_2024/products/EUR0OPSSNX_20233650000_07D_07D_SOL.SNX.gz",
                     STATIONS).stations
    rot = frames.enu_rotation_at(eur["MADR"].xyz)
    rows = []
    for f in files:
        s = read_sinex(f, STATIONS).stations
        for x in ("CEBR", "VILL", "YEBE"):
            if x in s and "MADR" in s:
                d = rot @ ((s[x].xyz - s["MADR"].xyz) - (eur[x].xyz - eur["MADR"].xyz)) * 1e3
                rows.append({"AC": f.name[:3], "station": x, "dE_mm": d[0], "dN_mm": d[1], "dU_mm": d[2]})
    return pd.DataFrame(rows)


def antex_comparison(log):
    dest = ROOT / "data/raw/audit_epn_ac/epn_20.atx"
    log.append(frd.fetch(f"{BKG}/EUREF/station/general/epn_20.atx", dest))
    out = {}
    for label, path in (("igs20_2287", ROOT / "data/raw/madrid_2024/products/igs20_2287.atx"),
                        ("epn_20", dest)):
        lines = path.read_text(errors="replace").splitlines()
        hits = [i for i, ln in enumerate(lines) if ln.startswith("SEPCHOKE_B3E6   NONE")
                and "TYPE / SERIAL NO" in ln]
        out[label] = [{"serial": lines[i][20:40].strip(), "method_line": lines[i + 1][:60].strip(),
                       "neu_G01_G02": [lines[j][:30].strip() for j in range(i, i + 40)
                                       if "NORTH / EAST / UP" in lines[j]][:2]} for i in hits]
    return out


def cebr_signal_test():
    pc = yaml.safe_load((ROOT / "configs/real_data/madrid_2024_processing.yaml").read_text())
    m = yaml.safe_load((ROOT / pc["dataset"]).read_text())
    raw, work = ROOT / "data/raw/madrid_2024", ROOT / "data/work/madrid_2024"
    snx = read_sinex(raw / "products/EUR0OPSSNX_20233650000_07D_07D_SOL.SNX.gz", STATIONS).stations
    params = {**pc["processing"], **pc["strategies"]["iflc_float"]}
    rows = []
    for doy in (1, 2, 3):
        obs, aux, excl = prn.day_inputs(m, raw, work, doy)
        for base, rover in (("MADR", "CEBR"), ("CEBR", "YEBE")):
            for sig in ("L2W", "L2L"):
                extra = {**params["extra"], "pos1-exclsats": " ".join(excl)}
                extra["misc-rnxopt1" if rover == "CEBR" else "misc-rnxopt2"] = f"-GL1C -G{sig}"
                cfg = ProcessingConfig(**{**params, "nav_systems": tuple(params["nav_systems"]),
                                          "antex_file": str(raw / "products/igs20_2287.atx"),
                                          "extra": extra})
                b = process_baseline(session_id=f"2024{doy:03d}", rover_obs=obs[rover], base_obs=obs[base],
                                     aux_files=aux, base_xyz=snx[base].xyz,
                                     workdir=work / "audit_signal" / sig / f"d{doy}", cfg=cfg,
                                     rover_station=rover, base_station=base, variant=pc["rtklib_variant"],
                                     selection="last")
                d = frames.vector_ecef_to_enu(b.vector_ecef - (snx[rover].xyz - snx[base].xyz),
                                              snx[base].xyz) * 1e3
                d = d if rover == "CEBR" else -d        # implied CEBR position error
                rows.append({"doy": doy, "baseline": f"{base}->{rover}", "cebr_l2_signal": sig,
                             "cebr_err_e_mm": d[0], "cebr_err_n_mm": d[1], "cebr_err_u_mm": d[2]})
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    log: list = []
    records = json.loads((ROOT / "results/real_network/processing/baselines_all.json").read_text())["records"]
    (OUT / "calibration_decomposition.json").write_text(json.dumps(calibration_decomposition(records), indent=1))
    sensitivity_summary().to_csv(OUT / "sensitivity_summary.csv", index=False)
    epn_ac_comparison(log).to_csv(OUT / "epn_ac_comparison.csv", index=False)
    (OUT / "antex_comparison.json").write_text(json.dumps(antex_comparison(log), indent=1))
    cebr_signal_test().to_csv(OUT / "cebr_l2_signal_test.csv", index=False)
    runs = [rie.experiment(net, net, 30, 4242, [10, 20, 30, 50], ["random"], [1])[0]
            for net in ("quiet_formal", "main_formal")]
    pd.concat(runs).to_csv(OUT / "injection_on_main_network.csv", index=False)
    (OUT / "audit_download_log.json").write_text(json.dumps(log, indent=1))
    print("audit outputs written to", OUT)


if __name__ == "__main__":
    main()
