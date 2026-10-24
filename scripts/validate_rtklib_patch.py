"""Validate the RTKLIB output-precision patch and the covariance reconstruction.

Uses the real GSI 0759/3040 test pair (tests/data/gsi_0759_3040) and both
binaries built by scripts/build_rtklib.sh. For several processing
configurations, with *every epoch* written (out-solstatic=all), it checks:

1. formatting-only: stock and patched (hp) outputs have identical epochs,
   solution status, satellite counts, ages and ratios, and every numeric
   field agrees to within half the stock print resolution (0.05 mm);
2. covariance reconstruction: the ECEF covariance rebuilt from the hp
   ``sdx..sdzx`` fields and rotated to the base ENU frame by
   ``gnss_adjust.gnss.frames`` reproduces RTKLIB's *own* ENU output
   (computed internally from the float32 covariance, hp-printed);
3. stock reconstruction error: the same reconstruction from the stock
   (0.1 mm) fields, compared with the hp ENU reference.

Writes results/rtklib_validation/{summary.json,epochs.csv,stock_vs_hp_enu.png}.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gnss_adjust.gnss import frames  # noqa: E402
from gnss_adjust.gnss.rtklib import (AntennaSetting, ProcessingConfig, binary_info,  # noqa: E402
                                     find_rnx2rtkp, run_rnx2rtkp)
from gnss_adjust.gnss.rtklib_pos import covariance_from_rtklib, parse_pos_file  # noqa: E402

DATA = ROOT / "tests" / "data" / "gsi_0759_3040"
BASE_XYZ = np.array([-3976219.5082, 3382372.5671, 3652512.9849])  # header approx of 0759
CONFIGS = {
    "static_combined": ProcessingConfig(elevation_mask_deg=15, solstatic="all"),
    "static_forward": ProcessingConfig(elevation_mask_deg=15, solstatic="all",
                                       solution_type="forward"),
    "kinematic_forward": ProcessingConfig(elevation_mask_deg=15, solstatic="all",
                                          pos_mode="kinematic", solution_type="forward"),
    "static_l1_forward": ProcessingConfig(elevation_mask_deg=15, solstatic="all",
                                          frequency="l1", solution_type="forward"),
}


def signed_sqrt(c: np.ndarray) -> np.ndarray:
    return np.sign(c) * np.sqrt(np.abs(c))


def enu_fields(q_enu: np.ndarray) -> np.ndarray:
    """RTKLIB ENU print fields (sde, sdn, sdu, sden, sdnu, sdue) of a covariance."""
    return np.array([np.sqrt(q_enu[0, 0]), np.sqrt(q_enu[1, 1]), np.sqrt(q_enu[2, 2]),
                     signed_sqrt(q_enu[0, 1]), signed_sqrt(q_enu[1, 2]),
                     signed_sqrt(q_enu[2, 0])])


def correlations(q: np.ndarray) -> np.ndarray:
    s = np.sqrt(np.diag(q))
    c = q / np.outer(s, s)
    return np.array([c[0, 1], c[1, 2], c[2, 0]])


def run(out_dir: Path) -> tuple[dict, pd.DataFrame]:
    work = ROOT / "data" / "work" / "rtklib_validation"
    ant = AntennaSetting(mode="explicit", antenna_type="TRM29659.00", delta_enu=(0, 0, 0))
    files = (DATA / "30400920.05o", DATA / "07590920.05o",
             [DATA / "07590920.05n", DATA / "30400920.05n"])
    rows = []
    summary: dict = {"binaries": {v: binary_info(find_rnx2rtkp(v)) for v in ("stock", "hp")},
                     "base_xyz": BASE_XYZ.tolist(), "configs": {}}
    for name, cfg in CONFIGS.items():
        parsed = {}
        for variant in ("stock", "hp"):
            for frame in ("xyz", "enu"):
                c = ProcessingConfig(**{**cfg.to_dict(), "extra": {"out-solformat": frame}})
                out = work / f"{name}_{variant}_{frame}.pos"
                run_rnx2rtkp(files[0], files[1], files[2], out, c, BASE_XYZ, ant, ant,
                             variant=variant)
                parsed[(variant, frame)] = parse_pos_file(out, frame=frame).records
        n = len(parsed[("hp", "xyz")])
        # RTKLIB b34 quirk: in combined mode the smoothed base position (rbs in
        # postpos.c) is left unset for equal-status epochs, so ENU output is
        # suppressed; XYZ output does not use it. ENU checks need forward runs.
        frames_used = ("xyz", "enu") if len(parsed[("hp", "enu")]) == n else ("xyz",)
        # 1. formatting-only check
        same_meta = all(
            len(parsed[(v, f)]) == n for v in ("stock", "hp") for f in frames_used
        ) and all(
            (a.time, a.q_status, a.n_sat, a.age, a.ratio) == (b.time, b.q_status, b.n_sat, b.age, b.ratio)
            for f in frames_used for a, b in zip(parsed[("stock", f)], parsed[("hp", f)])
        )
        max_field_dev = 0.0
        for f in frames_used:
            for a, b in zip(parsed[("stock", f)], parsed[("hp", f)]):
                va = np.r_[a.xyz, a.sd, a.sd_off]
                vb = np.r_[b.xyz, b.sd, b.sd_off]
                max_field_dev = max(max_field_dev, float(np.max(np.abs(va - vb))))
        cfg_sum = {"epochs": n, "enu_epochs": len(parsed[("hp", "enu")]),
                   "metadata_identical": bool(same_meta),
                   "max_abs_field_difference_m": max_field_dev,
                   "half_stock_resolution_m": 5e-5}
        # 2./3. reconstruction checks, epoch by epoch
        if "enu" not in frames_used:
            cfg_sum["note"] = ("ENU output suppressed by RTKLIB in combined mode "
                               "(unset smoothed base position); XYZ compared only")
            summary["configs"][name] = cfg_sum
            continue
        rec_err_hp = rec_err_vec = 0.0
        for k in range(n):
            hx, he = parsed[("hp", "xyz")][k], parsed[("hp", "enu")][k]
            sx = parsed[("stock", "xyz")][k]
            q_ref = covariance_from_rtklib(he.sd, he.sd_off)        # RTKLIB's own ENU
            q_hp = frames.covariance_ecef_to_enu(hx.covariance, BASE_XYZ)
            q_st = frames.covariance_ecef_to_enu(sx.covariance, BASE_XYZ)
            rec_err_hp = max(rec_err_hp, float(np.max(np.abs(enu_fields(q_hp) - np.r_[he.sd, he.sd_off]))))
            b_enu = frames.vector_ecef_to_enu(hx.xyz - BASE_XYZ, BASE_XYZ)
            rec_err_vec = max(rec_err_vec, float(np.max(np.abs(b_enu - he.xyz))))
            s_ref, s_st = np.sqrt(np.diag(q_ref)), np.sqrt(np.clip(np.diag(q_st), 0, None))
            c_ref, c_st = correlations(q_ref), correlations(q_st)
            rows.append({
                "config": name, "epoch": hx.time, "q_status": hx.q_status,
                "sd_e_ref_mm": s_ref[0] * 1e3, "sd_n_ref_mm": s_ref[1] * 1e3, "sd_u_ref_mm": s_ref[2] * 1e3,
                "rel_err_sd_e_stock": s_st[0] / s_ref[0] - 1, "rel_err_sd_n_stock": s_st[1] / s_ref[1] - 1,
                "rel_err_sd_u_stock": s_st[2] / s_ref[2] - 1,
                "corr_en_ref": c_ref[0], "corr_nu_ref": c_ref[1], "corr_ue_ref": c_ref[2],
                "corr_err_en_stock": c_st[0] - c_ref[0], "corr_err_nu_stock": c_st[1] - c_ref[1],
                "corr_err_ue_stock": c_st[2] - c_ref[2],
                "corr_sign_flip_stock": bool(np.any(np.sign(c_st) != np.sign(c_ref))),
                "stock_cov_positive_definite": bool(np.linalg.eigvalsh(sx.covariance)[0] > 0),
            })
        cfg_sum["max_abs_enu_field_error_hp_reconstruction_m"] = rec_err_hp
        cfg_sum["max_abs_enu_vector_error_hp_m"] = rec_err_vec
        summary["configs"][name] = cfg_sum
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "epochs.csv", index=False)
    rel = df[["rel_err_sd_e_stock", "rel_err_sd_n_stock", "rel_err_sd_u_stock"]].abs()
    cer = df[["corr_err_en_stock", "corr_err_nu_stock", "corr_err_ue_stock"]].abs()
    summary["stock_reconstruction"] = {
        "epochs": int(len(df)),
        "max_rel_sigma_error": float(rel.max().max()),
        "median_rel_sigma_error": float(np.median(rel.values)),
        "max_abs_correlation_error": float(cer.max().max()),
        "fraction_epochs_with_correlation_sign_flip": float(df["corr_sign_flip_stock"].mean()),
        "fraction_epochs_stock_cov_not_pd": float(1 - df["stock_cov_positive_definite"].mean()),
    }
    summary["verdict"] = {
        "patch_changes_formatting_only": all(
            c["metadata_identical"] and c["max_abs_field_difference_m"] <= 5e-5 + 1e-9
            for c in summary["configs"].values()),
        "hp_reconstruction_matches_rtklib_enu": all(
            c["max_abs_enu_field_error_hp_reconstruction_m"] < 1e-8
            for c in summary["configs"].values() if "note" not in c),
    }
    return summary, df


def plot(df: pd.DataFrame, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for col, lab in (("rel_err_sd_e_stock", "East"), ("rel_err_sd_n_stock", "North"),
                     ("rel_err_sd_u_stock", "Up")):
        axes[0].hist(100 * df[col], bins=40, alpha=0.6, label=lab)
    axes[0].set_xlabel("stock-file σ error vs high-precision reference [%]")
    axes[0].set_ylabel("epochs")
    axes[0].legend()
    for col, lab in (("corr_err_en_stock", "ρ_EN"), ("corr_err_nu_stock", "ρ_NU"),
                     ("corr_err_ue_stock", "ρ_UE")):
        axes[1].hist(df[col], bins=40, alpha=0.6, label=lab)
    axes[1].set_xlabel("stock-file correlation error [-]")
    axes[1].legend()
    fig.suptitle("Covariance reconstructed from 0.1 mm-rounded RTKLIB output\n"
                 "(real GSI 0759→3040 data, all epochs, 4 configurations)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=ROOT / "results" / "rtklib_validation")
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    summary, df = run(args.output)
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2))
    plot(df, args.output / "stock_vs_hp_enu.png")
    print(json.dumps({k: summary[k] for k in ("configs", "stock_reconstruction", "verdict")}, indent=2))


if __name__ == "__main__":
    main()
