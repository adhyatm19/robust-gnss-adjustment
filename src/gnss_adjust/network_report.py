"""Adjust a given baseline network with every method and tabulate the results.

Unlike :mod:`gnss_adjust.benchmark`, nothing here needs (or assumes) ground
truth: it is the entry point for real networks. Detection *metrics*
(precision/recall) are deliberately absent -- for unmodified real data there
are no outlier labels, and a statistical flag is not proof of a physical
error. Where independent reference coordinates are supplied, coordinate
differences are reported relative to the datum station (translation-free).

All residuals, directions and coordinate differences are additionally
expressed in the local East/North/Up frame (of the baseline's from-station,
or of the station itself) because ECEF components are hard to interpret.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import benchmark as bm
from .diagnostics import baseline_tests, global_model_test
from .gnss import frames
from .sd_outlier import SDConfig, detect_specific_direction, run_sd_snooping
from .types import AdjustmentResult, BaselineObservation, Station

logger = logging.getLogger(__name__)

ALL_REAL_METHODS = [*bm.ALL_METHODS, "SD-reject"]


@dataclass
class NetworkAdjustment:
    stations: list[Station]
    observations: list[BaselineObservation]
    outcomes: dict[str, bm.MethodOutcome]
    sd_local: Any                               # SpecificDirectionDetection on WLS
    sd_snooping: Any                            # SDSnoopingResult
    coordinates: pd.DataFrame
    baselines: pd.DataFrame
    summary: dict[str, Any]
    reference: dict[str, np.ndarray] | None = None
    notes: list[str] = field(default_factory=list)


def _final(detail: Any) -> AdjustmentResult:
    return detail.final if hasattr(detail, "final") else detail


def local_redundancy(result: AdjustmentResult) -> np.ndarray:
    """Per-baseline redundancy number tr(Qv_bb P_bb) (0 = uncontrolled, 3 = fully redundant)."""
    out = []
    for k in range(len(result.baseline_ids)):
        s = slice(3 * k, 3 * k + 3)
        out.append(float(np.trace(result.residual_covariance[s, s] @ result.P[s, s])))
    return np.array(out)


def adjust_network(
    stations: list[Station],
    observations: list[BaselineObservation],
    settings: bm.MethodSettings | None = None,
    sd_config: SDConfig | None = None,
    reference: dict[str, np.ndarray] | None = None,
    metadata: pd.DataFrame | None = None,
) -> NetworkAdjustment:
    """Run WLS, DIA (reject/inflate), Huber, Hampel and SD snooping on one network."""
    settings = settings or bm.MethodSettings()
    sd_cfg = sd_config or SDConfig(alpha_global=settings.dia.alpha_global,
                                   alpha_local=settings.dia.alpha_baseline,
                                   bonferroni=settings.dia.bonferroni,
                                   sigma0_prior_sq=settings.dia.sigma0_prior_sq)
    outcomes = bm.run_all_methods(stations, observations, settings)
    wls = outcomes["WLS"].detail
    sd_local = detect_specific_direction(wls, alpha=sd_cfg.alpha_local, bonferroni=sd_cfg.bonferroni)
    sd_snoop = run_sd_snooping(stations, observations, sd_cfg)
    outcomes["SD-reject"] = bm.MethodOutcome(
        method="SD-reject", coordinates=sd_snoop.final.coordinates,
        flagged=set(sd_snoop.flagged_baselines), weights=None,
        iterations=len(sd_snoop.iterations), converged=sd_snoop.status in ("clean", "adapted"),
        runtime_s=float("nan"), condition_number=sd_snoop.final.diagnostics.condition_number_N,
        variance_factor=sd_snoop.final.variance_factor, detail=sd_snoop)

    smap = {s.station_id: s for s in stations}
    fixed = [s.station_id for s in stations if s.is_fixed]
    datum = fixed[0]
    sigma0_sq = settings.dia.sigma0_prior_sq

    # ---------------- coordinates
    coord_rows = []
    for name, oc in outcomes.items():
        res = _final(oc.detail)
        for sid, xyz in oc.coordinates.items():
            row: dict[str, Any] = {"method": name, "station": sid, "fixed": sid in fixed,
                                   "x": xyz[0], "y": xyz[1], "z": xyz[2]}
            r_enu = frames.local_rotation(xyz)
            if sid in res.unknown_ids:
                k = res.unknown_ids.index(sid)
                qx = sigma0_sq * res.covariance[3 * k:3 * k + 3, 3 * k:3 * k + 3]
                qe = r_enu @ qx @ r_enu.T
                row.update(sigma_e_mm=1e3 * np.sqrt(qe[0, 0]), sigma_n_mm=1e3 * np.sqrt(qe[1, 1]),
                           sigma_u_mm=1e3 * np.sqrt(qe[2, 2]))
            if reference is not None and sid in reference and datum in reference:
                # translation-free comparison: (X_i - X_datum)_est - (X_i - X_datum)_ref
                d = (xyz - oc.coordinates[datum]) - (reference[sid] - reference[datum])
                de = r_enu @ d
                row.update(diff_e_mm=1e3 * de[0], diff_n_mm=1e3 * de[1], diff_u_mm=1e3 * de[2],
                           diff_3d_mm=1e3 * float(np.linalg.norm(d)))
            coord_rows.append(row)
    coords = pd.DataFrame(coord_rows)

    # ---------------- baselines
    meta = metadata.set_index("baseline_id") if metadata is not None else None
    base_rows = []
    sd_by_id = {t.baseline_id: t for t in sd_local.tests}
    for name, oc in outcomes.items():
        res = _final(oc.detail)
        tests = baseline_tests(res, alpha_local=settings.dia.alpha_local,
                               alpha_baseline=settings.dia.alpha_baseline,
                               bonferroni=settings.dia.bonferroni)
        rnum = local_redundancy(res)
        by_id = {o.baseline_id: o for o in observations}
        for k, t in enumerate(tests):
            o = by_id[t.baseline_id]
            r_enu = frames.local_rotation(smap[o.from_station].coords)
            v_enu = r_enu @ t.residual
            row = {
                "method": name, "baseline_id": t.baseline_id,
                "from": o.from_station, "to": o.to_station,
                "length_m": float(np.linalg.norm(o.vector)),
                "v_x_mm": 1e3 * t.residual[0], "v_y_mm": 1e3 * t.residual[1],
                "v_z_mm": 1e3 * t.residual[2],
                "v_e_mm": 1e3 * v_enu[0], "v_n_mm": 1e3 * v_enu[1], "v_u_mm": 1e3 * v_enu[2],
                "v_norm_mm": 1e3 * t.residual_norm,
                "max_abs_w": t.max_abs_component_stat, "w_critical": t.component_critical,
                "T_b": t.group_statistic, "T_b_critical": t.group_critical, "T_b_df": t.group_df,
                "flag_group": t.flagged_group, "local_redundancy": rnum[k],
                "flagged_by_method": t.baseline_id in oc.flagged,
                "robust_weight": (oc.weights or {}).get(t.baseline_id, np.nan),
            }
            if name == "WLS":
                sd = sd_by_id[t.baseline_id]
                row.update(sd_statistic=sd.sd_statistic, sd_critical=sd.critical_sd,
                           sd_significant=sd.significant, sd_testable=sd.testable,
                           sd_bias_norm_mm=1e3 * sd.outlier_magnitude)
                if sd.testable and np.all(np.isfinite(sd.direction)):
                    u = r_enu @ sd.direction
                    u = u if u[2] >= 0 else -u   # direction is defined up to sign
                    row.update(sd_dir_e=u[0], sd_dir_n=u[1], sd_dir_u=u[2],
                               sd_dir_elevation_deg=float(np.degrees(np.arcsin(np.clip(u[2], -1, 1)))))
                    # angular uncertainty of the bias direction from Cov(d_hat)=sigma0^2 Pbar_ii^-1
                    cov_d = sigma0_sq * np.linalg.inv(sd.reliability_block)
                    dvec = sd.outlier_vector
                    proj = np.eye(3) - np.outer(dvec, dvec) / (dvec @ dvec)
                    perp = np.linalg.eigvalsh(proj @ cov_d @ proj)[-1]
                    row["sd_dir_sigma_deg"] = float(np.degrees(np.arctan2(np.sqrt(perp),
                                                                          np.linalg.norm(dvec))))
            if meta is not None and t.baseline_id in meta.index:
                for col in meta.columns:
                    row[f"meta_{col}"] = meta.loc[t.baseline_id, col]
            base_rows.append(row)
    baselines = pd.DataFrame(base_rows)

    # ---------------- summary
    summary: dict[str, Any] = {
        "local_frame": ("ENU (ECEF rotated at each station)" if frames.looks_geocentric(
            smap[datum].coords) else "network's own Cartesian axes (not geocentric; no rotation)"),
        "datum_station": datum,
        "n_stations": len(stations), "n_baselines": len(observations),
        "n_observations": wls.diagnostics.n_observations,
        "n_parameters": wls.diagnostics.n_parameters,
        "rank_A": wls.diagnostics.rank_A,
        "degrees_of_freedom": wls.degrees_of_freedom,
        "condition_number_N": wls.diagnostics.condition_number_N,
        "connected": wls.diagnostics.connected,
        "sigma0_prior_sq": sigma0_sq,
        "settings": {"dia": vars(settings.dia), "huber": vars(settings.huber),
                     "hampel": vars(settings.hampel), "sd": vars(sd_cfg)},
        "methods": {},
    }
    for name, oc in outcomes.items():
        res = _final(oc.detail)
        g = global_model_test(res, alpha=settings.dia.alpha_global)
        entry = {"global_statistic": g.statistic, "global_critical": g.critical,
                 "global_dof": g.dof, "global_passed": g.passed,
                 "variance_factor": res.variance_factor,
                 "flagged_baselines": sorted(oc.flagged), "iterations": oc.iterations,
                 "converged": oc.converged}
        if hasattr(oc.detail, "status"):
            entry["status"] = oc.detail.status
        if oc.weights is not None:
            entry["weights"] = {str(k): v for k, v in sorted(oc.weights.items())}
        summary["methods"][name] = entry
    summary["sd_local_identified"] = sd_local.identified_baseline
    summary["sd_snooping_audit"] = [
        {"iteration": it.iteration, "global_statistic": it.global_statistic,
         "global_critical": it.global_critical, "identified": it.identified_baseline,
         "statistic": it.identified_statistic, "critical": it.local_critical, "action": it.action}
        for it in sd_snoop.iterations]
    summary["dia_audit"] = {
        name: [{"iteration": it.iteration, "global_statistic": it.global_statistic,
                "global_critical": it.global_critical, "identified": it.identified_baseline,
                "statistic": it.identified_statistic, "action": it.action}
               for it in outcomes[name].detail.iterations]
        for name in ("DIA-reject", "DIA-inflate") if name in outcomes}
    return NetworkAdjustment(stations=stations, observations=observations, outcomes=outcomes,
                             sd_local=sd_local, sd_snooping=sd_snoop, coordinates=coords,
                             baselines=baselines, summary=summary, reference=reference)


def read_reference_csv(path: str | Path) -> dict[str, np.ndarray]:
    """``station_id,x,y,z`` (extra columns such as epoch/frame are ignored)."""
    df = pd.read_csv(path)
    return {str(r.station_id): np.array([r.x, r.y, r.z], dtype=float) for r in df.itertuples()}


def write_outputs(adj: NetworkAdjustment, out_dir: Path, title: str = "Network adjustment",
                  preamble: list[str] | None = None) -> None:
    """CSV/JSON tables, figures and a Markdown report."""
    from . import io as gio
    from . import plotting as plots

    out_dir.mkdir(parents=True, exist_ok=True)
    adj.coordinates.to_csv(out_dir / "coordinates.csv", index=False)
    adj.baselines.to_csv(out_dir / "baselines.csv", index=False)
    gio.save_json(adj.summary, out_dir / "summary.json")

    wls_res = adj.outcomes["WLS"].detail
    wls_tests = baseline_tests(wls_res)
    flagged_any = set().union(*(oc.flagged for oc in adj.outcomes.values()))
    plots.plot_network_enu(adj.stations, adj.observations, flagged=flagged_any,
                           path=out_dir / "network_enu.png")
    plots.plot_group_statistics(wls_tests, path=out_dir / "group_statistics.png",
                                title="WLS 3D baseline test statistics")
    plots.plot_residual_norms({n: baseline_tests(_final(oc.detail))
                               for n, oc in adj.outcomes.items()},
                              path=out_dir / "residual_norms.png")
    weights = {n: oc.weights for n, oc in adj.outcomes.items() if oc.weights is not None}
    if weights:
        plots.plot_robust_weights(weights, path=out_dir / "robust_weights.png")

    s = adj.summary
    lines = [f"# {title}", ""]
    lines += preamble or []
    lines += [
        f"- stations: {s['n_stations']}, baselines: {s['n_baselines']}, "
        f"observations: {s['n_observations']}, parameters: {s['n_parameters']}",
        f"- datum: station `{s['datum_station']}` held fixed; rank(A) = {s['rank_A']}; "
        f"redundancy r = {s['degrees_of_freedom']}; cond(N) = {s['condition_number_N']:.3g}",
        f"- a-priori variance factor sigma0^2 = {s['sigma0_prior_sq']:g}",
        "",
        "## Global model test and method outcomes",
        "",
        "| method | T = v'Pv/sigma0^2 | critical | passed | sigma0_hat^2 | flagged | iterations |",
        "|---|---|---|---|---|---|---|",
    ]
    for name, m in s["methods"].items():
        lines.append(f"| {name} | {m['global_statistic']:.2f} | {m['global_critical']:.2f} | "
                     f"{m['global_passed']} | {m['variance_factor']:.3f} | "
                     f"{m['flagged_baselines'] or '-'} | {m['iterations']} |")
    lines += ["", "Statistical flags indicate inconsistency with the stochastic model; they do "
              "not by themselves prove a physical observation error.", ""]
    if "diff_3d_mm" in adj.coordinates:
        lines += ["## Coordinate differences to the reference (relative to the datum station)", ""]
        sub = adj.coordinates[~adj.coordinates.fixed][["method", "station", "diff_e_mm",
                                                        "diff_n_mm", "diff_u_mm", "sigma_e_mm",
                                                        "sigma_n_mm", "sigma_u_mm"]]
        lines += [sub.to_markdown(index=False, floatfmt=".2f"), ""]
    wls_b = adj.baselines[adj.baselines.method == "WLS"]
    cols = ["baseline_id", "from", "to", "length_m", "v_e_mm", "v_n_mm", "v_u_mm", "T_b",
            "T_b_critical", "local_redundancy", "sd_statistic", "sd_critical", "sd_dir_elevation_deg"]
    lines += ["## WLS residuals and local tests (ENU of from-station)", "",
              wls_b[[c for c in cols if c in wls_b]].to_markdown(index=False, floatfmt=".3g"), ""]
    (out_dir / "report.md").write_text("\n".join(lines))
