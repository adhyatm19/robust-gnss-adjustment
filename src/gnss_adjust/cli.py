"""Command-line interface.

Examples
--------
python -m gnss_adjust demo --seed 42 --output results/demo
python -m gnss_adjust demo --seed 42 --method all --outlier-baseline 3 \\
    --outlier-vector 0.035 -0.025 0.090 --output results/directional_outlier
python -m gnss_adjust benchmark --runs 200 --scenario single-outlier \\
    --output results/benchmark
python -m gnss_adjust adjust --network results/real_network/network \\
    --reference results/real_network/network/reference_coordinates.csv \\
    --output results/real_network/experiment_A
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import benchmark as bm
from . import io as gio
from . import plotting as plots
from .diagnostics import baseline_tests, global_model_test
from .simulation import NetworkConfig, NoiseModel, OutlierSpec, simulate_network

logger = logging.getLogger(__name__)

METHOD_CHOICES = ["all", *bm.ALL_METHODS]


def run_demo_experiment(cfg: NetworkConfig, methods: list[str], out_dir: Path,
                        settings: bm.MethodSettings | None = None) -> pd.DataFrame:
    """Simulate one network, run methods, write plots / CSV / JSON / report."""
    out_dir.mkdir(parents=True, exist_ok=True)
    settings = settings or bm.MethodSettings()
    stations, observations, truth = simulate_network(cfg)
    true_out = set(truth.outlier_baselines)
    outcomes = bm.run_all_methods(stations, observations, settings, methods,
                                  true_outliers=true_out)
    rows = bm.score_run(outcomes, truth, stations, observations, scenario="demo")
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "method_comparison.csv", index=False)
    gio.save_json(truth.config, out_dir / "experiment_config.json")
    gio.write_network_csv(stations, observations, out_dir / "network")

    # ---- summary JSON
    summary: dict[str, Any] = {
        "seed": truth.seed,
        "true_outlier_baselines": sorted(true_out),
        "injected_biases_m": {str(k): v for k, v in truth.outlier_vectors.items()},
        "methods": {},
    }
    for name, oc in outcomes.items():
        summary["methods"][name] = {
            "flagged_baselines": sorted(oc.flagged),
            "iterations": oc.iterations,
            "converged": oc.converged,
            "variance_factor": oc.variance_factor,
            "condition_number": oc.condition_number,
            "runtime_s": oc.runtime_s,
        }
    gio.save_json(summary, out_dir / "summary.json")

    # ---- plots
    detected = outcomes["DIA-reject"].flagged if "DIA-reject" in outcomes else set()
    downweighted: set[int] = set()
    weights_by_method: dict[str, dict[int, float]] = {}
    for name in ("Huber", "Hampel"):
        if name in outcomes and outcomes[name].weights is not None:
            weights_by_method[name] = outcomes[name].weights
            downweighted |= outcomes[name].flagged
    plots.plot_network(stations, observations, true_outliers=true_out,
                       detected=detected, downweighted=downweighted,
                       path=out_dir / "network.png")
    plots.plot_true_vs_estimated(
        truth.true_coordinates,
        {name: oc.coordinates for name, oc in outcomes.items()},
        fixed_station=next(s.station_id for s in stations if s.is_fixed),
        path=out_dir / "true_vs_estimated.png")
    fixed_ids = {s.station_id for s in stations if s.is_fixed}
    errors_by_method = {
        name: {sid: float(np.linalg.norm(oc.coordinates[sid] - truth.true_coordinates[sid]))
               for sid in oc.coordinates if sid not in fixed_ids}
        for name, oc in outcomes.items()
    }
    plots.plot_coordinate_errors(errors_by_method, path=out_dir / "coordinate_errors.png")

    tests_by_method = {}
    for name, oc in outcomes.items():
        res = oc.detail.final if hasattr(oc.detail, "final") else oc.detail
        tests_by_method[name] = baseline_tests(res, true_outliers=true_out)
    if "WLS" in outcomes:
        wls_tests = tests_by_method["WLS"]
        plots.plot_standardized_residuals(wls_tests, path=out_dir / "standardized_residuals.png",
                                          title="WLS component w-tests")
        plots.plot_group_statistics(wls_tests, path=out_dir / "group_statistics.png",
                                    title="WLS 3D baseline group statistics")
    plots.plot_residual_norms(tests_by_method, true_outliers=true_out,
                              path=out_dir / "residual_norms.png")
    if weights_by_method:
        plots.plot_robust_weights(weights_by_method, true_outliers=true_out,
                                  path=out_dir / "robust_weights.png")

    _write_demo_report(out_dir, cfg, truth, outcomes, df)
    return df


def _write_demo_report(out_dir: Path, cfg: NetworkConfig, truth: Any,
                       outcomes: dict[str, bm.MethodOutcome], df: pd.DataFrame) -> None:
    lines = [
        "# Demo experiment report",
        "",
        "> Synthetic baseline-level research prototype; not raw GNSS/RINEX processing.",
        "",
        f"- seed: `{truth.seed}`",
        f"- stations: {cfg.n_stations}, baselines: {cfg.n_baselines}",
        f"- true outlier baselines: {sorted(truth.outlier_baselines) or 'none'}",
        "",
        "## Method comparison",
        "",
        df.drop(columns=["run_id", "scenario", "seed"]).to_markdown(index=False, floatfmt=".5g"),
        "",
        "## DIA audit trails",
        "",
    ]
    for name in ("DIA-reject", "DIA-inflate"):
        if name not in outcomes:
            continue
        dia_res = outcomes[name].detail
        lines.append(f"### {name} (status: {dia_res.status})")
        lines.append("")
        lines.append("| iter | T_global | crit | passed | identified | T_b | action |")
        lines.append("|---|---|---|---|---|---|---|")
        for it in dia_res.iterations:
            lines.append(
                f"| {it.iteration} | {it.global_statistic:.2f} | {it.global_critical:.2f} "
                f"| {it.global_passed} | {it.identified_baseline} "
                f"| {'' if it.identified_statistic is None else f'{it.identified_statistic:.2f}'} "
                f"| {it.action} |")
        lines.append("")
    for name in ("Huber", "Hampel"):
        if name not in outcomes:
            continue
        r = outcomes[name].detail
        lines.append(f"### {name}: converged={r.converged} in {r.iterations} iterations; "
                     f"down-weighted baselines: {r.downweighted or 'none'}")
        lines.append("")
    (out_dir / "report.md").write_text("\n".join(lines))
    logger.info("wrote %s", out_dir / "report.md")


def _add_demo_parser(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("demo", help="run one synthetic experiment with plots and report")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--stations", type=int, default=8)
    p.add_argument("--baselines", type=int, default=16)
    p.add_argument("--method", choices=METHOD_CHOICES, default="all")
    p.add_argument("--a-mm", type=float, default=5.0, help="constant noise term [mm]")
    p.add_argument("--b-ppm", type=float, default=1.0, help="length-dependent noise [ppm]")
    p.add_argument("--student-t-dof", type=float, default=None,
                   help="use Student-t noise with this dof (default: Gaussian)")
    p.add_argument("--outlier-baseline", type=int, default=None)
    p.add_argument("--outlier-vector", type=float, nargs=3, default=None,
                   metavar=("DX", "DY", "DZ"), help="explicit outlier bias [m]")
    p.add_argument("--outlier-sigma", type=float, default=None,
                   help="outlier magnitude in sigma units (random direction)")
    p.add_argument("--outlier-m", type=float, default=None,
                   help="outlier magnitude in metres (random direction)")
    p.add_argument("--contamination", type=float, default=None,
                   help="fraction of baselines to corrupt at 20 sigma")
    p.add_argument("--alpha-global", type=float, default=0.05)
    p.add_argument("--alpha-local", type=float, default=0.001)
    p.add_argument("--huber-c", type=float, default=1.5)
    p.add_argument("--hampel-abc", type=float, nargs=3, default=(1.5, 3.5, 8.0),
                   metavar=("A", "B", "C"))
    p.add_argument("--output", type=Path, required=True)


def _add_benchmark_parser(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("benchmark", help="Monte Carlo benchmark over a scenario")
    p.add_argument("--scenario",
                   choices=[*bm.SCENARIOS.keys(), "magnitude-sweep", "contamination-sweep"],
                   default="single-outlier")
    p.add_argument("--runs", type=int, default=50)
    p.add_argument("--seed", type=int, default=1000)
    p.add_argument("--output", type=Path, required=True)


def _add_adjust_parser(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser("adjust", help="adjust a network given as CSV files (no simulation, "
                                      "no ground truth) with every method")
    p.add_argument("--network", type=Path, required=True,
                   help="directory with stations.csv, baselines.csv, covariance_blocks.csv "
                        "(optional baseline_metadata.csv)")
    p.add_argument("--reference", type=Path, default=None,
                   help="optional CSV station_id,x,y,z of independent reference coordinates")
    p.add_argument("--sigma0-prior-sq", type=float, default=1.0,
                   help="a-priori variance factor: true covariance = factor x supplied blocks")
    p.add_argument("--alpha-global", type=float, default=0.05)
    p.add_argument("--alpha-local", type=float, default=0.001)
    p.add_argument("--bonferroni", action="store_true")
    p.add_argument("--huber-c", type=float, default=1.5)
    p.add_argument("--hampel-abc", type=float, nargs=3, default=(1.5, 3.5, 8.0),
                   metavar=("A", "B", "C"))
    p.add_argument("--title", default="Network adjustment")
    p.add_argument("--output", type=Path, required=True)


def _settings_from_args(args: argparse.Namespace) -> bm.MethodSettings:
    settings = bm.MethodSettings()
    settings.dia.alpha_global = args.alpha_global
    settings.dia.alpha_local = args.alpha_local
    settings.dia.alpha_baseline = args.alpha_local
    settings.huber.huber_c = args.huber_c
    a, b, c = args.hampel_abc
    settings.hampel.hampel_a, settings.hampel.hampel_b, settings.hampel.hampel_c = a, b, c
    if hasattr(args, "bonferroni"):
        settings.dia.bonferroni = args.bonferroni
    if hasattr(args, "sigma0_prior_sq"):
        for cfg in (settings.dia, settings.huber, settings.hampel):
            cfg.sigma0_prior_sq = args.sigma0_prior_sq
    return settings


def _run_adjust(args: argparse.Namespace) -> int:
    from . import network_report as nr

    net_dir: Path = args.network
    stations, observations = gio.load_network_from_csv(
        net_dir / "stations.csv", net_dir / "baselines.csv", net_dir / "covariance_blocks.csv")
    meta_path = net_dir / "baseline_metadata.csv"
    metadata = pd.read_csv(meta_path) if meta_path.exists() else None
    reference = nr.read_reference_csv(args.reference) if args.reference else None
    adj = nr.adjust_network(stations, observations, _settings_from_args(args),
                            reference=reference, metadata=metadata)
    preamble = [f"- network: `{net_dir}`",
                f"- reference coordinates: `{args.reference}`" if args.reference else
                "- reference coordinates: none", ""]
    nr.write_outputs(adj, args.output, title=args.title, preamble=preamble)
    rows = [{"method": n, "T": m["global_statistic"], "crit": m["global_critical"],
             "sigma0_hat^2": m["variance_factor"], "flagged": m["flagged_baselines"]}
            for n, m in adj.summary["methods"].items()]
    print(pd.DataFrame(rows).to_string(index=False))
    print(f"\nresults written to {args.output}")
    return 0


def _demo_config_from_args(args: argparse.Namespace) -> NetworkConfig:
    outliers: list[OutlierSpec] = []
    if args.outlier_baseline is not None:
        outliers.append(OutlierSpec(
            baseline_id=args.outlier_baseline,
            vector=tuple(args.outlier_vector) if args.outlier_vector else None,
            magnitude_sigma=args.outlier_sigma,
            magnitude_m=args.outlier_m,
        ))
    elif args.outlier_vector or args.outlier_sigma or args.outlier_m:
        raise SystemExit("--outlier-vector/--outlier-sigma/--outlier-m need --outlier-baseline")
    return NetworkConfig(
        n_stations=args.stations,
        n_baselines=args.baselines,
        noise=NoiseModel(a_mm=args.a_mm, b_ppm=args.b_ppm, student_t_dof=args.student_t_dof),
        outliers=outliers,
        contamination_fraction=args.contamination,
        seed=args.seed,
    )


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    parser = argparse.ArgumentParser(
        prog="gnss_adjust",
        description="Robust GNSS baseline network adjustment "
                    "(baseline-level research prototype; synthetic benchmarks and "
                    "RTKLIB-derived real networks)")
    sub = parser.add_subparsers(dest="command", required=True)
    _add_demo_parser(sub)
    _add_benchmark_parser(sub)
    _add_adjust_parser(sub)
    args = parser.parse_args(argv)

    if args.command == "adjust":
        return _run_adjust(args)

    if args.command == "demo":
        cfg = _demo_config_from_args(args)
        methods = bm.ALL_METHODS if args.method == "all" else [args.method]
        settings = bm.MethodSettings()
        settings.dia.alpha_global = args.alpha_global
        settings.dia.alpha_local = args.alpha_local
        settings.dia.alpha_baseline = args.alpha_local
        settings.huber.huber_c = args.huber_c
        a, b, c = args.hampel_abc
        settings.hampel.hampel_a, settings.hampel.hampel_b, settings.hampel.hampel_c = a, b, c
        df = run_demo_experiment(cfg, methods, args.output, settings)
        print(df[["method", "rmse_m", "f1", "n_flagged", "iterations"]].to_string(index=False))
        print(f"\nresults written to {args.output}")
        return 0

    if args.command == "benchmark":
        out: Path = args.output
        if args.scenario == "magnitude-sweep":
            df = bm.run_magnitude_sweep([0, 3, 5, 8, 12, 20, 30, 50],
                                        runs_per_point=max(args.runs // 8, 1),
                                        base_seed=args.seed)
            agg = bm.export_benchmark(df, out, vars(args) | {"scenario": args.scenario},
                                      by=["outlier_magnitude_sigma", "method"])
            _sweep_plots(agg, "outlier_magnitude_sigma",
                         "outlier magnitude [sigma units]", out)
        elif args.scenario == "contamination-sweep":
            df = bm.run_contamination_sweep([0.0, 0.05, 0.1, 0.2, 0.3],
                                            runs_per_point=max(args.runs // 5, 1),
                                            base_seed=args.seed)
            agg = bm.export_benchmark(df, out, vars(args) | {"scenario": args.scenario},
                                      by=["contamination_fraction", "method"])
            _sweep_plots(agg, "contamination_fraction", "contamination fraction [-]", out)
        else:
            df = bm.run_monte_carlo(args.scenario, args.runs, base_seed=args.seed)
            agg = bm.export_benchmark(df, out, vars(args))
            errors = {m: df[df.method == m]["mean_pos_error_m"].tolist()
                      for m in df.method.unique()}
            plots.plot_monte_carlo_box(
                errors, path=out / "monte_carlo_box.png",
                title=f"Monte Carlo ({args.scenario}, {args.runs} runs)")
        print(agg.to_string(index=False))
        print(f"\nresults written to {out}")
        return 0
    return 1


def _sweep_plots(agg: pd.DataFrame, xcol: str, xlabel: str, out: Path) -> None:
    xs = sorted(agg[xcol].unique())
    rmse = {m: [float(agg[(agg.method == m) & (agg[xcol] == x)]["rmse_m_mean"].iloc[0])
                for x in xs] for m in agg.method.unique()}
    f1 = {m: [float(agg[(agg.method == m) & (agg[xcol] == x)]["f1_mean"].iloc[0])
              for x in xs] for m in agg.method.unique()}
    plots.plot_metric_vs_x(xs, rmse, xlabel, "coordinate RMSE [m]",
                           f"Coordinate RMSE vs {xlabel}", out / "rmse_vs_x.png", logy=True)
    plots.plot_metric_vs_x(xs, f1, xlabel, "detection F1 [-]",
                           f"Detection F1 vs {xlabel}", out / "f1_vs_x.png")
