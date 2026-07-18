"""Regenerate all demonstration results referenced in docs/EXPERIMENTS.md.

Usage:  python scripts/generate_demo_results.py [--fast]

Runs (deterministic seeds):
  1. clean-data experiment              -> results/clean
  2. directional-outlier experiment     -> results/directional_outlier
  3. Monte Carlo single-outlier         -> results/benchmark_single_outlier
  4. outlier-magnitude sweep            -> results/magnitude_sweep
  5. contamination sweep                -> results/contamination_sweep
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gnss_adjust import benchmark as bm                      # noqa: E402
from gnss_adjust import plotting as plots                    # noqa: E402
from gnss_adjust.cli import run_demo_experiment              # noqa: E402
from gnss_adjust.simulation import NetworkConfig, OutlierSpec  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fast", action="store_true", help="fewer Monte Carlo runs")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    results = ROOT / "results"

    mc_runs = 20 if args.fast else 100
    sweep_runs = 10 if args.fast else 30

    # 1. clean data
    run_demo_experiment(
        NetworkConfig(n_stations=8, n_baselines=16, seed=42),
        bm.ALL_METHODS, results / "clean")

    # 2. directional outlier
    run_demo_experiment(
        NetworkConfig(n_stations=8, n_baselines=16, seed=42,
                      outliers=[OutlierSpec(baseline_id=3, vector=(0.035, -0.025, 0.090))]),
        bm.ALL_METHODS, results / "directional_outlier")

    # 3. Monte Carlo, single 30-sigma outlier
    df = bm.run_monte_carlo("single-outlier", runs=mc_runs, base_seed=1000)
    bm.export_benchmark(df, results / "benchmark_single_outlier",
                        {"scenario": "single-outlier", "runs": mc_runs, "base_seed": 1000})
    errors = {m: df[df.method == m]["mean_pos_error_m"].tolist() for m in df.method.unique()}
    plots.plot_monte_carlo_box(errors, path=results / "benchmark_single_outlier" / "monte_carlo_box.png",
                               title=f"Monte Carlo (single-outlier, {mc_runs} runs)")

    # 4. magnitude sweep
    mags = [0, 3, 5, 8, 12, 20, 30, 50]
    df = bm.run_magnitude_sweep(mags, runs_per_point=sweep_runs, base_seed=5000)
    agg = bm.export_benchmark(df, results / "magnitude_sweep",
                              {"magnitudes_sigma": mags, "runs_per_point": sweep_runs,
                               "base_seed": 5000},
                              by=["outlier_magnitude_sigma", "method"])
    xs = sorted(agg["outlier_magnitude_sigma"].unique())
    rmse = {m: [float(agg[(agg.method == m) & (agg.outlier_magnitude_sigma == x)]["rmse_m_mean"].iloc[0])
                for x in xs] for m in agg.method.unique()}
    f1 = {m: [float(agg[(agg.method == m) & (agg.outlier_magnitude_sigma == x)]["f1_mean"].iloc[0])
              for x in xs] for m in agg.method.unique()}
    plots.plot_metric_vs_x(xs, rmse, "outlier magnitude [sigma]", "coordinate RMSE [m]",
                           "RMSE vs outlier magnitude", results / "magnitude_sweep" / "rmse_vs_magnitude.png",
                           logy=True)
    plots.plot_metric_vs_x(xs, f1, "outlier magnitude [sigma]", "detection F1 [-]",
                           "Detection F1 vs outlier magnitude",
                           results / "magnitude_sweep" / "f1_vs_magnitude.png")

    # 5. contamination sweep
    fracs = [0.0, 0.05, 0.1, 0.2, 0.3]
    df = bm.run_contamination_sweep(fracs, runs_per_point=sweep_runs, base_seed=9000)
    agg = bm.export_benchmark(df, results / "contamination_sweep",
                              {"fractions": fracs, "runs_per_point": sweep_runs,
                               "base_seed": 9000},
                              by=["contamination_fraction", "method"])
    xs = sorted(agg["contamination_fraction"].unique())
    rmse = {m: [float(agg[(agg.method == m) & (agg.contamination_fraction == x)]["rmse_m_mean"].iloc[0])
                for x in xs] for m in agg.method.unique()}
    plots.plot_metric_vs_x(xs, rmse, "contamination fraction [-]", "coordinate RMSE [m]",
                           "RMSE vs contamination fraction",
                           results / "contamination_sweep" / "rmse_vs_contamination.png",
                           logy=True)
    print("\nAll demonstration results regenerated under", results)


if __name__ == "__main__":
    main()
