"""End-to-end: CSV round-trip, demo pipeline, benchmark plumbing."""

import numpy as np
import pandas as pd
import pytest

from gnss_adjust import benchmark as bm
from gnss_adjust import io as gio
from gnss_adjust.cli import main as cli_main
from gnss_adjust.simulation import NetworkConfig, OutlierSpec, simulate_network
from gnss_adjust.wls import solve_wls


def test_csv_roundtrip_matches_in_memory(tmp_path, simulated_clean):
    stations, obs, _ = simulated_clean
    gio.write_network_csv(stations, obs, tmp_path)
    stations2, obs2 = gio.load_network_from_csv(
        tmp_path / "stations.csv", tmp_path / "baselines.csv",
        tmp_path / "covariance_blocks.csv")
    r1 = solve_wls(stations, obs)
    r2 = solve_wls(stations2, obs2)
    for sid in r1.coordinates:
        np.testing.assert_allclose(r2.coordinates[sid], r1.coordinates[sid], atol=1e-9)
    assert r2.vTPv == pytest.approx(r1.vTPv, rel=1e-9)


def test_demo_pipeline_produces_files(tmp_path):
    rc = cli_main([
        "demo", "--seed", "42", "--stations", "7", "--baselines", "13",
        "--outlier-baseline", "3", "--outlier-vector", "0.0", "0.0", "0.4",
        "--output", str(tmp_path / "demo"),
    ])
    assert rc == 0
    out = tmp_path / "demo"
    for f in ["method_comparison.csv", "experiment_config.json", "summary.json",
              "network.png", "true_vs_estimated.png", "coordinate_errors.png",
              "standardized_residuals.png", "group_statistics.png",
              "residual_norms.png", "robust_weights.png", "report.md",
              "network/stations.csv"]:
        assert (out / f).exists(), f"missing {f}"
    df = pd.read_csv(out / "method_comparison.csv")
    assert set(df.method) == set(bm.ALL_METHODS)
    # the planted 0.4 m outlier must be caught by DIA
    dia_row = df[df.method == "DIA-reject"].iloc[0]
    assert dia_row.tp == 1 and dia_row.fn == 0


def test_benchmark_monte_carlo_small(tmp_path):
    df = bm.run_monte_carlo("single-outlier", runs=3, base_seed=100)
    assert len(df) == 3 * len(bm.ALL_METHODS)
    agg = bm.export_benchmark(df, tmp_path, {"runs": 3})
    assert (tmp_path / "runs.csv").exists()
    assert (tmp_path / "aggregate.csv").exists()
    assert (tmp_path / "config.json").exists()
    assert set(agg.method) == set(bm.ALL_METHODS)
    # robust/DIA should beat plain WLS on average under a 30-sigma outlier
    rmse = agg.set_index("method")["rmse_m_mean"]
    assert rmse["DIA-reject"] < rmse["WLS"]
    assert rmse["Hampel"] < rmse["WLS"]


def test_student_t_scenario_runs():
    cfg = bm.SCENARIOS["student-t"](seed=3)
    stations, obs, truth = simulate_network(cfg)
    outcomes = bm.run_all_methods(stations, obs)
    assert set(outcomes) == set(bm.ALL_METHODS)
