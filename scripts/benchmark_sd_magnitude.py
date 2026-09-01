"""Sweep injected outlier magnitude and measure SD direction recovery."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gnss_adjust.sd_outlier import angular_error_degrees, detect_specific_direction  # noqa: E402
from gnss_adjust.simulation import NetworkConfig, OutlierSpec, simulate_network  # noqa: E402
from gnss_adjust.wls import solve_wls  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--magnitudes", type=float, nargs="+",
        default=(1, 2, 3, 5, 8, 10, 15, 20, 30, 50),
        help="outlier magnitudes in mean-component sigma units")
    parser.add_argument("--runs", type=int, default=50, help="runs per magnitude")
    parser.add_argument("--seed", type=int, default=5000)
    parser.add_argument("--baseline", type=int, default=3)
    parser.add_argument("--alpha", type=float, default=0.001)
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "sd_magnitude")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    reference = np.array([0.35, -0.25, 0.90], dtype=float)
    reference /= np.linalg.norm(reference)
    rows: list[dict[str, float | int | bool]] = []

    for magnitude in args.magnitudes:
        for run in range(args.runs):
            seed = args.seed + run
            cfg = NetworkConfig(
                n_stations=8,
                n_baselines=16,
                seed=seed,
                outliers=[OutlierSpec(
                    baseline_id=args.baseline,
                    direction=tuple(reference),
                    magnitude_sigma=magnitude,
                )],
            )
            stations, observations, truth = simulate_network(cfg)
            detection = detect_specific_direction(solve_wls(stations, observations), alpha=args.alpha)
            target = next(test for test in detection.tests if test.baseline_id == args.baseline)
            injected = truth.outlier_vectors[args.baseline]
            rows.append({
                "magnitude_sigma": magnitude,
                "run": run,
                "seed": seed,
                "detected_correctly": detection.identified_baseline == args.baseline,
                "target_significant": target.significant,
                "angular_error_deg": angular_error_degrees(target.direction, injected),
                "injected_magnitude_m": float(np.linalg.norm(injected)),
                "estimated_magnitude_m": target.outlier_magnitude,
                "sd_statistic": target.sd_statistic,
                "identity_error": target.equivalence_error,
            })

    frame = pd.DataFrame(rows)
    summary = frame.groupby("magnitude_sigma", as_index=False).agg(
        runs=("run", "count"),
        detection_rate=("detected_correctly", "mean"),
        significance_rate=("target_significant", "mean"),
        median_angular_error_deg=("angular_error_deg", "median"),
        p90_angular_error_deg=("angular_error_deg", lambda values: values.quantile(0.9)),
        median_estimated_magnitude_m=("estimated_magnitude_m", "median"),
        max_identity_error=("identity_error", "max"),
    )

    args.output.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output / "raw_runs.csv", index=False)
    summary.to_csv(args.output / "summary.csv", index=False)

    fig, ax_angle = plt.subplots(figsize=(8.0, 4.8))
    ax_rate = ax_angle.twinx()
    ax_angle.plot(
        summary["magnitude_sigma"], summary["median_angular_error_deg"],
        marker="o", color="#2457A7", label="median angular error")
    ax_angle.fill_between(
        summary["magnitude_sigma"],
        summary["median_angular_error_deg"],
        summary["p90_angular_error_deg"],
        color="#2457A7", alpha=0.15, label="median to 90th percentile")
    ax_rate.plot(
        summary["magnitude_sigma"], summary["detection_rate"],
        marker="s", color="#D1495B", label="correct detection rate")
    ax_angle.set(xlabel="injected magnitude [mean-component sigma]",
                 ylabel="direction error [degrees]",
                 title="Specific-direction recovery versus outlier magnitude")
    ax_rate.set_ylabel("correct detection rate")
    ax_rate.set_ylim(-0.02, 1.02)
    lines = ax_angle.lines + ax_rate.lines
    ax_angle.legend(lines, [line.get_label() for line in lines], loc="upper right")
    ax_angle.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(args.output / "direction_recovery_vs_magnitude.png", dpi=180)
    plt.close(fig)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()

