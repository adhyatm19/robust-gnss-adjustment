"""Sweep random sphere orientations at a fixed outlier magnitude."""

from __future__ import annotations

import argparse
import json
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
    parser.add_argument("--runs", type=int, default=200)
    parser.add_argument("--magnitude", type=float, default=20.0,
                        help="magnitude in mean-component sigma units")
    parser.add_argument("--seed", type=int, default=7000)
    parser.add_argument("--baseline", type=int, default=3)
    parser.add_argument("--alpha", type=float, default=0.001)
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "sd_direction")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rng = np.random.default_rng(args.seed + 99173)
    rows: list[dict[str, float | int | bool]] = []

    for run in range(args.runs):
        direction = rng.standard_normal(3)
        direction /= np.linalg.norm(direction)
        seed = args.seed + run
        cfg = NetworkConfig(
            n_stations=8,
            n_baselines=16,
            seed=seed,
            outliers=[OutlierSpec(
                baseline_id=args.baseline,
                direction=tuple(direction),
                magnitude_sigma=args.magnitude,
            )],
        )
        stations, observations, truth = simulate_network(cfg)
        detection = detect_specific_direction(solve_wls(stations, observations), alpha=args.alpha)
        target = next(test for test in detection.tests if test.baseline_id == args.baseline)
        injected = truth.outlier_vectors[args.baseline]
        rows.append({
            "run": run,
            "seed": seed,
            "direction_x": direction[0],
            "direction_y": direction[1],
            "direction_z": direction[2],
            "longitude_rad": float(np.arctan2(direction[1], direction[0])),
            "latitude_rad": float(np.arcsin(direction[2])),
            "angular_error_deg": angular_error_degrees(target.direction, injected),
            "detected_correctly": detection.identified_baseline == args.baseline,
            "target_significant": target.significant,
            "sd_statistic": target.sd_statistic,
            "identity_error": target.equivalence_error,
        })

    frame = pd.DataFrame(rows)
    summary = {
        "runs": args.runs,
        "magnitude_sigma": args.magnitude,
        "correct_detection_rate": float(frame["detected_correctly"].mean()),
        "target_significance_rate": float(frame["target_significant"].mean()),
        "median_angular_error_deg": float(frame["angular_error_deg"].median()),
        "p90_angular_error_deg": float(frame["angular_error_deg"].quantile(0.9)),
        "maximum_angular_error_deg": float(frame["angular_error_deg"].max()),
        "maximum_identity_error": float(frame["identity_error"].max()),
    }
    args.output.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output / "raw_runs.csv", index=False)
    (args.output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    fig = plt.figure(figsize=(9.0, 5.2))
    axis = fig.add_subplot(111, projection="mollweide")
    points = axis.scatter(
        frame["longitude_rad"], frame["latitude_rad"],
        c=frame["angular_error_deg"], cmap="viridis", s=28, alpha=0.85)
    missed = frame.loc[~frame["detected_correctly"]]
    if len(missed):
        axis.scatter(missed["longitude_rad"], missed["latitude_rad"],
                     marker="x", c="red", s=45, label="mislocalized")
        axis.legend(loc="lower left")
    axis.grid(alpha=0.3)
    axis.set_title("SD angular recovery across injected orientations", pad=18)
    colorbar = fig.colorbar(points, ax=axis, pad=0.08, shrink=0.85)
    colorbar.set_label("sign-invariant angular error [degrees]")
    fig.tight_layout()
    fig.savefig(args.output / "orientation_sweep.png", dpi=180)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

