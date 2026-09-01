"""Run the deterministic seed-42 specific-direction case study.

Usage:
    python scripts/run_sd_case_study.py
    python scripts/run_sd_case_study.py --output results/sd_case_study
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

from gnss_adjust.sd_outlier import (  # noqa: E402
    angular_error_degrees,
    detect_specific_direction,
    eliminate_along_specific_direction,
    eliminate_full_3d,
)
from gnss_adjust.simulation import NetworkConfig, OutlierSpec, simulate_network  # noqa: E402
from gnss_adjust.wls import solve_wls  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--baseline", type=int, default=3)
    parser.add_argument(
        "--vector", type=float, nargs=3, default=(0.035, -0.025, 0.090),
        metavar=("DX", "DY", "DZ"), help="additive baseline bias in metres")
    parser.add_argument("--alpha", type=float, default=0.001)
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "sd_case_study")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = NetworkConfig(
        n_stations=8,
        n_baselines=16,
        seed=args.seed,
        outliers=[OutlierSpec(baseline_id=args.baseline, vector=tuple(args.vector))],
    )
    stations, observations, truth = simulate_network(cfg)
    adjustment = solve_wls(stations, observations)
    detection = detect_specific_direction(adjustment, alpha=args.alpha)
    target = next(test for test in detection.tests if test.baseline_id == args.baseline)
    injected = truth.outlier_vectors[args.baseline]

    sd_adjusted = eliminate_along_specific_direction(stations, observations, target)
    full_3d_adjusted = eliminate_full_3d(stations, observations, args.baseline)
    coordinate_equivalence = max(
        float(np.linalg.norm(
            sd_adjusted.coordinates[station_id] - full_3d_adjusted.coordinates[station_id]
        ))
        for station_id in sd_adjusted.coordinates
    )

    args.output.mkdir(parents=True, exist_ok=True)
    test_rows = [{
        "baseline_id": test.baseline_id,
        "testable": test.testable,
        "significant": test.significant,
        "sd_statistic": test.sd_statistic,
        "statistic_3d": test.statistic_3d,
        "critical_sd": test.critical_sd,
        "critical_3d": test.critical_3d,
        "estimated_dx_m": test.outlier_vector[0],
        "estimated_dy_m": test.outlier_vector[1],
        "estimated_dz_m": test.outlier_vector[2],
        "estimated_magnitude_m": test.outlier_magnitude,
        "direction_x": test.direction[0],
        "direction_y": test.direction[1],
        "direction_z": test.direction[2],
        "reliability_rank": test.reliability_rank,
        "reliability_condition": test.reliability_condition,
        "equivalence_error": test.equivalence_error,
    } for test in detection.tests]
    pd.DataFrame(test_rows).to_csv(args.output / "baseline_tests.csv", index=False)

    summary = {
        "seed": args.seed,
        "planted_baseline": args.baseline,
        "detected_baseline": detection.identified_baseline,
        "injected_vector_m": injected.tolist(),
        "estimated_vector_m": target.outlier_vector.tolist(),
        "injected_magnitude_m": float(np.linalg.norm(injected)),
        "estimated_magnitude_m": target.outlier_magnitude,
        "angular_direction_error_deg": angular_error_degrees(target.direction, injected),
        "sd_statistic": target.sd_statistic,
        "critical_sd": target.critical_sd,
        "statistic_3d": target.statistic_3d,
        "critical_3d": target.critical_3d,
        "sd_squared_minus_3t": target.sd_statistic ** 2 - 3.0 * target.statistic_3d,
        "sd_vs_full_3d_max_coordinate_difference_m": coordinate_equivalence,
    }
    (args.output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    display = pd.DataFrame({"quantity": list(summary), "result": list(summary.values())})
    print(display.to_string(index=False))
    print(f"\nWrote {args.output / 'summary.json'}")


if __name__ == "__main__":
    main()

