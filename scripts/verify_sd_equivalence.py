"""Verify the paper's SD/3D identities across synthetic networks."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gnss_adjust.sd_outlier import (  # noqa: E402
    detect_specific_direction,
    eliminate_along_specific_direction,
    eliminate_full_3d,
)
from gnss_adjust.simulation import NetworkConfig, OutlierSpec, simulate_network  # noqa: E402
from gnss_adjust.wls import solve_wls  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--networks", type=int, default=100)
    parser.add_argument("--seed", type=int, default=9000)
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "sd_equivalence")
    parser.add_argument("--identity-tolerance", type=float, default=1e-10)
    parser.add_argument("--coordinate-tolerance-m", type=float, default=1e-8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    max_identity_error = 0.0
    max_parallel_error = 0.0
    max_coordinate_error = 0.0
    tested_blocks = 0
    untestable_blocks = 0
    coordinate_equivalence_cases = 0

    for run in range(args.networks):
        rng = np.random.default_rng(args.seed + run)
        direction = rng.standard_normal(3)
        direction /= np.linalg.norm(direction)
        cfg = NetworkConfig(
            n_stations=8,
            n_baselines=16,
            seed=args.seed + run,
            outliers=[OutlierSpec(
                baseline_id=3, direction=tuple(direction), magnitude_sigma=25.0)],
        )
        stations, observations, _ = simulate_network(cfg)
        detection = detect_specific_direction(solve_wls(stations, observations))
        for test in detection.tests:
            if not test.testable:
                untestable_blocks += 1
                continue
            tested_blocks += 1
            max_identity_error = max(max_identity_error, test.equivalence_error)
            if test.outlier_magnitude > 1e-12:
                max_parallel_error = max(
                    max_parallel_error,
                    float(np.linalg.norm(np.cross(test.direction, test.outlier_vector))),
                )

        # This algebraic equivalence does not depend on crossing a detection
        # threshold in a particular noisy realization.
        testable = [test for test in detection.tests if test.testable]
        if not testable:
            continue
        candidate = max(testable, key=lambda test: test.sd_statistic)
        sd_result = eliminate_along_specific_direction(stations, observations, candidate)
        full_result = eliminate_full_3d(stations, observations, candidate.baseline_id)
        coordinate_error = max(
            float(np.linalg.norm(
                sd_result.coordinates[station_id] - full_result.coordinates[station_id]
            ))
            for station_id in sd_result.coordinates
        )
        max_coordinate_error = max(max_coordinate_error, coordinate_error)
        coordinate_equivalence_cases += 1

    summary = {
        "networks": args.networks,
        "testable_baseline_blocks": tested_blocks,
        "untestable_baseline_blocks": untestable_blocks,
        "coordinate_equivalence_cases": coordinate_equivalence_cases,
        "max_abs_sd_squared_minus_3t": max_identity_error,
        "max_direction_parallel_cross_norm": max_parallel_error,
        "max_sd_vs_full_3d_coordinate_difference_m": max_coordinate_error,
        "identity_tolerance": args.identity_tolerance,
        "coordinate_tolerance_m": args.coordinate_tolerance_m,
        "passed": (
            max_identity_error <= args.identity_tolerance
            and max_coordinate_error <= args.coordinate_tolerance_m
        ),
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    if not summary["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
