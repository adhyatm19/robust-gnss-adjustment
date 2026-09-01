"""Preliminary shared-station antenna-height-style SD simulation.

The synthetic generator uses a local Cartesian box rather than physical ECEF
coordinates.  This experiment therefore injects a documented local-up proxy
into up to two baselines from one receiver/session, with the sign set by baseline
orientation.  Other sessions involving the receiver remain clean, following
the paper's N006/session-2 construction.  This is not a production antenna
calibration model.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gnss_adjust.sd_outlier import angular_error_degrees, detect_specific_direction  # noqa: E402
from gnss_adjust.simulation import NetworkConfig, simulate_network  # noqa: E402
from gnss_adjust.wls import solve_wls  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", type=int, default=100)
    parser.add_argument("--seed", type=int, default=11000)
    parser.add_argument("--station", default="S06")
    parser.add_argument("--height-error-m", type=float, default=0.08)
    parser.add_argument("--alpha", type=float, default=0.001)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "results" / "sd_antenna_height")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    local_up = np.array([0.35, 0.20, 0.915], dtype=float)
    local_up /= np.linalg.norm(local_up)
    baseline_rows: list[dict[str, float | int | bool]] = []
    run_rows: list[dict[str, float | int]] = []

    for run in range(args.runs):
        stations, observations, _ = simulate_network(NetworkConfig(
            n_stations=8, n_baselines=16, seed=args.seed + run))
        incident = [
            observation for observation in observations
            if args.station in (observation.from_station, observation.to_station)
        ]
        incident_ids = {observation.baseline_id for observation in incident}
        # The paper corrupts only baselines 3 and 11 from session 2 even though
        # N006 participates in other clean sessions.  Use two deterministic
        # incident baselines as the synthetic session counterpart.
        session_ids = {observation.baseline_id for observation in incident[:2]}
        corrupted = []
        signed_biases: dict[int, np.ndarray] = {}
        for observation in observations:
            sign = 0.0
            if observation.baseline_id in session_ids and observation.to_station == args.station:
                sign = 1.0
            elif observation.baseline_id in session_ids and observation.from_station == args.station:
                sign = -1.0
            bias = sign * args.height_error_m * local_up
            if sign:
                signed_biases[observation.baseline_id] = bias
                corrupted.append(replace(observation, vector=observation.vector + bias))
            else:
                corrupted.append(observation)

        detection = detect_specific_direction(solve_wls(stations, corrupted), alpha=args.alpha)
        incident_significant = 0
        off_station_significant = 0
        angles = []
        for test in detection.tests:
            affected = test.baseline_id in signed_biases
            angle = (
                angular_error_degrees(test.direction, signed_biases[test.baseline_id])
                if affected else float("nan")
            )
            if affected and test.significant:
                incident_significant += 1
            if not affected and test.significant:
                off_station_significant += 1
            if affected and np.isfinite(angle):
                angles.append(angle)
            baseline_rows.append({
                "run": run,
                "seed": args.seed + run,
                "baseline_id": test.baseline_id,
                "affected_session_baseline": affected,
                "incident_on_receiver": test.baseline_id in incident_ids,
                "significant": test.significant,
                "angular_error_deg": angle,
                "sd_statistic": test.sd_statistic,
            })
        run_rows.append({
            "run": run,
            "n_incident": len(signed_biases),
            "incident_significant": incident_significant,
            "off_station_significant": off_station_significant,
            "incident_recall": incident_significant / len(signed_biases),
            "median_incident_angular_error_deg": (
                float(np.median(angles)) if angles else float("nan")
            ),
        })

    baselines = pd.DataFrame(baseline_rows)
    runs = pd.DataFrame(run_rows)
    summary = {
        "runs": args.runs,
        "affected_station": args.station,
        "height_error_m": args.height_error_m,
        "local_up_proxy": local_up.tolist(),
        "mean_affected_session_baselines_per_run": float(runs["n_incident"].mean()),
        "mean_affected_session_recall": float(runs["incident_recall"].mean()),
        "mean_off_station_significant_per_run": float(runs["off_station_significant"].mean()),
        "median_incident_angular_error_deg": float(
            baselines.loc[baselines["affected_session_baseline"], "angular_error_deg"].median()),
        "scope_note": "two-baseline receiver/session error with a local-up proxy; not a physical ECEF antenna model",
    }
    args.output.mkdir(parents=True, exist_ok=True)
    baselines.to_csv(args.output / "baseline_tests.csv", index=False)
    runs.to_csv(args.output / "run_summary.csv", index=False)
    (args.output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    affected_angles = baselines.loc[
        baselines["affected_session_baseline"], "angular_error_deg"].dropna()
    fig, axis = plt.subplots(figsize=(7.5, 4.5))
    axis.hist(affected_angles, bins=20, color="#2457A7", alpha=0.85)
    axis.set(xlabel="sign-invariant direction error [degrees]", ylabel="baseline tests",
             title="Shared-station antenna-height-style direction recovery")
    axis.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(args.output / "incident_direction_error.png", dpi=180)
    plt.close(fig)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
