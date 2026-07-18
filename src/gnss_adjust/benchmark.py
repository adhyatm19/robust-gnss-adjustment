"""Benchmarking framework: run all methods on synthetic scenarios and score them.

Every number exported by this module is produced by actually running the
estimators — nothing is hard-coded or post-edited.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

from . import metrics as met
from .dia import DIAConfig, run_dia
from .robust import RobustConfig, robust_adjust
from .simulation import NetworkConfig, NoiseModel, OutlierSpec, simulate_network
from .types import BaselineObservation, SimulationTruth, Station
from .wls import solve_wls

logger = logging.getLogger(__name__)

ALL_METHODS = ["WLS", "DIA-reject", "DIA-inflate", "Huber", "Hampel"]


@dataclass
class MethodSettings:
    """Tuning shared by every benchmark run."""

    dia: DIAConfig = field(default_factory=DIAConfig)
    huber: RobustConfig = field(default_factory=lambda: RobustConfig(method="huber"))
    hampel: RobustConfig = field(default_factory=lambda: RobustConfig(method="hampel"))


@dataclass
class MethodOutcome:
    """Everything recorded about one method on one dataset."""

    method: str
    coordinates: dict[str, np.ndarray]
    flagged: set[int]                    # rejected / inflated / strongly down-weighted
    weights: dict[int, float] | None
    iterations: int
    converged: bool
    runtime_s: float
    condition_number: float
    variance_factor: float
    detail: Any                          # AdjustmentResult / DIAResult / RobustResult


def run_all_methods(
    stations: list[Station],
    observations: list[BaselineObservation],
    settings: MethodSettings | None = None,
    methods: list[str] | None = None,
    true_outliers: set[int] | None = None,
) -> dict[str, MethodOutcome]:
    """Run the selected estimation methods on one dataset."""
    st = settings or MethodSettings()
    methods = methods or ALL_METHODS
    out: dict[str, MethodOutcome] = {}

    for method in methods:
        t0 = time.perf_counter()
        if method == "WLS":
            r = solve_wls(stations, observations, sigma0_prior_sq=st.dia.sigma0_prior_sq)
            out[method] = MethodOutcome(
                method=method, coordinates=r.coordinates, flagged=set(), weights=None,
                iterations=1, converged=True, runtime_s=time.perf_counter() - t0,
                condition_number=r.diagnostics.condition_number_N,
                variance_factor=r.variance_factor, detail=r)
        elif method in ("DIA-reject", "DIA-inflate"):
            cfg_kwargs = {**st.dia.__dict__, "adaptation": "reject" if method == "DIA-reject" else "inflate"}
            r = run_dia(stations, observations, DIAConfig(**cfg_kwargs), true_outliers=true_outliers)
            out[method] = MethodOutcome(
                method=method, coordinates=r.final.coordinates,
                flagged=set(r.flagged_baselines), weights=None,
                iterations=len(r.iterations), converged=r.status in ("clean", "adapted"),
                runtime_s=time.perf_counter() - t0,
                condition_number=r.final.diagnostics.condition_number_N,
                variance_factor=r.final.variance_factor, detail=r)
        elif method in ("Huber", "Hampel"):
            cfg = st.huber if method == "Huber" else st.hampel
            r = robust_adjust(stations, observations, cfg)
            out[method] = MethodOutcome(
                method=method, coordinates=r.final.coordinates,
                flagged=set(r.downweighted), weights=r.weights,
                iterations=r.iterations, converged=r.converged,
                runtime_s=time.perf_counter() - t0,
                condition_number=r.final.diagnostics.condition_number_N,
                variance_factor=r.final.variance_factor, detail=r)
        else:
            raise ValueError(f"unknown method {method!r}")
    return out


def score_run(
    outcomes: dict[str, MethodOutcome],
    truth: SimulationTruth,
    stations: list[Station],
    observations: list[BaselineObservation],
    run_id: int = 0,
    scenario: str = "",
    extra: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """One CSV row per method with coordinate and detection metrics."""
    fixed = {s.station_id for s in stations if s.is_fixed}
    all_ids = {o.baseline_id for o in observations}
    true_out = set(truth.outlier_baselines)
    rows = []
    for method, oc in outcomes.items():
        cm = met.coordinate_metrics(truth.true_coordinates, oc.coordinates, exclude=fixed)
        dm = met.detection_metrics(true_out, oc.flagged, all_ids)
        rows.append({
            "run_id": run_id, "scenario": scenario, "method": method,
            "seed": truth.seed,
            "rmse_m": cm.rmse_m,
            "mean_pos_error_m": cm.mean_position_error_m,
            "max_pos_error_m": cm.max_position_error_m,
            "tp": dm.true_positives, "fp": dm.false_positives, "fn": dm.false_negatives,
            "precision": dm.precision, "recall": dm.recall, "f1": dm.f1,
            "fpr": dm.false_positive_rate,
            "localization_correct": dm.localization_correct,
            "iterations": oc.iterations, "converged": oc.converged,
            "runtime_s": oc.runtime_s,
            "condition_number": oc.condition_number,
            "variance_factor": oc.variance_factor,
            "n_flagged": len(oc.flagged),
            "n_true_outliers": len(true_out),
            **(extra or {}),
        })
    return rows


# ----------------------------------------------------------------------------
# Scenario registry
# ----------------------------------------------------------------------------

def _base_config(seed: int, **kw: Any) -> NetworkConfig:
    return NetworkConfig(n_stations=kw.pop("n_stations", 8),
                         n_baselines=kw.pop("n_baselines", 16),
                         seed=seed, **kw)


def scenario_clean(seed: int) -> NetworkConfig:
    return _base_config(seed)


def scenario_single_small(seed: int) -> NetworkConfig:
    return _base_config(seed, outliers=[OutlierSpec(baseline_id=3, magnitude_sigma=5.0)])


def scenario_single_large(seed: int) -> NetworkConfig:
    return _base_config(seed, outliers=[OutlierSpec(baseline_id=3, magnitude_sigma=30.0)])


def scenario_directional(seed: int) -> NetworkConfig:
    return _base_config(seed, outliers=[OutlierSpec(baseline_id=3, vector=(0.035, -0.025, 0.090))])


def scenario_multiple(seed: int) -> NetworkConfig:
    return _base_config(seed, outliers=[
        OutlierSpec(baseline_id=2, magnitude_sigma=25.0),
        OutlierSpec(baseline_id=7, magnitude_sigma=18.0),
        OutlierSpec(baseline_id=11, magnitude_sigma=35.0),
    ])


def scenario_student_t(seed: int) -> NetworkConfig:
    return _base_config(seed, noise=NoiseModel(student_t_dof=3.0))


def scenario_weak_geometry(seed: int) -> NetworkConfig:
    # minimal redundancy: spanning tree + 2 extra baselines
    return _base_config(seed, n_baselines=9,
                        outliers=[OutlierSpec(baseline_id=3, magnitude_sigma=30.0)])


def scenario_misscaled_cov(seed: int) -> NetworkConfig:
    # reported covariances 4x too optimistic
    return _base_config(seed, noise=NoiseModel(covariance_scale_error=0.25),
                        outliers=[OutlierSpec(baseline_id=3, magnitude_sigma=30.0)])


SCENARIOS: dict[str, Callable[[int], NetworkConfig]] = {
    "clean": scenario_clean,
    "single-small-outlier": scenario_single_small,
    "single-outlier": scenario_single_large,
    "directional-outlier": scenario_directional,
    "multiple-outliers": scenario_multiple,
    "student-t": scenario_student_t,
    "weak-geometry": scenario_weak_geometry,
    "misscaled-covariance": scenario_misscaled_cov,
}


# ----------------------------------------------------------------------------
# Monte Carlo and sweeps
# ----------------------------------------------------------------------------

def run_monte_carlo(
    scenario: str,
    runs: int,
    base_seed: int = 1000,
    settings: MethodSettings | None = None,
    methods: list[str] | None = None,
) -> pd.DataFrame:
    """Repeat one scenario ``runs`` times with consecutive seeds."""
    if scenario not in SCENARIOS:
        raise ValueError(f"unknown scenario {scenario!r}; choose from {sorted(SCENARIOS)}")
    factory = SCENARIOS[scenario]
    rows: list[dict[str, Any]] = []
    for r in range(runs):
        cfg = factory(base_seed + r)
        stations, observations, truth = simulate_network(cfg)
        outcomes = run_all_methods(stations, observations, settings, methods,
                                   true_outliers=set(truth.outlier_baselines))
        rows.extend(score_run(outcomes, truth, stations, observations,
                              run_id=r, scenario=scenario))
        if (r + 1) % 25 == 0:
            logger.info("%s: %d/%d runs done", scenario, r + 1, runs)
    return pd.DataFrame(rows)


def run_magnitude_sweep(
    magnitudes_sigma: list[float],
    runs_per_point: int,
    base_seed: int = 5000,
    settings: MethodSettings | None = None,
) -> pd.DataFrame:
    """Outlier magnitude sweep (in sigma units) on the single-outlier scenario."""
    rows: list[dict[str, Any]] = []
    for mag in magnitudes_sigma:
        for r in range(runs_per_point):
            cfg = _base_config(base_seed + r,
                               outliers=[OutlierSpec(baseline_id=3, magnitude_sigma=mag)]
                               if mag > 0 else [])
            stations, observations, truth = simulate_network(cfg)
            outcomes = run_all_methods(stations, observations, settings,
                                       true_outliers=set(truth.outlier_baselines))
            rows.extend(score_run(outcomes, truth, stations, observations, run_id=r,
                                  scenario="magnitude-sweep",
                                  extra={"outlier_magnitude_sigma": mag}))
    return pd.DataFrame(rows)


def run_contamination_sweep(
    fractions: list[float],
    runs_per_point: int,
    base_seed: int = 9000,
    settings: MethodSettings | None = None,
) -> pd.DataFrame:
    """Contamination-fraction sweep with 20-sigma outliers."""
    rows: list[dict[str, Any]] = []
    for frac in fractions:
        for r in range(runs_per_point):
            cfg = _base_config(base_seed + r,
                               contamination_fraction=frac if frac > 0 else None)
            stations, observations, truth = simulate_network(cfg)
            outcomes = run_all_methods(stations, observations, settings,
                                       true_outliers=set(truth.outlier_baselines))
            rows.extend(score_run(outcomes, truth, stations, observations, run_id=r,
                                  scenario="contamination-sweep",
                                  extra={"contamination_fraction": frac}))
    return pd.DataFrame(rows)


def aggregate(df: pd.DataFrame, by: list[str] | None = None) -> pd.DataFrame:
    """Aggregate per-run rows to summary statistics per method (and sweep value)."""
    by = by or ["scenario", "method"]
    numeric = ["rmse_m", "mean_pos_error_m", "max_pos_error_m", "precision",
               "recall", "f1", "fpr", "iterations", "runtime_s",
               "variance_factor", "n_flagged"]
    agg = df.groupby(by, sort=False)[numeric].agg(["mean", "median", "std"])
    agg.columns = ["_".join(c) for c in agg.columns]
    conv = df.groupby(by, sort=False)["converged"].mean().rename("convergence_rate")
    loc = df.groupby(by, sort=False)["localization_correct"].mean().rename("localization_rate")
    return pd.concat([agg, conv, loc], axis=1).reset_index()


def export_benchmark(df: pd.DataFrame, out_dir: str | Path, config: dict[str, Any],
                     by: list[str] | None = None) -> pd.DataFrame:
    """Write per-run CSV, aggregate CSV, and the JSON config; return the aggregate."""
    from . import io as gio
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "runs.csv", index=False)
    agg = aggregate(df, by=by)
    agg.to_csv(out / "aggregate.csv", index=False)
    gio.save_json(config, out / "config.json")
    logger.info("benchmark results written to %s", out)
    return agg
