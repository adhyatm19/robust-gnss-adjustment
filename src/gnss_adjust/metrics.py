"""Evaluation metrics: coordinate accuracy and outlier-detection quality."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class CoordinateMetrics:
    rmse_m: float                       # RMS over all coordinate components
    mean_position_error_m: float        # mean 3D station position error
    max_position_error_m: float
    per_station_error_m: dict[str, float] = field(default_factory=dict)
    normalized_rmse: float | None = None  # rmse / mean predicted coordinate sigma


@dataclass
class DetectionMetrics:
    true_positives: int
    false_positives: int
    false_negatives: int
    precision: float
    recall: float
    f1: float
    false_positive_rate: float
    localization_correct: bool          # flagged set exactly equals true set


def coordinate_metrics(
    true_coords: dict[str, np.ndarray],
    est_coords: dict[str, np.ndarray],
    exclude: set[str] | None = None,
    mean_coord_sigma_m: float | None = None,
) -> CoordinateMetrics:
    """Coordinate accuracy of estimated vs true station positions.

    ``exclude`` normally contains the fixed (datum) station, whose error is
    zero by construction and would dilute the averages.
    """
    exclude = exclude or set()
    ids = [sid for sid in est_coords if sid in true_coords and sid not in exclude]
    if not ids:
        raise ValueError("no comparable stations")
    diffs = np.array([est_coords[sid] - true_coords[sid] for sid in ids])  # (k, 3)
    per_station = {sid: float(np.linalg.norm(d)) for sid, d in zip(ids, diffs)}
    rmse = float(np.sqrt(np.mean(diffs ** 2)))
    pos_err = np.linalg.norm(diffs, axis=1)
    return CoordinateMetrics(
        rmse_m=rmse,
        mean_position_error_m=float(pos_err.mean()),
        max_position_error_m=float(pos_err.max()),
        per_station_error_m=per_station,
        normalized_rmse=(rmse / mean_coord_sigma_m) if mean_coord_sigma_m else None,
    )


def detection_metrics(
    true_outliers: set[int],
    flagged: set[int],
    all_baselines: set[int],
) -> DetectionMetrics:
    """Baseline-level detection quality against known ground truth."""
    tp = len(true_outliers & flagged)
    fp = len(flagged - true_outliers)
    fn = len(true_outliers - flagged)
    tn = len(all_baselines - true_outliers - flagged)
    precision = tp / (tp + fp) if (tp + fp) else (1.0 if not true_outliers else 0.0)
    recall = tp / (tp + fn) if (tp + fn) else 1.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) else 0.0
    fpr = fp / (fp + tn) if (fp + tn) else 0.0
    return DetectionMetrics(
        true_positives=tp, false_positives=fp, false_negatives=fn,
        precision=precision, recall=recall, f1=f1, false_positive_rate=fpr,
        localization_correct=(flagged == true_outliers),
    )
