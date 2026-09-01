"""Reliability matrix and 3x3 block helpers for the SD method.

The paper defines ``P_bar = P Q_v P`` and ``g = P_bar l``.  This project uses
``l + v = A x`` and hence ``v = -Q_v P l``.  Therefore ``g = -P v`` is the
same quantity, evaluated without multiplying kilometre-scale observations by
an annihilator and losing precision through cancellation.
"""

from __future__ import annotations

import numpy as np

from ..types import AdjustmentResult


def reliability_matrix(result: AdjustmentResult) -> np.ndarray:
    """Return the symmetric observation reliability matrix ``P Q_v P``."""
    p_bar = result.P @ result.residual_covariance @ result.P
    return 0.5 * (p_bar + p_bar.T)


def reliability_score(result: AdjustmentResult) -> np.ndarray:
    """Return stacked paper score vectors using ``g = -P v``."""
    return -(result.P @ result.residuals)


def block3(matrix: np.ndarray, row: int, column: int | None = None) -> np.ndarray:
    """Extract a 3x3 baseline block by zero-based positional index."""
    column = row if column is None else column
    rs = slice(3 * row, 3 * row + 3)
    cs = slice(3 * column, 3 * column + 3)
    return np.asarray(matrix[rs, cs], dtype=float)


def vector3(vector: np.ndarray, index: int) -> np.ndarray:
    """Extract one 3-vector by zero-based baseline position."""
    return np.asarray(vector[3 * index:3 * index + 3], dtype=float)


def baseline_reliability(
    result: AdjustmentResult,
    baseline_index: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return ``(P_bar_ii, g_i)`` for a baseline position."""
    if not 0 <= baseline_index < len(result.baseline_ids):
        raise IndexError(f"baseline index {baseline_index} is out of range")
    return (
        block3(reliability_matrix(result), baseline_index),
        vector3(reliability_score(result), baseline_index),
    )

