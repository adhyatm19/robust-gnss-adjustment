"""Construction and validation of 3x3 baseline covariance blocks.

Terminology (kept strict throughout the package):

* **covariance matrix** ``Q`` — actual second moments in m^2;
* **cofactor matrix** — covariance divided by an a-priori variance factor
  ``sigma0^2``; when ``sigma0^2 = 1`` (the package default) the two coincide;
* **weight matrix** ``P = Q^{-1}``.
"""

from __future__ import annotations

import logging

import numpy as np
import scipy.linalg

logger = logging.getLogger(__name__)

MM_PER_M = 1000.0
PPM = 1.0e-6


class CovarianceError(ValueError):
    """Raised when a covariance block fails validation."""


def baseline_sigma_m(a_mm: float, b_ppm: float, length_m: float) -> float:
    """Baseline standard deviation (metres) from the classic GNSS error model.

        sigma = sqrt( a_mm^2 + (b_ppm * L)^2 )

    ``a_mm`` is the constant part in millimetres; ``b_ppm`` scales with the
    baseline length ``length_m`` (metres); the result is returned in metres.
    """
    if a_mm < 0 or b_ppm < 0 or length_m < 0:
        raise ValueError("a_mm, b_ppm and length_m must be non-negative")
    a_m = a_mm / MM_PER_M
    b_m = b_ppm * PPM * length_m
    return float(np.hypot(a_m, b_m))


def component_sigmas(sigma_m: float, axis_scale: tuple[float, float, float]) -> np.ndarray:
    """Per-component standard deviations, allowing X/Y/Z precision to differ.

    ``axis_scale`` multiplies the base sigma per axis; e.g. ``(1, 1, 2)``
    makes the Z (height-dominated) component twice as noisy.
    """
    scale = np.asarray(axis_scale, dtype=float)
    if scale.shape != (3,) or np.any(scale <= 0):
        raise ValueError("axis_scale must be three positive numbers")
    return sigma_m * scale


def random_correlation(rng: np.random.Generator, max_abs_corr: float = 0.4) -> np.ndarray:
    """A random valid 3x3 correlation matrix with moderate off-diagonals.

    Correlations are drawn uniformly then the matrix is projected to positive
    definiteness by eigenvalue flooring and re-normalised to unit diagonal.
    """
    if not 0 <= max_abs_corr < 1:
        raise ValueError("max_abs_corr must be in [0, 1)")
    r_xy, r_xz, r_yz = rng.uniform(-max_abs_corr, max_abs_corr, size=3)
    R = np.array([[1.0, r_xy, r_xz], [r_xy, 1.0, r_yz], [r_xz, r_yz, 1.0]])
    w, V = np.linalg.eigh(R)
    if w.min() <= 1e-6:
        w = np.clip(w, 1e-6, None)
        R = V @ np.diag(w) @ V.T
        d = np.sqrt(np.diag(R))
        R = R / np.outer(d, d)
    return R


def build_covariance_block(sigmas_m: np.ndarray, correlation: np.ndarray) -> np.ndarray:
    """Assemble ``Q_b = D R D`` from component sigmas and a correlation matrix."""
    D = np.diag(np.asarray(sigmas_m, dtype=float))
    Q = D @ np.asarray(correlation, dtype=float) @ D
    return 0.5 * (Q + Q.T)  # exact symmetry


def validate_covariance_block(Q: np.ndarray, name: str = "Q", sym_tol: float = 1e-10) -> None:
    """Check symmetry and positive definiteness; raise ``CovarianceError`` otherwise."""
    Q = np.asarray(Q, dtype=float)
    if Q.shape != (3, 3):
        raise CovarianceError(f"{name}: expected shape (3, 3), got {Q.shape}")
    if not np.all(np.isfinite(Q)):
        raise CovarianceError(f"{name}: contains non-finite entries")
    scale = max(np.abs(Q).max(), 1e-300)
    if np.abs(Q - Q.T).max() > sym_tol * scale:
        raise CovarianceError(f"{name}: not symmetric within tolerance")
    try:
        np.linalg.cholesky(0.5 * (Q + Q.T))
    except np.linalg.LinAlgError as exc:
        raise CovarianceError(f"{name}: not positive definite") from exc


def weight_block(Q: np.ndarray) -> np.ndarray:
    """Weight matrix ``P = Q^{-1}`` via a Cholesky solve against the identity."""
    c, low = scipy.linalg.cho_factor(0.5 * (Q + Q.T))
    P = scipy.linalg.cho_solve((c, low), np.eye(3))
    return 0.5 * (P + P.T)


def block_diagonal(blocks: list[np.ndarray]) -> np.ndarray:
    """Dense block-diagonal matrix from a list of 3x3 blocks."""
    return scipy.linalg.block_diag(*blocks) if blocks else np.zeros((0, 0))
