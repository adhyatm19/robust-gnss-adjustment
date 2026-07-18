"""Core data structures for GNSS baseline network adjustment.

Conventions used throughout the package
---------------------------------------
* All coordinates and baseline vectors are in metres.
* All covariance matrices are in metres squared.
* A baseline observation from station ``i`` to station ``j`` observes

      l_ij = X_j - X_i + e_ij,      Cov(e_ij) = Q_ij  (full 3x3 block)

* The residual sign convention is

      v = A @ x_hat - l_reduced

  i.e. residual = (adjusted value) - (observed value). This convention is
  used consistently in code, tests, and documentation.
* ``Q`` denotes a covariance (or cofactor) matrix, ``P = Q^{-1}`` a weight
  matrix. When an a-priori variance factor ``sigma0^2`` differs from 1 the
  supplied blocks are treated as covariances scaled by that factor.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

Vector3 = np.ndarray  # shape (3,)
Matrix3 = np.ndarray  # shape (3, 3)


@dataclass(frozen=True)
class Station:
    """A network station with 3D Cartesian coordinates (metres)."""

    station_id: str
    coords: Vector3
    is_fixed: bool = False

    def __post_init__(self) -> None:
        coords = np.asarray(self.coords, dtype=float).reshape(3)
        object.__setattr__(self, "coords", coords)


@dataclass(frozen=True)
class BaselineObservation:
    """One observed 3D baseline vector with its full covariance block.

    ``vector`` is the observed (X_j - X_i) in metres; ``covariance`` is the
    full 3x3 covariance block in m^2.
    """

    baseline_id: int
    from_station: str
    to_station: str
    vector: Vector3
    covariance: Matrix3

    def __post_init__(self) -> None:
        object.__setattr__(self, "vector", np.asarray(self.vector, dtype=float).reshape(3))
        object.__setattr__(self, "covariance", np.asarray(self.covariance, dtype=float).reshape(3, 3))

    def reversed(self) -> "BaselineObservation":
        """Return the same physical baseline observed in the opposite direction."""
        return BaselineObservation(
            baseline_id=self.baseline_id,
            from_station=self.to_station,
            to_station=self.from_station,
            vector=-self.vector,
            covariance=self.covariance.copy(),
        )


@dataclass
class DesignSystem:
    """The linear system  l_reduced + v = A x  for the free parameters.

    ``unknown_ids`` gives the station ordering: station ``unknown_ids[k]``
    occupies parameter columns ``3k : 3k+3``.
    """

    A: np.ndarray                      # (3m, 3u) design matrix
    l_reduced: np.ndarray              # (3m,) observations minus fixed-station part
    unknown_ids: list[str]
    baseline_ids: list[int]
    fixed: dict[str, np.ndarray]       # fixed station coordinates


@dataclass
class NetworkDiagnostics:
    """Structural / numerical health report for an adjustment problem."""

    n_observations: int                # 3m
    n_parameters: int                  # 3u
    rank_A: int
    degrees_of_freedom: int
    condition_number_N: float
    connected: bool
    datum_station: str | None
    notes: list[str] = field(default_factory=list)


@dataclass
class AdjustmentResult:
    """Full output of one weighted least-squares adjustment."""

    coordinates: dict[str, np.ndarray]         # adjusted coords of ALL stations (fixed included)
    covariance: np.ndarray                     # Q_x for the unknown parameters (3u x 3u)
    unknown_ids: list[str]
    baseline_ids: list[int]
    residuals: np.ndarray                      # v = A x_hat - l_reduced, (3m,)
    residual_covariance: np.ndarray            # Q_v = Q_l - A Q_x A^T
    standardized_residuals: np.ndarray         # v_i / (sigma0_prior * sqrt(Qv_ii))
    adjusted_baselines: dict[int, np.ndarray]  # A x_hat mapped back per baseline
    variance_factor: float                     # a-posteriori sigma0^2 = v'Pv / r
    sigma0_prior_sq: float                     # a-priori variance factor used for tests
    vTPv: float
    degrees_of_freedom: int
    diagnostics: NetworkDiagnostics
    Q_l: np.ndarray                            # block-diagonal observation covariance used
    P: np.ndarray                              # weight matrix used

    def residual_block(self, k: int) -> np.ndarray:
        """Residual 3-vector of the k-th baseline (by position, not id)."""
        return self.residuals[3 * k: 3 * k + 3]

    def residual_cov_block(self, k: int) -> np.ndarray:
        """3x3 diagonal block of Q_v for the k-th baseline (by position)."""
        return self.residual_covariance[3 * k: 3 * k + 3, 3 * k: 3 * k + 3]


@dataclass
class BaselineTest:
    """Identification statistics for a single baseline."""

    baseline_id: int
    residual: np.ndarray               # (3,) metres
    residual_norm: float               # metres
    standardized_components: np.ndarray  # (3,) w-statistics
    max_abs_component_stat: float
    component_critical: float
    group_statistic: float             # v_b' pinv(Qv_bb) v_b / sigma0_prior^2
    group_critical: float
    group_df: int
    flagged_component: bool
    flagged_group: bool
    is_true_outlier: bool | None = None


@dataclass
class DIAIteration:
    """Audit-trail record of a single DIA iteration."""

    iteration: int
    global_statistic: float
    global_critical: float
    global_passed: bool
    identified_baseline: int | None
    identified_statistic: float | None
    local_critical: float | None
    action: str                        # "none" | "reject" | "inflate" | "stopped:<reason>"
    active_baselines: list[int]
    max_coordinate_change: float | None
    baseline_tests: list[BaselineTest] = field(default_factory=list)


@dataclass
class DIAResult:
    """Outcome of a full Detection-Identification-Adaptation run."""

    final: AdjustmentResult
    iterations: list[DIAIteration]
    flagged_baselines: list[int]       # baselines rejected or inflated, in order
    status: str                        # "clean" | "adapted" | "max_iterations" | "no_candidate" | "rank_guard"
    adaptation: str                    # "reject" | "inflate"


@dataclass
class RobustResult:
    """Outcome of a Huber/Hampel iteratively-reweighted adjustment."""

    final: AdjustmentResult
    method: str
    weights: dict[int, float]          # final robust multiplier per baseline
    t_values: dict[int, float]         # final standardized block residuals
    iterations: int
    converged: bool
    objective_history: list[float]
    max_coord_change_history: list[float]
    max_weight_change_history: list[float]
    scale_history: list[float]
    downweighted: list[int]            # baselines with final weight < 0.5


@dataclass
class SimulationTruth:
    """Ground truth of a synthetic experiment."""

    true_coordinates: dict[str, np.ndarray]
    outlier_baselines: list[int]
    outlier_vectors: dict[int, np.ndarray]     # injected bias per corrupted baseline
    seed: int
    config: dict[str, Any] = field(default_factory=dict)
