"""Weighted least-squares adjustment of a GNSS baseline network.

Solves  l_reduced + v = A x  with  P = Q_l^{-1}  via the normal equations

    N = A' P A,   u = A' P l_reduced,   x_hat = N^{-1} u

using a Cholesky factorisation of ``N`` (never an explicit inverse for the
parameter solve). ``Q_x = N^{-1}`` is obtained by solving against the
identity, which is required to report a covariance matrix.
"""

from __future__ import annotations

import logging

import numpy as np
import scipy.linalg

from . import covariance as cov
from . import network as net
from .types import AdjustmentResult, BaselineObservation, NetworkDiagnostics, Station

logger = logging.getLogger(__name__)


class AdjustmentError(RuntimeError):
    """Raised when the normal equations cannot be solved reliably."""


def solve_wls(
    stations: list[Station],
    observations: list[BaselineObservation],
    sigma0_prior_sq: float = 1.0,
    weight_multipliers: dict[int, float] | None = None,
    validate_covariances: bool = True,
) -> AdjustmentResult:
    """Run one weighted least-squares adjustment.

    Parameters
    ----------
    stations, observations:
        Network definition. Exactly the stations flagged ``is_fixed`` define
        the datum (their coordinates are held).
    sigma0_prior_sq:
        A-priori variance factor. The supplied covariance blocks are taken as
        ``sigma0_prior_sq`` times the cofactor blocks; with the default 1.0
        they are used directly as covariances.
    weight_multipliers:
        Optional robust multiplier ``w_b`` per baseline id. The effective
        weight block is ``w_b * P_b`` (equivalently ``Q_b / w_b``). Baselines
        with ``w_b == 0`` contribute nothing to the normal equations but are
        kept in the residual bookkeeping.
    """
    if validate_covariances:
        for obs in observations:
            cov.validate_covariance_block(obs.covariance, name=f"Q(baseline {obs.baseline_id})")

    system = net.build_design_system(stations, observations)
    A, l_red = system.A, system.l_reduced
    n, p = A.shape
    m = len(observations)

    w_mult = weight_multipliers or {}
    Q_blocks: list[np.ndarray] = []
    P_blocks: list[np.ndarray] = []
    for obs in observations:
        w = float(w_mult.get(obs.baseline_id, 1.0))
        if w < 0:
            raise ValueError(f"negative weight multiplier for baseline {obs.baseline_id}")
        P_b = w * cov.weight_block(obs.covariance)
        if w > 0:
            Q_b = obs.covariance / w
        else:  # excluded observation: infinite variance; store original for reporting
            Q_b = np.full((3, 3), np.inf)
            np.fill_diagonal(Q_b, np.inf)
        Q_blocks.append(Q_b)
        P_blocks.append(P_b)
    P = cov.block_diagonal(P_blocks)

    N = A.T @ P @ A
    u = A.T @ P @ l_red
    rank_A = int(np.linalg.matrix_rank(A))
    try:
        c, low = scipy.linalg.cho_factor(0.5 * (N + N.T))
    except np.linalg.LinAlgError as exc:
        raise AdjustmentError(
            "normal matrix is not positive definite: the network is likely "
            "disconnected, rank deficient, or too many observations have zero weight"
        ) from exc
    x_hat = scipy.linalg.cho_solve((c, low), u)
    Q_x = scipy.linalg.cho_solve((c, low), np.eye(p))  # solve, not explicit inverse
    Q_x = 0.5 * (Q_x + Q_x.T)

    v = A @ x_hat - l_red                      # residual convention: v = A x_hat - l
    vTPv = float(v @ P @ v)
    dof = n - rank_A
    variance_factor = vTPv / dof if dof > 0 else float("nan")

    # Residual covariance. Where an observation was excluded (w=0) its Q block
    # is infinite; report Q_v there as +inf on the diagonal block.
    AQxAT = A @ Q_x @ A.T
    Q_l = cov.block_diagonal([b if np.all(np.isfinite(b)) else np.zeros((3, 3)) for b in Q_blocks])
    Q_v = Q_l - AQxAT
    std_res = np.full(n, np.nan)
    sigma0_prior = float(np.sqrt(sigma0_prior_sq))
    for k, obs in enumerate(observations):
        sl = slice(3 * k, 3 * k + 3)
        if not np.all(np.isfinite(Q_blocks[k])):
            # excluded observation: predicted-residual covariance Q0 + A Qx A'
            Q_v[sl, sl] = obs.covariance + AQxAT[sl, sl]
        d = np.diag(Q_v[sl, sl]).copy()
        d[d < 1e-30] = np.nan  # fully constrained component: w-test undefined
        std_res[sl] = v[sl] / (sigma0_prior * np.sqrt(d))

    cond_N = float(np.linalg.cond(N))
    diag = NetworkDiagnostics(
        n_observations=n,
        n_parameters=p,
        rank_A=rank_A,
        degrees_of_freedom=dof,
        condition_number_N=cond_N,
        connected=net.is_connected(stations, observations),
        datum_station=net.datum_station(stations),
        notes=[],
    )
    if rank_A < p:
        diag.notes.append(f"rank deficiency: rank(A)={rank_A} < parameters={p}")
    if cond_N > 1e12:
        diag.notes.append(f"ill-conditioned normal matrix: cond(N)={cond_N:.3e}")

    coords = {s.station_id: s.coords.copy() for s in stations if s.is_fixed}
    for k, sid in enumerate(system.unknown_ids):
        coords[sid] = x_hat[3 * k: 3 * k + 3].copy()

    adjusted = {
        obs.baseline_id: coords[obs.to_station] - coords[obs.from_station]
        for obs in observations
    }

    return AdjustmentResult(
        coordinates=coords,
        covariance=Q_x,
        unknown_ids=system.unknown_ids,
        baseline_ids=system.baseline_ids,
        residuals=v,
        residual_covariance=Q_v,
        standardized_residuals=std_res,
        adjusted_baselines=adjusted,
        variance_factor=variance_factor,
        sigma0_prior_sq=sigma0_prior_sq,
        vTPv=vTPv,
        degrees_of_freedom=dof,
        diagnostics=diag,
        Q_l=Q_l,
        P=P,
    )
