"""Empirical calibration of formal baseline covariances from session repeatability.

For a station pair observed in ``D`` independent sessions with baseline
estimates ``b_d`` and formal covariances ``Q_d``, the weighted mean is
``b_bar = (sum P_d)^-1 sum P_d b_d`` and

    S = sum_d (b_d - b_bar)' P_d (b_d - b_bar)  ~  k^2 chi2(3(D-1))

if the true covariances are ``k^2 Q_d`` and the sessions are stochastically
independent. Pooling over pairs gives ``k_hat^2 = sum S / sum dof`` with an
exact chi-square confidence interval under those assumptions.

What this does and does not measure
-----------------------------------
Repeatability captures only error components that *vary between sessions*.
Errors constant over all sessions (antenna-model errors, reference-coordinate
or station-specific biases, day-repeating multipath) are invisible to it, so
``k^2`` calibrates *precision*, not *accuracy*. Positive correlation between
sessions (e.g. sidereally repeating multipath) biases ``k^2`` low.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import stats


@dataclass
class RepeatabilityResult:
    k2: float                          # pooled variance factor (true cov = k2 * formal)
    k2_ci95: tuple[float, float]
    dof: int
    n_pairs: int
    n_sessions: int
    per_pair: list[dict] = field(default_factory=list)
    enu_variance_ratio: dict[str, float] = field(default_factory=dict)  # empirical / formal


def repeatability_calibration(records: list[dict]) -> RepeatabilityResult:
    """Pooled variance factor from repeated pairs.

    ``records`` need ``from_station``, ``to_station``, ``session_id``,
    ``vector_ecef`` (3,), ``covariance_ecef`` (3x3), and ``vector_enu`` /
    ``covariance_enu`` for the component-wise ratios.
    """
    groups: dict[tuple[str, str], list[dict]] = {}
    for r in records:
        groups.setdefault((r["from_station"], r["to_station"]), []).append(r)
    s_tot, dof_tot = 0.0, 0
    per_pair = []
    num_enu = np.zeros(3)
    den_enu = np.zeros(3)
    for (a, b), rs in sorted(groups.items()):
        if len(rs) < 2:
            continue
        sessions = [r["session_id"] for r in rs]
        if len(set(sessions)) != len(sessions):
            raise ValueError(f"{a}->{b}: repeated session ids; sessions must be distinct")
        bs = [np.asarray(r["vector_ecef"], dtype=float) for r in rs]
        ps = [np.linalg.inv(np.asarray(r["covariance_ecef"], dtype=float)) for r in rs]
        bbar = np.linalg.solve(sum(ps), sum(p @ x for p, x in zip(ps, bs)))
        s = float(sum((x - bbar) @ p @ (x - bbar) for p, x in zip(ps, bs)))
        dof = 3 * (len(rs) - 1)
        s_tot += s
        dof_tot += dof
        # component-wise (ENU of the from-station) unweighted variance ratio
        e = np.array([r["vector_enu"] for r in rs])
        f = np.array([np.diag(r["covariance_enu"]) for r in rs])
        num_enu += ((e - e.mean(axis=0)) ** 2).sum(axis=0)
        den_enu += f.mean(axis=0) * (len(rs) - 1)
        per_pair.append({"from": a, "to": b, "sessions": len(rs), "S": s, "dof": dof,
                         "k2_pair": s / dof,
                         "empirical_sd_enu_mm": (1e3 * e.std(axis=0, ddof=1)).tolist(),
                         "formal_sd_enu_mm": (1e3 * np.sqrt(f.mean(axis=0))).tolist()})
    if dof_tot == 0:
        raise ValueError("no station pair observed in two or more sessions")
    k2 = s_tot / dof_tot
    lo = s_tot / stats.chi2.ppf(0.975, dof_tot)
    hi = s_tot / stats.chi2.ppf(0.025, dof_tot)
    return RepeatabilityResult(
        k2=k2, k2_ci95=(float(lo), float(hi)), dof=dof_tot, n_pairs=len(per_pair),
        n_sessions=len({r["session_id"] for r in records}), per_pair=per_pair,
        enu_variance_ratio={c: float(num_enu[i] / den_enu[i]) for i, c in enumerate("ENU")})
