"""Reproduce the experiment that drove the IRLS standardization design choice.

Three candidate references for the block-standardised residual
t_b = sqrt(v_b' M^{-1} v_b / 3) were tried during development:

  Scheme A: M = Qv_eff, the residual covariance of the *current weighted*
            solve (predicted covariance Q0 + A Qx A' when w = 0).
  Scheme B: M = Qv0 = Q0 - A Qx A', the original-weight residual covariance
            (eigenvalue-floored).
  Scheme C: M = Q0, i.e. the fixed original observation weight P0 = Q0^{-1}
            — the scheme shipped in gnss_adjust.robust.

Failure modes demonstrated below (deterministic, no randomness in case 1):

  A: the reference moves with the weights — down-weighting an observation
     inflates its effective covariance and shrinks its own statistic, so the
     weight oscillates. On the majority-bridge case this is a limit cycle;
     the gross baseline's weight visits 0 -> 0.08 -> 0.48 -> ... -> 1 -> 0
     and never settles.
  B: the reference *subtracts* the redundancy term, so it is smaller than
     the actual dispersion of a down-weighted observation's residual; t is
     overestimated, the weight drops further, the residual grows — positive
     feedback. On a redundant network with one 0.5 m outlier, ten clean
     baselines are dragged below w = 0.5 and Huber never converges.
  C: the reference is fixed and >= the residual dispersion of fully
     weighted observations, so feedback runs only through the coordinate
     solution: rejecting an outlier grows its own residual and shrinks the
     others' — monotone separation, clean convergence, no collateral
     down-weighting. Classic equivalent-weight practice (Danish method,
     IGG schemes), at the price of mild conservatism absorbed by the
     robust scale options.

Usage:  python scripts/compare_standardizations.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from gnss_adjust import network as net                      # noqa: E402
from gnss_adjust.covariance import weight_block             # noqa: E402
from gnss_adjust.robust import hampel_weight, huber_weight  # noqa: E402
from gnss_adjust.simulation import (                        # noqa: E402
    NetworkConfig, OutlierSpec, simulate_network)
from gnss_adjust.types import BaselineObservation, Station  # noqa: E402
from gnss_adjust.wls import solve_wls                       # noqa: E402


def t_scheme_A(res, obs_list, weights, P0):
    t = {}
    for k, o in enumerate(obs_list):
        v = res.residual_block(k)
        Qv = 0.5 * (res.residual_cov_block(k) + res.residual_cov_block(k).T)
        stat = float(v @ np.linalg.pinv(Qv, rcond=1e-10) @ v)
        t[o.baseline_id] = np.sqrt(max(stat, 0.0) / 3.0)
    return t


def t_scheme_B(res, obs_list, weights, P0):
    t = {}
    for k, o in enumerate(obs_list):
        v = res.residual_block(k)
        Qv_eff = 0.5 * (res.residual_cov_block(k) + res.residual_cov_block(k).T)
        w = weights[o.baseline_id]
        if w > 0:
            AQxAT = o.covariance / w - Qv_eff
            Qv0 = o.covariance - AQxAT
        else:
            Qv0 = Qv_eff
        Qv0 = 0.5 * (Qv0 + Qv0.T)
        ev, V = np.linalg.eigh(Qv0)
        ev = np.clip(ev, 1e-4 * np.trace(o.covariance) / 3.0, None)
        stat = float(v @ (V @ np.diag(1.0 / ev) @ V.T) @ v)
        t[o.baseline_id] = np.sqrt(max(stat, 0.0) / 3.0)
    return t


def t_scheme_C(res, obs_list, weights, P0):
    return {
        o.baseline_id: float(np.sqrt(max(
            res.residual_block(k) @ P0[o.baseline_id] @ res.residual_block(k),
            0.0) / 3.0))
        for k, o in enumerate(obs_list)
    }


T_FUNCS = {"A": t_scheme_A, "B": t_scheme_B, "C": t_scheme_C}


def run(scheme: str, stations, obs_list, wfun, iters: int = 25):
    """Bare IRLS loop with the production rank-guard floor, returning the
    weight trajectory."""
    P0 = {o.baseline_id: weight_block(o.covariance) for o in obs_list}
    weights = {o.baseline_id: 1.0 for o in obs_list}
    traj = []
    for _ in range(iters):
        solve_w = dict(weights)
        positive = [o for o in obs_list if solve_w[o.baseline_id] > 0]
        if len(positive) < len(obs_list) and not net.is_connected(stations, positive):
            solve_w = {b: max(x, 1e-6) for b, x in solve_w.items()}
        res = solve_wls(stations, obs_list, weight_multipliers=solve_w,
                        validate_covariances=False)
        t = T_FUNCS[scheme](res, obs_list, solve_w, P0)
        weights = {b: wfun(tv) for b, tv in t.items()}
        traj.append(dict(weights))
    return traj


def main() -> None:
    # Case 1: majority bridge — 2 clean + 1 gross (0.2 m) B->C observation.
    stations = [
        Station("A", np.zeros(3), is_fixed=True),
        Station("B", np.array([1000.0, 0.0, 0.0])),
        Station("C", np.array([2000.0, 0.0, 0.0])),
    ]
    q = np.eye(3) * 1e-4
    obs = [
        BaselineObservation(0, "A", "B", np.array([1000.0, 0.0, 0.0]), q),
        BaselineObservation(1, "A", "B", np.array([1000.005, 0.0, 0.0]), q),
        BaselineObservation(2, "B", "C", np.array([1000.0, 0.0, 0.0]), q),
        BaselineObservation(3, "B", "C", np.array([1000.003, 0.0, 0.0]), q),
        BaselineObservation(4, "B", "C", np.array([1000.0, 0.0, 0.2]), q),
    ]
    print("Case 1: majority bridge, Hampel — weight of the gross baseline (id 4)")
    for s in ("A", "C"):
        traj = run(s, stations, obs, hampel_weight)
        w4 = [round(w[4], 3) for w in traj]
        tail = w4[-8:]
        verdict = "LIMIT CYCLE" if len(set(tail)) > 1 else f"converged at {tail[-1]}"
        print(f"  scheme {s}: {w4[:12]} ... -> {verdict}")

    # Case 2: redundant network (8 stations, 16 baselines), 0.5 m Z outlier.
    cfg = NetworkConfig(n_stations=8, n_baselines=16, seed=7,
                        outliers=[OutlierSpec(baseline_id=3, vector=(0.0, 0.0, 0.5))])
    st2, obs2, _ = simulate_network(cfg)
    print("\nCase 2: redundant network + 0.5 m outlier on baseline 3, Huber")
    for s in ("B", "C"):
        final = run(s, st2, obs2, huber_weight, iters=40)[-1]
        collateral = sorted(b for b, w in final.items() if b != 3 and w < 0.5)
        print(f"  scheme {s}: w[outlier]={final[3]:.3f}; "
              f"clean baselines dragged below 0.5: {collateral or 'none'}")


if __name__ == "__main__":
    main()
