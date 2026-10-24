"""Independent verification of the SD statistic against the 3D outlier test.

The existing ``test_sd_equivalence`` checks ``|w*|^2 = 3 T`` where ``T`` is
*defined* in the SD code as ``|w*|^2 / 3`` -- an internal-consistency check.
This module verifies the substantive identities with code paths that share
nothing with ``solve_wls`` / ``sd_outlier`` / ``diagnostics``:

* a dense generalised-least-squares solver written here (Cholesky whitening
  + ``numpy.linalg.lstsq``, own design-matrix construction);
* the 3D outlier statistic defined from first principles as the
  likelihood-ratio reduction of the weighted sum of squares when baseline
  ``i`` receives a free 3D bias parameter (sigma0_prior^2 = 1):

      T_LR(i) = Omega_0 - Omega_a(i),   Omega = v' Q^{-1} v.

Expected (Teunissen 2000, Ch. 3-4; Nie et al. 2019):
  * for any positive-definite Q_l:  T_LR = (C'Pv)' (C'P Q_v P C)^{-1} (C'Pv),
    which is exactly the SD quantity ``g' Pbar_ii^{-1} g = |w*|^2``;
  * only for block-diagonal Q_l does this reduce to the DIA group statistic
    ``v_b' Q_v,bb^{-1} v_b``;
  * under H0 with a correctly scaled covariance, T ~ chi2(3).
"""

import numpy as np
import pytest
from scipy import stats

from gnss_adjust.diagnostics import baseline_tests
from gnss_adjust.sd_outlier import detect_specific_direction
from gnss_adjust.simulation import NetworkConfig, OutlierSpec, simulate_network
from gnss_adjust.wls import solve_wls


def _design(stations, observations):
    """Own design-matrix construction (independent of gnss_adjust.network)."""
    free = sorted(s.station_id for s in stations if not s.is_fixed)
    fixed = {s.station_id: s.coords for s in stations if s.is_fixed}
    col = {sid: 3 * k for k, sid in enumerate(free)}
    m = len(observations)
    A = np.zeros((3 * m, 3 * len(free)))
    l = np.zeros(3 * m)
    for k, o in enumerate(observations):
        rows = slice(3 * k, 3 * k + 3)
        l[rows] = o.vector
        for sid, sign in ((o.to_station, 1.0), (o.from_station, -1.0)):
            if sid in fixed:
                l[rows] -= sign * fixed[sid]
            else:
                A[rows, col[sid]:col[sid] + 3] = sign * np.eye(3)
    return A, l


def _gls(A, Q, l):
    """Dense GLS by whitening: returns (x_hat, v = A x - l, Omega)."""
    L = np.linalg.cholesky(Q)
    Aw = np.linalg.solve(L, A)
    lw = np.linalg.solve(L, l)
    x, *_ = np.linalg.lstsq(Aw, lw, rcond=None)
    v = A @ x - l
    return x, v, float(v @ np.linalg.solve(Q, v))


def _lr_statistic(A, Q, l, k):
    """Likelihood-ratio 3D outlier statistic and bias estimate for baseline k."""
    C = np.zeros((A.shape[0], 3))
    C[3 * k:3 * k + 3] = np.eye(3)
    _, _, om0 = _gls(A, Q, l)
    # model l = A x + C d + e  <=>  l + v = [A, C] [x; d] in this sign convention
    xa, _, oma = _gls(np.hstack([A, C]), Q, l)
    return om0 - oma, xa[-3:]


def _baarda_general(A, Q, l, k):
    P = np.linalg.inv(Q)
    _, v, _ = _gls(A, Q, l)
    Qv = Q - A @ np.linalg.solve(A.T @ P @ A, A.T)
    C = np.zeros((A.shape[0], 3))
    C[3 * k:3 * k + 3] = np.eye(3)
    g = C.T @ P @ v
    return float(g @ np.linalg.solve(C.T @ P @ Qv @ P @ C, g))


def _block_form(A, Q, l, k):
    P = np.linalg.inv(Q)
    _, v, _ = _gls(A, Q, l)
    Qv = Q - A @ np.linalg.solve(A.T @ P @ A, A.T)
    s = slice(3 * k, 3 * k + 3)
    return float(v[s] @ np.linalg.solve(Qv[s, s], v[s]))


@pytest.mark.parametrize("seed", [3, 7, 11])
@pytest.mark.parametrize("outlier", [False, True])
def test_lr_equals_sd_and_dia_for_block_diagonal_covariance(seed, outlier):
    outliers = [OutlierSpec(baseline_id=2, magnitude_sigma=12)] if outlier else []
    stations, obs, _ = simulate_network(NetworkConfig(n_stations=7, n_baselines=14, seed=seed,
                                                      outliers=outliers))
    A, l = _design(stations, obs)
    Q = np.zeros((A.shape[0], A.shape[0]))
    for k, o in enumerate(obs):
        Q[3 * k:3 * k + 3, 3 * k:3 * k + 3] = o.covariance
    res = solve_wls(stations, obs)
    sd = detect_specific_direction(res)
    dia = baseline_tests(res)
    checked = 0
    for k in range(len(obs)):
        if not sd.tests[k].testable:
            continue
        t_lr, d_lr = _lr_statistic(A, Q, l, k)
        assert sd.tests[k].sd_statistic ** 2 == pytest.approx(t_lr, rel=1e-7, abs=1e-9)
        assert dia[k].group_statistic == pytest.approx(t_lr, rel=1e-7, abs=1e-9)
        np.testing.assert_allclose(sd.tests[k].outlier_vector, d_lr, atol=1e-8)
        checked += 1
    assert checked >= 10


def test_general_correlated_covariance_baarda_vs_block_form():
    """With cross-baseline covariance the DIA block statistic is not the 3D test."""
    stations, obs, _ = simulate_network(NetworkConfig(n_stations=6, n_baselines=12, seed=5))
    A, l = _design(stations, obs)
    n = A.shape[0]
    rng = np.random.default_rng(0)
    blocks = np.zeros((n, n))
    for k, o in enumerate(obs):
        blocks[3 * k:3 * k + 3, 3 * k:3 * k + 3] = o.covariance
    # correlated model: D^(1/2) R D^(1/2) with R mixing all baselines (SPD)
    d = np.sqrt(np.diag(blocks))
    B = rng.normal(size=(n, n))
    R = 0.6 * np.corrcoef(B @ B.T) + 0.4 * np.eye(n)
    Q = np.outer(d, d) * R
    rel_gap = []
    for k in range(len(obs)):
        t_lr, _ = _lr_statistic(A, Q, l, k)
        assert _baarda_general(A, Q, l, k) == pytest.approx(t_lr, rel=1e-6)
        rel_gap.append(abs(_block_form(A, Q, l, k) - t_lr) / t_lr)
    assert max(rel_gap) > 0.05  # block form materially wrong under correlation


def test_null_distribution_is_chi2_3_and_global_chi2_r():
    """Normalisation check: with correct covariance, T_b ~ chi2(3), v'Pv ~ chi2(r)."""
    t_vals, g_vals, dofs = [], [], set()
    for seed in range(400):
        stations, obs, _ = simulate_network(NetworkConfig(n_stations=6, n_baselines=12,
                                                          seed=10_000 + seed))
        res = solve_wls(stations, obs)
        tests = [t for t in baseline_tests(res) if t.group_df == 3]
        t_vals.append(tests[seed % len(tests)].group_statistic)  # one per network: iid
        g_vals.append(res.vTPv)
        dofs.add(res.degrees_of_freedom)
    assert len(dofs) == 1
    r = dofs.pop()
    assert stats.kstest(t_vals, stats.chi2(3).cdf).pvalue > 1e-3
    assert np.mean(t_vals) == pytest.approx(3.0, abs=0.4)
    assert stats.kstest(g_vals, stats.chi2(r).cdf).pvalue > 1e-3
