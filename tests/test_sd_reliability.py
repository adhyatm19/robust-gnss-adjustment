"""Reliability-matrix construction for the specific-direction method."""

import numpy as np

from gnss_adjust.sd_outlier import baseline_reliability, reliability_matrix
from gnss_adjust.wls import solve_wls


def test_sd_reliability_is_symmetric_psd(simulated_clean):
    stations, observations, _ = simulated_clean
    result = solve_wls(stations, observations)
    p_bar = reliability_matrix(result)

    np.testing.assert_allclose(p_bar, p_bar.T, atol=1e-8)
    eigenvalues = np.linalg.eigvalsh(p_bar)
    assert eigenvalues.min() > -1e-8 * max(eigenvalues.max(), 1.0)


def test_sd_block_and_score_shapes(simulated_clean):
    stations, observations, _ = simulated_clean
    result = solve_wls(stations, observations)

    for index in range(len(observations)):
        p_ii, g_i = baseline_reliability(result, index)
        assert p_ii.shape == (3, 3)
        assert g_i.shape == (3,)
        np.testing.assert_allclose(p_ii, p_ii.T, atol=1e-8)


def test_score_matches_paper_reliability_product(simulated_clean):
    """Verify g=P_bar*l while computing it stably as -P*v."""
    from gnss_adjust import network as net
    from gnss_adjust.sd_outlier import reliability_score

    stations, observations, _ = simulated_clean
    result = solve_wls(stations, observations)
    system = net.build_design_system(stations, observations)
    expected = reliability_matrix(result) @ system.l_reduced

    np.testing.assert_allclose(reliability_score(result), expected, rtol=2e-8, atol=2e-8)

