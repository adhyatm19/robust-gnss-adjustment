"""Covariance construction and validation."""

import numpy as np
import pytest

from gnss_adjust import covariance as cov


def test_baseline_sigma_units():
    # 5 mm constant, 1 ppm at 10 km -> sqrt(5^2 + 10^2) mm = 11.18 mm
    s = cov.baseline_sigma_m(a_mm=5.0, b_ppm=1.0, length_m=10_000.0)
    assert s == pytest.approx(np.sqrt(0.005 ** 2 + 0.010 ** 2))


def test_generated_blocks_symmetric_positive_definite():
    rng = np.random.default_rng(123)
    for _ in range(50):
        sigma = cov.baseline_sigma_m(5.0, 1.0, rng.uniform(100, 20000))
        sigmas = cov.component_sigmas(sigma, (1.0, 1.0, 1.8))
        R = cov.random_correlation(rng, 0.4)
        Q = cov.build_covariance_block(sigmas, R)
        cov.validate_covariance_block(Q)          # raises on failure
        np.testing.assert_allclose(Q, Q.T, atol=1e-18)
        assert np.all(np.linalg.eigvalsh(Q) > 0)


def test_validate_rejects_asymmetric():
    Q = np.eye(3) * 1e-4
    Q[0, 1] = 1e-5
    with pytest.raises(cov.CovarianceError, match="symmetric"):
        cov.validate_covariance_block(Q)


def test_validate_rejects_indefinite():
    Q = np.diag([1e-4, 1e-4, -1e-4])
    with pytest.raises(cov.CovarianceError, match="positive definite"):
        cov.validate_covariance_block(Q)


def test_weight_block_is_inverse():
    rng = np.random.default_rng(5)
    R = cov.random_correlation(rng, 0.3)
    Q = cov.build_covariance_block(np.array([0.01, 0.01, 0.02]), R)
    P = cov.weight_block(Q)
    np.testing.assert_allclose(P @ Q, np.eye(3), atol=1e-10)
