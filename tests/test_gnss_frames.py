"""Coordinate-frame and covariance transformations (analytical checks)."""

import numpy as np
import pytest

from gnss_adjust.gnss import frames


@pytest.mark.parametrize("lat_deg,lon_deg,h", [
    (0.0, 0.0, 0.0), (35.5, 139.6, 120.0), (-33.9, 18.4, 25.0),
    (40.43, -4.25, 830.0), (89.9, 10.0, 3000.0), (26.5, 80.2, 126.0),
])
def test_geodetic_round_trip(lat_deg, lon_deg, h):
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    xyz = frames.geodetic_to_ecef(lat, lon, h)
    lat2, lon2, h2 = frames.ecef_to_geodetic(xyz)
    assert lat2 == pytest.approx(lat, abs=1e-11)
    assert lon2 == pytest.approx(lon, abs=1e-11)
    assert h2 == pytest.approx(h, abs=1e-6)


def test_enu_axes_at_equator_greenwich():
    r = frames.enu_rotation(0.0, 0.0)
    np.testing.assert_allclose(r[0], [0, 1, 0], atol=1e-15)   # east = +Y
    np.testing.assert_allclose(r[1], [0, 0, 1], atol=1e-15)   # north = +Z
    np.testing.assert_allclose(r[2], [1, 0, 0], atol=1e-15)   # up = +X


def test_up_vector_is_ellipsoid_normal():
    lat, lon = np.radians(40.43), np.radians(-4.25)
    p0 = frames.geodetic_to_ecef(lat, lon, 0.0)
    p1 = frames.geodetic_to_ecef(lat, lon, 1.0)
    np.testing.assert_allclose(frames.enu_rotation(lat, lon)[2], p1 - p0, atol=1e-9)


def test_rotation_orthonormal_and_round_trip():
    ref = frames.geodetic_to_ecef(np.radians(35.0), np.radians(139.0), 50.0)
    r = frames.enu_rotation_at(ref)
    np.testing.assert_allclose(r @ r.T, np.eye(3), atol=1e-15)
    assert np.linalg.det(r) == pytest.approx(1.0)   # proper rotation, no reflection
    d = np.array([-2022.77, 468.63, -2610.29])
    np.testing.assert_allclose(frames.vector_enu_to_ecef(frames.vector_ecef_to_enu(d, ref), ref),
                               d, atol=1e-9)


def test_covariance_rotation_preserves_invariants():
    rng = np.random.default_rng(1)
    a = rng.normal(size=(3, 3))
    q = a @ a.T * 1e-6
    ref = frames.geodetic_to_ecef(np.radians(-20.0), np.radians(60.0), 0.0)
    qe = frames.covariance_ecef_to_enu(q, ref)
    np.testing.assert_allclose(qe, qe.T, atol=0)
    np.testing.assert_allclose(np.linalg.eigvalsh(qe), np.linalg.eigvalsh(q), rtol=1e-12)
    d = rng.normal(size=3) * 1e-3
    de = frames.vector_ecef_to_enu(d, ref)
    # Mahalanobis form (hence every chi-square test) is frame invariant
    assert de @ np.linalg.solve(qe, de) == pytest.approx(d @ np.linalg.solve(q, d), rel=1e-12)
    np.testing.assert_allclose(frames.covariance_enu_to_ecef(qe, ref), q, atol=1e-20)


def test_pure_vertical_variance_maps_to_up():
    ref = frames.geodetic_to_ecef(np.radians(48.0), np.radians(11.0), 500.0)
    q_enu = np.diag([0.0, 0.0, 4e-6])
    q_ecef = frames.covariance_enu_to_ecef(q_enu, ref)
    up = frames.enu_rotation_at(ref)[2]
    np.testing.assert_allclose(q_ecef, 4e-6 * np.outer(up, up), atol=1e-20)
