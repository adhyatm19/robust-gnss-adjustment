"""Coordinate-frame and covariance transformations (ECEF <-> geodetic <-> local ENU).

Conventions
-----------
* ECEF Cartesian coordinates in metres on the WGS84 ellipsoid (GRS80 differs
  only in the flattening at the 1e-11 level; irrelevant for rotations).
* Geodetic latitude/longitude in radians, ellipsoidal height in metres.
* The ENU rotation matrix ``R`` has the local East, North and Up unit vectors
  (expressed in ECEF) as its **rows**, so a vector transforms as
  ``d_enu = R @ d_ecef`` and a covariance as ``Q_enu = R @ Q_ecef @ R.T``.
  ``R`` is orthonormal, hence the inverse transformation uses ``R.T``.
* A baseline vector and its covariance are rotated with the ENU frame of a
  chosen reference point (normally the from-station). Rotation preserves
  the quadratic form ``d' Q^{-1} d`` exactly, so statistical tests are
  frame-invariant; only their component-wise interpretation changes.
"""

from __future__ import annotations

import numpy as np

WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)


def geodetic_to_ecef(lat: float, lon: float, h: float) -> np.ndarray:
    """Geodetic (rad, rad, m) to ECEF (m)."""
    n = WGS84_A / np.sqrt(1.0 - WGS84_E2 * np.sin(lat) ** 2)
    return np.array([
        (n + h) * np.cos(lat) * np.cos(lon),
        (n + h) * np.cos(lat) * np.sin(lon),
        (n * (1.0 - WGS84_E2) + h) * np.sin(lat),
    ])


def ecef_to_geodetic(xyz: np.ndarray, tol: float = 1e-12, max_iter: int = 20) -> tuple[float, float, float]:
    """ECEF (m) to geodetic (rad, rad, m) by fixed-point iteration on latitude."""
    x, y, z = np.asarray(xyz, dtype=float).reshape(3)
    p = np.hypot(x, y)
    if p < 1e-9 and abs(z) < 1e-9:
        raise ValueError("geodetic coordinates undefined at the geocentre")
    lon = float(np.arctan2(y, x))
    lat = float(np.arctan2(z, p * (1.0 - WGS84_E2)))
    for _ in range(max_iter):
        n = WGS84_A / np.sqrt(1.0 - WGS84_E2 * np.sin(lat) ** 2)
        h = p / np.cos(lat) - n if abs(lat) < np.pi / 4 else z / np.sin(lat) - n * (1.0 - WGS84_E2)
        new_lat = float(np.arctan2(z, p * (1.0 - WGS84_E2 * n / (n + h))))
        if abs(new_lat - lat) < tol:
            lat = new_lat
            break
        lat = new_lat
    n = WGS84_A / np.sqrt(1.0 - WGS84_E2 * np.sin(lat) ** 2)
    h = p / np.cos(lat) - n if abs(lat) < np.pi / 4 else z / np.sin(lat) - n * (1.0 - WGS84_E2)
    return lat, lon, float(h)


def enu_rotation(lat: float, lon: float) -> np.ndarray:
    """Rotation matrix whose rows are the local East, North, Up unit vectors."""
    sl, cl = np.sin(lon), np.cos(lon)
    sp, cp = np.sin(lat), np.cos(lat)
    return np.array([
        [-sl, cl, 0.0],
        [-sp * cl, -sp * sl, cp],
        [cp * cl, cp * sl, sp],
    ])


def enu_rotation_at(xyz: np.ndarray) -> np.ndarray:
    """ENU rotation matrix at an ECEF point."""
    lat, lon, _ = ecef_to_geodetic(xyz)
    return enu_rotation(lat, lon)


def vector_ecef_to_enu(d_ecef: np.ndarray, ref_xyz: np.ndarray) -> np.ndarray:
    return enu_rotation_at(ref_xyz) @ np.asarray(d_ecef, dtype=float).reshape(3)


def vector_enu_to_ecef(d_enu: np.ndarray, ref_xyz: np.ndarray) -> np.ndarray:
    return enu_rotation_at(ref_xyz).T @ np.asarray(d_enu, dtype=float).reshape(3)


def covariance_ecef_to_enu(q_ecef: np.ndarray, ref_xyz: np.ndarray) -> np.ndarray:
    r = enu_rotation_at(ref_xyz)
    q = r @ np.asarray(q_ecef, dtype=float) @ r.T
    return 0.5 * (q + q.T)


def covariance_enu_to_ecef(q_enu: np.ndarray, ref_xyz: np.ndarray) -> np.ndarray:
    r = enu_rotation_at(ref_xyz)
    q = r.T @ np.asarray(q_enu, dtype=float) @ r
    return 0.5 * (q + q.T)


def looks_geocentric(xyz: np.ndarray) -> bool:
    """True if a point is plausibly ECEF on/near the Earth surface (|r| within 6.30-6.42 Mm).

    Used to decide whether an ENU rotation is meaningful; synthetic networks
    in a local Cartesian box are reported in their own axes instead.
    """
    r = float(np.linalg.norm(np.asarray(xyz, dtype=float)))
    return 6.30e6 < r < 6.42e6


def local_rotation(xyz: np.ndarray) -> np.ndarray:
    """ENU rotation for geocentric points, identity for local Cartesian networks."""
    return enu_rotation_at(xyz) if looks_geocentric(xyz) else np.eye(3)
