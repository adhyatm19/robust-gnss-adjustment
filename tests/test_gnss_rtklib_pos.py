"""RTKLIB .pos parsing and covariance reconstruction (no binary needed)."""

import numpy as np
import pytest

from gnss_adjust.gnss.rtklib_pos import (PosFormatError, check_covariance,
                                         covariance_from_rtklib, parse_pos_file,
                                         signed_square)

BASE = np.array([-3976219.5082, 3382372.5671, 3652512.9849])
HEADER = """% program   : RTKLIB ver.2.4.3
% inp file  : rover.obs
% inp file  : base.obs
% obs start : 2005/04/02 00:00:00.0 GPST (week1316 518400.0s)
% obs end   : 2005/04/02 00:59:30.0 GPST (week1316 521970.0s)
% pos mode  : Static
% ref pos   :  -3976219.5082000    3382372.5671000    3652512.9849000
%
% (x/y/z-ecef=WGS84,Q=1:fix,2:float,3:sbas,4:dgps,5:single,6:ppp,ns=# of satellites)
%  GPST                      x-ecef(m)      y-ecef(m)      z-ecef(m)   Q  ns   sdx(m)   sdy(m)   sdz(m)  sdxy(m)  sdyz(m)  sdzx(m) age(s)  ratio
"""


def _fields(q):
    sq = lambda c: np.sign(c) * np.sqrt(abs(c))  # RTKLIB sqvar
    return [np.sqrt(q[0, 0]), np.sqrt(q[1, 1]), np.sqrt(q[2, 2]),
            sq(q[0, 1]), sq(q[1, 2]), sq(q[2, 0])]


def _write(tmp_path, q, fmt):
    rover = BASE + np.array([-2022.7699041, 468.6279402, -2610.2896407])
    f = _fields(q)
    line = ("2005/04/02 00:00:00.000 " + " ".join(f"{v:.7f}" for v in rover)
            + "   1   5 " + " ".join(fmt.format(v) for v in f) + "  -0.01  673.1\n")
    p = tmp_path / "sol.pos"
    p.write_text(HEADER + line)
    return p


@pytest.fixture
def q_true():
    # strongly correlated, realistic ECEF block (real GSI magnitude, ~1-2 mm)
    sd = np.array([1.3e-3, 1.7e-3, 1.9e-3])
    r = np.array([[1, -0.76, -0.79], [-0.76, 1, 0.70], [-0.79, 0.70, 1]])
    return np.outer(sd, sd) * r


def test_signed_square_inverts_sqvar():
    c = np.array([-2.5e-6, 0.0, 9e-7])
    s = np.sign(c) * np.sqrt(np.abs(c))
    np.testing.assert_allclose(signed_square(s), c, rtol=1e-15)


def test_covariance_field_order(q_true):
    np.testing.assert_allclose(covariance_from_rtklib(_fields(q_true)[:3], _fields(q_true)[3:]),
                               q_true, rtol=1e-12)


def test_parse_high_precision_reconstructs_exactly(tmp_path, q_true):
    pf = parse_pos_file(_write(tmp_path, q_true, "{:.10f}"))
    assert pf.frame == "xyz" and len(pf.records) == 1
    np.testing.assert_allclose(pf.ref_pos, BASE, atol=1e-7)
    rec = pf.records[0]
    assert rec.q_status == 1 and rec.status_name == "fix" and rec.n_sat == 5
    assert rec.resolution_m == pytest.approx(1e-10)
    np.testing.assert_allclose(rec.covariance, q_true, rtol=1e-6, atol=1e-15)
    chk = check_covariance(rec, BASE, n_samples=300)
    assert chk.positive_definite and not chk.rounding_sensitive


def test_stock_rounding_is_detected(tmp_path, q_true):
    rec = parse_pos_file(_write(tmp_path, q_true, "{:8.4f}")).records[0]
    assert rec.resolution_m == pytest.approx(1e-4)
    chk = check_covariance(rec, BASE, n_samples=500)
    assert chk.positive_definite
    assert chk.rounding_sensitive
    assert chk.rounding_max_rel_sigma_change_enu > 0.05


def test_wrong_format_rejected(tmp_path):
    p = tmp_path / "llh.pos"
    p.write_text("%  GPST   latitude(deg) longitude(deg)  height(m)   Q  ns\n"
                 "2005/04/02 00:00:00.000 35.0 139.0 50.0 1 5\n")
    with pytest.raises(PosFormatError):
        parse_pos_file(p)


def test_enu_header_detected(tmp_path):
    p = tmp_path / "enu.pos"
    p.write_text("%  GPST                  e-baseline(m)  n-baseline(m)  u-baseline(m)   Q  ns   "
                 "sde(m)   sdn(m)   sdu(m)  sden(m)  sdnu(m)  sdue(m) age(s)  ratio\n"
                 "2005/04/02 00:57:00.000 953.6738 -3196.1393 4.6482 1 5 0.0007 0.0009 0.0025 "
                 "-0.0002 0.0009 -0.0010 -0.01 673.1\n")
    rec = parse_pos_file(p, frame="enu").records[0]
    np.testing.assert_allclose(rec.xyz, [953.6738, -3196.1393, 4.6482])
    assert rec.covariance[0, 1] == pytest.approx(-4e-8)
    with pytest.raises(PosFormatError):
        parse_pos_file(p, frame="xyz")
