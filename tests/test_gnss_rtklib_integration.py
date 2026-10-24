"""End-to-end RINEX -> RTKLIB -> baseline on real GSI data (needs the pinned build).

Skipped when scripts/build_rtklib.sh has not been run.
"""

from pathlib import Path

import numpy as np
import pytest

from gnss_adjust.gnss import frames
from gnss_adjust.gnss.baseline import process_baseline
from gnss_adjust.gnss.rtklib import AntennaSetting, ProcessingConfig, RTKLIBError, find_rnx2rtkp

DATA = Path(__file__).parent / "data" / "gsi_0759_3040"
BASE = np.array([-3976219.5082, 3382372.5671, 3652512.9849])
ANT = "TRM29659.00"

try:
    find_rnx2rtkp("hp")
    find_rnx2rtkp("stock")
    HAVE_RTKLIB = True
except RTKLIBError:
    HAVE_RTKLIB = False

pytestmark = pytest.mark.skipif(not HAVE_RTKLIB, reason="pinned RTKLIB not built")


def _solve(tmp_path, tag, rover=DATA / "30400920.05o", variant="hp", **kw):
    return process_baseline(
        session_id=tag, rover_obs=rover, base_obs=DATA / "07590920.05o",
        aux_files=[DATA / "07590920.05n", DATA / "30400920.05n"], base_xyz=BASE,
        workdir=tmp_path / tag, cfg=ProcessingConfig(elevation_mask_deg=15),
        variant=variant, **kw)


@pytest.fixture(scope="module")
def module_tmp(tmp_path_factory):
    return tmp_path_factory.mktemp("rtk")


@pytest.fixture(scope="module")
def reference(module_tmp):
    return _solve(module_tmp, "ref")


def test_real_baseline_fixed_and_plausible(reference):
    b = reference
    assert (b.from_station, b.to_station) == ("0759", "3040")
    assert b.status == "fix" and b.ratio > 3.0
    assert b.length_m == pytest.approx(3335.389, abs=0.01)
    assert b.covariance_source == "rtklib_formal_float32"
    assert b.covariance_check.positive_definite
    assert not b.covariance_check.rounding_sensitive
    sig_enu = np.sqrt(np.diag(b.covariance_enu()))
    assert np.all(sig_enu > 1e-4) and np.all(sig_enu < 1e-2)  # formal: 0.1-10 mm
    assert sig_enu[2] > sig_enu[0] and sig_enu[2] > sig_enu[1]  # up is weakest


def test_command_uses_config_only(reference):
    cmd = reference.run.command
    assert cmd[1] == "-k" and cmd[3] == "-o"
    assert not any(c.startswith("-") for c in cmd[5:])  # no overriding CLI options


def test_stock_and_hp_agree_to_print_resolution(module_tmp, reference):
    stock = _solve(module_tmp, "stock", variant="stock")
    np.testing.assert_allclose(stock.vector_ecef, reference.vector_ecef, atol=5e-5 + 1e-9)
    assert stock.status == reference.status and stock.ratio == reference.ratio
    assert stock.covariance_check.rounding_sensitive  # 0.1 mm rounding detected


def _header_copy(tmp_path, h):
    lines = (DATA / "30400920.05o").read_text().splitlines(True)
    out = [(f"{h:14.4f}{0.0:14.4f}{0.0:14.4f}".ljust(60) + "ANTENNA: DELTA H/E/N\n")
           if "ANTENNA: DELTA H/E/N" in ln else ln for ln in lines]
    d = tmp_path / "edited"
    d.mkdir(exist_ok=True)
    (d / "30400920.05o").write_text("".join(out))
    return d / "30400920.05o"


def test_antenna_height_perturbation_is_applied(module_tmp, reference):
    up_rover = frames.enu_rotation_at(reference.rover_xyz)[2]
    up_base = frames.enu_rotation_at(BASE)[2]
    rov = _solve(module_tmp, "rov10", rover_antenna=AntennaSetting("explicit", ANT, (0, 0, 0.10)))
    bas = _solve(module_tmp, "bas10", base_antenna=AntennaSetting("explicit", ANT, (0, 0, 0.10)))
    # marker = ARP - h*up: a too-large rover height lowers the rover marker
    np.testing.assert_allclose(rov.vector_ecef - reference.vector_ecef, -0.10 * up_rover, atol=2e-4)
    np.testing.assert_allclose(bas.vector_ecef - reference.vector_ecef, +0.10 * up_base, atol=2e-4)
    assert rov.rover_antenna.delta_enu == (0, 0, 0.10)


def test_antenna_precedence_rinex_star_overrides_config(module_tmp, reference):
    star = _solve(module_tmp, "star", rover_antenna=AntennaSetting("rinex", ANT, (0, 0, 0.10)))
    np.testing.assert_allclose(star.vector_ecef, reference.vector_ecef, atol=1e-9)


def test_antenna_precedence_header_used_only_with_star(module_tmp, reference):
    edited = _header_copy(module_tmp, 0.10)
    star = _solve(module_tmp, "hdr_star", rover=edited, rover_antenna=AntennaSetting("rinex"))
    expl = _solve(module_tmp, "hdr_expl", rover=edited,
                  rover_antenna=AntennaSetting("explicit", ANT, (0, 0, 0)))
    rov = _solve(module_tmp, "rov10b", rover_antenna=AntennaSetting("explicit", ANT, (0, 0, 0.10)))
    np.testing.assert_allclose(star.vector_ecef, rov.vector_ecef, atol=1e-9)        # header read
    np.testing.assert_allclose(expl.vector_ecef, reference.vector_ecef, atol=1e-9)  # header ignored
    # and antenna_from_header (the default) reproduces the edited header
    default = _solve(module_tmp, "hdr_default", rover=edited)
    np.testing.assert_allclose(default.vector_ecef, rov.vector_ecef, atol=1e-9)
