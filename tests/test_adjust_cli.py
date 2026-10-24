"""The truth-free ``adjust`` CLI and network report, plus datum invariance."""

import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from gnss_adjust import io as gio
from gnss_adjust import network_report as nr
from gnss_adjust.cli import main
from gnss_adjust.diagnostics import baseline_tests
from gnss_adjust.sd_outlier import detect_specific_direction
from gnss_adjust.simulation import NetworkConfig, OutlierSpec, simulate_network
from gnss_adjust.types import Station
from gnss_adjust.wls import solve_wls

ROOT = Path(__file__).resolve().parents[1]
ECEF_OFFSET = np.array([4849202.0, -360329.0, 4114913.0])  # ~Madrid on the ellipsoid


def _ecef_network(seed=7, outliers=()):
    stations, obs, truth = simulate_network(NetworkConfig(n_stations=6, n_baselines=12, seed=seed,
                                                          outliers=list(outliers)))
    stations = [replace(s, coords=s.coords + ECEF_OFFSET) for s in stations]
    ref = {k: v + ECEF_OFFSET for k, v in truth.true_coordinates.items()}
    return stations, obs, ref, truth


def test_adjust_cli_on_examples(tmp_path):
    out = tmp_path / "adj"
    assert main(["adjust", "--network", str(ROOT / "examples"), "--output", str(out)]) == 0
    for f in ("coordinates.csv", "baselines.csv", "summary.json", "report.md",
              "network_enu.png", "group_statistics.png", "residual_norms.png"):
        assert (out / f).exists(), f
    s = json.loads((out / "summary.json").read_text())
    assert set(s["methods"]) == {"WLS", "DIA-reject", "DIA-inflate", "Huber", "Hampel", "SD-reject"}
    assert s["degrees_of_freedom"] == s["n_observations"] - s["rank_A"]
    assert "not geocentric" in s["local_frame"]


def test_adjust_ecef_network_with_reference_and_outlier(tmp_path):
    stations, obs, ref, truth = _ecef_network(outliers=[OutlierSpec(baseline_id=4,
                                                                    vector=(0.0, 0.03, 0.12))])
    net_dir = tmp_path / "net"
    gio.write_network_csv(stations, obs, net_dir)
    pd.DataFrame([{"station_id": k, "x": v[0], "y": v[1], "z": v[2]} for k, v in ref.items()]
                 ).to_csv(net_dir / "reference.csv", index=False)
    out = tmp_path / "adj"
    assert main(["adjust", "--network", str(net_dir), "--reference", str(net_dir / "reference.csv"),
                 "--output", str(out)]) == 0
    s = json.loads((out / "summary.json").read_text())
    assert s["local_frame"].startswith("ENU")
    assert s["methods"]["DIA-reject"]["flagged_baselines"] == [4]
    assert s["methods"]["SD-reject"]["flagged_baselines"] == [4]
    coords = pd.read_csv(out / "coordinates.csv")
    dia = coords[(coords.method == "DIA-reject") & ~coords.fixed]
    assert dia.diff_3d_mm.max() < 60.0          # back near truth after rejection
    b = pd.read_csv(out / "baselines.csv")
    assert {"v_e_mm", "v_n_mm", "v_u_mm", "sd_dir_u", "local_redundancy"} <= set(b.columns)
    wls4 = b[(b.method == "WLS") & (b.baseline_id == 4)].iloc[0]
    assert wls4.sd_significant and wls4.sd_dir_elevation_deg > 30  # bias mostly vertical-ish
    # redundancy numbers sum to the total redundancy
    wls = b[b.method == "WLS"]
    assert wls.local_redundancy.sum() == pytest.approx(s["degrees_of_freedom"], rel=1e-9)


@pytest.mark.parametrize("seed", [2, 9])
def test_datum_invariance_of_residuals_and_tests(seed):
    """Changing the fixed station (a minimal constraint) must not change v, Qv or any test."""
    stations, obs, _, _ = _ecef_network(seed=seed)
    alt = [Station(s.station_id, s.coords, is_fixed=(s.station_id == stations[3].station_id))
           for s in stations]
    r1, r2 = solve_wls(stations, obs), solve_wls(alt, obs)
    # ECEF-magnitude cancellation leaves ~1e-8 m differences (double precision)
    np.testing.assert_allclose(r1.residuals, r2.residuals, atol=5e-8)
    np.testing.assert_allclose(r1.residual_covariance, r2.residual_covariance, atol=1e-12)
    assert r1.vTPv == pytest.approx(r2.vTPv, rel=1e-7)
    t1, t2 = baseline_tests(r1), baseline_tests(r2)
    for a, b in zip(t1, t2):
        # ~1e-6 relative: absolute ECEF parameterisation, negligible for test decisions
        assert a.group_statistic == pytest.approx(b.group_statistic, rel=1e-5)
    s1, s2 = detect_specific_direction(r1), detect_specific_direction(r2)
    for a, b in zip(s1.tests, s2.tests):
        if a.testable:
            np.testing.assert_allclose(a.outlier_vector, b.outlier_vector, atol=5e-8)


def test_wls_precision_at_ecef_magnitude():
    """Absolute ECEF coordinates (~6.4e6 m) lose no meaningful precision in WLS."""
    stations, obs, _, _ = _ecef_network(seed=3)
    local = [replace(s, coords=s.coords - ECEF_OFFSET) for s in stations]
    r_ecef, r_loc = solve_wls(stations, obs), solve_wls(local, obs)
    for sid in r_loc.coordinates:
        np.testing.assert_allclose(r_ecef.coordinates[sid] - ECEF_OFFSET, r_loc.coordinates[sid],
                                   atol=1e-7)
    np.testing.assert_allclose(r_ecef.residuals, r_loc.residuals, atol=1e-7)
