"""SINEX reading, satellite screening, repeatability calibration and network assembly."""

import numpy as np
import pytest
from scipy import stats

from gnss_adjust.gnss.calibration import repeatability_calibration
from gnss_adjust.gnss.network_build import (NetworkDesignError, build_network,
                                            check_session_independence)
from gnss_adjust.gnss.satscreen import rinex2_gps_unhealthy, sp3_bad_clock_satellites
from gnss_adjust.gnss.sinex import read_sinex, sinex_epoch
from gnss_adjust.wls import solve_wls

SINEX = """%=SNX 2.01 TST 24:001:00000 TST 23:365:00000 24:007:00000 P 00003 1 S
+SITE/ANTENNA
 MADR  A    1 P 23:365:00000 24:007:00000 AOAD/M_T        NONE 389
-SITE/ANTENNA
+SITE/ECCENTRICITY
 MADR  A    1 P 23:365:00000 24:007:00000 UNE   0.0254   0.0000   0.0000
-SITE/ECCENTRICITY
+SOLUTION/ESTIMATE
*INDEX TYPE__ CODE PT SOLN _REF_EPOCH__ UNIT S __ESTIMATED VALUE____ _STD_DEV___
     1 STAX   MADR  A    1 24:003:43185 m    2  4.84920218874871E+06 1.74238E-04
     2 STAY   MADR  A    1 24:003:43185 m    2 -3.60328594652260E+05 5.22905E-05
     3 STAZ   MADR  A    1 24:003:43185 m    2  4.11491343122203E+06 1.40638E-04
-SOLUTION/ESTIMATE
"""


def test_sinex_reader(tmp_path):
    p = tmp_path / "t.snx"
    p.write_text(SINEX)
    sol = read_sinex(p)
    s = sol.stations["MADR"]
    np.testing.assert_allclose(s.xyz, [4849202.18874871, -360328.59465226, 4114913.43122203])
    np.testing.assert_allclose(s.sigma_xyz, [1.74238e-4, 5.22905e-5, 1.40638e-4])
    assert s.antenna_type == "AOAD/M_T        NONE"
    np.testing.assert_allclose(s.eccentricity_une, [0.0254, 0, 0])
    assert s.ref_epoch == sinex_epoch("24:003:43185")
    assert read_sinex(p, ["CEBR"]).stations == {}


def test_satellite_screening(tmp_path):
    sp3 = tmp_path / "x.sp3"
    sp3.write_text(
        "*  2024  1  1  0  0  0.00000000\n"
        "PG01  -1901.698416 -15399.006346 -21780.734592    254.900827\n"
        "PG27  -1901.698416 -15399.006346 -21780.734592 999999.999999\n"
        "*  2024  1  2  0  0  0.00000000\n"
        "PG08  -1901.698416 -15399.006346 -21780.734592 999999.999999\n")
    assert sp3_bad_clock_satellites(sp3, (2024, 1, 1)) == {"G27"}
    assert sp3_bad_clock_satellites(sp3) == {"G27", "G08"}
    from pathlib import Path
    gsi_nav = Path(__file__).parent / "data" / "gsi_0759_3040" / "07590920.05n"
    assert rinex2_gps_unhealthy(gsi_nav) == set()          # real file: all healthy
    lines = gsi_nav.read_text().splitlines(True)
    h = next(i for i, ln in enumerate(lines) if "END OF HEADER" in ln) + 1
    rec6 = lines[h + 6]
    lines[h + 6] = rec6[:22] + f"{63.0:19.12E}".replace("E", "D") + rec6[41:]
    bad = tmp_path / "bad.05n"
    bad.write_text("".join(lines))
    assert rinex2_gps_unhealthy(bad) == {f"G{int(lines[h][:2]):02d}"}


def _records(rng, k2=1.0, n_sessions=4, pairs=(("A", "B"), ("A", "C"), ("B", "C")), bias=None):
    truth = {"A": np.array([4849202.0, -360328.0, 4114913.0])}
    truth["B"] = truth["A"] + np.array([1000.0, 2000.0, -500.0])
    truth["C"] = truth["A"] + np.array([-3000.0, 1500.0, 800.0])
    out = []
    for s in range(n_sessions):
        for a, b in pairs:
            q = np.diag([4e-6, 1e-6, 9e-6])
            v = truth[b] - truth[a] + rng.multivariate_normal(np.zeros(3), k2 * q)
            if bias and (f"S{s}", a, b) in bias:
                v = v + bias[(f"S{s}", a, b)]
            out.append({"session_id": f"S{s}", "from_station": a, "to_station": b,
                        "vector_ecef": v, "covariance_ecef": q, "vector_enu": v,
                        "covariance_enu": q, "status": "float", "fix_fraction": 0.0,
                        "length_m": float(np.linalg.norm(v)),
                        "covariance_source": "test", "covariance_check": {"rounding_sensitive": False}})
    return out, truth


def test_repeatability_calibration_recovers_k2():
    rng = np.random.default_rng(3)
    est = []
    for _ in range(200):
        recs, _ = _records(rng, k2=4.0)
        est.append(repeatability_calibration(recs).k2)
    assert np.mean(est) == pytest.approx(4.0, rel=0.05)      # unbiased
    recs, _ = _records(np.random.default_rng(9), k2=4.0, n_sessions=40)
    r = repeatability_calibration(recs)
    assert r.k2_ci95[0] < 4.0 < r.k2_ci95[1]
    assert r.dof == 3 * 3 * 39 and r.n_pairs == 3


def test_repeatability_requires_repeats():
    recs, _ = _records(np.random.default_rng(0), n_sessions=1)
    with pytest.raises(ValueError):
        repeatability_calibration(recs)


def test_trivial_baselines_rejected_and_allowed():
    edges = [("S1", "A", "B", 1), ("S1", "B", "C", 2), ("S1", "A", "C", 3)]
    with pytest.raises(NetworkDesignError, match="trivial"):
        check_session_independence(edges)
    assert check_session_independence(edges, allow_trivial=True) == {"S1": [1, 2, 3]}
    # the same three pairs from three sessions are algebraically independent
    ok = [("S1", "A", "B", 1), ("S2", "B", "C", 2), ("S3", "A", "C", 3)]
    assert check_session_independence(ok) == {"S1": [1], "S2": [2], "S3": [3]}


def test_build_network_loops_and_redundancy():
    rng = np.random.default_rng(1)
    bias = {("S1", "B", "C"): np.array([0.0, 0.0, 0.05])}
    recs, truth = _records(rng, k2=1e-12, n_sessions=3, bias=bias)  # ~noiseless
    sel = [("S0", "A", "B"), ("S0", "A", "C"), ("S1", "B", "C"), ("S1", "A", "B"),
           ("S2", "A", "C"), ("S2", "B", "C")]
    approx = {k: v + 1.0 for k, v in truth.items()}
    d = build_network(recs, sel, "A", approx, truth["A"])
    assert d.redundancy == 3 * 6 - 3 * 2
    assert len(d.loops) == 6 - 2 and not any(lp.trivial for lp in d.loops)
    for lp in d.loops:
        has_bias = 3 in lp.baselines        # baseline id 3 = (S1, B, C)
        assert (np.linalg.norm(lp.misclosure) > 0.04) == has_bias
    res = solve_wls(d.stations, d.observations)
    assert res.degrees_of_freedom == d.redundancy
    with pytest.raises(NetworkDesignError):
        build_network(recs, [("S0", "A", "B"), ("S0", "B", "C"), ("S0", "A", "C")], "A",
                      approx, truth["A"])


def test_loop_chi2_distribution_under_h0():
    """Cross-session loop statistics follow chi2(3) when the block model holds."""
    vals = []
    for seed in range(150):
        recs, truth = _records(np.random.default_rng(100 + seed), n_sessions=3)
        sel = [("S0", "A", "B"), ("S1", "B", "C"), ("S2", "A", "C")]
        d = build_network(recs, sel, "A", truth, truth["A"])
        vals += [lp.chi2 for lp in d.loops]
    assert stats.kstest(vals, stats.chi2(3).cdf).pvalue > 1e-3
