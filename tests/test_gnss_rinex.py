"""RINEX metadata, decompression and overlap checks on the real GSI test files."""

import datetime as dt
import gzip
from pathlib import Path

import numpy as np
import pytest

from gnss_adjust.gnss.rinex import (RinexError, _decompressed_name, check_overlap,
                                    prepare_rinex, read_obs_meta, sha256_file)

DATA = Path(__file__).parent / "data" / "gsi_0759_3040"


def test_fixture_checksums_unchanged():
    for line in (DATA / "SHA256SUMS").read_text().splitlines():
        digest, name = line.split()
        assert sha256_file(DATA / name) == digest, f"{name} was modified"


def test_header_and_epochs_rinex2():
    m = read_obs_meta(DATA / "07590920.05o")
    assert m.version == pytest.approx(2.10)
    assert m.station_id == "0759"
    assert m.receiver_type == "TRIMBLE 5700"
    assert m.antenna_type == "TRM29659.00"
    np.testing.assert_allclose(m.approx_xyz, [-3976219.5082, 3382372.5671, 3652512.9849])
    np.testing.assert_allclose(m.antenna_delta_hen, [0.0, 0.0, 0.0])
    assert m.obs_types["*"] == ["L1", "C1", "L2", "P2"]
    assert m.n_epochs == 120
    assert m.median_interval_s == pytest.approx(30.0, abs=0.01)
    assert m.epoch_first == dt.datetime(2005, 4, 2, 0, 0, 0)
    assert m.systems_in_data == ["G"]
    assert m.warnings == []


def test_overlap_ok_and_failure_modes():
    a = read_obs_meta(DATA / "07590920.05o")
    b = read_obs_meta(DATA / "30400920.05o")
    ov = check_overlap(a, b)
    assert ov.ok and ov.common_epochs == 120 and ov.overlap_s == pytest.approx(3570, abs=1)
    assert not check_overlap(a, b, min_overlap_s=7200).ok
    short = read_obs_meta(DATA / "30400920.05o")
    short.epochs = short.epochs[:10]
    assert not check_overlap(a, short).ok


def test_prepare_never_modifies_source_and_decompresses_gzip(tmp_path):
    src_dir = tmp_path / "archive"
    src_dir.mkdir()
    gz = src_dir / "30400920.05o.gz"
    gz.write_bytes(gzip.compress((DATA / "30400920.05o").read_bytes()))
    before = sha256_file(gz)
    out = prepare_rinex(gz, tmp_path / "work")
    assert out.parent == tmp_path / "work" and out.name == "30400920.05o"
    assert out.read_bytes() == (DATA / "30400920.05o").read_bytes()
    assert sha256_file(gz) == before
    assert sorted(p.name for p in src_dir.iterdir()) == ["30400920.05o.gz"]


@pytest.mark.parametrize("name,expected", [
    ("MADR00ESP_R_20240010000_01D_30S_MO.crx.gz", "MADR00ESP_R_20240010000_01D_30S_MO.rnx"),
    ("madr0010.24d.Z", "madr0010.24o"),
    ("07590920.05o", "07590920.05o"),
    ("brdc0010.24n.gz", "brdc0010.24n"),
])
def test_decompressed_names(name, expected):
    assert _decompressed_name(name) == expected


def test_non_observation_file_rejected():
    with pytest.raises(RinexError):
        read_obs_meta(DATA / "07590920.05n")


def _prepare_worker(args):
    src, work = args
    out = prepare_rinex(src, work)
    return read_obs_meta(out).n_epochs


def test_prepare_is_safe_under_concurrency(tmp_path):
    """Regression: parallel workers decompressing one file must never read a partial copy."""
    import concurrent.futures as cf
    src_dir = tmp_path / "archive"
    src_dir.mkdir()
    # ~15 MB, like a real 24 h file, so that a non-atomic write is observable
    text = (DATA / "30400920.05o").read_text()
    head, body = text.split("END OF HEADER\n", 1)
    big = head + "END OF HEADER\n" + body * 200
    gz = src_dir / "30400920.05o.gz"
    gz.write_bytes(gzip.compress(big.encode(), compresslevel=1))
    work = tmp_path / "work"
    with cf.ProcessPoolExecutor(max_workers=8) as ex:
        n = list(ex.map(_prepare_worker, [(gz, work)] * 32))
    assert n == [120 * 200] * 32
    assert sorted(p.name for p in work.iterdir()) == ["30400920.05o"]  # no leftovers
