"""Download a real GNSS dataset described by a YAML manifest (see configs/real_data/).

Files go to data/raw/<name>/{obs,nav,orbit,products} (gitignored). Existing
files are never overwritten. Every file is logged with its source URL, size
and SHA-256 in data/raw/<name>/download_log.json, and a copy of the log is
written to results/real_network/<name>_download_log.json for the record.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def gps_week(year: int, doy: int) -> int:
    d = dt.date(year, 1, 1) + dt.timedelta(days=doy - 1)
    return (d - dt.date(1980, 1, 6)).days // 7


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path, retries: int = 4) -> dict:
    """Download url to dest (unless present); retry with exponential backoff."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        tmp = dest.with_suffix(dest.suffix + ".part")
        for attempt in range(retries + 1):
            try:
                with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as fh:
                    while chunk := r.read(1 << 20):
                        fh.write(chunk)
                tmp.rename(dest)
                break
            except (urllib.error.URLError, ConnectionError, TimeoutError) as exc:
                if attempt == retries:
                    raise RuntimeError(f"failed to download {url}: {exc}") from exc
                time.sleep(2 ** (attempt + 1))
    return {"url": url, "file": str(dest.relative_to(ROOT)), "bytes": dest.stat().st_size,
            "sha256": sha256(dest)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("manifest", type=Path)
    args = ap.parse_args()
    m = yaml.safe_load(args.manifest.read_text())
    base = m["archive"].rstrip("/")
    out = ROOT / "data" / "raw" / m["name"]
    log: list[dict] = []

    for doy in m["days"]:
        year, yy = m["year"], m["year"] % 100
        for sta in m["stations"]:
            name = m["obs_template"].format(long=m["station_long_names"][sta], year=year, doy=doy)
            last_err = None
            for d in m["obs_dirs"]:
                try:
                    log.append(fetch(f"{base}/{d.format(year=year, doy=doy)}/{name}",
                                     out / "obs" / name) | {"kind": "obs", "station": sta, "doy": doy})
                    break
                except RuntimeError as exc:
                    last_err = exc
            else:
                raise RuntimeError(f"{sta} day {doy}: not found in any obs dir ({last_err})")
        nav = m["nav_template"].format(doy=doy, yy=yy)
        log.append(fetch(f"{base}/{m['nav_dir'].format(year=year, doy=doy)}/{nav}", out / "nav" / nav)
                   | {"kind": "nav", "doy": doy})
    for year, doy in m["orbit_days"]:
        orb = m["orbit_template"].format(year=year, doy=doy)
        d = m["orbit_dir"].format(gpsweek=gps_week(year, doy))
        log.append(fetch(f"{base}/{d}/{orb}", out / "orbit" / orb) | {"kind": "orbit", "doy": doy,
                                                                       "year": year})
    for key in ("antex", "reference_sinex"):
        url = m[key]["url"]
        log.append(fetch(f"{base}/{url}", out / "products" / Path(url).name) | {"kind": key})

    record = {"dataset": m["name"], "manifest": str(args.manifest), "retrieved_utc":
              dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
              "total_bytes": sum(e["bytes"] for e in log), "files": log}
    (out / "download_log.json").write_text(json.dumps(record, indent=2))
    res = ROOT / "results" / "real_network"
    res.mkdir(parents=True, exist_ok=True)
    (res / f"{m['name']}_download_log.json").write_text(json.dumps(record, indent=2))
    print(f"{len(log)} files, {record['total_bytes'] / 1e6:.1f} MB in {out}")


if __name__ == "__main__":
    sys.exit(main())
