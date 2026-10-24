"""Satellite screening from orbit/clock products and broadcast health.

Why this exists: RTKLIB 2.4.3 b34 with precise ephemerides falls back to the
*broadcast* clock when an SP3 clock is missing (``999999.999999``) and in
that path ignores the broadcast health flag. A satellite that the IGS has
flagged (bad clock) and the control segment has set unhealthy is then used
anyway; on the Madrid data (2024-01-01, G27) its ~100 m pseudorange error
failed RTKLIB's single-point chi-square test at every epoch while it was
visible, silently discarding six hours of observations.

Screening rule (applied per processing day, logged, never silent):
a GPS satellite is excluded if, during the day, its SP3 clock is flagged
missing at any epoch, or any of its broadcast ephemerides has a non-zero
health word.
"""

from __future__ import annotations

from pathlib import Path

SP3_BAD_CLOCK = 999999.0


def sp3_bad_clock_satellites(sp3_path: str | Path, day: tuple[int, int, int] | None = None) -> set[str]:
    """Satellites (e.g. ``"G27"``) with a missing/bad clock in an SP3 file.

    If ``day=(year, month, day)`` is given only epochs on that date count.
    """
    bad: set[str] = set()
    in_day = day is None
    with open(sp3_path, encoding="ascii", errors="replace") as fh:
        for line in fh:
            if line.startswith("*"):
                parts = line[1:].split()
                in_day = day is None or (int(parts[0]), int(parts[1]), int(parts[2])) == day
            elif line.startswith("P") and in_day and len(line) >= 60:
                try:
                    clk = float(line[46:60])
                except ValueError:
                    continue
                if clk >= SP3_BAD_CLOCK:
                    bad.add(line[1:4].replace(" ", "0"))
    return bad


def rinex2_gps_unhealthy(nav_path: str | Path) -> set[str]:
    """GPS satellites with a non-zero SV health word in a RINEX 2 GPS navigation file."""
    unhealthy: set[str] = set()
    with open(nav_path, encoding="ascii", errors="replace") as fh:
        lines = fh.readlines()
    i = next(k for k, ln in enumerate(lines) if "END OF HEADER" in ln) + 1
    while i + 7 < len(lines):
        prn = int(lines[i][0:2])
        # broadcast orbit 6 (8th line of the record), second field: SV health
        health = float(lines[i + 6][22:41].replace("D", "E"))
        if health != 0.0:
            unhealthy.add(f"G{prn:02d}")
        i += 8
    return unhealthy
