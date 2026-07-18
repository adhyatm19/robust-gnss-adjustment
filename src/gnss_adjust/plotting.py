"""Publication-quality matplotlib figures (no seaborn).

Every function returns a ``matplotlib.figure.Figure`` and optionally saves a
PNG. All axes are labelled with SI units.
"""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .types import AdjustmentResult, BaselineObservation, BaselineTest, Station

logger = logging.getLogger(__name__)

METHOD_COLORS = {
    "WLS": "#4053d3",
    "DIA-reject": "#00b25d",
    "DIA-inflate": "#00beff",
    "Huber": "#ddb310",
    "Hampel": "#b51d14",
}


def _save(fig: plt.Figure, path: str | Path | None) -> plt.Figure:
    fig.tight_layout()
    if path is not None:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(p, dpi=150)
        logger.info("wrote %s", p)
    return fig


def plot_network(
    stations: list[Station],
    observations: list[BaselineObservation],
    true_outliers: set[int] | None = None,
    detected: set[int] | None = None,
    downweighted: set[int] | None = None,
    path: str | Path | None = None,
    title: str = "GNSS baseline network (XY projection)",
) -> plt.Figure:
    """XY projection with fixed/unknown stations and baseline classes."""
    true_outliers = true_outliers or set()
    detected = detected or set()
    downweighted = downweighted or set()
    coords = {s.station_id: s.coords for s in stations}
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    for o in observations:
        p, q = coords[o.from_station], coords[o.to_station]
        is_out = o.baseline_id in true_outliers
        is_det = o.baseline_id in detected
        is_dw = o.baseline_id in downweighted
        color, lw, ls = "0.6", 1.0, "-"
        if is_out:
            color, lw = METHOD_COLORS["Hampel"], 2.2
        if is_det:
            ls = "--"
            lw = max(lw, 2.2)
            if not is_out:
                color = METHOD_COLORS["DIA-reject"]
        elif is_dw and not is_out:
            color, ls, lw = METHOD_COLORS["Huber"], ":", 2.0
        ax.plot([p[0], q[0]], [p[1], q[1]], color=color, lw=lw, ls=ls, zorder=1)
        mid = 0.5 * (p + q)
        ax.annotate(str(o.baseline_id), (mid[0], mid[1]), fontsize=7, color="0.35")
    for s in stations:
        if s.is_fixed:
            ax.plot(s.coords[0], s.coords[1], "s", color="k", ms=10, zorder=3)
        else:
            ax.plot(s.coords[0], s.coords[1], "o", color=METHOD_COLORS["WLS"], ms=7, zorder=3)
        ax.annotate(s.station_id, (s.coords[0], s.coords[1]),
                    textcoords="offset points", xytext=(6, 6), fontsize=8)
    handles = [
        plt.Line2D([], [], marker="s", color="k", ls="", label="fixed station"),
        plt.Line2D([], [], marker="o", color=METHOD_COLORS["WLS"], ls="", label="unknown station"),
        plt.Line2D([], [], color="0.6", label="clean baseline"),
        plt.Line2D([], [], color=METHOD_COLORS["Hampel"], lw=2.2, label="true outlier"),
        plt.Line2D([], [], color=METHOD_COLORS["DIA-reject"], lw=2.2, ls="--", label="detected"),
        plt.Line2D([], [], color=METHOD_COLORS["Huber"], lw=2, ls=":", label="down-weighted"),
    ]
    ax.legend(handles=handles, loc="best", fontsize=8)
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="datalim")
    return _save(fig, path)


def plot_true_vs_estimated(
    true_coords: dict[str, np.ndarray],
    estimates: dict[str, dict[str, np.ndarray]],   # method -> {station: coords}
    fixed_station: str,
    path: str | Path | None = None,
    scale: float = 200.0,
) -> plt.Figure:
    """True positions with per-method error vectors (exaggerated by ``scale``)."""
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    for sid, xyz in true_coords.items():
        marker = "s" if sid == fixed_station else "o"
        ax.plot(xyz[0], xyz[1], marker, color="k", ms=7 if sid == fixed_station else 5)
        ax.annotate(sid, (xyz[0], xyz[1]), textcoords="offset points", xytext=(5, 5), fontsize=8)
    for method, coords in estimates.items():
        color = METHOD_COLORS.get(method, "0.4")
        for sid, est in coords.items():
            if sid == fixed_station or sid not in true_coords:
                continue
            d = (est - true_coords[sid]) * scale
            ax.annotate("", xy=(true_coords[sid][0] + d[0], true_coords[sid][1] + d[1]),
                        xytext=(true_coords[sid][0], true_coords[sid][1]),
                        arrowprops=dict(arrowstyle="->", color=color, lw=1.5))
        ax.plot([], [], color=color, label=method)
    ax.legend(fontsize=8)
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_title(f"True vs estimated positions (error vectors x{scale:g})")
    ax.set_aspect("equal", adjustable="datalim")
    return _save(fig, path)


def plot_coordinate_errors(
    errors_by_method: dict[str, dict[str, float]],   # method -> {station: error m}
    path: str | Path | None = None,
) -> plt.Figure:
    """Grouped bar chart of per-station 3D position error by method."""
    methods = list(errors_by_method)
    stations = sorted({s for e in errors_by_method.values() for s in e})
    x = np.arange(len(stations))
    width = 0.8 / max(len(methods), 1)
    fig, ax = plt.subplots(figsize=(9, 4.5))
    for k, m in enumerate(methods):
        vals = [errors_by_method[m].get(s, np.nan) * 1000 for s in stations]
        ax.bar(x + (k - (len(methods) - 1) / 2) * width, vals, width,
               label=m, color=METHOD_COLORS.get(m, None))
    ax.set_xticks(x, stations)
    ax.set_ylabel("3D position error [mm]")
    ax.set_xlabel("station")
    ax.set_title("Station position error by method")
    ax.legend(fontsize=8)
    return _save(fig, path)


def plot_residual_norms(
    tests_by_method: dict[str, list[BaselineTest]],
    true_outliers: set[int] | None = None,
    path: str | Path | None = None,
) -> plt.Figure:
    """Baseline residual norm per baseline, one line per method."""
    fig, ax = plt.subplots(figsize=(9, 4.5))
    for method, tests in tests_by_method.items():
        bids = [t.baseline_id for t in tests]
        norms = [t.residual_norm * 1000 for t in tests]
        ax.plot(bids, norms, "o-", ms=4, label=method,
                color=METHOD_COLORS.get(method, None))
    for b in (true_outliers or set()):
        ax.axvline(b, color="0.85", zorder=0)
    ax.set_xlabel("baseline id")
    ax.set_ylabel("residual norm [mm]")
    ax.set_title("Baseline residual norms (grey lines: true outliers)")
    ax.legend(fontsize=8)
    return _save(fig, path)


def plot_standardized_residuals(
    tests: list[BaselineTest],
    path: str | Path | None = None,
    title: str = "Component standardized residuals (w-tests)",
) -> plt.Figure:
    """|w| statistics for each baseline component with the critical line."""
    fig, ax = plt.subplots(figsize=(9, 4.5))
    comp_names = ["X", "Y", "Z"]
    markers = ["o", "^", "s"]
    for c in range(3):
        bids = [t.baseline_id for t in tests]
        vals = [abs(t.standardized_components[c]) for t in tests]
        ax.plot(bids, vals, markers[c], ms=5, label=f"{comp_names[c]} component")
    if tests:
        ax.axhline(tests[0].component_critical, color="r", ls="--",
                   label=f"critical |w| = {tests[0].component_critical:.2f}")
    ax.set_xlabel("baseline id")
    ax.set_ylabel("|w| statistic [-]")
    ax.set_title(title)
    ax.legend(fontsize=8)
    return _save(fig, path)


def plot_group_statistics(
    tests: list[BaselineTest],
    path: str | Path | None = None,
    title: str = "3D baseline group test statistics",
) -> plt.Figure:
    """Per-baseline chi-square group statistic against its critical value."""
    fig, ax = plt.subplots(figsize=(9, 4.5))
    bids = [t.baseline_id for t in tests]
    stats_ = [t.group_statistic for t in tests]
    colors = ["#b51d14" if t.flagged_group else "#4053d3" for t in tests]
    ax.bar(bids, stats_, color=colors)
    if tests:
        ax.axhline(tests[0].group_critical, color="r", ls="--",
                   label=f"chi-square critical = {tests[0].group_critical:.1f}")
    ax.set_yscale("log")
    ax.set_xlabel("baseline id")
    ax.set_ylabel("group statistic $T_b$ [-] (log scale)")
    ax.set_title(title)
    ax.legend(fontsize=8)
    return _save(fig, path)


def plot_robust_weights(
    weights_by_method: dict[str, dict[int, float]],
    true_outliers: set[int] | None = None,
    path: str | Path | None = None,
) -> plt.Figure:
    """Final robust weight multiplier per baseline."""
    fig, ax = plt.subplots(figsize=(9, 4.5))
    for method, w in weights_by_method.items():
        bids = sorted(w)
        ax.plot(bids, [w[b] for b in bids], "o-", ms=4, label=method,
                color=METHOD_COLORS.get(method, None))
    for b in (true_outliers or set()):
        ax.axvline(b, color="0.85", zorder=0)
    ax.set_ylim(-0.05, 1.1)
    ax.set_xlabel("baseline id")
    ax.set_ylabel("robust weight multiplier [-]")
    ax.set_title("Final robust weights (grey lines: true outliers)")
    ax.legend(fontsize=8)
    return _save(fig, path)


def plot_metric_vs_x(
    x_values: list[float],
    series: dict[str, list[float]],
    xlabel: str,
    ylabel: str,
    title: str,
    path: str | Path | None = None,
    logy: bool = False,
) -> plt.Figure:
    """Generic curve plot used for RMSE / F1 vs magnitude or contamination."""
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for name, ys in series.items():
        ax.plot(x_values, ys, "o-", label=name, color=METHOD_COLORS.get(name, None))
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    return _save(fig, path)


def plot_monte_carlo_box(
    errors_by_method: dict[str, list[float]],
    path: str | Path | None = None,
    title: str = "Monte Carlo distribution of mean station position error",
) -> plt.Figure:
    """Box plot of per-run mean position error for each method."""
    methods = list(errors_by_method)
    data = [np.asarray(errors_by_method[m]) * 1000 for m in methods]
    fig, ax = plt.subplots(figsize=(8, 5))
    bp = ax.boxplot(data, tick_labels=methods, showfliers=True, patch_artist=True)
    for patch, m in zip(bp["boxes"], methods):
        patch.set_facecolor(METHOD_COLORS.get(m, "0.7"))
        patch.set_alpha(0.6)
    ax.set_ylabel("mean station position error [mm]")
    ax.set_title(title)
    ax.grid(axis="y", alpha=0.3)
    return _save(fig, path)
