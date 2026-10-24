"""Research figures for the real-data validation (results/real_network/figures/).

Reads only files written by the real-data scripts; computes nothing new except
plotting summaries. Method colours are a validated categorical palette in a
fixed order (WLS, DIA-reject, DIA-inflate, Huber, Hampel, SD-reject); every
series also has its own marker/line style so identity never relies on colour.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results" / "real_network"
OUT = RES / "figures"
METHODS = ["WLS", "DIA-reject", "DIA-inflate", "Huber", "Hampel", "SD-reject"]
COLOR = dict(zip(METHODS, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]))
STYLE = {"WLS": ("o", "-"), "DIA-reject": ("s", "-"), "DIA-inflate": ("^", "-"),
         "Huber": ("D", "-"), "Hampel": ("v", "-"), "SD-reject": ("x", "--")}
INK, MUTED = "#0b0b0b", "#52514e"
plt.rcParams.update({"axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.grid": True, "grid.color": "#e4e3df",
                     "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
                     "font.size": 9, "lines.linewidth": 2, "lines.markersize": 6,
                     "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb"})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / name, dpi=150)
    plt.close(fig)


def fig_repeatability():
    d = pd.read_csv(RES / "processing/baselines_all.csv")
    d = d[d.strategy == "iflc_float"].copy()
    d["pair"] = d["from"] + "→" + d["to"]
    d["session"] = d.session.astype(str)
    days = sorted(d.session.unique())
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=False)
    for ax, c, lab in zip(axes, "enu", ("East", "North", "Up")):
        for i, (pair, g) in enumerate(d.groupby("pair")):
            dev = g.set_index("session")[f"diff_ref_{c}_mm"]
            dev = dev - dev.mean()
            ax.plot([days.index(s) + 1 for s in dev.index], dev.values, marker="os^Dvx"[i % 6],
                    lw=1.2, ms=5, label=pair, color=list(COLOR.values())[i % 6])
        ax.axvspan(0.5, 3.5, color="#e4e3df", alpha=0.5, lw=0)
        ax.set_xticks(range(1, len(days) + 1), [s[-3:] for s in days])
        ax.set_xlabel("day of year 2024 (shaded: calibration days)")
        ax.set_ylabel(f"{lab} deviation from pair mean [mm]")
        ax.set_title(f"Day-to-day repeatability, {lab}")
    axes[0].legend(fontsize=7, title="baseline", title_fontsize=7)
    fig.suptitle("Real Madrid baselines (RTKLIB iono-free float, 24 h sessions)", color=INK)
    save(fig, "repeatability_by_day.png")


def fig_formal_vs_empirical():
    cal = json.loads((RES / "calibration.json").read_text())
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    for ax, strat, title in zip(axes, ("iflc_float", "l12_fixed_short"),
                                ("iono-free float (all pairs)", "L1/L2 fixed (short pairs)")):
        rows = cal[strat]["all_days"]["per_pair"]
        labels = [f"{r['from']}→{r['to']}" for r in rows]
        x = np.arange(len(rows))
        for j, (comp, mk) in enumerate(zip("ENU", "os^")):
            emp = [r["empirical_sd_enu_mm"][j] for r in rows]
            frm = [r["formal_sd_enu_mm"][j] for r in rows]
            ax.scatter(x + 0.22 * (j - 1), np.array(emp) / np.array(frm), marker=mk, s=40,
                       color=["#2a78d6", "#eb6834", "#1baf7a"][j], label=comp, zorder=3)
        ax.axhline(1.0, color=MUTED, ls="--", lw=1)
        ax.set_yscale("log")
        ax.set_xticks(x, labels, rotation=30, fontsize=8)
        ax.set_ylabel("empirical σ / formal σ (6 days) [-]")
        ax.set_title(title)
        ax.legend(title="component", fontsize=8)
    fig.suptitle("Calibration of RTKLIB formal baseline uncertainty against day-to-day repeatability",
                 color=INK)
    save(fig, "formal_vs_empirical_sigma.png")


def fig_loops():
    nets = [("quiet_formal", "quiet days, cross-session loops"),
            ("main_formal", "evaluation days, cross-session loops"),
            ("trivial_day004", "same-session (trivial) loops")]
    fig, ax = plt.subplots(figsize=(7.5, 4))
    for i, (net, lab) in enumerate(nets):
        des = json.loads((RES / "networks" / net / "network_design.json").read_text())
        v = [lp["chi2_3"] for lp in des["loops"]]
        ax.scatter(np.full(len(v), i) + np.linspace(-0.12, 0.12, len(v)), v,
                   color=["#2a78d6", "#eb6834", "#1baf7a"][i], marker="os^"[i], s=45, zorder=3)
        ax.text(i + 0.18, np.mean(v), f"mean {np.mean(v):.1f}", va="center", color=MUTED, fontsize=8)
    for q, ls in ((0.5, ":"), (0.95, "--")):
        ax.axhline(stats.chi2.ppf(q, 3), color=MUTED, ls=ls, lw=1,
                   label=f"χ²(3) {int(q * 100)}% quantile = {stats.chi2.ppf(q, 3):.2f}")
    ax.set_xticks(range(len(nets)), [lab for _, lab in nets], fontsize=8)
    ax.set_ylabel("loop misclosure statistic w'C⁻¹w [-] (block-diagonal C)")
    ax.set_title("Loop misclosures, formal RTKLIB covariance")
    ax.legend(fontsize=8)
    save(fig, "loop_misclosures.png")


def fig_experiment_a():
    b = pd.read_csv(RES / "experiment_A/main_formal/baselines.csv")
    w = b[b.method == "WLS"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    ax = axes[0]
    colors = ["#eb6834" if f else "#2a78d6" for f in w.flag_group]
    ax.bar(w.baseline_id, w.T_b, color=colors, width=0.6)
    ax.axhline(w.T_b_critical.iloc[0], color=MUTED, ls="--", lw=1,
               label=f"critical χ²(3), α=0.001: {w.T_b_critical.iloc[0]:.1f}")
    ax.set_xlabel("baseline id (days 004-006)")
    ax.set_ylabel("WLS 3D test statistic T_b = |w*|² [-]")
    ax.set_title("Experiment A: unmodified real data")
    ax.legend(fontsize=8)
    ax = axes[1]
    ax.bar(w.baseline_id, w.sd_dir_elevation_deg, color="#6250d6", width=0.6)
    ax.set_ylim(0, 90)
    ax.set_xlabel("baseline id")
    ax.set_ylabel("elevation of SD bias direction [deg]")
    ax.set_title("SD direction (90° = local vertical)")
    ax = axes[2]
    for m in ("Huber", "Hampel"):
        s = b[b.method == m]
        mk, ls = STYLE[m]
        ax.plot(s.baseline_id, s.robust_weight, marker=mk, ls=ls, color=COLOR[m], label=m)
    ax.set_ylim(-0.05, 1.1)
    ax.set_xlabel("baseline id")
    ax.set_ylabel("final robust weight multiplier [-]")
    ax.set_title("Robust weights")
    ax.legend(fontsize=8)
    save(fig, "experimentA_statistics.png")


def fig_coordinates():
    c = pd.read_csv(RES / "experiment_A/main_formal/coordinates.csv")
    c = c[~c.fixed]
    stations = sorted(c.station.unique())
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), sharey=False)
    for ax, comp, lab in zip(axes, "enu", ("East", "North", "Up")):
        for j, m in enumerate(METHODS):
            s = c[c.method == m].set_index("station").loc[stations]
            x = np.arange(len(stations)) + 0.12 * (j - 2.5)
            mk, _ = STYLE[m]
            ax.errorbar(x, s[f"diff_{comp}_mm"], yerr=2 * s[f"sigma_{comp}_mm"], fmt=mk,
                        color=COLOR[m], ms=5, capsize=2, lw=1, label=m)
        ax.axhline(0, color=MUTED, lw=1)
        ax.set_xticks(range(len(stations)), stations)
        ax.set_ylabel(f"{lab}: adjusted − EPN [mm] (±2σ formal)")
        ax.set_title(lab)
    axes[0].legend(fontsize=7)
    fig.suptitle("Coordinates relative to datum MADR vs EPN weekly combined solution (days 004-006)",
                 color=INK)
    save(fig, "coordinate_comparison.png")


def fig_injection():
    a = pd.read_csv(RES / "experiment_BC/aggregate.csv")
    geoms = sorted(a.geometry.unique(), reverse=True)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), sharex=True)
    for i, n in enumerate((1, 2)):
        for j, g in enumerate(geoms):
            ax = axes[i, j]
            s = a[(a.n_outliers == n) & (a.geometry == g) & (a.direction_mode == "random")]
            for m in METHODS:
                t = s[s.method == m].sort_values("magnitude_mm")
                mk, ls = STYLE[m]
                ax.plot(t.magnitude_mm, t.f1, marker=mk, ls=ls, color=COLOR[m], label=m)
            ax.set_title(f"{g}, {n} injected baseline{'s' if n > 1 else ''}")
            ax.set_ylabel("detection F1 [-]")
            ax.set_ylim(-0.03, 1.05)
            ax.set_xscale("log")
            ticks = sorted(s.magnitude_mm.unique())
            ax.set_xticks(ticks, [str(t) for t in ticks])
            ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    for ax in axes[1]:
        ax.set_xlabel("injected bias magnitude [mm] (random direction)")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("Experiment B: CONTROLLED injections into real-data-derived baselines", color=INK)
    save(fig, "experimentB_f1.png")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, col, lab in ((axes[0], "rmse_vs_clean_mm", "coordinate RMSE vs unmodified-data solution [mm]"),
                         (axes[1], "fpr", "false-positive rate [-]")):
        s = a[(a.n_outliers == 1) & (a.geometry == geoms[0]) & (a.direction_mode == "random")]
        for m in METHODS:
            t = s[s.method == m].sort_values("magnitude_mm")
            mk, ls = STYLE[m]
            ax.plot(t.magnitude_mm, t[col], marker=mk, ls=ls, color=COLOR[m], label=m)
        ax.set_xscale("log")
        ax.set_xlabel("injected bias magnitude [mm]")
        ax.set_ylabel(lab)
    axes[0].legend(fontsize=7)
    fig.suptitle(f"Experiment B ({geoms[0]}, one injected baseline)", color=INK)
    save(fig, "experimentB_rmse_fpr.png")


def fig_sd():
    d = pd.read_csv(RES / "experiment_BC/sd_directions.csv")
    d = d[d.n_outliers == 1]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    ax = axes[0]
    lam = np.logspace(-0.5, 2.5, 200)
    c3, c1 = stats.chi2.ppf(0.999, 3), stats.norm.ppf(0.9995) ** 2
    ax.plot(lam, stats.ncx2.sf(c3, 3, lam), color="#2a78d6", label="theory: 3D / SD test, χ²(3)")
    ax.plot(lam, stats.ncx2.sf(c1, 1, lam), color="#eb6834", ls="--",
            label="theory: prespecified true direction, χ²(1)")
    bins = np.logspace(-0.5, 2.5, 13)
    d["bin"] = pd.cut(d.noncentrality, bins)
    g = d.groupby("bin", observed=True)
    mid = np.array([iv.mid for iv in g.groups.keys()])
    ax.plot(mid, g.sd_significant.mean(), "o", color="#2a78d6", label="empirical: SD (= 3D) test")
    ax.plot(mid, g.detect_1d_true.mean(), "s", color="#eb6834", label="empirical: 1D true direction")
    ax.plot(mid, g.detect_1d_up.mean(), "^", color="#1baf7a", label="empirical: 1D local-up (often wrong)")
    ax.set_xscale("log")
    ax.set_xlabel("non-centrality λ = dᵀ P̄_ii d [-]")
    ax.set_ylabel("detection rate, α = 0.001 [-]")
    ax.set_title("Optimised vs prespecified direction")
    ax.legend(fontsize=7)
    ax = axes[1]
    for k, (geo, gd) in enumerate(d[d.direction_mode == "random"].groupby("geometry")):
        q = gd.groupby("magnitude_mm").angle_deg
        ax.plot(q.median().index, q.median().values, marker="os"[k], color=["#2a78d6", "#eb6834"][k],
                label=f"{geo}: median")
        ax.plot(q.quantile(0.9).index, q.quantile(0.9).values, marker="os"[k], ls=":",
                color=["#2a78d6", "#eb6834"][k], label=f"{geo}: 90th pct")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ticks = sorted(d.magnitude_mm.unique())
    ax.set_xticks(ticks, [str(t) for t in ticks])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("injected bias magnitude [mm]")
    ax.set_ylabel("angular error of SD direction [deg]")
    ax.set_title("Experiment C: direction recovery")
    ax.legend(fontsize=7)
    ax = axes[2]
    r = (d.angle_deg / d.angle_sigma_pred_deg).clip(upper=6)
    ax.hist(r, bins=30, color="#6250d6", rwidth=0.9)
    ax.axvline(np.median(r), color=INK, ls="--", lw=1, label=f"median {np.median(r):.2f}")
    ax.set_xlabel("angular error / predicted angular σ [-]")
    ax.set_ylabel("injections")
    ax.set_title("Predicted direction uncertainty")
    ax.legend(fontsize=8)
    fig.suptitle("Experiment C: specific-direction validation (controlled injections, real network)",
                 color=INK)
    save(fig, "experimentC_sd_direction.png")


def fig_antenna():
    d = pd.read_csv(RES / "experiment_D/summary.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for k, ((case, base), g) in enumerate(d.groupby(["case", "baseline"])):
        g = g.sort_values("dh_m")
        c = ["#2a78d6", "#eb6834", "#1baf7a"][k]
        axes[0].plot(100 * g.dh_m, g.sd_statistic, marker="o", color=c, label=f"{case} {base}: |w*| (SD/3D)")
        axes[0].plot(100 * g.dh_m, g.w_up_1d.abs(), marker="s", ls="--", color=c,
                     label=f"{case} {base}: |w| vertical 1D")
        axes[1].plot(100 * g.dh_m, g.sd_angle_to_vertical_deg, marker="o", color=c, label=f"{case} {base}")
    axes[0].axhline(d.sd_critical.iloc[0], color=MUTED, ls="-", lw=1, label="critical |w*| (χ²(3))")
    axes[0].axhline(d.w_up_critical.iloc[0], color=MUTED, ls=":", lw=1, label="critical |w| (N(0,1))")
    axes[0].set_xlabel("antenna height error entered in RTKLIB [cm]")
    axes[0].set_ylabel("test statistic [-]")
    axes[0].legend(fontsize=7)
    axes[0].set_title("Detection of a reprocessed antenna-height error")
    axes[1].set_xlabel("antenna height error [cm]")
    axes[1].set_ylabel("SD direction: angle to local vertical [deg]")
    axes[1].set_title("SD direction recovers the vertical")
    axes[1].legend(fontsize=7)
    fig.suptitle("Experiment D: RINEX reprocessed with perturbed antenna height (D1 YEBE, D2 VILL)",
                 color=INK)
    save(fig, "experimentD_antenna.png")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for f in (fig_repeatability, fig_formal_vs_empirical, fig_loops, fig_experiment_a,
              fig_coordinates, fig_injection, fig_sd, fig_antenna):
        f()
        print("ok", f.__name__)


if __name__ == "__main__":
    sys.exit(main())
