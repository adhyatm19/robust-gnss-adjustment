"""Streamlit demonstration app for robust GNSS baseline network adjustment.

Run with:  streamlit run app.py
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import streamlit as st

from gnss_adjust import benchmark as bm
from gnss_adjust import plotting as plots
from gnss_adjust.dia import DIAConfig
from gnss_adjust.diagnostics import baseline_tests, global_model_test
from gnss_adjust.metrics import coordinate_metrics, detection_metrics
from gnss_adjust.robust import RobustConfig
from gnss_adjust.simulation import NetworkConfig, NoiseModel, OutlierSpec, simulate_network

st.set_page_config(page_title="Robust GNSS Network Adjustment", layout="wide")

st.title("Robust GNSS Baseline Network Adjustment")
st.info("**Synthetic baseline-level research prototype; not raw GNSS/RINEX processing.** "
        "Baseline vectors and covariances are simulated; ground truth is known, so "
        "detection quality can be scored exactly.")

# ----------------------------------------------------------------- sidebar
with st.sidebar:
    st.header("Network")
    seed = st.number_input("Random seed", value=42, step=1)
    n_stations = st.slider("Stations", 4, 20, 8)
    max_b = n_stations * (n_stations - 1) // 2
    n_baselines = st.slider("Baselines", n_stations - 1, max_b,
                            min(2 * (n_stations - 1), max_b))
    a_mm = st.slider("Noise a [mm]", 1.0, 20.0, 5.0)
    b_ppm = st.slider("Noise b [ppm]", 0.0, 5.0, 1.0)
    heavy = st.checkbox("Student-t heavy-tailed noise (dof=3)")

    st.header("Outlier")
    inject = st.checkbox("Inject outlier", value=True)
    outlier_baseline = st.number_input("Corrupted baseline id", 0, n_baselines - 1, 3)
    mode = st.radio("Outlier specification", ["vector", "sigma magnitude"])
    if mode == "vector":
        c1, c2, c3 = st.columns(3)
        ox = c1.number_input("dX [m]", value=0.035, format="%.3f")
        oy = c2.number_input("dY [m]", value=-0.025, format="%.3f")
        oz = c3.number_input("dZ [m]", value=0.090, format="%.3f")
    else:
        mag_sigma = st.slider("Magnitude [sigma]", 1.0, 60.0, 25.0)

    st.header("Methods")
    methods = st.multiselect("Run methods", bm.ALL_METHODS, default=bm.ALL_METHODS)
    alpha_global = st.select_slider("alpha (global test)", [0.001, 0.01, 0.05, 0.1], value=0.05)
    alpha_local = st.select_slider("alpha (local tests)", [0.0001, 0.001, 0.01], value=0.001)
    huber_c = st.slider("Huber c", 1.0, 3.0, 1.5, 0.1)
    hampel_a = st.slider("Hampel a", 1.0, 3.0, 1.5, 0.1)
    hampel_b = st.slider("Hampel b", 2.0, 6.0, 3.5, 0.1)
    hampel_c = st.slider("Hampel c", 4.0, 15.0, 8.0, 0.5)

# ----------------------------------------------------------------- simulate
outliers = []
if inject:
    if mode == "vector":
        outliers = [OutlierSpec(baseline_id=int(outlier_baseline), vector=(ox, oy, oz))]
    else:
        outliers = [OutlierSpec(baseline_id=int(outlier_baseline), magnitude_sigma=mag_sigma)]

cfg = NetworkConfig(
    n_stations=int(n_stations), n_baselines=int(n_baselines), seed=int(seed),
    noise=NoiseModel(a_mm=a_mm, b_ppm=b_ppm, student_t_dof=3.0 if heavy else None),
    outliers=outliers,
)
stations, observations, truth = simulate_network(cfg)
true_out = set(truth.outlier_baselines)

if hampel_a >= hampel_b or hampel_b >= hampel_c:
    st.error("Hampel thresholds must satisfy a < b < c.")
    st.stop()

settings = bm.MethodSettings(
    dia=DIAConfig(alpha_global=alpha_global, alpha_local=alpha_local,
                  alpha_baseline=alpha_local),
    huber=RobustConfig(method="huber", huber_c=huber_c),
    hampel=RobustConfig(method="hampel", hampel_a=hampel_a, hampel_b=hampel_b,
                        hampel_c=hampel_c),
)

if not methods:
    st.warning("Select at least one method in the sidebar.")
    st.stop()

outcomes = bm.run_all_methods(stations, observations, settings, methods,
                              true_outliers=true_out)
rows = bm.score_run(outcomes, truth, stations, observations, scenario="app")
df = pd.DataFrame(rows).drop(columns=["run_id", "scenario", "seed"])

# ----------------------------------------------------------------- layout
left, right = st.columns([1.1, 1])

with left:
    st.subheader("Network")
    detected = outcomes["DIA-reject"].flagged if "DIA-reject" in outcomes else set()
    downweighted = set()
    for name in ("Huber", "Hampel"):
        if name in outcomes:
            downweighted |= outcomes[name].flagged
    fig = plots.plot_network(stations, observations, true_outliers=true_out,
                             detected=detected, downweighted=downweighted)
    st.pyplot(fig)
    if true_out:
        st.caption(f"True corrupted baseline(s): {sorted(true_out)}; injected bias "
                   + ", ".join(f"{k}: {np.round(v, 4).tolist()} m"
                               for k, v in truth.outlier_vectors.items()))

with right:
    st.subheader("Method comparison")
    show = df[["method", "rmse_m", "mean_pos_error_m", "max_pos_error_m",
               "precision", "recall", "f1", "n_flagged", "iterations",
               "converged", "variance_factor"]]
    st.dataframe(show.style.format({
        "rmse_m": "{:.4f}", "mean_pos_error_m": "{:.4f}", "max_pos_error_m": "{:.4f}",
        "precision": "{:.2f}", "recall": "{:.2f}", "f1": "{:.2f}",
        "variance_factor": "{:.2f}"}), width="stretch")
    st.caption("RMSE and position errors in metres against known ground truth; "
               "precision/recall/F1 score baseline-level outlier flags.")

    st.subheader("True vs estimated coordinates")
    fixed_id = next(s.station_id for s in stations if s.is_fixed)
    fig2 = plots.plot_true_vs_estimated(
        truth.true_coordinates,
        {name: oc.coordinates for name, oc in outcomes.items()},
        fixed_station=fixed_id)
    st.pyplot(fig2)

st.divider()
c1, c2 = st.columns(2)

with c1:
    st.subheader("Residual diagnostics (WLS)")
    base = outcomes.get("WLS") or next(iter(outcomes.values()))
    res = base.detail.final if hasattr(base.detail, "final") else base.detail
    gt = global_model_test(res, alpha=alpha_global)
    st.metric("Global chi-square statistic",
              f"{gt.statistic:.1f}",
              f"critical {gt.critical:.1f} ({'PASS' if gt.passed else 'FAIL'})",
              delta_color="normal" if gt.passed else "inverse")
    tests = baseline_tests(res, alpha_local=alpha_local, alpha_baseline=alpha_local,
                           true_outliers=true_out)
    st.pyplot(plots.plot_group_statistics(tests))
    st.pyplot(plots.plot_standardized_residuals(tests))

with c2:
    st.subheader("Robust weights")
    weights_by_method = {n: outcomes[n].weights for n in ("Huber", "Hampel")
                         if n in outcomes and outcomes[n].weights is not None}
    if weights_by_method:
        st.pyplot(plots.plot_robust_weights(weights_by_method, true_outliers=true_out))
        for n, oc in outcomes.items():
            if oc.weights is not None:
                r = oc.detail
                st.caption(f"**{n}**: converged={r.converged} in {r.iterations} "
                           f"iterations; down-weighted: {r.downweighted or 'none'}")
    else:
        st.caption("Run Huber or Hampel to see robust weights.")

    st.subheader("DIA audit trail")
    shown = False
    for name in ("DIA-reject", "DIA-inflate"):
        if name not in outcomes:
            continue
        shown = True
        dia_res = outcomes[name].detail
        st.markdown(f"**{name}** — status: `{dia_res.status}`, "
                    f"flagged: {dia_res.flagged_baselines or 'none'}")
        audit = pd.DataFrame([{
            "iter": it.iteration,
            "T_global": round(it.global_statistic, 2),
            "critical": round(it.global_critical, 2),
            "passed": it.global_passed,
            "identified": it.identified_baseline,
            "T_baseline": None if it.identified_statistic is None
                          else round(it.identified_statistic, 2),
            "action": it.action,
        } for it in dia_res.iterations])
        st.dataframe(audit, width="stretch", hide_index=True)
    if not shown:
        st.caption("Run a DIA method to see the audit trail.")
