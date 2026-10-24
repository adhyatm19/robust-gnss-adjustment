"""HYPOTHETICAL cross-baseline correlation sensitivity (not measured covariance).

Model per session s: baseline error e_ij = eta_j - eta_i + eps_ij, with
station-session errors eta ~ N(0, (rho/2) Qbar_s) and eps ~ N(0, (1-rho) Q_ij),
Qbar_s = mean formal block of session s. Marginal baseline variance ~ Q_ij;
two baselines sharing a station in a session get correlation of order rho/2.
"""
import json
from pathlib import Path
import numpy as np, pandas as pd
ROOT = Path(__file__).resolve().parents[2]
from scipy import stats

def load(net):
    d = f"{ROOT}/results/real_network/networks/{net}/"
    st = pd.read_csv(d + "stations.csv"); bl = pd.read_csv(d + "baselines.csv")
    cv = pd.read_csv(d + "covariance_blocks.csv"); md = pd.read_csv(d + "baseline_metadata.csv")
    o = st[st.is_fixed == 1].iloc[0][["x", "y", "z"]].values.astype(float)
    free = sorted(st[st.is_fixed == 0].station_id); col = {s: 3*i for i, s in enumerate(free)}
    fixed = {r.station_id: np.array([r.x, r.y, r.z]) - o for r in st.itertuples() if r.is_fixed}
    m = len(bl); A = np.zeros((3*m, 3*len(free))); l = np.zeros(3*m); Qb = []
    for k, r in bl.iterrows():
        s = slice(3*k, 3*k+3); l[s] = [r.dx, r.dy, r.dz]
        c = cv[cv.baseline_id == r.baseline_id].iloc[0]
        Qb.append(np.array([[c.q_xx, c.q_xy, c.q_xz], [c.q_xy, c.q_yy, c.q_yz], [c.q_xz, c.q_yz, c.q_zz]]))
        for sid, sg in ((r.to_station, 1), (r.from_station, -1)):
            if sid in fixed: l[s] -= sg * fixed[sid]
            else: A[s, col[sid]:col[sid]+3] = sg * np.eye(3)
    sess = md.set_index("baseline_id").loc[bl.baseline_id, "session"].astype(str).values
    return A, l, Qb, sess, bl

def corr_cov(Qb, sess, bl, rho):
    m = len(Qb); Q = np.zeros((3*m, 3*m))
    for k in range(m): Q[3*k:3*k+3, 3*k:3*k+3] = (1 - rho) * Qb[k]
    for s in set(sess):
        ks = [k for k in range(m) if sess[k] == s]
        Qbar = np.mean([Qb[k] for k in ks], axis=0) * rho / 2
        stations = set(bl.from_station[ks]) | set(bl.to_station[ks])
        for sta in stations:
            g = np.zeros((3*m, 3))
            for k in ks:
                if bl.to_station[k] == sta: g[3*k:3*k+3] = np.eye(3)
                if bl.from_station[k] == sta: g[3*k:3*k+3] = -np.eye(3)
            Q += g @ Qbar @ g.T
    return 0.5 * (Q + Q.T)

def tests(A, Q, l):
    P = np.linalg.inv(Q); N = A.T @ P @ A; x = np.linalg.solve(N, A.T @ P @ l); v = A @ x - l
    Qv = Q - A @ np.linalg.solve(N, A.T); om = float(v @ P @ v); T = []
    for k in range(A.shape[0] // 3):
        C = np.zeros((A.shape[0], 3)); C[3*k:3*k+3] = np.eye(3)
        g = C.T @ P @ v; T.append(float(g @ np.linalg.solve(C.T @ P @ Qv @ P @ C, g)))
    return om, np.array(T)

out = {}
A, l, Qb, sess, bl = load("main_formal")
for rho in (0.0, 0.5, 1.0):
    Q = corr_cov(Qb, sess, bl, rho)
    om, T = tests(A, Q, l)
    ev = np.linalg.eigvalsh(Q)[0]
    out[f"expA_rho{rho}"] = {"vTPv": om, "T_b": T.round(2).tolist(),
                             "flag_order": [int(bl.baseline_id[i]) for i in np.argsort(-T)[:3]],
                             "min_eig": float(ev)}
    print(f"rho={rho}: v'Pv={om:.1f} (crit {stats.chi2.ppf(.95,18):.1f}); T_b={np.round(T,1)}")
# size of block-diagonal tests if correlation is real (resampled)
rng = np.random.default_rng(5)
A, l, Qb, sess, bl = load("quiet_formal")
Qbd = corr_cov(Qb, sess, bl, 0.0)
for rho in (0.5, 1.0):
    Qt = corr_cov(Qb, sess, bl, rho); L = np.linalg.cholesky(Qt)
    g = loc = 0; n = 4000; crit = stats.chi2.ppf(.999, 3)
    for _ in range(n):
        e = L @ rng.normal(size=len(l))
        om, T = tests(A, Qbd, A @ np.zeros(A.shape[1]) + e)
        g += om > stats.chi2.ppf(.95, 18); loc += np.mean(T > crit)
    out[f"size_rho{rho}"] = {"global_reject": g / n, "local_exceedance": loc / n}
    print(f"true rho={rho}, block-diag tests: global reject {g/n:.3f} (nominal 0.05), local {loc/n:.4f} (nominal 0.001)")
(ROOT / "results/real_network/audit/hypothetical_cross_covariance.json").write_text(json.dumps(out, indent=1))
