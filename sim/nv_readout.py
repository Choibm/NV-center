"""
NV-center charge-state readout: reproduction of Fig. 5 of
D'Anjou, Kuret, Childress, Burkard, Phys. Rev. X 6, 011017 (2016).

Three readout methods are compared (error rate eps vs average time T):
  1. Photon counting           — total photons over [0, tf] vs threshold nu
  2. Nonadaptive MLE           — log-likelihood ratio at fixed tf, sign decision
  3. Adaptive MLE              — sequential update with stopping at posterior

The 2x2 update matrix M(dn) is the dn-th Taylor coefficient (in z) of the
matrix exponential exp((L - K + zK) dt), giving the joint conditional
probability of the charge state and dn detected photons in time dt.
We extract the Taylor coefficients via FFT-based contour integration on
the unit circle in the complex plane.
"""

import sys
import time
import numpy as np
from scipy.linalg import expm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Force UTF-8 stdout (Windows consoles often default to cp949/cp1252)
try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass


# ---------------------------------------------------------------------------
# Physical parameters (Table I of D'Anjou 2016)
# ---------------------------------------------------------------------------
GAMMA_P = 720.0   # Hz, photon detection rate in NV-  (|+>)
GAMMA_M = 50.0    # Hz, photon detection rate in NV0  (|->)
CAP_GAMMA_P = 3.6   # Hz, ionization        NV- -> NV0
CAP_GAMMA_M = 0.98  # Hz, recombination     NV0 -> NV-

# Generators
L = np.array([[-CAP_GAMMA_P, CAP_GAMMA_M],
              [ CAP_GAMMA_P, -CAP_GAMMA_M]])
K = np.array([[GAMMA_P, 0.0],
              [0.0,     GAMMA_M]])

# Time discretisation
DT      = 0.1e-3       # 0.1 ms bin
T_MAX   = 25e-3        # 25 ms
N_BINS  = int(round(T_MAX / DT))   # 250
MAX_DN  = 5            # P(dn > 5) ~ 1e-10 per bin


# ---------------------------------------------------------------------------
# Build the M(dn) update matrices via contour integration
# ---------------------------------------------------------------------------
def build_M_matrices(L, K, dt, max_dn=5, n_contour=128, radius=1.0):
    """Return list of M(0)..M(max_dn).

    M(n) = (1/(2 pi i)) ∮ exp((L-K+zK) dt) z^{-n-1} dz
         = (1/N) Σ_k f(r e^{i 2π k/N}) e^{-i 2π k n / N} / r^n
    where f(z) = expm((L-K+zK) dt).  The unit circle keeps the
    truncation error tiny because the radius of convergence is infinite
    (matrix exponential is entire in z).
    """
    A = L - K
    angles = np.arange(n_contour) * 2.0 * np.pi / n_contour
    f_values = np.empty((n_contour, 2, 2), dtype=complex)
    for k, theta in enumerate(angles):
        z = radius * np.exp(1j * theta)
        f_values[k] = expm((A + z * K) * dt)
    M_list = []
    for n in range(max_dn + 1):
        weights = np.exp(-1j * n * angles) / radius**n
        M_n = np.tensordot(weights, f_values, axes=1) / n_contour
        M_list.append(np.real(M_n))
    return M_list


M_LIST  = build_M_matrices(L, K, DT, MAX_DN)
M_ARRAY = np.stack(M_LIST)  # (MAX_DN+1, 2, 2)

# Sanity checks
sumM = M_ARRAY.sum(axis=0)
expL = expm(L * DT)
assert np.allclose(sumM, expL, atol=1e-10), "Sum_n M(n) must equal exp(L dt)"
# Each row of (1,1).M(n) gives the photon-count distribution conditional on
# starting state.  Test against a Poisson(gamma * dt) (when charge dynamics
# is frozen) — should be close at this short timescale.
ph_dist_minus = np.array([(np.array([1.0, 1.0]) @ Mn @ np.array([1.0, 0.0]))
                          for Mn in M_LIST])
ph_dist_zero  = np.array([(np.array([1.0, 1.0]) @ Mn @ np.array([0.0, 1.0]))
                          for Mn in M_LIST])
assert np.isclose(ph_dist_minus.sum(), 1.0, atol=1e-9)
assert np.isclose(ph_dist_zero.sum(),  1.0, atol=1e-9)


# ---------------------------------------------------------------------------
# Combined trajectory generator + log-likelihood tracker
# ---------------------------------------------------------------------------
def simulate(M_array, initial_state, n_traj, n_bins, rng):
    """Generate trajectories and compute λ_t at every bin in one pass.

    initial_state: 0 -> NV- (|+>), 1 -> NV0 (|->)
    Returns
    -------
    photon_counts : (n_traj, n_bins) int8     - dn observed per bin
    log_lambdas   : (n_traj, n_bins+1) float32 - λ_t = ln P(record|+) / P(record|-)
    """
    photon_counts = np.empty((n_traj, n_bins), dtype=np.int8)
    log_lambdas   = np.zeros((n_traj, n_bins + 1), dtype=np.float32)

    # State posterior given the *true* initial state (used for sampling)
    rho = np.zeros((n_traj, 2)); rho[:, initial_state] = 1.0
    # Hypothesis-conditional unnormalised vectors (renormalised each step)
    v_p = np.zeros((n_traj, 2)); v_p[:, 0] = 1.0    # start in NV-
    v_m = np.zeros((n_traj, 2)); v_m[:, 1] = 1.0    # start in NV0
    log_p = np.zeros(n_traj)
    log_m = np.zeros(n_traj)

    n_photon = M_array.shape[0]
    eps_floor = 1e-300

    for k in range(n_bins):
        # Candidate updates for every dn in 0..MAX_DN.
        # M_array (6,2,2) acting on rho (n_traj,2):
        # candidates[n, j, i] = sum_l M_array[n, i, l] rho[j, l]
        candidates = np.matmul(M_array, rho.T).transpose(0, 2, 1)  # (6, n_traj, 2)
        probs = candidates.sum(axis=2)                              # (6, n_traj)

        # Sample dn for each trajectory using the per-trajectory cdf.
        probs_T = probs.T                                           # (n_traj, 6)
        np.clip(probs_T, 0.0, None, out=probs_T)
        probs_T /= probs_T.sum(axis=1, keepdims=True)
        cdf = probs_T.cumsum(axis=1)
        u = rng.random(n_traj)
        deltas = (cdf < u[:, None]).sum(axis=1).astype(np.int8)
        np.clip(deltas, 0, n_photon - 1, out=deltas)
        photon_counts[:, k] = deltas

        idx = np.arange(n_traj)
        rho = candidates[deltas, idx, :] \
              / np.maximum(probs[deltas, idx], eps_floor)[:, None]

        # Apply the same M(dn) to the two hypothesis vectors.
        M_batch = M_array[deltas]                                   # (n_traj, 2, 2)
        v_p_un = np.matmul(M_batch, v_p[..., None])[..., 0]
        v_m_un = np.matmul(M_batch, v_m[..., None])[..., 0]
        np_norm = v_p_un.sum(axis=1)
        nm_norm = v_m_un.sum(axis=1)
        v_p = v_p_un / np.maximum(np_norm, eps_floor)[:, None]
        v_m = v_m_un / np.maximum(nm_norm, eps_floor)[:, None]
        log_p += np.log(np.maximum(np_norm, eps_floor))
        log_m += np.log(np.maximum(nm_norm, eps_floor))

        log_lambdas[:, k + 1] = (log_p - log_m).astype(np.float32)

    return photon_counts, log_lambdas


# ---------------------------------------------------------------------------
# Run Monte Carlo
# ---------------------------------------------------------------------------
N_TRAJ = 100_000
rng = np.random.default_rng(20160117)

t0 = time.time()
print(f"Generating {N_TRAJ} trajectories from each initial state...")
photons_p, lam_p = simulate(M_ARRAY, 0, N_TRAJ, N_BINS, rng)
photons_m, lam_m = simulate(M_ARRAY, 1, N_TRAJ, N_BINS, rng)
print(f"  done in {time.time() - t0:.1f} s")


# ---------------------------------------------------------------------------
# Method 1: photon counting (optimal integer threshold per tf)
# ---------------------------------------------------------------------------
print("Photon counting analysis...")
cum_p = photons_p.cumsum(axis=1).astype(np.int32)
cum_m = photons_m.cumsum(axis=1).astype(np.int32)

tf_indices = np.arange(1, N_BINS + 1)
T_pc = tf_indices * DT

eps_pc = np.empty(N_BINS)
nu_pc  = np.empty(N_BINS, dtype=int)
for j, tf_idx in enumerate(tf_indices):
    np_arr = cum_p[:, tf_idx - 1]
    nm_arr = cum_m[:, tf_idx - 1]
    n_max = max(np_arr.max(), nm_arr.max())
    nus = np.arange(n_max + 1)
    # vectorised over candidate thresholds
    eps_p_cand = (np_arr[:, None] <= nus[None, :]).mean(axis=0)
    eps_m_cand = (nm_arr[:, None] >  nus[None, :]).mean(axis=0)
    eps_cand = 0.5 * (eps_p_cand + eps_m_cand)
    best = int(np.argmin(eps_cand))
    eps_pc[j] = eps_cand[best]
    nu_pc[j]  = best


# ---------------------------------------------------------------------------
# Method 2: nonadaptive MLE
# ---------------------------------------------------------------------------
print("Nonadaptive MLE analysis...")
T_nmle = T_pc.copy()
eps_nmle = np.empty(N_BINS)
for j, tf_idx in enumerate(tf_indices):
    eps_p_val = (lam_p[:, tf_idx] <= 0.0).mean()
    eps_m_val = (lam_m[:, tf_idx] >  0.0).mean()
    eps_nmle[j] = 0.5 * (eps_p_val + eps_m_val)


# ---------------------------------------------------------------------------
# Method 3: adaptive MLE — Pareto sweep over (p_+, p_-)
# ---------------------------------------------------------------------------
print("Adaptive MLE Pareto sweep...")

def adaptive(lam, theta_p, theta_m, dt, n_bins):
    """Return per-trajectory (decision, stop_time)."""
    above = lam >= theta_p
    below = lam <= theta_m
    crossed = above | below
    crossed[:, 0] = False  # never decide before first observation

    first = np.argmax(crossed, axis=1)            # 0 if none
    has   = crossed.any(axis=1)
    times = first.astype(np.float64) * dt
    times[~has] = n_bins * dt

    n = lam.shape[0]
    rng_ = np.arange(n)
    crossed_value = lam[rng_, first]
    decisions = np.where(has,
                         np.where(crossed_value >= theta_p, 1, -1),
                         np.where(lam[:, -1] >= 0.0, 1, -1))
    return decisions, times


# Symmetric thresholds parametrised by a single posterior cutoff p, plus
# a small asymmetric sweep to capture any prior bias.
p_sym = 1.0 - np.logspace(-0.3, -8, 80)          # 0.5 ... 1 - 1e-8
asym_factors = np.array([0.5, 0.7, 1.0, 1.4, 2.0])

results = []
t1 = time.time()
for p_plus in p_sym:
    eps_plus = 1.0 - p_plus                        # i.e. the "1-p_+" complement
    for f in asym_factors:
        eps_minus = eps_plus * f
        if eps_minus >= 0.5:
            continue
        p_minus = eps_minus
        theta_p = np.log(p_plus  / (1.0 - p_plus))
        theta_m = np.log(p_minus / (1.0 - p_minus))
        if theta_p <= theta_m:
            continue
        d_p, t_p = adaptive(lam_p, theta_p, theta_m, DT, N_BINS)
        d_m, t_m = adaptive(lam_m, theta_p, theta_m, DT, N_BINS)
        eps_p_val = (d_p == -1).mean()
        eps_m_val = (d_m ==  1).mean()
        eps  = 0.5 * (eps_p_val + eps_m_val)
        Tavg = 0.5 * (t_p.mean() + t_m.mean())
        results.append((Tavg, eps, p_plus, p_minus))
print(f"  swept {len(results)} threshold pairs in {time.time() - t1:.1f} s")
results = np.array(results)


def pareto_frontier(points):
    """Lower-left Pareto frontier of (T, eps): minimise both."""
    order = np.argsort(points[:, 0])
    pts = points[order]
    front = []
    best_eps = np.inf
    for p in pts:
        if p[1] < best_eps:
            front.append(p)
            best_eps = p[1]
    return np.array(front)


pareto = pareto_frontier(results[:, :2])


# ---------------------------------------------------------------------------
# Speedup at a reference error rate
# ---------------------------------------------------------------------------
ref_eps = 0.019  # paper's photon-counting floor

# nonadaptive: smallest tf at which eps_nmle <= ref_eps
nmle_ok = np.where(eps_nmle <= ref_eps)[0]
adapt_ok = np.where(pareto[:, 1] <= ref_eps)[0]
tf_ref = T_nmle[nmle_ok[0]] * 1e3 if len(nmle_ok) else np.nan
T_ref  = pareto[adapt_ok[0], 0] * 1e3 if len(adapt_ok) else np.nan
speedup = tf_ref / T_ref if T_ref and not np.isnan(T_ref) else np.nan


# ---------------------------------------------------------------------------
# Plot Figure 5
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8.5, 6.0))

ax.plot(T_pc * 1e3,  eps_pc,  color="tab:blue",   linestyle="-",
        linewidth=2.0, label="Photon counting")
ax.plot(T_nmle * 1e3, eps_nmle, color="tab:orange", linestyle="--",
        linewidth=2.0, label="Nonadaptive MLE")
ax.plot(pareto[:, 0] * 1e3, pareto[:, 1], color="tab:green", linestyle=":",
        linewidth=2.5, label="Adaptive MLE")

eps_min_pc   = float(eps_pc.min())
eps_min_nmle = float(eps_nmle.min())
eps_min_amle = float(pareto[:, 1].min())

ax.axhline(eps_min_pc,   color="tab:blue",   linestyle=":", alpha=0.45)
ax.axhline(eps_min_nmle, color="tab:orange", linestyle=":", alpha=0.45)
ax.axhline(eps_min_amle, color="tab:green",  linestyle=":", alpha=0.45)

if not np.isnan(tf_ref):
    ax.axvline(tf_ref, color="tab:orange", linestyle=":", alpha=0.45)
if not np.isnan(T_ref):
    ax.axvline(T_ref,  color="tab:green",  linestyle=":", alpha=0.45)

ax.set_xscale("linear")
ax.set_yscale("log")
ax.set_xlim(0.0, 25.0)
ax.set_ylim(8e-3, 0.55)
ax.set_xlabel("Average readout time  T  (ms)", fontsize=12)
ax.set_ylabel("Error rate  ε", fontsize=12)
ax.set_title("NV charge-state readout (D'Anjou et al. PRX 6, 011017)", fontsize=12)
ax.grid(True, which="both", alpha=0.3)
ax.legend(fontsize=11, loc="upper right")

txt = (f"min ε:   PC = {eps_min_pc*100:.2f}%   NMLE = {eps_min_nmle*100:.2f}%"
       f"   AMLE = {eps_min_amle*100:.2f}%")
if not np.isnan(speedup):
    txt += f"\nspeedup at ε={ref_eps:.3f}:  tf/T = {speedup:.2f}"
ax.text(0.02, 0.04, txt, transform=ax.transAxes, fontsize=10,
        bbox=dict(facecolor="white", alpha=0.85, edgecolor="0.7"))

plt.tight_layout()
import os
script_dir = os.path.dirname(os.path.abspath(__file__))
out_path = os.path.join(script_dir, "figure5.png")
plt.savefig(out_path, dpi=160)
plt.close(fig)

# Cache the curve data so the plot can be re-rendered without rerunning the MC
np.savez(os.path.join(script_dir, "curves.npz"),
         T_pc=T_pc, eps_pc=eps_pc, nu_pc=nu_pc,
         T_nmle=T_nmle, eps_nmle=eps_nmle,
         pareto=pareto, results=results,
         eps_min_pc=eps_min_pc, eps_min_nmle=eps_min_nmle,
         eps_min_amle=eps_min_amle,
         tf_ref=tf_ref, T_ref=T_ref, speedup=speedup, ref_eps=ref_eps,
         M_array=M_ARRAY)


# ---------------------------------------------------------------------------
# Diagnostic figure: 4×2 panel of related quantities
# ---------------------------------------------------------------------------
print("Building diagnostic figure...")

tf_pc_opt_idx   = int(np.argmin(eps_pc))
tf_nmle_opt_idx = int(np.argmin(eps_nmle))
tf_pc_opt       = T_pc[tf_pc_opt_idx]
tf_nmle_opt     = T_nmle[tf_nmle_opt_idx]

# Operating point for the adaptive visualisations (panels C, D)
p_plus_show  = 0.99
p_minus_show = 0.01
theta_p_show = np.log(p_plus_show  / (1.0 - p_plus_show))
theta_m_show = np.log(p_minus_show / (1.0 - p_minus_show))

fig2, axes = plt.subplots(4, 2, figsize=(13, 16))

# --- [A] Photon-count histograms at the photon-counting optimum --------------
axA = axes[0, 0]
n_p_at = cum_p[:, tf_pc_opt_idx - 1]
n_m_at = cum_m[:, tf_pc_opt_idx - 1]
n_max = max(int(n_p_at.max()), int(n_m_at.max()))
bins = np.arange(n_max + 2) - 0.5
axA.hist(n_p_at, bins=bins, alpha=0.55, color="tab:blue",
         density=True, label=r"started in NV$^-$")
axA.hist(n_m_at, bins=bins, alpha=0.55, color="tab:red",
         density=True, label=r"started in NV$^0$")
axA.axvline(nu_pc[tf_pc_opt_idx] + 0.5, color="k", linestyle="--",
            label=fr"optimal $\nu$ = {nu_pc[tf_pc_opt_idx]}")
axA.set_xlabel("Total photons  n")
axA.set_ylabel("Probability")
axA.set_title(fr"[A] Photon-count distribution at $t_f$ = {tf_pc_opt*1e3:.1f} ms")
axA.legend()
axA.grid(True, alpha=0.3)

# --- [B] λ_tf histograms at the nonadaptive MLE optimum ---------------------
axB = axes[0, 1]
lp = lam_p[:, tf_nmle_opt_idx].astype(np.float64)
lm = lam_m[:, tf_nmle_opt_idx].astype(np.float64)
lo, hi = min(lp.min(), lm.min()), max(lp.max(), lm.max())
edges = np.linspace(lo, hi, 90)
axB.hist(lp, bins=edges, alpha=0.55, color="tab:blue",
         density=True, label=r"started in NV$^-$")
axB.hist(lm, bins=edges, alpha=0.55, color="tab:red",
         density=True, label=r"started in NV$^0$")
axB.axvline(0.0, color="k", linestyle="--", label=r"$\lambda$ = 0")
axB.set_xlabel(r"Log-likelihood ratio $\lambda_{t_f}$")
axB.set_ylabel("Probability density")
axB.set_title(fr"[B] MLE evidence at $t_f$ = {tf_nmle_opt*1e3:.1f} ms")
axB.legend()
axB.grid(True, alpha=0.3)

# --- [C] Sample posterior trajectories with adaptive thresholds ------------
axC = axes[1, 0]
n_show = 12
times_grid = np.arange(N_BINS + 1) * DT * 1e3
sigmoid = lambda x: 1.0 / (1.0 + np.exp(-x))


def first_stop(lam_traj, theta_p, theta_m):
    above = lam_traj >= theta_p
    below = lam_traj <= theta_m
    crossed = above | below
    crossed[0] = False
    if crossed.any():
        return int(np.argmax(crossed))
    return len(lam_traj) - 1


for i in range(n_show):
    pt = sigmoid(lam_p[i].astype(np.float64))
    si = first_stop(lam_p[i].astype(np.float64), theta_p_show, theta_m_show)
    axC.plot(times_grid[:si + 1], pt[:si + 1], color="tab:blue",
             alpha=0.55, linewidth=1.0)
    axC.plot(times_grid[si], pt[si], "o", color="tab:blue", markersize=4)
for i in range(n_show):
    pt = sigmoid(lam_m[i].astype(np.float64))
    si = first_stop(lam_m[i].astype(np.float64), theta_p_show, theta_m_show)
    axC.plot(times_grid[:si + 1], pt[:si + 1], color="tab:red",
             alpha=0.55, linewidth=1.0)
    axC.plot(times_grid[si], pt[si], "o", color="tab:red", markersize=4)
axC.axhline(p_plus_show,  color="k", linestyle="--", alpha=0.6,
            label=fr"$p_+$ = {p_plus_show}")
axC.axhline(p_minus_show, color="k", linestyle=":",  alpha=0.6,
            label=fr"$p_-$ = {p_minus_show}")
axC.set_xlabel("Time  t  (ms)")
axC.set_ylabel(r"Posterior  $p_t$ = $P($NV$^-|$record$)$")
axC.set_title(f"[C] Sample posterior trajectories ({n_show} from each state, dots = stop)")
axC.set_xlim(0.0, T_MAX * 1e3)
axC.set_ylim(0.0, 1.0)
axC.grid(True, alpha=0.3)
axC.legend(loc="center right", fontsize=9)

# --- [D] Stopping-time distribution at the chosen adaptive threshold -------
axD = axes[1, 1]
d_p, t_p = adaptive(lam_p, theta_p_show, theta_m_show, DT, N_BINS)
d_m, t_m = adaptive(lam_m, theta_p_show, theta_m_show, DT, N_BINS)
bins_t = np.linspace(0.0, T_MAX * 1e3, 60)
axD.hist(t_p * 1e3, bins=bins_t, alpha=0.55, color="tab:blue",
         density=True, label=fr"started in NV$^-$, $\langle T\rangle$={t_p.mean()*1e3:.2f} ms")
axD.hist(t_m * 1e3, bins=bins_t, alpha=0.55, color="tab:red",
         density=True, label=fr"started in NV$^0$, $\langle T\rangle$={t_m.mean()*1e3:.2f} ms")
axD.set_xlabel("Stopping time  T  (ms)")
axD.set_ylabel("Probability density")
axD.set_title(fr"[D] Adaptive stop-time distribution ($p_+={p_plus_show}, p_-={p_minus_show}$)")
axD.legend(fontsize=9)
axD.grid(True, alpha=0.3)

# --- [E] Per-class error rates --------------------------------------------
axE = axes[2, 0]
eps_pc_p   = np.empty(N_BINS); eps_pc_m   = np.empty(N_BINS)
for j, tf_idx in enumerate(tf_indices):
    nu = int(nu_pc[j])
    eps_pc_p[j] = (cum_p[:, tf_idx - 1] <= nu).mean()
    eps_pc_m[j] = (cum_m[:, tf_idx - 1] >  nu).mean()
eps_nmle_p = (lam_p[:, 1:] <= 0).mean(axis=0)
eps_nmle_m = (lam_m[:, 1:] >  0).mean(axis=0)

axE.plot(T_pc * 1e3, eps_pc_p,   color="tab:blue",   linestyle="-",
         label=r"PC, $\varepsilon^{(+)}$ (true NV$^-\!\to$NV$^0$)")
axE.plot(T_pc * 1e3, eps_pc_m,   color="tab:blue",   linestyle="--",
         label=r"PC, $\varepsilon^{(-)}$ (true NV$^0\!\to$NV$^-$)")
axE.plot(T_nmle * 1e3, eps_nmle_p, color="tab:orange", linestyle="-",
         label=r"NMLE, $\varepsilon^{(+)}$")
axE.plot(T_nmle * 1e3, eps_nmle_m, color="tab:orange", linestyle="--",
         label=r"NMLE, $\varepsilon^{(-)}$")
axE.set_yscale("log")
axE.set_xlabel("Time  T  (ms)")
axE.set_ylabel("Per-class error rate")
axE.set_title("[E] Per-class errors  (asymmetry from γ$_+\\neq$γ$_-$, Γ$_+\\neq$Γ$_-$)")
axE.set_xlim(0.0, T_MAX * 1e3)
axE.set_ylim(5e-3, 0.6)
axE.legend(fontsize=8, loc="upper right")
axE.grid(True, which="both", alpha=0.3)

# --- [F] Optimal photon-counting threshold ν vs tf -------------------------
axF = axes[2, 1]
axF.step(T_pc * 1e3, nu_pc, where="post", color="tab:blue", linewidth=2.0)
axF.set_xlabel(r"$t_f$  (ms)")
axF.set_ylabel(r"Optimal threshold  $\nu$")
axF.set_title(r"[F] PC: integer threshold that minimises $\varepsilon$ at each $t_f$")
axF.set_xlim(0.0, T_MAX * 1e3)
axF.grid(True, alpha=0.3)

# --- [G] Mean evidence ⟨λ_t⟩ vs t with ±1σ band ---------------------------
axG = axes[3, 0]
mean_p = lam_p.mean(axis=0)
mean_m = lam_m.mean(axis=0)
std_p  = lam_p.std(axis=0)
std_m  = lam_m.std(axis=0)
times_full = np.arange(N_BINS + 1) * DT * 1e3
axG.plot(times_full, mean_p, color="tab:blue",
         label=r"$\langle\lambda_t\rangle$ | NV$^-$")
axG.fill_between(times_full, mean_p - std_p, mean_p + std_p,
                 color="tab:blue", alpha=0.18)
axG.plot(times_full, mean_m, color="tab:red",
         label=r"$\langle\lambda_t\rangle$ | NV$^0$")
axG.fill_between(times_full, mean_m - std_m, mean_m + std_m,
                 color="tab:red", alpha=0.18)
axG.axhline(0.0, color="k", linewidth=0.6)
# Theoretical separation rate (Chernoff-like, ignoring charge dynamics):
# d/dt ⟨λ | +⟩ = γ_+ ln(γ_+/γ_-) - (γ_+ - γ_-)
slope_p = GAMMA_P * np.log(GAMMA_P / GAMMA_M) - (GAMMA_P - GAMMA_M)
slope_m = GAMMA_M * np.log(GAMMA_P / GAMMA_M) - (GAMMA_P - GAMMA_M)
axG.plot(times_full, slope_p * times_full * 1e-3, "k:", linewidth=1.0,
         label=fr"slope $\gamma_+\ln(\gamma_+/\gamma_-)-(\gamma_+-\gamma_-)$ = {slope_p:.1f} /s")
axG.plot(times_full, slope_m * times_full * 1e-3, "k:", linewidth=1.0)
axG.set_xlabel("Time  t  (ms)")
axG.set_ylabel(r"Log-likelihood ratio  $\lambda_t$")
axG.set_title(r"[G] Evidence growth: $\langle\lambda_t\rangle$  $\pm$ 1$\sigma$")
axG.set_xlim(0.0, T_MAX * 1e3)
axG.legend(fontsize=9, loc="lower left")
axG.grid(True, alpha=0.3)

# --- [H] Pareto sweep: T vs ε scatter, colour = log10(odds asymmetry) -----
axH = axes[3, 1]
log_p_plus  = np.log10(1.0 - results[:, 2])    # log10(1 - p_+)
log_p_minus = np.log10(results[:, 3])          # log10(p_-)
sc = axH.scatter(results[:, 0] * 1e3, results[:, 1],
                 c=log_p_minus - log_p_plus,
                 cmap="coolwarm", s=12, alpha=0.7,
                 vmin=-2, vmax=2)
axH.plot(pareto[:, 0] * 1e3, pareto[:, 1], "k-", linewidth=1.5,
         label="Pareto frontier")
axH.set_xscale("linear")
axH.set_yscale("log")
axH.set_xlabel("Average readout time  T  (ms)")
axH.set_ylabel("Error rate  ε")
axH.set_title(r"[H] Adaptive threshold sweep, colour = $\log_{10}\!\frac{p_-}{1-p_+}$")
cbar = plt.colorbar(sc, ax=axH, pad=0.02)
cbar.set_label("asymmetry of cutoffs", fontsize=9)
axH.legend(loc="upper right", fontsize=9)
axH.grid(True, which="both", alpha=0.3)
axH.set_xlim(0.0, T_MAX * 1e3)
axH.set_ylim(8e-3, 0.55)

plt.tight_layout()
diag_path = os.path.join(script_dir, "figure_diagnostics.png")
plt.savefig(diag_path, dpi=150)
plt.close(fig2)
print(f"  saved {diag_path}")


# ---------------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------------
print()
print("=" * 60)
print("PARAMETERS")
print("=" * 60)
print(f"  γ+  = {GAMMA_P:6.1f} Hz   (NV- photon rate)")
print(f"  γ-  = {GAMMA_M:6.1f} Hz   (NV0 photon rate)")
print(f"  Γ+  = {CAP_GAMMA_P:6.2f} Hz   (NV- -> NV0 ionisation)")
print(f"  Γ-  = {CAP_GAMMA_M:6.2f} Hz   (NV0 -> NV- recombination)")
print(f"  δt  = {DT*1e3:.3f} ms,  t_M = {T_MAX*1e3:.1f} ms,  N = {N_BINS}")
print(f"  N_traj per state = {N_TRAJ}")
print()
print("=" * 60)
print("UPDATE MATRICES (verification)")
print("=" * 60)
np.set_printoptions(precision=6, suppress=True)
for n in range(MAX_DN + 1):
    print(f"  M({n}) =")
    for row in M_LIST[n]:
        print(f"    [{row[0]:+.6e}  {row[1]:+.6e}]")
print(f"  Σ_n M(n) =")
for row in sumM:
    print(f"    [{row[0]:+.6e}  {row[1]:+.6e}]")
print(f"  exp(L δt) =")
for row in expL:
    print(f"    [{row[0]:+.6e}  {row[1]:+.6e}]")
print(f"  max |Σ M − exp(L δt)| = {np.abs(sumM - expL).max():.2e}")
print()
print("=" * 60)
print("RESULTS")
print("=" * 60)
print(f"  Photon counting   min ε = {eps_min_pc*100:6.3f}%   "
      f"at tf = {T_pc[np.argmin(eps_pc)]*1e3:.2f} ms  "
      f"(ν = {nu_pc[np.argmin(eps_pc)]})")
print(f"  Nonadaptive MLE   min ε = {eps_min_nmle*100:6.3f}%   "
      f"at tf = {T_nmle[np.argmin(eps_nmle)]*1e3:.2f} ms")
print(f"  Adaptive MLE      min ε = {eps_min_amle*100:6.3f}%   "
      f"at T  = {pareto[np.argmin(pareto[:, 1]), 0]*1e3:.2f} ms")
if not np.isnan(speedup):
    print()
    print(f"  At ε = {ref_eps:.3f}:")
    print(f"    nonadaptive tf = {tf_ref:.2f} ms")
    print(f"    adaptive     T = {T_ref:.2f} ms")
    print(f"    speedup tf/T  = {speedup:.2f}")
print()
print(f"Figure saved to {out_path}")
