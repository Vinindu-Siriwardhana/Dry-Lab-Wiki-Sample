"""
sa_run.py -- global sensitivity analysis: Morris screening and Sobol indices.

HKU iGEM 2026, Dry Lab.  TEAM REFERENCE COPY.

Self-contained implementations of the Saltelli estimator and Morris elementary
effects, so the package needs nothing beyond numpy (SALib gives the same
answers; this is here so the method is readable and the bundle installs
anywhere).

IMPORTANT -- what is actually evaluated
---------------------------------------
A Sobol analysis at N = 4096 over 10 factors needs ~49,000 model evaluations.
Running the full PDE that many times is not practical, so the two outputs are
evaluated with ANALYTIC SURROGATES that preserve the dependence structure of
the full model:

  Output A, scaled dose psi delivered to the macrophage plane
      All of the loaded dose escapes (Da << 1, verified in F3), so the mass
      delivered per unit area is eps_g L_gel C_V0.  It spreads into tissue with
      the exponential profile set by the penetration depth Lambda, and the
      depth-average over the target depth delta is

          C_bar = (M / (eps_t delta)) (1 - exp(-delta/Lambda)),  psi = C_bar/K_D,P

  Output B, closure time t_c
      t_c = (R_w / c) * 1.10 with c = 2 sqrt(D_n r_eff), the 10 % being the
      front-formation transient measured against the full radial solver in F4.
      Samples with r_eff <= 0 never close; they are counted separately as the
      STALL FRACTION and capped at 365 d so the variance decomposition stays
      well defined.

Both surrogates were checked against the full solvers at the central parameter
set.  A ranking computed this way is a statement about the CURRENT PRIORS as
much as about the physics -- the widest priors inevitably rank high, which is
said out loud on the wiki rather than left for a reader to notice.
"""

import numpy as np
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v3model as m
import params as P

DAY = 86400.0
UM = 1e-9
CAP_DAYS = 365.0


# =============================================================================
# priors
# =============================================================================
# (name, low, high)  -- all log-uniform, as in the provenance appendix
PRIORS_A = [
    ("$K_{D,P}$",        1.0 * UM,   100.0 * UM),    # 2 decades, T4
    ("$C_{V,0}$",        20.0 * UM,  800.0 * UM),    # 40x, design variable
    (r"$\delta$",        0.03,       0.15),          # modelling choice
    ("$k_{prot,V}$",     1e-6,       1e-4),          # 2 decades, T4
    ("$L_{gel}$",        300e-4,     2000e-4),       # design variable
    ("$k_{cl,V}$",       1e-7,       1e-5),          # T4
    ("$D_{0,V}$",        2.5e-6,     4.5e-6),        # T2
    (r"$\varepsilon_g$", 0.88,       0.93),          # T3
    (r"$\xi$",           5.0,        50.0),          # T3
]

PRIORS_B = [
    ("$r_0$",            0.3 / DAY,  3.0 / DAY),
    ("$C_e$",            0.1,        10.0),          # C_e / EC50 ratio basis
    ("$K_I^{r}$",        0.05,       1.0),
    ("$EC_{50,F}$",      0.1,        10.0),
    ("$d_0$",            0.05 / DAY, 1.0 / DAY),
    ("$I_e$",            0.05,       1.0),
    (r"$\beta$",         0.1,        0.95),
    (r"$\Delta d$",      0.05 / DAY, 1.5 / DAY),
    ("$E_{max}$",        0.2,        3.0),
    ("$D_n$",            1e-10,      1e-8),
]


def _scale(U, priors):
    """Map a unit hypercube sample to log-uniform priors."""
    lo = np.array([p[1] for p in priors])
    hi = np.array([p[2] for p in priors])
    return np.exp(np.log(lo) + U * (np.log(hi) - np.log(lo)))


# =============================================================================
# models (surrogates -- see module docstring)
# =============================================================================

def model_A(X):
    """Scaled dose psi delivered to the macrophage plane."""
    KD, CV0, delta, kprot, Lgel, kcl, D0V, epsg, xi = X.T
    Deff_t = m.D_eff(D0V, P.A_V14_NM, xi, epsg) * 0.4      # tissue is tighter
    lam = np.sqrt(Deff_t / (kprot + kcl))
    M = epsg * Lgel * CV0                                   # mol per cm^2
    C_bar = M / (P.EPS_TIS * delta) * (1 - np.exp(-delta / lam)) * (lam / delta)
    return C_bar / KD


def model_B(X):
    """Closure time in days; capped at CAP_DAYS when the front stalls."""
    r0, Ce, KIr, EC50, d0, Ie, beta, dd, Emax, Dn = X.T
    rp = r0 * (1 + Emax * Ce / (EC50 + Ce)) * (1 - beta * Ie / (KIr + Ie))
    d = d0 + dd * Ie ** 2 / (0.3 ** 2 + Ie ** 2)
    reff = rp - d
    c = 2 * np.sqrt(Dn * np.maximum(reff, 1e-30))
    tc = np.where(reff > 0, P.R_W_DEFAULT / np.maximum(c, 1e-30) / DAY * 1.10,
                  np.inf)
    return np.minimum(tc, CAP_DAYS), reff <= 0


# =============================================================================
# Saltelli sampling and Sobol estimators
# =============================================================================

def saltelli_sample(N, k, seed=0):
    """A, B and the k cross-sampled AB matrices.  Total N(k+2) rows."""
    rng = np.random.default_rng(seed)
    A = rng.random((N, k))
    B = rng.random((N, k))
    ABs = []
    for i in range(k):
        AB = A.copy()
        AB[:, i] = B[:, i]
        ABs.append(AB)
    return A, B, ABs


def sobol_indices(fA, fB, fABs):
    """
    First-order (Saltelli 2010) and total-order (Jansen 1999) estimators.

        S1_i = mean( fB (fAB_i - fA) ) / Var
        ST_i = mean( (fA - fAB_i)^2 ) / (2 Var)
    """
    var = np.var(np.concatenate([fA, fB]), ddof=1)
    S1, ST = [], []
    for fAB in fABs:
        S1.append(np.mean(fB * (fAB - fA)) / var)
        ST.append(np.mean((fA - fAB) ** 2) / (2 * var))
    return np.array(S1), np.array(ST)


def morris(model, priors, trajectories=200, levels=8, seed=1, scalar=True):
    """Elementary-effect screening; returns mu* for each factor."""
    rng = np.random.default_rng(seed)
    k = len(priors)
    delta = levels / (2 * (levels - 1))
    EE = np.zeros((trajectories, k))
    for t in range(trajectories):
        base = rng.integers(0, levels // 2, k) / (levels - 1)
        order = rng.permutation(k)
        pts = [base.copy()]
        cur = base.copy()
        for i in order:
            cur = cur.copy()
            cur[i] = min(cur[i] + delta, 1.0)
            pts.append(cur)
        Xu = np.array(pts)
        Y = model(_scale(Xu, priors))
        Y = Y[0] if isinstance(Y, tuple) else Y
        for j, i in enumerate(order):
            step = Xu[j + 1, i] - Xu[j, i]
            EE[t, i] = (Y[j + 1] - Y[j]) / (step if step else np.nan)
    return np.nanmean(np.abs(EE), axis=0)


def run(N=4096, seed=0, verbose=True):
    """Run both analyses and return a results dict."""
    out = {}
    for tag, priors, model in (("A", PRIORS_A, model_A),
                               ("B", PRIORS_B, model_B)):
        k = len(priors)
        A, B, ABs = saltelli_sample(N, k, seed)
        ev = lambda U: model(_scale(U, priors))
        rA, rB = ev(A), ev(B)
        stall = None
        if isinstance(rA, tuple):
            fA, sA = rA; fB, sB = rB
            stall = np.concatenate([sA, sB]).mean()
            fABs = [ev(M_)[0] for M_ in ABs]
        else:
            fA, fB = rA, rB
            fABs = [ev(M_) for M_ in ABs]
        # Analyse the RAW output, not its logarithm.  This is the standard
        # convention and it matters: taking logs of a near-multiplicative model
        # makes it near-additive, which would push sum(S1) to ~1 and hide
        # exactly the interaction structure the analysis exists to measure.
        S1, ST = sobol_indices(fA, fB, fABs)
        mu = morris(model, priors, trajectories=200, levels=8, seed=seed + 1)
        if tag == "B":
            mu = np.nan_to_num(mu, nan=0.0, posinf=0.0)
        out[tag] = dict(names=[p[0] for p in priors], S1=S1, ST=ST, mu=mu,
                        stall=stall, n_eval=N * (k + 2), f=np.concatenate([fA, fB]))
        if verbose:
            print(f"\noutput {tag}:  {N*(k+2)} evaluations, sum S1 = {S1.sum():.3f}")
            idx = np.argsort(-ST)
            for i in idx:
                print(f"   {out[tag]['names'][i]:>18}  mu*={mu[i]:10.3g} "
                      f"S1={S1[i]:6.3f}  ST={ST[i]:6.3f}")
            if stall is not None:
                print(f"   stall fraction (r_eff <= 0): {stall:.1%}")
    return out


if __name__ == "__main__":
    run()
