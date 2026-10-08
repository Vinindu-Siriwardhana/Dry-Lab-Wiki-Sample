"""
v3model.py -- reference implementation of Mathematical Model V3.

HKU iGEM 2026, Dry Lab.  TEAM REFERENCE COPY: this module exists so that every
number and every curve on the Modelling wiki page can be regenerated from
source by anyone on the team.  It is deliberately dependency-light
(numpy + scipy only) and deliberately un-clever: readability beats speed here,
because the point is that a teammate can check the maths against the equation
reference document line by line.

Contents
--------
Transport coefficients   amsden_H, D_eff
Model 1 (release)        crank_profile, crank_cumulative, release_time,
                         ReleaseSolver (method of lines, perfect sink)
Model 3-4 (signalling)   theta_inhib, theta_lps, p65, psi_50, dose_required,
                         dose_required_v2
Model 5 (NO readout)     inos_response, nitrite_reduced
Model 6 (closure)        FisherKPP (radial + Cartesian), front_speed,
                         coupling_rp, coupling_d
Model 7 (design)         loading_required

Units are CGS throughout: cm, s, mol cm^-3.  Concentrations are INTERSTITIAL
(mol per cm^3 of pore fluid), per Section 1.3 of the equation reference --
porosity never appears inside a PDE.

Convenience: 1 uM = 1e-9 mol/cm^3.
"""

import numpy as np
from scipy.sparse import diags, identity, csc_matrix
from scipy.sparse.linalg import splu

UM = 1e-9          # 1 micromolar, in mol cm^-3
DAY = 86400.0      # s
MIN = 60.0


# =============================================================================
# Transport coefficients  (equation reference, Section 1)
# =============================================================================

def amsden_H(a_solute_nm, xi_nm, a_fibre_nm=0.6):
    """
    Amsden obstruction-scaling steric hindrance factor.

        H = exp[ -pi ( (a + a_f) / (xi + 2 a_f) )^2 ],   0 < H <= 1

    a_solute_nm : solute hydrodynamic radius  (0.87 nm V14, 2.34 nm FGF2-G3)
    xi_nm       : hydrogel mesh size
    a_fibre_nm  : gelatin chain radius, ~0.6 nm
    """
    a = np.asarray(a_solute_nm, dtype=float)
    xi = np.asarray(xi_nm, dtype=float)
    return np.exp(-np.pi * ((a + a_fibre_nm) / (xi + 2 * a_fibre_nm)) ** 2)


def D_eff(D0, a_solute_nm, xi_nm, porosity=0.91, a_fibre_nm=0.6):
    """
    Pore (effective) diffusivity, D_eff = D0 * H / tau, with the Bruggeman
    closure tau = porosity^(-1/2).

    NOTE the porosity does *not* multiply the result again: it cancels out of
    the transient term (Eq. basis).  Applying it twice was the V2 error.
    """
    H = amsden_H(a_solute_nm, xi_nm, a_fibre_nm)
    tau = porosity ** -0.5
    return D0 * H / tau


# =============================================================================
# Model 1 -- release from the hydrogel
# =============================================================================

def crank_profile(x, t, C0, D, L, k_deg=0.0, n_terms=400):
    """
    Closed form for Fickian desorption from a plane sheet sealed on one face,
    perfect sink at x = 0 (Crank, Mathematics of Diffusion, 4.3.2):

        C(x,t) = C0 e^{-kt} sum_m 4/((2m+1)pi)
                 sin[(2m+1)pi x /2L] exp[-D(2m+1)^2 pi^2 t /4L^2]

    This is the designated regression test for the numerical solver.
    """
    x = np.atleast_1d(np.asarray(x, float))
    t = np.atleast_1d(np.asarray(t, float))
    m = np.arange(n_terms)
    lam = (2 * m + 1) * np.pi / (2 * L)
    coef = 4.0 / ((2 * m + 1) * np.pi)
    # shape (nt, nx)
    S = np.einsum("m,xm,tm->tx",
                  coef,
                  np.sin(np.outer(x, lam)),
                  np.exp(-D * np.outer(t, lam ** 2)))
    return C0 * np.exp(-k_deg * t)[:, None] * S


def crank_cumulative(t, D, L, n_terms=2000):
    """Cumulative fractional release M_t/M_inf (no degradation)."""
    t = np.atleast_1d(np.asarray(t, float))
    m = np.arange(n_terms)
    a = (2 * m + 1) ** 2 * np.pi ** 2
    coef = 8.0 / ((2 * m + 1) ** 2 * np.pi ** 2)
    S = np.einsum("m,tm->t", coef, np.exp(-D * np.outer(t, a) / (4 * L ** 2)))
    return 1.0 - S


def release_time(fraction, D, L, hi=None):
    """Time at which cumulative fractional release reaches `fraction`."""
    from scipy.optimize import brentq
    tau_late = 4 * L ** 2 / (np.pi ** 2 * D)
    hi = hi or tau_late * 50
    f = lambda t: crank_cumulative(np.array([t]), D, L)[0] - fraction
    return brentq(f, 1e-6, hi)


def tau_late(D, L):
    """Late-time decay constant, 4L^2 / (pi^2 D)."""
    return 4 * L ** 2 / (np.pi ** 2 * D)


def damkohler(k_deg, L, D):
    """Gel Damkohler number; << 1 means release beats in-gel degradation."""
    return k_deg * L ** 2 / D


def penetration_depth(D_tissue, k_loss):
    """Tissue penetration depth Lambda = sqrt(D_t / k)."""
    return np.sqrt(D_tissue / k_loss)


class ReleaseSolver:
    """
    Method-of-lines solver for the single-domain gel release problem with a
    perfect sink, used to verify second-order spatial convergence against
    crank_profile().

        dC/dt = D d2C/dx2 - k C,   C(0,t)=0,  dC/dx|_L = 0,  C(x,0)=C0

    Time integration is Crank-Nicolson with a fixed step, so the observed error
    is dominated by the spatial discretisation (that is the point: we are
    measuring the spatial order).
    """

    def __init__(self, D, L, nx, k_deg=0.0):
        self.D, self.L, self.nx, self.k = D, L, nx, k_deg
        self.h = L / nx
        self.x = np.linspace(0, L, nx + 1)
        # unknowns: nodes 1..nx  (node 0 pinned at 0 by the sink)
        n = nx
        main = -2.0 * np.ones(n)
        main[-1] = -2.0                  # mirror node handles the zero-flux face
        lower = np.ones(n - 1)
        upper = np.ones(n - 1)
        A = diags([lower, main, upper], [-1, 0, 1], format="lil")
        A[n - 1, n - 2] = 2.0            # ghost-node reflection at x = L
        A = csc_matrix(A) * (D / self.h ** 2) - k_deg * identity(n, format="csc")
        self.A = csc_matrix(A)

    def solve(self, C0, times, dt, rannacher=4):
        """
        Return C at the requested `times`, shape (len(times), nx+1).

        The initial data are discontinuous at the corner (x,t)=(0,0) -- loaded
        gel against a perfect sink -- so dC/dx|_0 diverges as t^-1/2 and plain
        Crank-Nicolson rings badly on fine grids.  We therefore use RANNACHER
        STARTUP: the first `rannacher` steps are taken with fully implicit
        (backward) Euler at half the step, which damps the high-frequency modes,
        after which CN resumes and recovers second order.  Without this the
        nx=1600 error is ~1e-3 instead of ~1e-7.
        """
        n = self.nx
        I = identity(n, format="csc")
        lu_cn = splu(csc_matrix(I - 0.5 * dt * self.A))
        Rm = csc_matrix(I + 0.5 * dt * self.A)
        lu_be = splu(csc_matrix(I - 0.5 * dt * self.A))   # BE at half step

        C = np.full(n, float(C0))
        out, t, done = [], 0.0, 0
        for tt in np.atleast_1d(times):
            nsteps = int(round((tt - t) / dt))
            for _ in range(nsteps):
                if done < rannacher:
                    C = lu_be.solve(lu_be.solve(C))   # two BE half-steps
                else:
                    C = lu_cn.solve(Rm @ C)
                done += 1
            t += nsteps * dt
            out.append(np.concatenate([[0.0], C]))
        return np.array(out)

    def retained_mass(self, C):
        """Trapezoidal integral of C over the gel, per unit area."""
        return np.trapezoid(C, self.x) if hasattr(np, "trapezoid") \
            else np.trapz(C, self.x)


# =============================================================================
# Models 3 and 4 -- competitive occupancy, signalling, dose inversion
# =============================================================================

def theta_inhib(C_V, KD_P, L_lps, KD_L):
    """Fraction of MD2 occupied by V14. C_V sits OUTSIDE the KD_P bracket."""
    return C_V / (KD_P * (1 + L_lps / KD_L) + C_V)


def theta_lps(C_V, KD_P, L_lps, KD_L):
    """Fraction of MD2 occupied by LPS -- the pro-inflammatory species."""
    a = L_lps / KD_L
    return a / (1 + a + C_V / KD_P)


def theta_lps_scaled(psi, a):
    """Non-dimensional form: only psi = C_V/KD_P and a = [L]/KD_L appear."""
    return a / (1 + a + psi)


def suppression_ratio(psi, a):
    """theta_LPS(psi)/theta_LPS(0) = (1+a)/(1+a+psi).  No fitted constants."""
    return (1 + a) / (1 + a + psi)


def psi_50(a):
    """Exact scaled dose for 50% suppression of active receptor: psi50 = 1 + a."""
    return 1.0 + a


def p65_scaled(psi, a, gamma, h):
    """Normalised nuclear p65, pi = (a gamma)^h / ((a gamma)^h + (1+a+psi)^h)."""
    return (a * gamma) ** h / ((a * gamma) ** h + (1 + a + psi) ** h)


def dose_required(KD_P, L_lps, KD_L, TLR_tot, TLR_tgt):
    """
    V3 corrected required-dose inversion:

        [P]_req = KD_P [ ([L]/KD_L)(TLR_tot/TLR_tgt - 1) - 1 ]

    Asymptotically LINEAR in [L].  Returns NaN where the feasibility condition
    fails (a negative dose is not a dose).
    """
    a = L_lps / KD_L
    val = KD_P * (a * (TLR_tot / TLR_tgt - 1.0) - 1.0)
    return np.where(val > 0, val, np.nan)


def dose_required_v2(KD_P, L_lps, KD_L, TLR_tot, TLR_tgt):
    """
    The V2 equation, retained only so the error can be plotted:

        [P]_req^V2 = KD_P (1 + [L]/KD_L) ( TLR_tot [L] / (KD_L TLR_tgt) - 1 )

    It factors out KD_app, which is not legitimate when inverting for the dose,
    and is QUADRATIC in [L].
    """
    a = L_lps / KD_L
    return KD_P * (1 + a) * (TLR_tot * a / TLR_tgt - 1.0)


# =============================================================================
# Model 5 -- iNOS induction and the NO readout
# =============================================================================

def inos_response(t, p65_rel, N_cell, q_max, K_iNOS_rel=0.5, p=2.0,
                  k_turn=1 / (4 * 3600), k_scav=1 / 8.0, phi=0.8):
    """
    Integrate the iNOS turnover ODE and accumulate nitrite.

        dR/dt   = k_turn [ N q_max p65^p/(K^p + p65^p) - R ]
        [NO]    ~= R / k_scav                       (quasi-steady, free NO)
        d[NO2-]/dt = phi R                          (what Griess measures)

    p65_rel is the NORMALISED nuclear p65 (0..1), so K_iNOS_rel is on the same
    scale.  Returns (R_rel, free_NO, nitrite) on the grid t.
    """
    drive = N_cell * q_max * p65_rel ** p / (K_iNOS_rel ** p + p65_rel ** p)
    R = np.zeros_like(t, dtype=float)
    NO2 = np.zeros_like(t, dtype=float)
    for i in range(1, len(t)):
        dt = t[i] - t[i - 1]
        R[i] = R[i - 1] + dt * k_turn * (drive - R[i - 1])          # explicit
        NO2[i] = NO2[i - 1] + dt * phi * 0.5 * (R[i] + R[i - 1])    # trapezoid
    return R, R / k_scav, NO2


def nitrite_reduced(C_V, IC50, n_app):
    """
    The two-parameter reduced form that IS identifiable from the Griess curve:

        [NO2-](C)/[NO2-](0) = 1 / (1 + (C/IC50)^n)

    The mechanistic chain carries >= 14 parameters against a ~3-d.o.f. sigmoid,
    so only these two may be reported as fitted.
    """
    return 1.0 / (1.0 + (np.asarray(C_V, float) / IC50) ** n_app)


# =============================================================================
# Model 6 -- fibroblast-driven closure (Fisher-KPP)
# =============================================================================

def coupling_rp(r0, C_e, EC50_F, E_max, I_e, beta, K_I_r):
    """r_p = r0 [1 + Emax C/(EC50+C)] [1 - beta I/(K_I + I)]   (Tier-4 closure)"""
    return r0 * (1 + E_max * C_e / (EC50_F + C_e)) * (1 - beta * I_e / (K_I_r + I_e))


def coupling_d(d0, dd, I_e, K_I_d, s=2.0):
    """d = d0 + dd I^s/(K_I^s + I^s)   (Tier-4 closure)"""
    return d0 + dd * I_e ** s / (K_I_d ** s + I_e ** s)


def front_speed(D_n, r_eff):
    """Pulled-front speed c = 2 sqrt(D_n r_eff); zero if r_eff <= 0."""
    r_eff = np.asarray(r_eff, float)
    return np.where(r_eff > 0, 2 * np.sqrt(D_n * np.maximum(r_eff, 0)), 0.0)


def front_width(D_n, r_eff):
    """Front width w_f ~ sqrt(D_n / r_eff).  c*w_f ~ D_n, c/w_f ~ r_eff."""
    return np.sqrt(D_n / np.maximum(r_eff, 1e-30))


class FisherKPP:
    """
    Radial (or Cartesian) Fisher-KPP front with a death term:

        dn/dt = (1/r) d/dr ( r D_n dn/dr ) + r_p n (1 - n/K) - d n

    Explicit in time (the diffusion number is small: D_n ~ 1e-9 cm^2/s against
    dr ~ 10 um gives dr^2/D ~ 1000 s, so dt = 200 s is comfortably stable and
    exactly reproduces the analytic front speed to a few per cent).

    r_p and d may be constants or callables of t (seconds) -- the latter is how
    a transient therapeutic pulse is driven in Section "durability".
    """

    def __init__(self, D_n, K=1.0, R_inf=1.0, nr=1000, geometry="radial"):
        self.D_n, self.K, self.R_inf, self.nr = D_n, K, R_inf, nr
        self.geometry = geometry
        self.r = np.linspace(0.0, R_inf, nr + 1)
        self.dr = self.r[1] - self.r[0]

    def initial(self, R_w, w=75e-4):
        """Smoothed wound void: n = K/2 [1 + tanh((r - R_w)/w)]."""
        return 0.5 * self.K * (1 + np.tanh((self.r - R_w) / w))

    def _lap(self, n):
        d2 = np.zeros_like(n)
        d2[1:-1] = (n[2:] - 2 * n[1:-1] + n[:-2]) / self.dr ** 2
        if self.geometry == "radial":
            # curvature term (D/r) dn/dr; at r=0 the symmetric limit is 2*d2
            d2[1:-1] += (n[2:] - n[:-2]) / (2 * self.dr * self.r[1:-1])
            d2[0] = 2 * (n[1] - n[0]) / self.dr ** 2 * 2
        else:
            d2[0] = 2 * (n[1] - n[0]) / self.dr ** 2
        d2[-1] = 0.0
        return self.D_n * d2

    def run(self, R_w, r_p, d, t_end, dt=None, w=75e-4, record_every=50,
            fixed_threshold=None):
        """
        Integrate to t_end.  Returns dict with t, r_front, n_final, closed, t_c.

        The front is tracked at n = 0.5 n_ss with n_ss = K(1 - d/r_p), the
        steady state BEHIND the front -- not 0.5K, which is unreachable once
        d/r_p >= 0.5, i.e. exactly in the chronic regime the death term exists
        to represent.

        `fixed_threshold` pins the front-tracking level instead of recomputing
        n_ss each step.  Use it whenever r_p and d vary in time: once the drug
        effect fades and r_eff goes negative, n_ss collapses toward zero and a
        moving threshold would report a spurious front.  For the durability runs
        we pin it at 0.5 n_ss evaluated at the THERAPEUTIC rates, i.e. the state
        the tissue was in while it was healing.
        """
        if dt is None:
            # explicit stability: dt <= dr^2/(2 D).  A safety factor of 0.2 also
            # keeps the reaction step accurate when r_p is large.
            dt = 0.2 * self.dr ** 2 / self.D_n
        n = self.initial(R_w, w)
        rp_f = r_p if callable(r_p) else (lambda t: r_p)
        d_f = d if callable(d) else (lambda t: d)

        ts, fronts = [], []
        t, step, t_c = 0.0, 0, np.inf
        nsteps = int(t_end / dt)
        while step < nsteps:
            rp, dd = rp_f(t), d_f(t)
            n = n + dt * (self._lap(n) + rp * n * (1 - n / self.K) - dd * n)
            np.clip(n, 0.0, None, out=n)
            n[-1] = self.K
            t += dt
            step += 1
            if step % record_every == 0 or step == nsteps:
                if fixed_threshold is not None:
                    thr = fixed_threshold
                else:
                    nss = self.K * max(1 - dd / rp, 1e-6) if rp > 0 else 1e-6
                    thr = 0.5 * nss
                rf = self._front(n, thr)
                ts.append(t); fronts.append(rf)
                if rf <= 0 and not np.isfinite(t_c):
                    t_c = t
                if rf <= 0 and t_c == np.inf:
                    t_c = t
        ts = np.array(ts); fronts = np.array(fronts)
        closed = bool(np.any(fronts <= 1e-9))
        t_c = ts[np.argmax(fronts <= 1e-9)] if closed else np.inf
        return dict(t=ts, r_front=fronts, n=n, closed=closed, t_c=t_c, dt=dt)

    def _front(self, n, thresh):
        """Linear interpolation for the outermost crossing of `thresh`."""
        above = n >= thresh
        idx = np.argmax(above)       # first index where n >= thresh
        if above[0]:
            return 0.0
        if idx == 0:
            return self.R_inf
        n0, n1 = n[idx - 1], n[idx]
        if n1 == n0:
            return self.r[idx]
        frac = (thresh - n0) / (n1 - n0)
        return self.r[idx - 1] + frac * self.dr


# =============================================================================
# Model 7 -- design inversion
# =============================================================================

def loading_required(C_target, delta, L_gel, f_rel=0.5,
                     eps_gel=0.9, eps_tis=1.0):
    """
    Lower bound on the FREE loading needed to hold C_target at depth delta:

        C_V0 >~ C_target eps_t delta / (f_rel eps_g L_gel)

    Neglects proteolysis and clearance, so the real requirement is higher.
    With sorptive retention the TOTAL that must dissolve is R_V times this
    (see twodomain.free_from_total).
    """
    return C_target * eps_tis * delta / (f_rel * eps_gel * L_gel)
