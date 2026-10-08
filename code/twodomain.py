"""
twodomain.py -- monolithic gel + tissue delivery solver.

HKU iGEM 2026, Dry Lab.  TEAM REFERENCE COPY.

This is the solver behind the retention-envelope and tissue-profile results.
It solves the two coupled domains together, with the partition and flux-matching
conditions imposed as INTERNAL interface conditions rather than as prescribed
boundary fluxes:

    gel      0 <= x <= L_gel     R_V dC/dt = D_g d2C/dx2 - k_deg C
    tissue  -L_tis <= x <= 0         dC/dt = D_t d2C/dx2 - k_prot C

    dC_g/dx|_{L_gel} = 0                             impermeable backing
    C_g(0) = K_p C_t(0)                              partition
    eps_g D_g dC_g/dx|_{0+} = eps_t D_t dC_t/dx|_{0-}  flux matching (per total area)
    C_t(-L_tis) = 0                                  systemic sink at dermal base

Two modelling points that are easy to get wrong and are enforced here:

1.  The deep boundary is a Dirichlet SINK and there is NO distributed clearance
    term as well.  Both represent drug leaving into the circulation; V2 imposed
    both and suppressed the tissue concentration twice.

2.  With sorptive retention the loaded concentration is the TOTAL (free+bound),
    so the free initial concentration is C_total/R_V.  Treating the loading as
    free silently loads R_V times more drug than solubility allows and
    over-predicts sustained release by orders of magnitude.  Use
    free_from_total() -- it exists to make that error hard to make.

Interface treatment: the interface is massless, so C_t(0) is determined
algebraically from its two neighbours,

    C_t0 = (a C_g1 + b C_t1) / (b + a K_p),
    a = eps_g D_g / h_g,   b = eps_t D_t / h_t,

and C_g0 = K_p C_t0.  Both eliminated nodes are substituted into the equations
for the first interior node on each side, which keeps the system sparse.
"""

import numpy as np
from scipy.sparse import lil_matrix, csc_matrix, identity
from scipy.sparse.linalg import splu

UM = 1e-9     # 1 micromolar in mol cm^-3
DAY = 86400.0


def free_from_total(C_total, R_V):
    """
    Free (diffusible) loaded concentration given the TOTAL dissolved loading.

        C_free,0 = C_total,0 / R_V

    Retention is paid for in driving concentration, which is why R_V has an
    optimum rather than a maximum.
    """
    return C_total / R_V


class TwoDomain:
    def __init__(self, D_gel, D_tis, L_gel, L_tis, nx_gel=200, nx_tis=400,
                 eps_gel=0.91, eps_tis=0.20, K_p=1.0,
                 k_deg_gel=0.0, k_prot=1e-5, R_V=1.0):
        self.__dict__.update(locals()); del self.self

        self.hg = L_gel / nx_gel
        self.ht = L_tis / nx_tis
        self.xg = np.linspace(0, L_gel, nx_gel + 1)       # gel nodes, x>=0
        self.xt = -np.linspace(0, L_tis, nx_tis + 1)      # tissue nodes, x<=0

        # unknowns: gel 1..nx_gel  (ng values), tissue 1..nx_tis-1 (nt values)
        self.ng = nx_gel
        self.nt = nx_tis - 1
        self.N = self.ng + self.nt

        a = eps_gel * D_gel / self.hg
        b = eps_tis * D_tis / self.ht
        self.denom = b + a * K_p
        self.ca = a / self.denom      # C_t0 = ca*C_g1 + cb*C_t1
        self.cb = b / self.denom
        self.A = self._operator()

    # -- indexing helpers ----------------------------------------------------
    def ig(self, j):  return j - 1              # gel node j = 1..ng
    def it(self, i):  return self.ng + i - 1    # tissue node i = 1..nt

    def _operator(self):
        ng, nt, hg, ht = self.ng, self.nt, self.hg, self.ht
        Dg, Dt, Kp = self.D_gel, self.D_tis, self.K_p
        A = lil_matrix((self.N, self.N))

        # ---- gel interior nodes -------------------------------------------
        for j in range(1, ng + 1):
            r = self.ig(j)
            c = Dg / hg ** 2
            A[r, r] += -2 * c - self.k_deg_gel
            if j == 1:
                # neighbour at j=0 is C_g0 = K_p*C_t0 = K_p*(ca*C_g1 + cb*C_t1)
                A[r, self.ig(1)] += c * Kp * self.ca
                A[r, self.it(1)] += c * Kp * self.cb
            else:
                A[r, self.ig(j - 1)] += c
            if j == ng:
                A[r, self.ig(ng - 1)] += c      # ghost reflection, zero flux
            else:
                A[r, self.ig(j + 1)] += c
            if j == ng:
                A[r, r] += 0.0
        # retardation: apparent diffusivity D/R_V in the gel
        for j in range(1, ng + 1):
            r = self.ig(j)
            A[r, :] = A[r, :] / self.R_V

        # ---- tissue interior nodes ----------------------------------------
        for i in range(1, nt + 1):
            r = self.it(i)
            c = Dt / ht ** 2
            A[r, r] += -2 * c - self.k_prot
            if i == 1:
                A[r, self.ig(1)] += c * self.ca     # neighbour is C_t0
                A[r, self.it(1)] += c * self.cb
            else:
                A[r, self.it(i - 1)] += c
            if i < nt:
                A[r, self.it(i + 1)] += c
            # i == nt: neighbour is the Dirichlet node C_t(-L_tis)=0, contributes 0
        return csc_matrix(A)

    # -- time integration ----------------------------------------------------
    def run(self, C_free_0, t_end, segments=None, record_every=1):
        """
        Integrate from a uniformly loaded gel into clean tissue.

        `segments` is a list of (dt, n_steps) pairs; the default ramps the step
        from 0.5 s through 600 s so the fast release transient (minutes) and the
        slow tissue phase (days) are both resolved without a stiff solve.

        Returns dict(t, C_gel, C_tis, C_t0, depth_avg) with depth_avg the free
        concentration averaged over the top `delta` of tissue (set by
        set_depth_window, default 1 mm).
        """
        if segments is None:
            segments = self._default_segments(t_end)

        C = np.zeros(self.N)
        C[:self.ng] = C_free_0                      # gel loaded, tissue clean

        I = identity(self.N, format="csc")
        ts, recs = [0.0], [C.copy()]
        t = 0.0
        for dt, nsteps in segments:
            lu = splu(csc_matrix(I - dt * self.A))  # implicit Euler, A fixed
            for s in range(nsteps):
                C = lu.solve(C)
                t += dt
                if (s + 1) % record_every == 0 or s == nsteps - 1:
                    ts.append(t); recs.append(C.copy())
                if t >= t_end:
                    break
            if t >= t_end:
                break

        R = np.array(recs); ts = np.array(ts)
        Cg1 = R[:, self.ig(1)]
        Ct1 = R[:, self.it(1)]
        Ct0 = self.ca * Cg1 + self.cb * Ct1
        return dict(t=ts, C=R, C_t0=Ct0,
                    C_gel=np.column_stack([self.K_p * Ct0, R[:, :self.ng]]),
                    C_tis=np.column_stack([Ct0, R[:, self.ng:],
                                           np.zeros(len(ts))]))

    def _default_segments(self, t_end):
        """Geometric ramp: resolve minutes first, then coast in 10-min steps."""
        segs = [(0.5, 120), (5.0, 108), (60.0, 54), (600.0, 0)]
        used = sum(dt * n for dt, n in segs)
        remaining = max(t_end - used, 0.0)
        segs[-1] = (600.0, int(np.ceil(remaining / 600.0)))
        return segs

    def depth_average(self, res, delta):
        """
        Mean FREE concentration over the top `delta` cm of tissue.

        This is the pharmacologically conservative metric: the macrophage target
        sits in a thin surface band where the concentration is always ABOVE this
        average, so days-above-target computed from it is a lower bound.
        """
        x = np.concatenate([[0.0], -self.xt[1:]])          # depths, 0..L_tis
        prof = res["C_tis"]                                 # (nt_times, nodes)
        mask = x <= delta
        xm = x[mask]
        return np.trapezoid(prof[:, mask], xm, axis=1) / delta \
            if hasattr(np, "trapezoid") else \
            np.trapz(prof[:, mask], xm, axis=1) / delta

    def days_above(self, res, delta, target):
        """Days for which the depth-averaged free concentration exceeds target."""
        avg = self.depth_average(res, delta)
        t = res["t"]
        above = avg >= target
        if not np.any(above):
            return 0.0, float(avg.max())
        # total time above threshold, trapezoidal on the recorded grid
        dur = np.sum(np.diff(t) * (above[1:] | above[:-1])) / DAY
        return float(dur), float(avg.max())


# =============================================================================
# verification -- run as a script
# =============================================================================

def verify(verbose=True):
    """
    Four checks, all of which must pass before any result is taken from this
    solver.  Returns a dict of the measured errors.
    """
    out = {}

    # 1. tissue operator vs analytic steady state, Dirichlet at both ends.
    #    Built standalone so the test cannot be passed vacuously by everything
    #    decaying to zero (an earlier version of this check did exactly that).
    D, L, k = 1.0e-6, 0.4, 1e-4
    lam = np.sqrt(D / k)
    errs = []
    for nx in (200, 400, 800):
        h = L / nx
        n = nx - 1                       # interior nodes
        main = (-2 * D / h ** 2 - k) * np.ones(n)
        off = (D / h ** 2) * np.ones(n - 1)
        A = np.diag(main) + np.diag(off, 1) + np.diag(off, -1)
        b = np.zeros(n); b[0] = -D / h ** 2 * 1.0      # C(0) = 1
        num = np.linalg.solve(A, b)
        x = np.linspace(0, L, nx + 1)[1:-1]
        ana = np.sinh((L - x) / lam) / np.sinh(L / lam)
        errs.append(np.max(np.abs(num - ana)))
    out["steady_state_err"] = errs
    out["steady_order"] = (np.log2(errs[0] / errs[1]), np.log2(errs[1] / errs[2]))

    # 2. partition recovery: no decay, sealed both ends -> C_gel/C_tis -> K_p
    Kp = 2.5
    td = TwoDomain(D_gel=3e-6, D_tis=1e-6, L_gel=500e-4, L_tis=500e-4,
                   nx_gel=100, nx_tis=100, eps_gel=1.0, eps_tis=1.0,
                   K_p=Kp, k_prot=0.0, k_deg_gel=0.0)
    td.A = td.A.tolil()
    td.A[td.it(td.nt), :] = 0.0                     # seal the deep end
    r = td.it(td.nt); c = td.D_tis / td.ht ** 2
    td.A[r, r] = -c; td.A[r, td.it(td.nt - 1)] = c
    td.A = csc_matrix(td.A)
    res = td.run(1.0, t_end=40 * 3600, segments=[(20.0, 7200)])
    ratio = res["C_gel"][-1].mean() / res["C_tis"][-1][:-1].mean()
    out["partition_err"] = abs(ratio - Kp) / Kp

    # 3. pure gel-side first-order decay matches the analytic exponential
    td = TwoDomain(D_gel=1e-12, D_tis=1e-12, L_gel=500e-4, L_tis=500e-4,
                   nx_gel=20, nx_tis=20, K_p=1.0, k_prot=0.0, k_deg_gel=1e-4)
    res = td.run(1.0, t_end=3600, segments=[(0.5, 7200)])
    out["decay_err"] = abs(res["C_gel"][-1][5] - np.exp(-1e-4 * 3600))

    # 4. mass budget: with no decay anywhere, everything loaded must leave
    #    through the deep sink.  Checks the interface flux is conservative.
    td = TwoDomain(D_gel=3.25e-6, D_tis=1.0e-6, L_gel=500e-4, L_tis=0.4,
                   nx_gel=100, nx_tis=200, eps_gel=0.91, eps_tis=0.20,
                   K_p=1.0, k_prot=0.0, k_deg_gel=0.0)
    res = td.run(1.0, t_end=20 * DAY)
    xg = td.xg; xt = -td.xt
    m0 = td.eps_gel * td.L_gel * 1.0
    mg = td.eps_gel * np.trapezoid(res["C_gel"][-1], xg)
    mt = td.eps_tis * np.trapezoid(res["C_tis"][-1], xt)
    out["mass_remaining_frac"] = (mg + mt) / m0

    if verbose:
        print("two-domain solver verification")
        print(f"  steady state   err {errs[-1]:.2e}  observed order "
              f"{out['steady_order'][1]:.2f}")
        print(f"  partition K_p  rel err {out['partition_err']:.2%}")
        print(f"  gel decay      abs err {out['decay_err']:.2e}")
        print(f"  mass budget    fraction still in domain after 20 d: "
              f"{out['mass_remaining_frac']:.3%}")
    return out


if __name__ == "__main__":
    verify()
