"""
F4 -- Solver verification.  Nothing downstream is worth reading if this panel
is not right, so it appears above the fold rather than in an appendix.

(a) spatial convergence of the method-of-lines release solver against the
    closed-form Crank series, with a slope-2 reference triangle.
(b) numerical vs analytic Fisher-KPP pulled-front speed, parity plot.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import style as S, params as P, v3model as m

S.apply()
DAY = 86400.0
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.9))

# ---- (a) spatial convergence ----------------------------------------------
D = m.D_eff(P.D0_V14, P.A_V14_NM, P.XI_NM, P.EPS_GEL)
L = P.L_GEL_ALT          # 500 um, the benchmark geometry of verify.py
times = np.array([30., 60., 120., 300., 600., 900.])
nxs = np.array([100, 200, 400, 800, 1600])
errs = []
for nx in nxs:
    sol = m.ReleaseSolver(D, L, nx)
    Cn = sol.solve(1.0, times, dt=0.4 * (100 / nx))
    Ce = m.crank_profile(sol.x, times, 1.0, D, L)
    errs.append(np.max(np.abs(Cn - Ce)))
errs = np.array(errs)

ax1.loglog(nxs, errs, "o-", color=S.C_V14, ms=6)
ref = errs[0] * (nxs[0] / nxs) ** 2
ax1.loglog(nxs, ref, ls=":", color=S.C_GREY, lw=1.4)
ax1.plot([400, 800, 800, 400], [ref[2], ref[2], ref[3], ref[2]],
         color=S.C_DERIVED, lw=1.1)
ax1.text(830, np.sqrt(ref[2] * ref[3]), "2", color=S.C_DERIVED, fontsize=9,
         weight="bold", va="center")
ax1.text(430, ref[2] * 1.25, "1", color=S.C_DERIVED, fontsize=9, weight="bold")
for nx, e in zip(nxs, errs):
    pass
ax1.text(0.03, 0.06, f"observed order {np.log2(errs[-2] / errs[-1]):.2f}\nfinal error {errs[-1]:.2e}",
         transform=ax1.transAxes, fontsize=9, weight="bold", color=S.C_V14)
ax1.set_xlabel("grid points $n_x$")
ax1.set_ylabel(r"max $|C_{num}-C_{exact}|\,/\,C_0$")
ax1.set_title("(a)  release solver vs closed form", loc="left")

# ---- (b) front speed parity ------------------------------------------------
cases = [(1e-9, 1.0, 0.0, 60), (1e-8, 1.0, 0.0, 20),
         (1e-9, 1.0, 0.3, 70), (1e-9, 2.0, 0.5, 50)]
an, nu, labs = [], [], []
for Dn, r0, d, tend in cases:
    reff = (r0 - d) / DAY
    c_an = m.front_speed(Dn, reff)
    sol = m.FisherKPP(Dn, R_inf=1.6, nr=1600, geometry="cartesian")
    res = sol.run(1.4, r0 / DAY, d / DAY, t_end=tend * DAY, record_every=100)
    t, rf = res["t"], res["r_front"]
    msk = (rf > 0.35) & (rf < 1.05)
    c_nu = -np.polyfit(t[msk], rf[msk], 1)[0]
    an.append(c_an * DAY * 1e4); nu.append(c_nu * DAY * 1e4)
    labs.append(f"$D_n$={Dn:.0e}, $r_0$={r0}, $d$={d}")

an, nu = np.array(an), np.array(nu)
lim = [0, 1.12 * max(an.max(), nu.max())]
ax2.plot(lim, lim, color=S.C_GREY, ls=":", lw=1.4, label="1 : 1")
ax2.plot(an, nu, "o", color=S.C_FGF, ms=7)
for a, n_, lb in zip(an, nu, labs):
    ax2.annotate(f"{abs(n_-a)/a:.1%}", (a, n_), textcoords="offset points",
                 xytext=(8, -11), fontsize=8, color=S.C_DERIVED)
ax2.set_xlim(lim); ax2.set_ylim(lim)
ax2.set_xlabel(r"analytic  $c = 2\sqrt{D_n r_{eff}}$   (µm d$^{-1}$)")
ax2.set_ylabel("numerical front speed  (µm d$^{-1}$)")
ax2.set_title("(b)  Fisher–KPP front speed", loc="left")
ax2.legend(loc="upper left")
rel = abs(nu - an) / an
ax2.text(0.97, 0.05, f"errors {rel.min():.1%}–{rel.max():.1%}\nall below 1:1, as expected\n"
         f"for a discretised pulled front", transform=ax2.transAxes, ha="right",
         fontsize=8.5, weight="bold", color=S.C_FGF, linespacing=1.4)

S.save(fig, "fig04_verification")
print("  convergence errors:", " ".join(f"{e:.2e}" for e in errs))
print("  front-speed errors:", " ".join(f"{abs(n_-a)/a:.2%}" for a, n_ in zip(an, nu)))
