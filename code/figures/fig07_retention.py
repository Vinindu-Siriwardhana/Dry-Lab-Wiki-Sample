"""
F7 -- Retention has an OPTIMUM, not a maximum.  The headline design result.

Holding the TOTAL loading at the solubility ceiling, raising the retardation
factor R_V buys residence time and spends driving concentration in exact
proportion (C_free,0 = C_total,0 / R_V).  The result is a sharp optimum near
R_V ~ 10 and a collapse to ZERO coverage beyond R_V ~ 30.

The duration is plotted as a BAND, because it is sensitive to two parameters
the equation reference quotes only as ranges: the depth to the systemic sink
(L_tis = 2-4 mm) and the tissue diffusivity (D_t = 0.7-1.8e-6 cm2/s).  The
LOCATION of the optimum and the collapse are robust to both; the absolute
number of days is not, and is not claimed to be.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import style as S, params as P, v3model as m
from twodomain import TwoDomain, free_from_total, UM, DAY

S.apply()
Dg = m.D_eff(P.D0_V14, P.A_V14_NM, P.XI_NM, P.EPS_GEL)


def run(R_V, k_prot, L_tis, D_tis):
    td = TwoDomain(D_gel=Dg, D_tis=D_tis, L_gel=P.L_GEL, L_tis=L_tis,
                   nx_gel=60, nx_tis=200, eps_gel=P.EPS_GEL,
                   eps_tis=P.EPS_TIS, K_p=P.K_P, k_prot=k_prot, R_V=R_V)
    res = td.run(free_from_total(P.C_SOLUBILITY, R_V), t_end=32 * DAY)
    return td.days_above(res, P.DELTA, P.C_TARGET_V14)


RVs = np.logspace(0, 2, 13)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.6, 6.2), sharex=True,
                               gridspec_kw=dict(hspace=0.12, height_ratios=[1.35, 1]))

summary = {}
for k_prot, col in zip(P.K_PROT_RANGE, (S.C_V14_L, S.C_V14, "#0B3C5D")):
    cen, lo, hi, pk = [], [], [], []
    for R in RVs:
        d0, p0 = run(R, k_prot, P.L_TIS, P.D_TIS_V14)
        dl, _ = run(R, k_prot, 0.2, P.D_TIS_V14_HI)     # fastest drain
        dh, _ = run(R, k_prot, 0.4, P.D_TIS_V14_LO)     # slowest drain
        cen.append(d0); lo.append(dl); hi.append(dh); pk.append(p0 / UM)
    cen, lo, hi, pk = map(np.array, (cen, lo, hi, pk))
    lab = f"$k_{{prot}}$ = {k_prot:.0e} s$^{{-1}}$"
    ax1.fill_between(RVs, lo, hi, color=col, alpha=0.16, lw=0)
    ax1.plot(RVs, cen, color=col, label=lab)
    ax2.plot(RVs, pk, color=col, label=lab)
    summary[k_prot] = (RVs, cen, lo, hi, pk)

# mark the optimum on the slowest-proteolysis curve
RVo, ceno = summary[1e-6][0], summary[1e-6][1]
io = int(np.argmax(ceno))
ax1.plot([RVo[io]], [ceno[io]], "*", ms=15, color=S.C_THRESH, zorder=6)
ax1.annotate(f"optimum\n$R_V$ ≈ {RVo[io]:.0f},  {ceno[io]:.1f} d",
             (RVo[io], ceno[io]), textcoords="offset points", xytext=(14, -6),
             fontsize=9, weight="bold", color=S.C_THRESH)

first_zero = RVs[np.argmax(summary[1e-6][1] <= 0)] if np.any(summary[1e-6][1] <= 0) else None
if first_zero:
    for ax in (ax1, ax2):
        ax.axvspan(first_zero, 100, color=S.C_BAD, alpha=0.09, lw=0)
    ax1.text(np.sqrt(first_zero * 100), ax1.get_ylim()[1] * 0.62,
             "coverage\ncollapses\nto ZERO", ha="center", fontsize=9.5,
             weight="bold", color=S.C_BAD)

ax1.annotate("the gel as built\n(no heparin, $R_V$ = 1)", (1.0, summary[1e-5][1][0]),
             textcoords="offset points", xytext=(10, 26), fontsize=8.5,
             color=S.C_DERIVED,
             arrowprops=dict(arrowstyle="->", color=S.C_DERIVED, lw=1))

ax2.axhline(P.C_TARGET_V14 / UM, color=S.C_THRESH, ls=":", lw=1.6)
ax2.text(1.05, 28, "25 µM therapeutic target", fontsize=8.4, color=S.C_THRESH)
ax2.set_yscale("log")
ax2.set_xscale("log")
ax2.set_xlabel("retardation factor  $R_V$   (retention strength)")
ax2.set_ylabel("peak tissue conc.  (µM)")
ax1.set_ylabel("days above the 25 µM target")
ax1.set_title("More retention is not better: $R_V$ has a sharp optimum")
ax1.legend(loc="upper left")
ax1.set_ylim(bottom=0)
for ax in (ax1, ax2):
    ax.set_xticks([1, 3, 10, 30, 100], ["1", "3", "10", "30", "100"])
S.tier4_note(ax1, "band = L$_{tis}$ 2–4 mm, D$_t$ 0.7–1.8e-6 cm²/s", (0.98, 0.52))
S.save(fig, "fig07_retention")

print(f"{'R_V':>6} " + " ".join(f"{k:>22}" for k in P.K_PROT_RANGE))
for i, R in enumerate(RVs):
    if R in RVs[[0, 3, 6, 9, 12]]:
        row = " ".join(f"{summary[k][1][i]:6.1f} d [{summary[k][2][i]:4.1f}-"
                       f"{summary[k][3][i]:4.1f}] {summary[k][4][i]:5.0f}uM"
                       for k in P.K_PROT_RANGE)
        print(f"{R:6.1f} {row}")
