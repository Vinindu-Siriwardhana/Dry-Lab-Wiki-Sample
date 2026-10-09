"""
F7 -- What the as-built gel supplies, and why adding retention would not rescue it.

The gel is heparin-free, so R_V = 1 is LOCKED (constants register) and the
retention lever is out of the design space.  This figure shows the
counterfactual anyway, because it is the reason the lever was dropped:
holding the TOTAL loading at the solubility ceiling, raising R_V buys residence
time and spends driving concentration in exact proportion
(C_free,0 = C_total,0 / R_V).  With a 150 um gel the reservoir is so small that
the peak tissue concentration falls below the 25 uM target after only a few-fold
retention, and coverage collapses to ZERO.

The duration is plotted as a BAND, because it is sensitive to two parameters
known only as ranges: the depth to the systemic sink (L_tis = 1-2 mm) and the
tissue diffusivity (D_t = 0.7-1.8e-6 cm2/s).

Writes _coverage.npy (hours above target at R_V = 1, central/low/high, for the
register k_prot) -- fig10 reads it.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import style as S, params as P, v3model as m
from twodomain import TwoDomain, free_from_total, UM, DAY

S.apply()
Dg = m.D_eff(P.D0_V14, P.A_V14_NM, P.XI_NM, P.EPS_GEL)
SEG = [(0.5, 120), (5.0, 108), (30.0, 2800), (300.0, 600)]   # ~24 h at 30 s, then 2 d


def run(R_V, k_prot, L_tis, D_tis):
    td = TwoDomain(D_gel=Dg, D_tis=D_tis, L_gel=P.L_GEL, L_tis=L_tis,
                   nx_gel=60, nx_tis=200, eps_gel=P.EPS_GEL,
                   eps_tis=P.EPS_TIS, K_p=P.K_P, k_prot=k_prot,
                   k_deg_gel=P.K_DEG_GEL, R_V=R_V)
    res = td.run(free_from_total(P.C_SOLUBILITY, R_V), t_end=3 * DAY, segments=SEG)
    d, pk = td.days_above(res, P.DELTA, P.C_TARGET_V14)
    return d * 24.0, pk


RVs = np.logspace(0, np.log10(30), 13)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.6, 6.2), sharex=True,
                               gridspec_kw=dict(hspace=0.12, height_ratios=[1.35, 1]))

KPS = (1e-4, P.K_PROT_V14, 1e-6)
summary = {}
for k_prot, col in zip(KPS, (S.C_V14_L, S.C_V14, "#0B3C5D")):
    cen, lo, hi, pk = [], [], [], []
    for R in RVs:
        h0, p0 = run(R, k_prot, P.L_TIS, P.D_TIS_V14)
        hl, _ = run(R, k_prot, P.L_TIS_RANGE[0], P.D_TIS_V14_HI)   # fastest drain
        hh, _ = run(R, k_prot, P.L_TIS_RANGE[1], P.D_TIS_V14_LO)   # slowest drain
        cen.append(h0); lo.append(hl); hi.append(hh); pk.append(p0 / UM)
    cen, lo, hi, pk = map(np.array, (cen, lo, hi, pk))
    lab = f"$k_{{prot}}$ = {k_prot:.2g} s$^{{-1}}$" + \
        ("  (register value)" if k_prot == P.K_PROT_V14 else "")
    ax1.fill_between(RVs, lo, hi, color=col, alpha=0.16, lw=0)
    ax1.plot(RVs, cen, color=col, label=lab)
    ax2.plot(RVs, pk, color=col, label=lab)
    summary[k_prot] = (RVs, cen, lo, hi, pk)

# the device as built
c0 = summary[P.K_PROT_V14]
ax1.plot([1.0], [c0[1][0]], "*", ms=15, color=S.C_THRESH, zorder=6)
ax1.annotate(f"the gel as built\n(no heparin, $R_V$ = 1):  {c0[1][0]:.1f} h\n"
             f"(band {c0[2][0]:.1f}–{c0[3][0]:.1f} h)",
             (1.0, c0[1][0]), textcoords="offset points", xytext=(16, 62),
             fontsize=9, weight="bold", color=S.C_THRESH,
             bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.85),
             arrowprops=dict(arrowstyle="->", color=S.C_THRESH, lw=1))

z = summary[1e-6][1]
first_zero = RVs[np.argmax(z <= 0)] if np.any(z <= 0) else None
if first_zero:
    for ax in (ax1, ax2):
        ax.axvspan(first_zero, RVs[-1], color=S.C_BAD, alpha=0.09, lw=0)
    ax1.text(np.sqrt(first_zero * RVs[-1]), max(c0[3].max(), 1) * 0.55,
             "coverage\ncollapses\nto ZERO", ha="center", fontsize=9.5,
             weight="bold", color=S.C_BAD)

ax2.axhline(P.C_TARGET_V14 / UM, color=S.C_THRESH, ls=":", lw=1.6)
ax2.text(1.05, 28, "25 µM therapeutic target", fontsize=8.4, color=S.C_THRESH)
ax2.set_yscale("log")
ax2.set_xscale("log")
ax2.set_xlabel("retardation factor  $R_V$   (counterfactual: locked at 1 in the device)")
ax2.set_ylabel("peak tissue conc.  (µM)")
ax1.set_ylabel("hours above the 25 µM target")
ax1.set_title("One patch covers hours, not days, and retention cannot fix it")
ax1.legend(loc="upper right", fontsize=8.2)
ax1.set_ylim(bottom=0)
for ax in (ax1, ax2):
    ax.set_xticks([1, 2, 3, 5, 10, 30], ["1", "2", "3", "5", "10", "30"])
    ax.set_xlim(1, RVs[-1])
S.tier4_note(ax1, "band = L$_{tis}$ 1–2 mm, D$_t$ 0.7–1.8e-6 cm²/s", (0.70, 0.02))
S.save(fig, "fig07_retention")

np.save(os.path.join(os.path.dirname(__file__), "_coverage.npy"),
        np.array([c0[1][0], c0[2][0], c0[3][0]]) / 24.0)      # days
print(f"{'R_V':>6} " + " ".join(f"{k:>24.2g}" for k in KPS))
for i, R in enumerate(RVs):
    row = " ".join(f"{summary[k][1][i]:5.2f} h [{summary[k][2][i]:4.2f}-"
                   f"{summary[k][3][i]:5.2f}] {summary[k][4][i]:5.0f}uM"
                   for k in KPS)
    print(f"{R:6.2f} {row}")
