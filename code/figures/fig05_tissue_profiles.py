"""
F5 -- V14 concentration profiles in tissue, from the two-domain solver.

Shows where the drug actually is over time, the penetration depth Lambda, and
how little of the window sits above the 25 uM therapeutic target.  Run for the
gel AS BUILT: 150 um, heparin-free (R_V = 1), loaded at the solubility ceiling.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import style as S, params as P, v3model as m
from twodomain import TwoDomain, free_from_total, UM, DAY

S.apply()
Dg = m.D_eff(P.D0_V14, P.A_V14_NM, P.XI_NM, P.EPS_GEL)
fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.9), sharey=True)
SEG = [(0.5, 120), (5.0, 108), (30.0, 1400)]          # resolve to ~12 h

peak = 0.0
for ax, kprot in zip(axes, (P.K_PROT_V14, 1e-6)):
    td = TwoDomain(D_gel=Dg, D_tis=P.D_TIS_V14, L_gel=P.L_GEL, L_tis=P.L_TIS,
                   nx_gel=60, nx_tis=200, eps_gel=P.EPS_GEL, eps_tis=P.EPS_TIS,
                   K_p=P.K_P, k_prot=kprot, k_deg_gel=P.K_DEG_GEL, R_V=P.R_V)
    res = td.run(free_from_total(P.C_SOLUBILITY, P.R_V), t_end=12 * 3600,
                 segments=SEG)
    depth = np.concatenate([[0.0], -td.xt[1:]]) * 10          # mm
    peak = max(peak, res["C_tis"].max() / UM)

    want = [1, 5, 15, 60, 120, 240, 480]                        # minutes
    cmap = plt.cm.viridis(np.linspace(0.08, 0.88, len(want)))
    for mn, col in zip(want, cmap):
        i = int(np.argmin(np.abs(res["t"] - mn * 60)))
        ax.plot(depth, res["C_tis"][i] / UM, color=col, lw=1.7,
                label=f"{mn:g} min" if mn < 60 else f"{mn/60:g} h")

    lam = m.penetration_depth(P.D_TIS_V14, kprot) * 10
    if lam < P.L_TIS * 10:
        ax.axvline(lam, color=S.C_DERIVED, ls="-.", lw=1.1)
        ax.text(lam + 0.03, 330, f"$\\Lambda$ = {lam:.2f} mm", fontsize=8.4,
                color=S.C_DERIVED, rotation=90, va="top")
    else:
        ax.text(0.97, 0.6, f"$\\Lambda$ = {lam:.0f} mm\n(beyond the sink)",
                transform=ax.transAxes, ha="right", fontsize=8.4,
                color=S.C_DERIVED)
    ax.axvline(P.DELTA * 10, color=S.C_GREY, ls="--", lw=1.0)
    ax.axhline(25, color=S.C_THRESH, ls=":", lw=1.6)
    ax.set_xlabel("depth into tissue  (mm)")
    ax.set_title(f"$k_{{prot}}$ = {kprot:.2g} s$^{{-1}}$", loc="left")
    ax.set_xlim(0, P.L_TIS * 10)

axes[0].set_ylabel("free V14 in tissue  (µM)")
axes[0].text(P.L_TIS * 10 * 0.62, 30, "25 µM target", fontsize=8.4, color=S.C_THRESH)
axes[0].text(P.DELTA * 10 + 0.02, 230, "target depth\n$\\delta$ = %.0f µm" % (P.DELTA * 1e4),
             fontsize=7.8, color="#777777")
axes[0].set_ylim(0, 360)
axes[1].legend(ncol=2, loc="upper right", title="time after application",
               title_fontsize=8)
fig.suptitle("The as-built gel empties in minutes; the dermal sink drains it in hours",
             fontsize=11, weight="bold", y=1.0)
S.save(fig, "fig05_tissue_profiles")
print(f"  peak free V14 in tissue {peak:.0f} uM; "
      f"Lambda at k_prot={P.K_PROT_V14:.2e}: "
      f"{m.penetration_depth(P.D_TIS_V14, P.K_PROT_V14)*1e4:.0f} um")
