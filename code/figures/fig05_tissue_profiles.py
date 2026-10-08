"""
F5 -- V14 concentration profiles in tissue, from the two-domain solver.

Shows where the drug actually is over time, the penetration depth Lambda, and
how little of the window sits above the 25 uM therapeutic target.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import style as S, params as P, v3model as m
from twodomain import TwoDomain, free_from_total, UM, DAY

S.apply()
Dg = m.D_eff(P.D0_V14, P.A_V14_NM, P.XI_NM, P.EPS_GEL)
fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.9), sharey=True)

for ax, kprot in zip(axes, (1e-4, 1e-6)):
    td = TwoDomain(D_gel=Dg, D_tis=P.D_TIS_V14, L_gel=P.L_GEL, L_tis=P.L_TIS,
                   nx_gel=60, nx_tis=200, eps_gel=P.EPS_GEL, eps_tis=P.EPS_TIS,
                   K_p=P.K_P, k_prot=kprot, R_V=10.0)
    res = td.run(free_from_total(P.C_SOLUBILITY, 10.0), t_end=12 * DAY)
    depth = np.concatenate([[0.0], -td.xt[1:]]) * 10          # mm

    want = [0.25, 1, 4, 24, 72, 168, 288]                      # hours
    cmap = plt.cm.viridis(np.linspace(0.08, 0.88, len(want)))
    for h, col in zip(want, cmap):
        i = int(np.argmin(np.abs(res["t"] - h * 3600)))
        ax.plot(depth, res["C_tis"][i] / UM, color=col, lw=1.7,
                label=f"{h:g} h" if h < 24 else f"{h/24:g} d")

    lam = m.penetration_depth(P.D_TIS_V14, kprot) * 10
    ax.axvline(lam, color=S.C_DERIVED, ls="-.", lw=1.1)
    ax.text(lam + 0.05, 320, f"$\\Lambda$ = {lam:.1f} mm", fontsize=8.4,
            color=S.C_DERIVED, rotation=90, va="top")
    ax.axhline(25, color=S.C_THRESH, ls=":", lw=1.6)
    ax.set_xlabel("depth into tissue  (mm)")
    ax.set_title(f"$k_{{prot}}$ = {kprot:.0e} s$^{{-1}}$", loc="left")
    ax.set_xlim(0, P.L_TIS * 10)

axes[0].set_ylabel("free V14 in tissue  (µM)")
axes[0].text(0.55, 27, "25 µM target", fontsize=8.4, color=S.C_THRESH)
axes[0].set_ylim(0, 340)
axes[1].legend(ncol=2, loc="upper right", title="time after application",
               title_fontsize=8)
fig.suptitle("Protease resistance, not gel chemistry, decides how long the dose lasts",
             fontsize=11, weight="bold", y=1.0)
S.save(fig, "fig05_tissue_profiles")
