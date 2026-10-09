"""F6 -- Steric hindrance: the transport coefficients were derived, not assumed."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import style as S, params as P, v3model as m

S.apply()
xi = np.linspace(3, 50, 400)
fig, ax = plt.subplots(figsize=(5.0, 3.4))
for name, a, D0, col in ((f"V14  ({P.A_V14_NM:g} nm)", P.A_V14_NM, P.D0_V14, S.C_V14),
                         (f"FGF2-G3  ({P.A_FGF_NM:g} nm)", P.A_FGF_NM, P.D0_FGF, S.C_FGF)):
    ax.plot(xi, m.D_eff(D0, a, xi, P.EPS_GEL) / D0, color=col, label=name)
ax.axvline(P.XI_NM, color=S.C_GREY, ls=":", lw=1.4)
ax.text(P.XI_NM + 0.7, 0.30, f"$\\xi$ = {P.XI_NM:.0f} nm\nused throughout\n(MEASURE THIS)",
        fontsize=8, color="#666666")
ax.set_xlabel("hydrogel mesh size  $\\xi$  (nm)")
ax.set_ylabel("$D_{eff}/D_0$")
ax.set_ylim(0, 1.02)
ax.set_title("The gel barely slows either molecule", fontsize=10)
ax.legend(loc="lower right")
S.save(fig, "fig06_hindrance")
print(f"  at xi=10 nm: V14 {m.D_eff(P.D0_V14,P.A_V14_NM,10,P.EPS_GEL)/P.D0_V14:.2f}, "
      f"FGF {m.D_eff(P.D0_FGF,P.A_FGF_NM,10,P.EPS_GEL)/P.D0_FGF:.2f}")
