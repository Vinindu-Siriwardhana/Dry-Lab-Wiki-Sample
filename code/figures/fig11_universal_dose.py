"""
F11 -- The V14 arm runs before K_D,P is ever measured.

In the scaled variables psi = C_V/K_D,P and a = [L]/K_D,L the whole of Models
3-4 collapses onto four numbers, and the suppression curve carries NO fitted
constants at all.  Half-maximal suppression sits exactly at psi_50 = 1 + a.

Measuring K_D,P converts the abscissa to a real concentration.  It calibrates
the model; it does not unblock it.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import style as S, params as P, v3model as m

S.apply()
psi = np.logspace(-2, 4.3, 700)
fig, ax = plt.subplots(figsize=(6.6, 4.2))

cols = plt.cm.viridis(np.linspace(0.05, 0.85, 5))
for a, c in zip((0.1, 1, 10, 100, 1000), cols):
    y = m.suppression_ratio(psi, a)
    ax.plot(psi, y, color=c, label=f"$a = [L]/K_{{D,L}}$ = {a:g}")
    p50 = m.psi_50(a)
    ax.plot([p50], [0.5], "o", ms=6, color=c, zorder=5)

ax.axhline(0.5, color=S.C_GREY, ls=":", lw=1.3)
ax.annotate(r"$\psi_{50} = 1 + a$   exactly" "\n" r"no fitted constant anywhere",
            (3e3, 0.5), textcoords="offset points", xytext=(-12, 34),
            ha="right", fontsize=10, weight="bold", color=S.C_DERIVED)
ax.set_xscale("log")
ax.set_xlabel(r"scaled dose  $\psi = C^{t}_{V}\,/\,K_{D,P}$")
ax.set_ylabel(r"active receptor,  $\theta_{LPS}(\psi)/\theta_{LPS}(0)$")
ax.set_title("A dose–response with no fitted parameters")
ax.set_ylim(0, 1.03)
ax.legend(loc="lower left", fontsize=8.3)
S.save(fig, "fig11_universal_dose")
