"""
F12 -- V3 corrects a SCALING LAW in V2, not just a prefactor.

Inverting for the required dose needs psi isolated ADDITIVELY; V2 factored out
K_D,app, which is legitimate in the occupancy expression and in Cheng-Prusoff
but not here.  The consequence is that V2 was quadratic in endotoxin load where
the mechanism is linear -- so it mispredicted exactly the quantity it claimed
to deliver.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import style as S, params as P, v3model as m

S.apply()
a = np.logspace(-0.5, 3.3, 400)
KD_P, KD_L, T_ratio = 1.0, 1.0, 5.0
v3 = m.dose_required(KD_P, a * KD_L, KD_L, T_ratio, 1.0)
v2 = m.dose_required_v2(KD_P, a * KD_L, KD_L, T_ratio, 1.0)

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.4, 4.0),
                              gridspec_kw=dict(width_ratios=[1.3, 1]))
ax.plot(a, v3, color=S.C_V14, label="V3, corrected  —  linear in $[L]$")
ax.plot(a, v2, color=S.C_GREY, ls="--", label="V2  —  quadratic in $[L]$")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"scaled endotoxin load  $a = [L]/K_{D,L}$")
ax.set_ylabel(r"required dose  $[P]_{req}\,/\,K_{D,P}$")
ax.set_title("The error grows with wound severity", loc="left", fontsize=10)
ax.legend(loc="upper left")

for aa in (1, 10, 100, 1000):
    r = m.dose_required_v2(KD_P, aa, KD_L, T_ratio, 1.0) / \
        m.dose_required(KD_P, aa, KD_L, T_ratio, 1.0)
    ax.annotate(f"{r:.0f}×", (aa, m.dose_required_v2(KD_P, aa, KD_L, T_ratio, 1.0)),
                textcoords="offset points", xytext=(4, 6), fontsize=8.4,
                color=S.C_DERIVED, weight="bold")

ratio = v2 / v3
ax2.plot(a, ratio, color=S.C_FGF)
ax2.set_xscale("log"); ax2.set_yscale("log")
ax2.set_xlabel(r"$a = [L]/K_{D,L}$")
ax2.set_ylabel("V2 over-estimate  (×)")
ax2.set_title("Over-estimate factor", loc="left", fontsize=10)
ax2.axhline(1, color=S.C_GREY, ls=":", lw=1.2)
ax2.text(0.5, 0.92, f"{ratio.min():.1f}× to {ratio.max():.0f}×",
         transform=ax2.transAxes, ha="center", fontsize=11, weight="bold",
         color=S.C_FGF)
S.save(fig, "fig12_dose_correction")
print(f"  over-estimate spans {ratio.min():.1f}x to {ratio.max():.0f}x")
