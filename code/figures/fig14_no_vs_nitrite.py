"""
F14 -- Free NO and accumulated nitrite are not the same observable.

Free NO is quasi-steady within minutes (1/k_scav = 100 s) and sits at nanomolar
levels; nitrite accumulates all day and reaches micromolar.  The Griess assay
measures the second one.  Fitting the model's instantaneous free-NO expression
to a Griess curve compares two quantities that differ by orders of magnitude.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import style as S, params as P, v3model as m

S.apply()
t = np.linspace(0, 24 * 3600, 4000)
N_cell = 1e6            # cells cm^-3  (2e5 cells in 200 uL)
q_max = P.Q_CHRONIC     # mol cell^-1 s^-1, sustained iNOS rate (register)
p65 = 0.85              # normalised, LPS-stimulated

R, freeNO, NO2 = m.inos_response(t, p65, N_cell, q_max, K_iNOS_rel=P.K_INOS_REL,
                                 p=P.P_INOS, k_turn=P.K_TURN, k_scav=P.K_SCAV,
                                 phi=P.PHI_NO2)

fig, ax = plt.subplots(figsize=(6.4, 3.8))
ax2 = ax.twinx()
l1, = ax.plot(t / 3600, NO2 / P.UM, color=S.C_V14, lw=2.2)
l2, = ax2.plot(t / 3600, freeNO / P.UM * 1000, color=S.C_FGF, lw=2.0, ls="--")

ax.set_xlabel("time after stimulation  (h)")
ax.set_ylabel("accumulated nitrite  (µM)", color=S.C_V14)
ax2.set_ylabel("free NO  (nM)", color=S.C_FGF)
ax.tick_params(axis="y", colors=S.C_V14)
ax2.tick_params(axis="y", colors=S.C_FGF)
ax2.grid(False)
ax.set_title("What Griess measures (blue) vs what the balance gives (red)",
             fontsize=10)
ax.legend([l1, l2], ["$[NO_2^-]$ — the measured observable",
                     "free $[NO]$ — quasi-steady, never measured"],
          loc="upper left", fontsize=8.5)
ax.text(0.97, 0.06, f"ratio at 24 h ≈ {NO2[-1]/freeNO[-1]:,.0f}×",
        transform=ax.transAxes, ha="right", fontsize=9.5, weight="bold",
        color=S.C_DERIVED)
S.save(fig, "fig14_no_vs_nitrite")
print(f"  nitrite at 24 h {NO2[-1]/P.UM:.1f} uM; free NO {freeNO[-1]/P.UM*1000:.2f} nM; "
      f"ratio {NO2[-1]/freeNO[-1]:,.0f}x")
