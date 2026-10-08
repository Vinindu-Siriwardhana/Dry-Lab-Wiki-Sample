"""
F8 -- Solubility, not potency, is the binding constraint.

Eq. (dosing) bounds the FREE loading.  With sorptive retention the TOTAL that
must actually dissolve in the casting buffer is R_V times larger, because
C_free,0 = C_total,0 / R_V.  Plotting that against the saturation ceiling
recovers the R_V ~ 3-10 optimum of F7 from a completely independent direction --
two routes, same answer.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import style as S, params as P, v3model as m

S.apply()
UM = P.UM
C_free_req = m.loading_required(P.C_TARGET_V14, P.DELTA, P.L_GEL,
                                P.F_REL, P.EPS_GEL, eps_tis=1.0)
R = np.logspace(0, 2, 300)
C_tot = R * C_free_req

fig, ax = plt.subplots(figsize=(6.3, 4.0))
ax.plot(R, C_tot / UM, color=S.C_V14, lw=2.2,
        label="total V14 that must dissolve  $= R_V \\times C_{free,0}$")
ax.axhline(P.C_SOLUBILITY / UM, color=S.C_FGF, lw=2.0)
ax.text(1.06, P.C_SOLUBILITY / UM * 1.1,
        f"saturation ceiling ≈ {P.C_SOLUBILITY/UM:.0f} µM  (1 mg mL$^{{-1}}$)",
        color=S.C_FGF, fontsize=8.8, weight="bold")

Rmax = P.C_SOLUBILITY / C_free_req
ax.fill_between(R, 1, np.minimum(C_tot / UM, P.C_SOLUBILITY / UM),
                where=(R <= Rmax), color=S.C_GOOD, alpha=0.13, lw=0)
ax.axvline(Rmax, color=S.C_DERIVED, ls="--", lw=1.2)
ax.text(Rmax * 1.05, 2.3e1, f"formulable only\nbelow $R_V$ = {Rmax:.0f}",
        fontsize=8.6, color=S.C_DERIVED, weight="bold")

for Rv in (1, 3, 10, 30):
    c = Rv * C_free_req / UM
    ok = c <= P.C_SOLUBILITY / UM
    ax.plot([Rv], [c], "o", ms=7, color=S.C_V14 if ok else S.C_BAD, zorder=5)
    ax.annotate(f"$R_V$={Rv}\n{c:.0f} µM" + ("" if ok else "\nIMPOSSIBLE"),
                (Rv, c), textcoords="offset points", xytext=(-6, 12),
                fontsize=8.2, ha="right",
                color=S.C_V14 if ok else S.C_BAD,
                weight="normal" if ok else "bold")

ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(1, 100); ax.set_ylim(20, 4000)
ax.set_xticks([1, 3, 10, 30, 100], ["1", "3", "10", "30", "100"])
ax.set_xlabel("retardation factor  $R_V$")
ax.set_ylabel("required total loading  (µM)")
ax.set_title("Retention is paid for in solubility")
ax.legend(loc="lower right")
S.save(fig, "fig08_solubility")
print(f"  free loading required: {C_free_req/UM:.1f} uM "
      f"({C_free_req/UM*P.MW_V14/1000:.0f} ug/mL)")
for Rv in (1, 3, 10, 30):
    print(f"   R_V={Rv:3d} -> total {Rv*C_free_req/UM:6.0f} uM "
          f"({Rv*C_free_req/P.C_SOLUBILITY:.0%} of ceiling)")
