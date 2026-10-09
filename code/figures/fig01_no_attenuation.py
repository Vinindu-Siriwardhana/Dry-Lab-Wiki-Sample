"""
F1 -- Why NO-responsive sensing was geometrically infeasible.

Steady penetration depth for diffusion with first-order consumption,
L_d = sqrt(D_eff / k_decay), k_decay = ln2 / t_half.  Attenuation across a
stack of thickness L is exp(-L/L_d).  Plotted against the actual 0.5-3 mm
bacterial-cellulose + hydrogel stack of the V2 device.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import style as S, params as P, v3model as m

S.apply()

def Ld(inv_tau, t_half):
    D = P.D0_NO * m.amsden_H(0.15, P.XI_NM) * inv_tau   # H ~ 1 for NO
    return np.sqrt(D / (np.log(2) / t_half))

fig, ax = plt.subplots(figsize=(6.4, 4.0))
L = np.linspace(0, 0.32, 600)          # cm

lds = [Ld(it, th) for it in P.INV_TAU for th in P.NO_HALFLIFE]
ld_lo, ld_hi = min(lds), max(lds)

ax.fill_between(L * 10, np.exp(-L / ld_lo), np.exp(-L / ld_hi),
                color=S.C_V14_L, alpha=0.45, lw=0,
                label=f"$L_d$ = {ld_lo*1e4:.0f}–{ld_hi*1e4:.0f} µm (swept)")
ax.plot(L * 10, np.exp(-L / Ld(0.5, P.NO_HALFLIFE[1])), color=S.C_V14,
        label="central case, $L_d$ = %.0f µm" % (Ld(0.5, P.NO_HALFLIFE[1]) * 1e4))

ax.axvspan(*P.STACK_MM, color=S.C_GREY, alpha=0.35, lw=0)
ax.text(1.75, 1e-6, "actual device stack\n0.5 – 3 mm", ha="center", va="center",
        fontsize=9, color="#555555", weight="bold")

ax.axhline(1e-2, color=S.C_THRESH, ls=":", lw=1.4)
ax.text(0.04, 1.3e-2, "most optimistic NorR activation threshold",
        fontsize=7.8, color=S.C_THRESH)

for Lx, lab in ((0.1, None), (0.2, None)):
    pass
Ls = P.STACK_MM[0] / 10                    # thinnest stack, cm
a_hi, a_lo = np.exp(-Ls / ld_hi), max(np.exp(-Ls / ld_lo), 1e-12)
ax.annotate("", xy=(Ls * 10, a_lo), xytext=(Ls * 10, a_hi),
            arrowprops=dict(arrowstyle="<->", color=S.C_DERIVED, lw=1.3))
ax.text(Ls * 10 + 0.08, a_hi * 1e-2,
        "at %.1f mm:\n$\\leq 10^{%d}$" % (Ls * 10, int(np.ceil(np.log10(a_hi)))),
        fontsize=8.5, va="center", color=S.C_DERIVED, weight="bold")

ax.set_yscale("log")
ax.set_ylim(1e-12, 2)
ax.set_xlim(0, 3.2)
ax.set_xlabel("depth into the patch stack, $L$  (mm)")
ax.set_ylabel("fraction of wound NO reaching the chassis,  $e^{-L/L_d}$")
ax.set_title("NO could never have reached the sensor")
ax.legend(loc="upper right")
S.save(fig, "fig01_no_attenuation")
print(f"  L_d range {ld_lo*1e4:.0f}-{ld_hi*1e4:.0f} um; "
      f"attenuation at {Ls*1e4:.0f} um: {np.exp(-Ls/ld_hi):.1e} to {np.exp(-Ls/ld_lo):.1e}")
