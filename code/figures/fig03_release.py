"""
F3 -- Release is a bolus, not a sustained delivery, and the two molecules do
not arrive together.

Cumulative fractional release from the closed-form Crank series, plotted on a
log time axis against the 7-21 day wound-closure window.  The empty gap between
the curves and the shaded window IS the 60-5000x timescale separation.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import style as S, params as P, v3model as m

S.apply()
Dv = m.D_eff(P.D0_V14, P.A_V14_NM, P.XI_NM, P.EPS_GEL, P.A_FIBRE_NM)
Df = m.D_eff(P.D0_FGF, P.A_FGF_NM, P.XI_NM, P.EPS_GEL, P.A_FIBRE_NM)

fig, ax = plt.subplots(figsize=(7.0, 4.2))
t = np.logspace(0, np.log10(25 * 86400), 900)

cases = [("V14",      Dv, P.L_GEL_THIN, S.C_V14, "-",  2.1),
         ("FGF2-G3",  Df, P.L_GEL_THIN, S.C_FGF, "-",  2.1),
         ("V14",      Dv, P.L_GEL,      S.C_V14, "--", 1.4),
         ("FGF2-G3",  Df, P.L_GEL,      S.C_FGF, "--", 1.4)]

t95 = {}
for name, D, L, col, ls, lw in cases:
    f = m.crank_cumulative(t, D, L)
    lab = f"{name}, $L_{{gel}}$ = {L*1e4:.0f} µm"
    ax.plot(t / 60, f, color=col, ls=ls, lw=lw, label=lab)
    tt = m.release_time(0.95, D, L); t95[(name, L)] = tt
    ax.plot([tt / 60], [0.95], "o", color=col, ms=5, zorder=5)

ax.axvspan(7 * 1440, 21 * 1440, color=S.C_GOOD, alpha=0.13, lw=0)
ax.text(12.5 * 1440, 0.46, "wound closure\nhappens here\n(7 – 21 days)",
        ha="center", va="center", fontsize=9.5, color="#36691a", weight="bold")

a, b = t95[("V14", P.L_GEL_THIN)] / 60, t95[("FGF2-G3", P.L_GEL_THIN)] / 60
ax.annotate("", xy=(a, 0.80), xytext=(b, 0.80),
            arrowprops=dict(arrowstyle="<->", color=S.C_DERIVED, lw=1.5))
ax.text(np.sqrt(a * b), 0.845, f"{b/a:.1f}× arrival separation\n(pure size effect)",
        ha="center", fontsize=8.6, color=S.C_DERIVED, weight="bold")

ax.axhline(0.95, color=S.C_GREY, lw=0.9, ls=":")
ax.text(1.15, 0.965, "95 % released", fontsize=7.8, color="#777777")

ax.set_xscale("log")
ax.set_xlim(1, 25 * 1440)
ax.set_ylim(0, 1.04)
ticks = [1, 10, 60, 240, 1440, 10080, 10080 * 3]
ax.set_xticks(ticks, ["1 min", "10 min", "1 h", "4 h", "1 d", "1 wk", "3 wk"])
ax.set_xlabel("time after application")
ax.set_ylabel("cumulative fractional release  $f_i(t)$")
ax.set_title("Both actives are gone before healing begins")
ax.legend(loc="lower right", bbox_to_anchor=(0.995, 0.03))
S.save(fig, "fig03_release")
for k, v in t95.items():
    print(f"  95% release {k[0]:8s} L={k[1]*1e4:4.0f}um : {v/60:6.1f} min")
print(f"  separation factor {b/a:.2f}x")
