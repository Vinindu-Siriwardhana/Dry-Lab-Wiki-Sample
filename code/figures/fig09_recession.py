"""
F9 -- An insufficient pulse is worse than it looks: the wound REOPENS.

Rates come from the register closures, not from hand-picked numbers:

  chronic baseline   C_e = 0, I_e = theta_LPS(0) = 0.94 (1 uM LPS, K_D,L 65 nM)
                     r_p = 0.50/d, d = 0.51/d  -> r_eff ~ -0.015/d, never closes
  therapeutic pulse  FGF2-G3 at EC50, V14 at psi_50 (I_e halved to 0.47)
                     r_p = 1.11/d, d = 0.45/d  -> r_eff = +0.66/d

The pulse lasts T days, after which the drug effect relaxes with the
effect-compartment time constant 1/k_e0 = 18 h.

When r_eff returns negative the front does not merely stall.  The stable state
behind it, n_ss = K(1 - d/r_p), becomes zero, so the repopulated tissue dies
back and the front RECEDES past where it started.  A sub-threshold treatment
therefore produces transient apparent healing followed by complete relapse --
a recognisable clinical picture in chronic wounds, and a falsifiable prediction.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import style as S, params as P, v3model as m

S.apply()
DAY = P.DAY


def pulse(T_days):
    """Therapeutic rates for T days, then exponential relaxation to chronic."""
    T = T_days * DAY
    def s(t):
        return 1.0 if t <= T else np.exp(-(t - T) / P.TAU_RELAX)
    rp = lambda t: P.R_CHRONIC + (P.R_PULSE - P.R_CHRONIC) * s(t)
    dd = lambda t: P.D_CHRONIC + (P.D_PULSE - P.D_CHRONIC) * s(t)
    return rp, dd


THR = 0.5 * P.K_CARRY * (1 - P.D_PULSE / P.R_PULSE)     # pinned front level


def closes(R_w_cm, T_days, t_end_days=110, nr=700):
    rp, dd = pulse(T_days)
    sol = m.FisherKPP(P.D_N, K=P.K_CARRY, R_inf=R_w_cm * 2.2, nr=nr)
    res = sol.run(R_w_cm, rp, dd, t_end=t_end_days * DAY, record_every=40,
                  fixed_threshold=THR, w=P.W_SMOOTH)
    return res


def min_pulse(R_w_cm, lo=1.0, hi=60.0, tol=0.4):
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if closes(R_w_cm, mid)["closed"]:
            hi = mid
        else:
            lo = mid
    return hi


# ---- panel: trajectories at R_w = 2 mm -------------------------------------
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.1),
                              gridspec_kw=dict(width_ratios=[1.45, 1]))

R_w = 0.2
for T, col in zip((5, 10, 14, 20), (S.C_BAD, "#C97BB0", S.C_V14_L, S.C_V14)):
    res = closes(R_w, T, t_end_days=60)
    rf = res["r_front"].copy()
    # The model describes a CHRONIC WOUND BED.  Once the front reaches the
    # centre the wound is closed and the bed no longer exists, so the equations
    # are not run past that point -- t_c is defined as the first time r_f = 0.
    # We therefore hold the trace at zero after closure rather than continuing
    # to apply chronic rates to healed tissue, which the model does not describe.
    z = np.argmax(rf <= 1e-9) if np.any(rf <= 1e-9) else None
    if z is not None:
        rf[z:] = 0.0
    ax.plot(res["t"] / DAY, np.clip(rf, 0, 0.26) * 10, color=col,
            label=f"{T} d pulse" + ("  ✓ closes" if res["closed"] else "  ✗ relapses"))
    ax.axvline(T, color=col, ls=":", lw=0.8, alpha=0.55)

ax.axhline(R_w * 10, color=S.C_GREY, ls="--", lw=1.1)
ax.text(40, R_w * 10 + 0.06, "original wound edge", fontsize=8, color="#777777")
ax.annotate("front advances,\nthen RECEDES", (14, 2.35), fontsize=9.5,
            weight="bold", color=S.C_BAD, ha="center")
ax.set_xlabel("time  (days)")
ax.set_ylabel("front position  $r_f$  (mm)")
ax.set_title("Wound radius 2 mm, chronic baseline", loc="left", fontsize=10)
ax.legend(loc="lower right", fontsize=8.4)
ax.set_ylim(0, 2.75)
ax.text(28, 2.63, "trace leaves the domain: complete die-back",
        fontsize=7.6, color="#888888", style="italic")

# ---- panel: minimum pulse vs wound radius ----------------------------------
radii_mm = np.array([1.0, 2.0, 3.0, 5.0])
mins = [min_pulse(r / 10) for r in radii_mm]
ax2.plot(radii_mm, mins, "o-", color=S.C_FGF, ms=7)
for r, t in zip(radii_mm, mins):
    ax2.annotate(f"{t:.0f} d", (r, t), textcoords="offset points",
                 xytext=(7, -12), fontsize=8.6, color=S.C_FGF, weight="bold")
ax2.set_xlabel("wound radius  $R_w$  (mm)")
ax2.set_ylabel("minimum pulse for permanent closure  (days)")
ax2.set_title("Below this line, the wound reopens", loc="left", fontsize=10)
ax2.set_xlim(0.5, 5.6); ax2.set_ylim(0, max(mins) * 1.2)
S.tier4_note(ax2, "every rate here is a Tier-4 closure")
S.save(fig, "fig09_recession")

np.save(os.path.join(os.path.dirname(__file__), "_minpulse.npy"),
        np.vstack([radii_mm, mins]))
for r, t in zip(radii_mm, mins):
    print(f"  R_w = {r:.0f} mm -> minimum pulse {t:.1f} d")
