"""
F10 -- The device specification: how often the patch must be changed.

Joins the two halves of the answer.

(a) The rising line is what Model 6 demands: the CONTINUOUS therapeutic pulse
    needed for permanent closure at a given wound radius (F9).  The horizontal
    band is what Models 1-2 supply: hours above the 25 uM target from one
    as-built patch (500 um alginate, R_V = 1; F7).  The gap is about two orders
    of magnitude, so no wound in the 1-5 mm range closes on one patch.

(b) Repeat dosing, run through the register's effect compartment
    (k_e0 = 1/18 h).  A patch applied every `interval` hours drives the effect
    state s(t) toward 1 for the hours it keeps tissue above target and lets it
    relax toward 0 otherwise:

        ds/dt = k_e0 (u(t) - s),   u = 1 while the latest patch is above target

    r_p and d interpolate between the chronic and therapeutic rates of F9 with
    s, and the front is tracked at the same pinned threshold.  This gives the
    longest change interval that still closes the wound -- the number the
    protocol actually needs.

Everything on the demand side rests on Tier-4 closures with no planned
experiment.  The ordering and the existence of a maximum interval are robust;
the absolute numbers are not.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import style as S, params as P, v3model as m

S.apply()
DAY = P.DAY
HERE = os.path.dirname(__file__)
mp = np.load(os.path.join(HERE, "_minpulse.npy"))
radii, pulses = mp[0], mp[1]
COV_CEN, COV_LO, COV_HI = np.load(os.path.join(HERE, "_coverage.npy"))   # days

THR = 0.5 * P.K_CARRY * (1 - P.D_PULSE / P.R_PULSE)      # as in F9


def closure_time(interval_h, R_w_cm, cov_d=COV_CEN, t_end_days=120, dt=60.0):
    """Closure time (days) when a fresh patch is applied every interval_h."""
    t = np.arange(0.0, t_end_days * DAY + dt, dt)
    u = ((t / 3600.0) % interval_h) < cov_d * 24.0
    s = np.zeros_like(t)
    a = np.exp(-dt * P.KE0)
    for i in range(1, len(t)):                 # exact step of the linear ODE
        s[i] = u[i - 1] + (s[i - 1] - u[i - 1]) * a
    rp = lambda tt: P.R_CHRONIC + (P.R_PULSE - P.R_CHRONIC) * s[min(int(tt / dt), len(s) - 1)]
    dd = lambda tt: P.D_CHRONIC + (P.D_PULSE - P.D_CHRONIC) * s[min(int(tt / dt), len(s) - 1)]
    sol = m.FisherKPP(P.D_N, K=P.K_CARRY, R_inf=R_w_cm * 2.2, nr=500)
    res = sol.run(R_w_cm, rp, dd, t_end=t_end_days * DAY, record_every=40,
                  fixed_threshold=THR, w=P.W_SMOOTH)
    return res["t_c"] / DAY


def max_interval(R_w_cm, lo=None, hi=24.0, tol=0.25):
    lo = lo or COV_CEN * 24.0
    if np.isfinite(closure_time(hi, R_w_cm)):
        return hi
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if np.isfinite(closure_time(mid, R_w_cm)):
            lo = mid
        else:
            hi = mid
    return lo


fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.0, 4.3),
                              gridspec_kw=dict(width_ratios=[1, 1.15]))

# ---- (a) supply vs demand ---------------------------------------------------
rr = np.linspace(radii.min(), radii.max(), 200)
dem = np.interp(rr, radii, pulses)
ax.fill_between(rr, COV_LO, COV_HI, color=S.C_V14, alpha=0.16, lw=0,
                label=f"one patch: {COV_LO*24:.1f}–{COV_HI*24:.1f} h above target")
ax.axhline(COV_CEN, color=S.C_V14, lw=2.0,
           label=f"one patch, central: {COV_CEN*24:.1f} h")
ax.plot(rr, dem, color=S.C_FGF, lw=2.2, ls="--",
        label="continuous pulse needed to close (Tier-4)")
for r, p in zip(radii, pulses):
    ax.annotate(f"{p/COV_CEN:.0f}×", (r, p), textcoords="offset points",
                xytext=(2, 8), fontsize=8.2, color=S.C_DERIVED, ha="left" if r < 3 else "right",
                weight="bold")
ax.set_yscale("log")
ax.set_ylim(0.02, 100)
ax.set_xlim(radii.min() - 0.2, radii.max() + 0.2)
ax.set_xlabel("wound radius  $R_w$  (mm)")
ax.set_ylabel("days of therapeutic cover")
ax.set_title("(a)  one patch: hours; closure: weeks", loc="left",
             fontsize=10)
ax.legend(loc="center right", fontsize=7.8)

# ---- (b) closure time vs change interval -----------------------------------
intervals = np.array([COV_CEN * 24, 8, 10, 12, 15, 18, 20, 22, 24])
summary = {}
for R_mm, col in ((1.0, S.C_V14_L), (2.0, S.C_V14), (3.0, "#0B3C5D")):
    tc = np.array([closure_time(h, R_mm / 10) for h in intervals])
    hmax = max_interval(R_mm / 10)
    summary[R_mm] = (tc, hmax)
    fin = np.isfinite(tc)
    ax2.plot(intervals[fin], tc[fin], "o-", color=col, ms=5,
             label=f"$R_w$ = {R_mm:g} mm")
    ax2.plot(intervals[~fin], np.full((~fin).sum(), 118), "x", color=col, ms=7)
    ax2.axvline(hmax, color=col, ls=":", lw=1.0)

daily_closes = all(np.isfinite(closure_time(P.REAPPLY_D * 24, R / 10))
                   for R in (1.0, 2.0, 3.0))
ax2.axvline(P.REAPPLY_D * 24, color=S.C_BAD, lw=1.4)
ax2.text(P.REAPPLY_D * 24 * 0.87, 104,
         "daily dressing\nchange:\n" + ("closes" if daily_closes else "never closes"),
         ha="right", va="top", fontsize=8.4, color=S.C_BAD, weight="bold")
hm = summary[2.0][1]
ax2.text(hm * 0.97, 58, f"longest interval\nthat closes:\n{hm:.1f} h", ha="right",
         fontsize=8.4, color=S.C_DERIVED, weight="bold")
ax2.text(0.02, 0.97, "× = never closes (within 120 d)", transform=ax2.transAxes,
         va="top", fontsize=7.8, color="#777777")
ax2.set_xscale("log")
ax2.set_xticks([6, 8, 12, 18, 24], ["6", "8", "12", "18", "24"])
ax2.xaxis.set_minor_formatter(plt.NullFormatter())
ax2.set_xlim(5.5, 27)
ax2.set_ylim(0, 125)
ax2.set_xlabel("patch change interval  (h)")
ax2.set_ylabel("closure time  $t_c$  (days)")
ax2.set_title("(b)  the protocol the model asks for", loc="left", fontsize=10)
ax2.legend(loc="upper left", bbox_to_anchor=(0.0, 0.9), fontsize=8)
S.tier4_note(ax2, "every closure rate here is Tier-4", (0.70, 0.02))
S.save(fig, "fig10_envelope")

for r, p in zip(radii, pulses):
    print(f"   R_w={r:.0f} mm: needs {p:.1f} d continuous = {p/COV_CEN:.0f} patch-lifetimes")
for R_mm, (tc, hmax) in summary.items():
    row = " ".join(f"{h:.1f}h:{t:.1f}d" for h, t in zip(intervals, tc))
    print(f"   R_w={R_mm:g} mm: max interval {hmax:.1f} h | {row}")
