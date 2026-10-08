"""
F10 -- The device specification: which wounds this patch can close.

Joins the two halves of the answer.  The RISING line is what Model 6 demands:
the pulse duration needed for permanent closure at a given wound radius (F9).
The HORIZONTAL band is what Models 1-2 can supply: days above the 25 uM target
from a single patch at the optimal R_V ~ 10 (F7).

Where the demand curve crosses the supply band is the largest wound one patch
can close.  Everything larger needs periodic replacement.

The crossing is reported as a RANGE, not a point.  The supply band carries the
L_tis and D_t uncertainty of F7; the demand line rests entirely on Tier-4
coupling closures with no planned experiment.  The ORDERING and the existence of
a crossover are robust; the absolute radius is not.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import style as S, params as P

S.apply()
mp = np.load(os.path.join(os.path.dirname(__file__), "_minpulse.npy"))
radii, pulses = mp[0], mp[1]

# single-patch coverage at the R_V ~ 10 optimum (central, and band) -- from F7
COV_CEN, COV_LO, COV_HI = 11.2, 4.5, 25.9

fig, ax = plt.subplots(figsize=(7.0, 4.4))
rr = np.linspace(1.0, 5.0, 200)   # only where F9 actually computed
dem = np.interp(rr, radii, pulses)

ax.fill_between(rr, COV_LO, COV_HI, color=S.C_V14, alpha=0.14, lw=0,
                label=f"one patch supplies {COV_LO:.0f}–{COV_HI:.0f} d above target")
ax.axhline(COV_CEN, color=S.C_V14, lw=2.0,
           label=f"one patch, central estimate: {COV_CEN:.1f} d")
ax.plot(rr, dem, color=S.C_FGF, lw=2.2, ls="--",
        label="pulse required for permanent closure (Tier-4)")

cross = lambda cov: float(np.interp(cov, dem, rr))
x_cen, x_lo, x_hi = cross(COV_CEN), cross(COV_LO), cross(COV_HI)
x_hi = min(x_hi, 5.0)   # F9 was only computed out to 5 mm

ax.fill_between(rr, 0, dem, where=(rr <= x_cen), color=S.C_GOOD, alpha=0.12, lw=0)
ax.axvline(x_cen, color=S.C_DERIVED, lw=1.3)
ax.plot([x_lo, x_hi], [COV_LO, COV_HI], "|", color=S.C_DERIVED, ms=11, mew=1.6)
ax.plot([x_cen], [COV_CEN], "*", ms=16, color=S.C_THRESH, zorder=6)

ax.annotate(f"one patch closes wounds up to\n"
            f"$R_w$ ≈ {x_cen:.1f} mm  (range {x_lo:.1f}–{x_hi:.1f} mm)",
            (x_cen, COV_CEN), textcoords="offset points", xytext=(14, 18),
            fontsize=9.5, weight="bold", color=S.C_DERIVED)

for n, lo_r in ((2, x_cen), (3, 3.6)):
    pass
ax.text(1.5, 2.0, "ONE PATCH", fontsize=11, weight="bold", color="#3b7a1c",
        ha="center")
ax.text(4.0, 5.0, "REPLACEMENT REQUIRED\nevery ~%.0f days" % COV_CEN,
        fontsize=10, weight="bold", color=S.C_FGF, ha="center")

ax.set_xlabel("wound radius  $R_w$  (mm)")
ax.set_ylabel("days of therapeutic cover required / supplied")
ax.set_title("What this device can and cannot close")
ax.set_xlim(1.0, 5.0); ax.set_ylim(0, 28)
ax.legend(loc="upper left", fontsize=8.4)
S.save(fig, "fig10_envelope")
print(f"  single-patch crossover: central {x_cen:.2f} mm, range {x_lo:.2f}-{x_hi:.2f} mm")
for r in (1, 2, 3, 5):
    need = np.interp(r, radii, pulses)
    print(f"   R_w={r} mm: needs {need:.1f} d -> "
          f"{max(1, int(np.ceil(need/COV_CEN)))} patch(es) at the central estimate")
