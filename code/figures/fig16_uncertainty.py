"""
F16 -- Closure under the full joint prior, not a single best guess.

Monte Carlo over the Model 6 priors of sa_run.PRIORS_B.  The ribbon is the
5th-95th percentile of the area-closure fraction; the flat traces along the
bottom are the samples with r_eff <= 0, which never close at all.

That stall fraction is the strongest available argument for keeping the death
term: a model without it would have returned a finite closure time for every
one of those samples.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import style as S, params as P, v3model as m
import sa_run

S.apply()
DAY = P.DAY
rng = np.random.default_rng(3)
N = 3000
U = rng.random((N, len(sa_run.PRIORS_B)))
X = sa_run._scale(U, sa_run.PRIORS_B)
tc, stalled = sa_run.model_B(X)

r0, Ce, KIr, EC50, d0, Ie, beta, dd, Emax, Dn = X.T
rp = r0 * (1 + Emax * Ce / (EC50 + Ce)) * (1 - beta * Ie / (KIr + Ie))
d = d0 + dd * Ie ** 2 / (0.3 ** 2 + Ie ** 2)
reff = rp - d
c = m.front_speed(Dn, reff)

t = np.linspace(0, 60 * DAY, 400)
R_w = P.R_W_DEFAULT
rf = np.clip(R_w - np.outer(c, t) / 1.10, 0, None)     # 1.10 = transient factor
Ac = 1 - (rf / R_w) ** 2
Ac[stalled] = 0.0

lo, med, hi = np.percentile(Ac, [5, 50, 95], axis=0)

fig, ax = plt.subplots(figsize=(6.6, 4.1))
ax.fill_between(t / DAY, lo, hi, color=S.C_FGF, alpha=0.18, lw=0,
                label="5th – 95th percentile of the joint prior")
ax.plot(t / DAY, med, color=S.C_FGF, lw=2.2, label="median")
idx = np.where(stalled)[0][:60]
for i in idx:
    ax.plot(t / DAY, Ac[i], color=S.C_BAD, lw=0.6, alpha=0.25)
ax.plot([], [], color=S.C_BAD, lw=1.2, alpha=0.6,
        label=f"never close ($r_p \\leq d$):  {stalled.mean():.1%} of samples")

ax.set_xlabel("time  (days)")
ax.set_ylabel("area closure  $\\mathcal{A}_c(t)$")
ax.set_ylim(0, 1.02); ax.set_xlim(0, 60)
ax.set_title(f"{stalled.mean():.0%} of the prior never heals at all  ($R_w$ = {R_w*10:.0f} mm)")
ax.legend(loc="lower right", fontsize=8.4)
S.tier4_note(ax, "every Model 6 parameter is Tier-4", (0.98, 0.52))
S.save(fig, "fig16_uncertainty")
print(f"  stall fraction {stalled.mean():.1%}; median closure "
      f"{np.median(tc[~stalled]):.1f} d")
