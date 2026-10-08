"""
F18 -- Where every parameter comes from, module by module.

T1 measured in house / T2 mechanistic correlation / T3 named literature
analogue / T4 order-of-magnitude prior.  The hatched overlay counts parameters
with NO planned experiment at all.

Models 1-2 are predictive now.  Model 6 is a structured hypothesis, and the bar
says so without anyone having to read the appendix.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import style as S, params as P

S.apply()
mods = list(P.PROVENANCE.keys())[::-1]
tiers = ["T1", "T2", "T3", "T4"]
cols = ["#1A5E2A", "#4393C3", "#F4A582", "#B2182B"]
labels = {"T1": "T1  measured in house", "T2": "T2  mechanistic correlation",
          "T3": "T3  literature analogue", "T4": "T4  order-of-magnitude prior"}

fig, ax = plt.subplots(figsize=(7.4, 3.6))
y = np.arange(len(mods))
left = np.zeros(len(mods))
for t, c in zip(tiers, cols):
    v = np.array([P.PROVENANCE[mm][t] for mm in mods])
    ax.barh(y, v, left=left, color=c, height=0.62, label=labels[t])
    for i, (vv, ll) in enumerate(zip(v, left)):
        if vv:
            ax.text(ll + vv / 2, i, str(vv), ha="center", va="center",
                    fontsize=8.2, color="white", weight="bold")
    left += v

unp = np.array([P.PROVENANCE[mm]["unplanned"] for mm in mods])
tot = left
ax.barh(y, unp, left=tot - unp, height=0.62, facecolor="none",
        edgecolor="#222222", hatch="///", lw=1.0,
        label="no planned experiment")

ax.set_yticks(y, mods)
ax.set_xlabel("number of parameters")
ax.set_title("Parameter provenance by module", fontsize=10.5)
ax.legend(loc="lower right", fontsize=8, ncol=1)
ax.set_xlim(0, tot.max() * 1.45)
S.save(fig, "fig18_provenance")
