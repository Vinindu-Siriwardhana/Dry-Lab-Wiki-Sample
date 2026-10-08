"""
F17 -- Why the interactive slider is NOT the uncertainty analysis.

A one-parameter-at-a-time sweep is exactly what a slider performs, and it can
only ever reach the first-order variance.  Everything in the interactions is
invisible to it.  This figure sits next to the slider on purpose: the limitation
of our own showpiece is worth volunteering.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import style as S

S.apply()
sA, sB, _ = np.load(os.path.join(os.path.dirname(__file__), "_sobol.npy"))

fig, ax = plt.subplots(figsize=(6.0, 2.5))
labels = ["A: delivered dose $\\psi$", "B: closure time $t_c$"]
vals = [sA, sB]
y = np.arange(2)
ax.barh(y, vals, color=S.C_V14, height=0.52, label="reachable by a slider ($\\sum S_1$)")
ax.barh(y, 1 - np.array(vals), left=vals, color=S.C_BAD, alpha=0.32,
        height=0.52, label="interactions a slider cannot reach")
for i, v in enumerate(vals):
    ax.text(v / 2, i, f"{v:.0%}", ha="center", va="center", color="white",
            weight="bold", fontsize=10)
    ax.text(v + (1 - v) / 2, i, f"{1-v:.0%}", ha="center", va="center",
            color="#4a1a55", weight="bold", fontsize=10)
ax.set_yticks(y, labels)
ax.set_xlim(0, 1); ax.set_xticks([0, .25, .5, .75, 1], ["0", "25 %", "50 %", "75 %", "100 %"])
ax.set_xlabel("share of output variance")
ax.set_title("A slider is presentation, not analysis", fontsize=10)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.42), ncol=2, fontsize=8.2)
ax.grid(False)
S.save(fig, "fig17_slider")
