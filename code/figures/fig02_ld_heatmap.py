"""
F2 -- Penetration depth across the whole swept parameter range.

Shows the Act-1 conclusion does not depend on one parameter choice: L_d stays
far below the stack thickness everywhere in the (1/tau, t_half) box.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import style as S, params as P, v3model as m

S.apply()
inv_tau = np.array(P.INV_TAU)
t_half = np.array(P.NO_HALFLIFE)
H = m.amsden_H(0.15, P.XI_NM)
G = np.array([[np.sqrt(P.D0_NO * H * it / (np.log(2) / th)) * 1e4
               for th in t_half] for it in inv_tau])

fig, ax = plt.subplots(figsize=(4.3, 3.1))
im = ax.imshow(G, cmap="Blues", vmin=0, vmax=G.max() * 1.2, aspect="auto")
for i in range(len(inv_tau)):
    for j in range(len(t_half)):
        ax.text(j, i, f"{G[i,j]:.0f}", ha="center", va="center",
                color="white" if G[i, j] > 0.6 * G.max() else "#1a1a1a",
                fontsize=10, weight="bold")
ax.set_xticks(range(len(t_half)), [f"{t:g} s" for t in t_half])
ax.set_yticks(range(len(inv_tau)), [f"{v:.1f}" for v in inv_tau])
ax.set_xlabel("NO biological half-life $t_{1/2}$")
ax.set_ylabel(r"$1/\tau$  (tortuosity group)")
ax.set_title("Penetration depth $L_d$  (µm)", fontsize=10)
ax.grid(False)
cb = fig.colorbar(im, ax=ax, fraction=0.046)
cb.outline.set_visible(False)
S.save(fig, "fig02_ld_heatmap")
print(f"  L_d spans {G.min():.0f}-{G.max():.0f} um; stack is 500-3000 um")
