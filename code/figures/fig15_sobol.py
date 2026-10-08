"""
F15 -- Global sensitivity: the analysis refutes two of our own predictions.

Paired first-order (S1) and total-order (ST) Sobol indices for the two outputs
that matter.  The gap between S1 and ST for a factor IS its interaction
strength.

Two factors are marked because an earlier draft predicted they would dominate
and the computed indices say otherwise: k_prot,V for delivered dose and D_n for
closure.  Both are refuted.  (k_prot,V is still the most useful parameter to
ACT on -- see the design section -- because variance and leverage are different
questions.)
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
import style as S
import sa_run

S.apply()
R = sa_run.run(N=4096, verbose=False)
np.save(os.path.join(os.path.dirname(__file__), "_sobol.npy"),
        np.array([R["A"]["S1"].sum(), R["B"]["S1"].sum(), R["B"]["stall"]]))

REFUTED = {"A": "$k_{prot,V}$", "B": "$D_n$"}
fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.6))

for ax, tag, title, col in zip(axes, ("A", "B"),
                               ("A:  scaled dose $\\psi$ delivered",
                                "B:  closure time $t_c$"),
                               (S.C_V14, S.C_FGF)):
    d = R[tag]
    order = np.argsort(d["ST"])
    names = [d["names"][i] for i in order]
    S1, ST = d["S1"][order], d["ST"][order]
    y = np.arange(len(names))

    ax.barh(y + 0.19, ST, height=0.36, color=col, alpha=0.42,
            label="$S_T$  total order")
    ax.barh(y - 0.19, np.clip(S1, 0, None), height=0.36, color=col,
            label="$S_1$  first order")
    ax.set_yticks(y, names)
    for i, nm in enumerate(names):
        if nm == REFUTED[tag]:
            ax.get_yticklabels()[i].set_color(S.C_BAD)
            ax.get_yticklabels()[i].set_fontweight("bold")
            ax.annotate("predicted to dominate — REFUTED",
                        (ST[i] + 0.02, i), fontsize=8.4, color=S.C_BAD,
                        weight="bold", va="center")
    ax.set_xlabel("Sobol index")
    ax.set_title(title, loc="left", fontsize=10)
    ax.set_xlim(0, max(ST.max() * 1.45, 0.9))
    ax.legend(loc="upper right", fontsize=8.4)
    ax.text(0.985, 0.74, f"$\\sum S_1$ = {d['S1'].sum():.3f}\n"
            f"{1-d['S1'].sum():.0%} of the variance\nlives in interactions",
            transform=ax.transAxes, ha="right", va="top", fontsize=8.6,
            weight="bold", color="#444444", linespacing=1.4)

fig.suptitle("What to measure next — and what turned out not to matter",
             fontsize=11.5, weight="bold")
S.save(fig, "fig15_sobol")
print(f"  A: sum S1 = {R['A']['S1'].sum():.3f}")
print(f"  B: sum S1 = {R['B']['S1'].sum():.3f}, stall fraction {R['B']['stall']:.1%}")
