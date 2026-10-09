"""
style.py — the single visual system for every V3 model figure.

Import this first in every figure script. It fixes the colour semantics so a
reader who learns the convention on one figure reads every other figure in two
seconds:

    V14 arm      blue      #2166AC
    FGF2-G3 arm  red       #B2182B
    derived      charcoal  #3A3A3A
    V2/deleted   grey      #BEBEBE
    threshold    orange    #E08214

    solid line   quantity is constrained by data, a derivation, or geometry
    dashed line  quantity rests on a Tier-4 order-of-magnitude prior

Blue/red is safe under deuteranopia and protanopia; nothing in these figures is
distinguished by colour alone -- line style, marker and direct labels always
carry the same information.

HKU iGEM 2026 Dry Lab.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cycler import cycler
import os

# --- semantic colours -------------------------------------------------------
C_V14      = "#2166AC"   # V14 arm
C_V14_L    = "#92C5DE"   # V14, light (bands, secondary curves)
C_FGF      = "#B2182B"   # FGF2-G3 arm
C_FGF_L    = "#F4A582"   # FGF, light
C_DERIVED  = "#3A3A3A"   # derived / neutral quantities
C_GREY     = "#BEBEBE"   # V2 material, deleted branches
C_THRESH   = "#E08214"   # thresholds, targets, ceilings
C_GOOD     = "#4D9221"   # feasible region
C_BAD      = "#762A83"   # infeasible / stalled
C_BG       = "#FFFFFF"

# line style semantics
LS_CONSTRAINED = "-"
LS_TIER4       = "--"

FIGDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "assets")


def apply():
    """Apply the shared rcParams. Call once at the top of a figure script."""
    plt.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "savefig.facecolor": C_BG,
        "figure.facecolor": C_BG,
        "axes.facecolor": C_BG,
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "axes.labelsize": 9.5,
        "axes.edgecolor": "#444444",
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": "#DDDDDD",
        "grid.linewidth": 0.6,
        "grid.alpha": 0.9,
        "legend.frameon": False,
        "legend.fontsize": 8.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "xtick.color": "#333333",
        "ytick.color": "#333333",
        "lines.linewidth": 1.9,
        "axes.prop_cycle": cycler(color=[C_V14, C_FGF, C_DERIVED,
                                         C_THRESH, C_GOOD, C_BAD]),
        "svg.fonttype": "path",      # embed glyphs: identical look on any machine
    })


def save(fig, name):
    """Write <assets>/<name>.svg and .png, and echo the path."""
    os.makedirs(FIGDIR, exist_ok=True)
    for ext in ("svg", "png"):
        p = os.path.join(FIGDIR, f"{name}.{ext}")
        fig.savefig(p, format=ext)
        if ext == "svg":
            _strip_dtd(p)
    print(f"  wrote {name}.svg / .png")
    plt.close(fig)


def _strip_dtd(path):
    """
    Remove the XML prolog and DOCTYPE that matplotlib writes.

    Wiki and artifact hosting both reject XML supporting files carrying DTD
    machinery, and nothing in these figures needs it -- the SVG renders
    identically without it.
    """
    with open(path, encoding="utf-8") as fh:
        txt = fh.read()
    i = txt.find("<svg")
    if i > 0:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(txt[i:])


def tier4_note(ax, text="dashed = Tier-4 prior", loc=(0.98, 0.02)):
    ax.text(*loc, text, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=7.5, color="#777777", style="italic")
