"""
F13 -- The one place the model meets data already in hand.

The Griess assay measures ACCUMULATED NITRITE, so the quantity to fit is the
two-parameter reduced form

    [NO2-](C)/[NO2-](0) = 1 / (1 + (C/IC50)^n)

and nothing else.  The mechanistic chain from C_V to nitrite carries at least
fourteen parameters against a sigmoid with roughly three degrees of freedom, so
the internal constants are NOT identifiable from this dataset and fitted values
for them must not be reported.

HOW TO USE THIS SCRIPT
----------------------
Put the measured curve in  code/data/nitrite_dose_response.csv
(see nitrite_dose_response_TEMPLATE.csv for the columns) and re-run.  Until that
file exists the script renders a clearly-marked SYNTHETIC placeholder so the
page layout can be finished -- the placeholder must NOT be published.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import style as S, params as P, v3model as m

S.apply()
HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "data", "nitrite_dose_response.csv")

SYNTHETIC = not os.path.exists(CSV)
if SYNTHETIC:
    rng = np.random.default_rng(7)
    C = np.repeat([0.0, 0.5, 1.0, 5.0, 10.0, 25.0, 50.0], 3)
    true_ic50, true_n = 7.5, 1.25
    y = m.nitrite_reduced(np.maximum(C, 1e-6), true_ic50, true_n)
    y = np.where(C == 0, 1.0, y) * (1 + rng.normal(0, 0.055, C.size))
else:
    import csv as _csv
    rows = [r for r in _csv.DictReader(
        [ln for ln in open(CSV) if not ln.startswith("#")])]
    C = np.array([float(r["C_V14_uM"]) for r in rows])
    raw = np.array([float(r["nitrite_uM"]) for r in rows])
    y = raw / raw[C == 0].mean()

# ---- fit the reduced form --------------------------------------------------
mask = C > 0
f = lambda c, ic50, n: m.nitrite_reduced(c, ic50, n)
p0 = [np.median(C[mask]), 1.0]
popt, pcov = curve_fit(f, C[mask], y[mask], p0=p0,
                       bounds=([1e-3, 0.3], [1e4, 6.0]), maxfev=20000)
perr = np.sqrt(np.diag(pcov))
ic50, n_app = popt

# ---- plot ------------------------------------------------------------------
fig, (ax, axr) = plt.subplots(2, 1, figsize=(6.4, 5.3), sharex=True,
                              gridspec_kw=dict(height_ratios=[3, 1], hspace=0.08))
cc = np.logspace(np.log10(max(C[mask].min() / 3, 1e-3)),
                 np.log10(C.max() * 3), 400)

# group replicates for display
uC = np.unique(C[mask])
mu = np.array([y[C == c].mean() for c in uC])
sd = np.array([y[C == c].std(ddof=1) if (C == c).sum() > 1 else 0 for c in uC])

ax.plot(cc, f(cc, *popt), color=S.C_V14, lw=2.2,
        label="reduced form, 2 parameters")
res_lo = f(cc, ic50 - perr[0], n_app)
res_hi = f(cc, ic50 + perr[0], n_app)
ax.fill_between(cc, res_lo, res_hi, color=S.C_V14, alpha=0.15, lw=0)
ax.errorbar(uC, mu, yerr=sd, fmt="o", color=S.C_DERIVED, ms=6, capsize=3,
            lw=1.2, label="Griess, mean ± s.d.")
ax.axhline(0.5, color=S.C_GREY, ls=":", lw=1.2)
ax.axvline(ic50, color=S.C_THRESH, ls="--", lw=1.4)
ax.annotate(f"$IC_{{50}}^{{NO}}$ = {ic50:.2f} ± {perr[0]:.2f} µM\n"
            f"$n_{{app}}$ = {n_app:.2f} ± {perr[1]:.2f}",
            (ic50, 0.5), textcoords="offset points", xytext=(12, 26),
            fontsize=10, weight="bold", color=S.C_THRESH)
ax.set_xscale("log")
ax.set_ylabel("nitrite, normalised to untreated")
ax.set_ylim(0, 1.25)
ax.legend(loc="lower left", fontsize=8.5)
ax.set_title("V14 suppresses NO production in RAW 264.7")

resid = y[mask] - f(C[mask], *popt)
axr.axhline(0, color=S.C_GREY, lw=1)
axr.plot(C[mask], resid, "o", color=S.C_DERIVED, ms=5)
axr.set_ylabel("resid.")
axr.set_xlabel("V14 concentration applied  (µM)")
axr.set_ylim(-max(0.12, np.abs(resid).max() * 1.4),
             max(0.12, np.abs(resid).max() * 1.4))

if SYNTHETIC:
    for a_ in (ax, axr):
        a_.patch.set_alpha(0)
    fig.text(0.5, 0.55, "SYNTHETIC\nPLACEHOLDER", fontsize=42, color="#D33",
             alpha=0.16, ha="center", va="center", rotation=22, weight="bold")
    ax.text(0.5, 1.14, "⚠  NOT REAL DATA — drop the Griess CSV into "
            "code/data/ and re-run", transform=ax.transAxes, ha="center",
            fontsize=9.5, weight="bold", color="#C00")

S.save(fig, "fig13_nitrite_fit")
print(f"  {'SYNTHETIC' if SYNTHETIC else 'MEASURED'} data")
print(f"  IC50 = {ic50:.3f} +/- {perr[0]:.3f} uM,  n_app = {n_app:.3f} +/- {perr[1]:.3f}")
print(f"  RMS residual = {np.sqrt((resid**2).mean()):.4f}")
