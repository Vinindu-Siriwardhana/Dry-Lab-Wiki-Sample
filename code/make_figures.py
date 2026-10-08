#!/usr/bin/env python3
"""
make_figures.py -- regenerate every figure on the Modelling wiki page.

    python3 make_figures.py            # all figures
    python3 make_figures.py 07 13      # just those two

Writes SVG (used by the page) and PNG (for slides and the poster) into
../assets/.  Every figure on the wiki comes from exactly one script in
figures/, and that script is linked from the figure's caption.

Run verify.py first if you have changed anything in v3model.py or twodomain.py.
"""
import os, sys, subprocess, time

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")

ORDER = [
    ("01", "fig01_no_attenuation.py",  "NO attenuation -- why sensing failed"),
    ("02", "fig02_ld_heatmap.py",      "penetration depth across the swept range"),
    ("03", "fig03_release.py",         "release is a bolus; 2.8x arrival separation"),
    ("04", "fig04_verification.py",    "solver verification (run this one first)"),
    ("05", "fig05_tissue_profiles.py", "V14 depth profiles in tissue"),
    ("06", "fig06_hindrance.py",       "steric hindrance vs mesh size"),
    ("07", "fig07_retention.py",       "THE retention optimum"),
    ("08", "fig08_solubility.py",      "solubility ceiling"),
    ("09", "fig09_recession.py",       "front recession + minimum pulse (SLOW, ~4 min)"),
    ("10", "fig10_envelope.py",        "device envelope (needs 09 first)"),
    ("11", "fig11_universal_dose.py",  "universal scaled dose-response"),
    ("12", "fig12_dose_correction.py", "V2 dose-equation correction"),
    ("13", "fig13_nitrite_fit.py",     "nitrite fit -- NEEDS THE WET-LAB CSV"),
    ("14", "fig14_no_vs_nitrite.py",   "free NO vs cumulative nitrite"),
    ("15", "fig15_sobol.py",           "Sobol tornado (SLOW, ~1 min)"),
    ("16", "fig16_uncertainty.py",     "closure under the joint prior"),
    ("17", "fig17_slider.py",          "slider deficit (needs 15 first)"),
    ("18", "fig18_provenance.py",      "parameter provenance"),
]


def main(which=None):
    t0 = time.time()
    fails = []
    for num, script, desc in ORDER:
        if which and num not in which:
            continue
        print(f"\n[{num}] {desc}\n     {script}")
        r = subprocess.run([sys.executable, os.path.join(FIGDIR, script)],
                           capture_output=True, text=True)
        print(r.stdout.rstrip())
        if r.returncode != 0:
            fails.append(script)
            print(r.stderr.rstrip()[-1500:])
    print(f"\ndone in {time.time()-t0:.0f} s")
    if fails:
        print("FAILED:", ", ".join(fails))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or None))
