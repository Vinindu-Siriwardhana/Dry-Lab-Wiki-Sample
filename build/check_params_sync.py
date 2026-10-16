#!/usr/bin/env python3
"""
check_params_sync.py -- fail loudly if any hand-maintained copy of a model
parameter has drifted from the master, New Graphs/model/params.py.

The model's numbers live in more than one place because the browser cannot
import Python. Every copy is checked here:

  1. _remote_check/code/*.py      shipped copy of the model package
                                  (must be byte-identical to New Graphs/model)
  2. build/js/sim.js  `var P`     the interactives' parameter block
  3. build/js/i3.js               Interactive 1's signalling / NO constants
  4. constants_register.xlsx      both copies: identical to each other, and the
                                  design target equal to params.C_TARGET_V14
  5. New Graphs/WL_kdp_bound_25uM.py  panel dose == sim.js PANEL_DOSE_UM

build_page.py runs this first and refuses to build on a mismatch.
Run it directly after changing any parameter:

    python build/check_params_sync.py
"""
import importlib.util
import io
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WIKI = os.path.dirname(HERE)                    # _remote_check
ROOT = os.path.dirname(WIKI)                    # project folder
MASTER = os.path.join(ROOT, "New Graphs", "model")
SHIPPED = os.path.join(WIKI, "code")
MODULES = ("params.py", "v3model.py", "twodomain.py", "sa_run.py", "style.py")

problems = []


def bad(msg):
    problems.append(msg)


def load_params():
    spec = importlib.util.spec_from_file_location("params_master", os.path.join(MASTER, "params.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def same(a, b, rel=1e-6):
    return math.isclose(a, b, rel_tol=rel, abs_tol=1e-30)


# ---- 1. shipped model package ----------------------------------------------
for f in MODULES:
    a, b = os.path.join(MASTER, f), os.path.join(SHIPPED, f)
    if not os.path.exists(b):
        bad("code/%s is missing (copy it from New Graphs/model/)" % f)
        continue
    ta = io.open(a, encoding="utf-8").read().replace("\r\n", "\n")
    tb = io.open(b, encoding="utf-8").read().replace("\r\n", "\n")
    if ta != tb:
        bad("code/%s differs from New Graphs/model/%s -- copy the master over it" % (f, f))

P = load_params()
UM = P.UM


# ---- 2. sim.js parameter block ---------------------------------------------
def js_block(path, start):
    t = io.open(path, encoding="utf-8").read()
    i = t.index(start) + len(start)
    j = t.index("};", i)
    body = re.sub(r"//[^\n]*", "", t[i:j])
    out = {}
    for part in body.split(","):
        if ":" not in part:
            continue
        k, v = part.split(":", 1)
        out[k.strip()] = v.strip()
    return out


def js_eval(expr):
    env = {"HOUR": 3600.0, "DAY": 86400.0, "LN2": math.log(2)}
    expr = expr.replace("Math.LN2", "LN2")
    if not re.fullmatch(r"[0-9eE.+\-*/() A-Z_0-9]+", expr):
        raise ValueError("unparseable expression %r" % expr)
    return float(eval(expr, {"__builtins__": {}}, env))


SIM = os.path.join(HERE, "js", "sim.js")
sim = {k: js_eval(v) for k, v in js_block(SIM, "var P = {").items()}
SIM_MAP = {                     # sim.js key : value expected from params.py
    "D0_V14": P.D0_V14, "A_V14": P.A_V14_NM, "A_FIBRE": P.A_FIBRE_NM, "XI": P.XI_NM,
    "EPS_GEL": P.EPS_GEL, "EPS_TIS": P.EPS_TIS, "D_TIS": P.D_TIS_V14, "L_GEL": P.L_GEL,
    "L_TIS": P.L_TIS, "K_P": P.K_P, "T_HALF_V14": P.T_HALF_V14, "C_SOL_UM": P.C_SOLUBILITY / UM,
    "DELTA": P.DELTA, "KD_P_UM": P.KD_P_LOWER / UM, "TARGET_UM": P.C_TARGET_V14 / UM,
    "D_N": P.D_N, "R0": P.R0, "D0": P.D0_DEATH, "DD": P.DELTA_D, "E_MAX": P.E_MAX,
    "EC50": P.EC50_F, "BETA": P.BETA_R, "K_IR": P.K_I_R, "K_ID": P.K_I_D, "S_I": P.S_I,
    "KE0": P.KE0, "W": P.W_SMOOTH, "KD_L_UM": P.KD_L / UM, "L_UM": P.L_LPS / UM,
    "C_LOAD_V14_UM": P.C_LOAD_V14 / UM, "F_REL_V14": P.F_REL_V14,
    "C_LOAD_FGF_NGML": P.C_LOAD_FGF_NGML, "F_REL_FGF": P.F_REL_FGF,
}
for k, want in SIM_MAP.items():
    if k not in sim:
        bad("sim.js P is missing %s (params.py value %g)" % (k, want))
    elif not same(sim[k], want):
        bad("sim.js P.%s = %g but params.py gives %g" % (k, sim[k], want))

# ---- 3. i3.js constants ----------------------------------------------------
I3 = io.open(os.path.join(HERE, "js", "i3.js"), encoding="utf-8").read()
m = re.search(r"var KD_L = ([0-9.eE+-]+), GAMMA = ([0-9.eE+-]+), H = ([0-9.eE+-]+);", I3)
if not m:
    bad("i3.js: could not find the `var KD_L = .., GAMMA = .., H = ..;` line")
else:
    for name, got, want in (("KD_L", m.group(1), P.KD_L / UM), ("GAMMA", m.group(2), P.GAMMA),
                            ("H", m.group(3), P.HILL_H)):
        if not same(float(got), want):
            bad("i3.js %s = %s but params.py gives %g" % (name, got, want))
no = {k: js_eval(v) for k, v in js_block(os.path.join(HERE, "js", "i3.js"), "var NO = {").items()}
for k, want in (("q", P.Q_CHRONIC), ("K", P.K_INOS_REL), ("p", P.P_INOS),
                ("kturn", P.K_TURN), ("phi", P.PHI_NO2)):
    if k not in no or not same(no[k], want):
        bad("i3.js NO.%s = %s but params.py gives %g" % (k, no.get(k), want))

# ---- 4. constants register -------------------------------------------------
try:
    import openpyxl
    regs = [os.path.join(ROOT, "constants_register (1).xlsx"), os.path.join(SHIPPED, "constants_register.xlsx")]
    sheets = []
    for r in regs:
        if not os.path.exists(r):
            bad("register missing: %s" % os.path.relpath(r, ROOT))
            continue
        wb = openpyxl.load_workbook(r, data_only=True)
        sheets.append({ws.title: [[c for c in row] for row in ws.iter_rows(values_only=True)] for ws in wb})
        tgt = wb["8_model7_design"]["C2"].value
        if not same(float(tgt), P.C_TARGET_V14 / UM):
            bad("%s: 8_model7_design!C2 (C_V_target) = %s but params.py gives %g"
                % (os.path.relpath(r, ROOT), tgt, P.C_TARGET_V14 / UM))
    if len(sheets) == 2 and sheets[0] != sheets[1]:
        bad("the two constants_register copies differ -- edit one, then copy it over the other")
except ImportError:
    print("  (openpyxl not installed: register not checked)")

# ---- 5. flow-panel dose ----------------------------------------------------
wl = io.open(os.path.join(ROOT, "New Graphs", "WL_kdp_bound_25uM.py"), encoding="utf-8").read()
mm = re.search(r"^C_V = ([0-9.]+)", wl, re.M)
if mm and "PANEL_DOSE_UM" in sim and not same(float(mm.group(1)), sim["PANEL_DOSE_UM"]):
    bad("WL_kdp_bound_25uM.py C_V = %s but sim.js PANEL_DOSE_UM = %g" % (mm.group(1), sim["PANEL_DOSE_UM"]))

# ---- report ----------------------------------------------------------------
if problems:
    print("PARAMETER COPIES OUT OF SYNC (%d):" % len(problems))
    for p in problems:
        print("  - " + p)
    print("Master is New Graphs/model/params.py. Fix the copies, then re-run.")
    sys.exit(1)
print("parameter copies in sync: code/ package, sim.js (%d keys), i3.js, register x2, panel dose" % len(SIM_MAP))
