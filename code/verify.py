#!/usr/bin/env python3
"""
verify.py -- run this before trusting any number out of this package.

Four checks on the release solver and four on the two-domain solver.  Nothing
in the results section is worth reading if these do not pass.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v3model as m
import twodomain as td
import params as P

DAY = 86400.0
ok = True


def check(name, got, want, tol, fmt="{:.3e}"):
    global ok
    good = abs(got - want) <= tol
    ok &= good
    print(f"  [{'PASS' if good else 'FAIL'}] {name:<46} "
          f"{fmt.format(got)}  (expect {fmt.format(want)} ± {fmt.format(tol)})")


# Every expected value below is for the gel AS BUILT -- 1% w/v calcium
# alginate, 2.5 mm (was 500 um until 2026-10-16), xi = 14 nm, eps = 0.99 --
# and the inputs come from
# params.py rather than being repeated here, so that changing the register
# changes these checks too.  An earlier version of this file hard-coded the
# superseded 5% GelMA gel (xi = 10 nm, eps = 0.91) and so printed ALL CHECKS
# PASSED while asserting release times that disagreed with the documents.
print("\n1. transport coefficients  (alginate as built: xi = %g nm, eps = %g)"
      % (P.XI_NM, P.EPS_GEL))
Dv = m.D_eff(P.D0_V14, P.A_V14_NM, P.XI_NM, P.EPS_GEL, P.A_FIBRE_NM)
Df = m.D_eff(P.D0_FGF, P.A_FGF_NM, P.XI_NM, P.EPS_GEL, P.A_FIBRE_NM)
check("D_eff V14 (cm^2/s)", Dv, 3.478e-6, 5e-9)
check("D_eff FGF2-G3 (cm^2/s)", Df, 1.327e-6, 5e-9)

print("\n2. release timescales (closed form)")
check("95% release, V14, 2.5 mm (min)",
      m.release_time(0.95, Dv, P.L_GEL) / 60, 338.1, 0.5, "{:.2f}")
check("95% release, FGF2-G3, 2.5 mm (min)",
      m.release_time(0.95, Df, P.L_GEL) / 60, 886.3, 0.5, "{:.2f}")
check("arrival separation factor",
      m.release_time(0.95, Df, P.L_GEL) / m.release_time(0.95, Dv, P.L_GEL),
      2.62, 0.02, "{:.2f}")

print("\n3. release solver vs the Crank series (spatial order)")
times = np.array([30., 60., 120., 300., 600., 900.])
errs = []
for nx in (100, 200, 400, 800):
    s = m.ReleaseSolver(Dv, 500e-4, nx)
    Cn = s.solve(1.0, times, dt=0.4 * (100 / nx))
    Ce = m.crank_profile(s.x, times, 1.0, Dv, 500e-4)
    errs.append(np.max(np.abs(Cn - Ce)))
order = np.log2(errs[-2] / errs[-1])
check("observed spatial order", order, 2.0, 0.08, "{:.2f}")
check("max error at nx=800", errs[-1], 8.08e-7, 5e-8)

print("\n4. Fisher-KPP front speed vs 2 sqrt(D r_eff)")
sol = m.FisherKPP(1e-9, R_inf=1.6, nr=1600, geometry="cartesian")
res = sol.run(1.4, 1.0 / DAY, 0.0, t_end=60 * DAY, record_every=100)
t, rf = res["t"], res["r_front"]
msk = (rf > 0.35) & (rf < 1.05)
c_nu = -np.polyfit(t[msk], rf[msk], 1)[0]
c_an = m.front_speed(1e-9, 1.0 / DAY)
# Across the four cases reported in the documents the error spans 1.8-6.2%, the
# numerical speed always sitting just below the analytic one as expected for a
# discretised pulled front.  This case is the 2.05% one.
check("front speed relative error", abs(c_nu - c_an) / c_an, 0.0, 0.07, "{:.2%}")

print("\n5. front failure is reproduced (chronic regime)")
sol = m.FisherKPP(1e-9, R_inf=1.1, nr=700)
stall = sol.run(0.5, 0.8 / DAY, 1.0 / DAY, t_end=60 * DAY, record_every=50)
print(f"  [{'PASS' if not stall['closed'] else 'FAIL'}] r_p < d  =>  wound never closes"
      f"                 t_c = {stall['t_c']}")
ok &= not stall["closed"]

print("\n6. two-domain solver")
v = td.verify(verbose=False)
check("tissue steady state, order", v["steady_order"][1], 2.0, 0.1, "{:.2f}")
check("partition K_p recovery (rel)", v["partition_err"], 0.0, 5e-3, "{:.2%}")
check("gel first-order decay (abs)", v["decay_err"], 0.0, 1e-4)
check("mass left in domain after 20 d, no decay",
      v["mass_remaining_frac"], 0.0, 5e-3, "{:.2%}")

print("\n" + ("ALL CHECKS PASSED" if ok else "SOMETHING FAILED -- do not report results"))
sys.exit(0 if ok else 1)
