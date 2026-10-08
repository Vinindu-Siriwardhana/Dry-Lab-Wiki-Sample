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

DAY = 86400.0
ok = True


def check(name, got, want, tol, fmt="{:.3e}"):
    global ok
    good = abs(got - want) <= tol
    ok &= good
    print(f"  [{'PASS' if good else 'FAIL'}] {name:<46} "
          f"{fmt.format(got)}  (expect {fmt.format(want)} ± {fmt.format(tol)})")


print("\n1. transport coefficients")
check("D_eff V14 at xi=10 nm (cm^2/s)",
      m.D_eff(3.6e-6, 0.87, 10, 0.91), 3.25e-6, 5e-8)
check("D_eff FGF2-G3 at xi=10 nm (cm^2/s)",
      m.D_eff(1.5e-6, 2.34, 10, 0.91), 1.15e-6, 5e-8)

print("\n2. release timescales (closed form)")
Dv = m.D_eff(3.6e-6, 0.87, 10, 0.91)
Df = m.D_eff(1.5e-6, 2.34, 10, 0.91)
check("95% release, V14, 500 um (min)",
      m.release_time(0.95, Dv, 500e-4) / 60, 14.5, 0.3, "{:.1f}")
check("95% release, FGF2-G3, 500 um (min)",
      m.release_time(0.95, Df, 500e-4) / 60, 40.9, 0.4, "{:.1f}")
check("arrival separation factor",
      m.release_time(0.95, Df, 500e-4) / m.release_time(0.95, Dv, 500e-4),
      2.82, 0.05, "{:.2f}")

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
check("max error at nx=800", errs[-1], 8.2e-7, 3e-7)

print("\n4. Fisher-KPP front speed vs 2 sqrt(D r_eff)")
sol = m.FisherKPP(1e-9, R_inf=1.6, nr=1600, geometry="cartesian")
res = sol.run(1.4, 1.0 / DAY, 0.0, t_end=60 * DAY, record_every=100)
t, rf = res["t"], res["r_front"]
msk = (rf > 0.35) & (rf < 1.05)
c_nu = -np.polyfit(t[msk], rf[msk], 1)[0]
c_an = m.front_speed(1e-9, 1.0 / DAY)
check("front speed relative error", abs(c_nu - c_an) / c_an, 0.0, 0.05, "{:.2%}")

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
