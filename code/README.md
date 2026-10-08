# Mathematical Model V3 — figure and solver code

**HKU iGEM 2026 · Dry Lab · TEAM REFERENCE COPY — not a wiki deliverable**

This bundle regenerates every figure on the Modelling wiki page from source.
It exists so that anyone on the team (and anyone judging us) can check a number
rather than take it on trust.

---

## Run it

```bash
pip install numpy scipy matplotlib        # that is the entire dependency list
python3 verify.py                         # ~40 s — do this first
python3 make_figures.py                   # ~8 min — writes ../assets/*.svg + *.png
python3 make_figures.py 07 13             # or just the figures you want
```

`verify.py` must print **ALL CHECKS PASSED** before any number here is quoted
anywhere. It checks the release solver against the closed-form Crank series
(second-order convergence), the Fisher–KPP solver against the analytic front
speed, and the two-domain solver four ways (steady state, partition recovery,
first-order decay, mass budget).

---

## Layout

| file | what it is |
|---|---|
| `params.py` | **the single source of truth.** No value is hard-coded in a solver. Change a number here and every dependent figure changes. |
| `v3model.py` | transport coefficients, Model 1 release (closed form + method of lines), Models 3–4 occupancy and dose inversion, Model 5 NO/nitrite, Model 6 Fisher–KPP |
| `twodomain.py` | monolithic gel + tissue solver with partition, flux matching and the Dirichlet systemic sink |
| `sa_run.py` | Morris screening and Sobol indices, implemented directly (no SALib needed) |
| `verify.py` | the regression suite |
| `make_figures.py` | regenerates everything |
| `figures/figNN_*.py` | one script per figure, linked from that figure's caption on the wiki |
| `data/` | where the wet-lab CSVs go |

---

## Things that are easy to get wrong

These are in the code as comments too, but they are worth saying once here.

1. **Porosity cancels out of the PDE.** It appears in exactly two places:
   interfacial flux matching, and converting loaded mass to concentration.
   Applying it inside the diffusivity *and* in the concentration basis
   (as V2 did) under-predicts spreading.

2. **With retention, the loading is the TOTAL, not the free concentration.**
   `C_free,0 = C_total,0 / R_V`. Use `twodomain.free_from_total()`. Treating
   the loading as free silently loads `R_V`× more drug than solubility allows
   and over-predicts sustained release by orders of magnitude.

3. **Clearance is distributed *or* at the deep boundary, never both.** Both
   represent drug leaving into the circulation.

4. **Track the closure front at `0.5 n_ss`, not `0.5 K`.** With a death term the
   state behind the front is `n_ss = K(1 − d/r_p) < K`, so a threshold fixed at
   `0.5 K` is unreachable exactly in the chronic regime the death term exists to
   represent. When `r_p` and `d` vary in time, pin the threshold
   (`fixed_threshold=`) or you will report a spurious front.

5. **The release corner is singular.** The initial data are discontinuous at
   `(x,t) = (0,0)`, so `∂C/∂x|₀` diverges as `t^(−1/2)`. Crank–Nicolson rings
   on fine grids; `ReleaseSolver` uses Rannacher startup (a few backward-Euler
   half-steps) to damp it. Sample the release flux on a log grid, or use the
   retained-mass form, or you will over-count early release several-fold.

6. **Sobol indices are computed on the raw output, not its logarithm.** Taking
   logs of a near-multiplicative model makes it near-additive, pushes ΣS₁ to ≈1
   and hides the interaction structure the analysis exists to measure.

---

## ⚠ Open items — reconcile these before the wiki freeze

This code was written to reproduce the V3 equation-reference document. Most
numbers match it exactly. **Four do not, and someone needs to decide which
version is right.**

| quantity | equation-reference doc | this code | why they differ |
|---|---|---|---|
| days above target at `R_V = 10`, `k_prot = 1e-6` | 9.5 d | 11.2 d (band 4.5–25.9 d) | The duration is very sensitive to `L_tis` (2–4 mm) and `D_t` (0.7–1.8e-6), which the doc quotes only as ranges. Over those ranges the answer spans 4.5–25.9 d. **The doc should state which values produced 9.5 d.** The *location* of the optimum and the collapse beyond `R_V ≈ 30` are robust to both. |
| minimum pulse for closure | 5 / 12 / 20 / 30 d | 4.4 / 8.9 / 13.4 / 22.0 d | Same ordering and scale, ~20 % lower. Most likely the front-threshold convention and the closure criterion. Worth one hour with both scripts side by side. |
| largest single-patch wound | 1.5 mm radius | 2.5 mm (range 1.0–5.0) | Follows directly from the two rows above — our coverage is longer *and* our pulse requirement is shorter, so the crossover moves out. |
| Sobol `δ` and `k_prot,V` | S_T = 0.433 / 0.145 | 0.228 / 0.003 | **Our SA uses analytic surrogates, not the full PDE** (49 k PDE solves is not practical). The surrogates are documented in `sa_run.py`. Both analyses agree on the two conclusions that matter — `K_D,P` dominates, and `k_prot,V` and `D_n` are *not* dominant — and on ΣS₁ (0.374 vs 0.389). **If the team's own `sa_run.py` + SALib results exist, use those on the wiki and delete ours.** |

Everything else — `D_eff`, all four release times, the 2.8× separation, the
Damköhler numbers, `L_d` = 65–217 µm, the attenuation range, the 55.6 µM free
loading, 556 µM at `R_V = 10` (75 % of the ceiling), the peak tissue
concentrations, radial/Cartesian closure ratios 0.953 and 0.896 — reproduces the
document to within rounding.

---

## Figure 13 needs real data

`figures/fig13_nitrite_fit.py` is the only place the model meets data already in
hand, and it is currently rendering a **synthetic placeholder with a warning
watermark**. To make it real:

1. Put the measured V14-dose vs nitrite curve in
   `data/nitrite_dose_response.csv` (see `nitrite_dose_response_TEMPLATE.csv`).
2. `python3 make_figures.py 13`
3. The watermark disappears by itself and the fitted `IC₅₀^NO` and `n_app`
   appear with their standard errors.

Record the LPS concentration used. An IC₅₀ is only interpretable at a stated
`[L]`, and without an LPS titration it is `K_D,app`, not `K_D,P`.

**Do not publish the placeholder.**
