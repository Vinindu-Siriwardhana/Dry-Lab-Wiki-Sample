# Dry-Lab-Wiki-Sample

HKU iGEM 2026 — Dry Lab Modelling wiki page (V3, "without heparin" scenario).

## Live site

**https://vinindu-siriwardhana.github.io/Dry-Lab-Wiki-Sample/**

## What's here

- `index.html` — the single-page Modelling wiki (self-contained, theme-aware).
- `assets/` — all 18 figures (SVG + PNG) and the downloadable model-code bundle.
- `code/` — the figure + solver source code (`verify.py`, `make_figures.py`, per-figure scripts, params, data template). See `code/README.md`.
- `build/` — the page generator (`build_page.py`, `page_parts.py`, `body.html`, `page.css`) used to assemble `index.html`.

## Regenerate the page

```bash
pip install numpy scipy matplotlib
cd code && python3 verify.py          # must print ALL CHECKS PASSED
cd code && python3 make_figures.py    # rewrites ../assets/*.svg + *.png
cd build && python3 build_page.py     # rewrites ../index.html
```
