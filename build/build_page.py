#!/usr/bin/env python3
"""
build_page.py -- assemble the Modelling wiki page.

Reads page.css + css/*.css + body.html + widgets/*.html + js/*.js, injects the
model map and the live code listings (straight out of code/, so the page can
never drift from what actually ran), and writes ../index.html.
"""
import os, html, datetime, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CODE = os.path.join(ROOT, "code")
sys.path.insert(0, HERE)
from page_parts import MODEL_MAP, LEGEND

FILES = [
    ("params.py",       "Every parameter, with its provenance tier. The single source of truth."),
    ("v3model.py",      "Transport coefficients, Crank series, release solver, occupancy and dose inversion, NO readout, Fisher&#8211;KPP."),
    ("twodomain.py",    "Monolithic gel + tissue solver: partition, flux matching, Dirichlet systemic sink."),
    ("sa_run.py",       "Morris screening and Sobol indices, implemented directly."),
    ("verify.py",       "The regression suite. Run this first."),
    ("make_figures.py", "Regenerates every figure into assets/."),
    ("style.py",        "The shared visual system &#8212; colour semantics and line-style meanings."),
]
FIGS = [(f, "") for f in sorted(os.listdir(os.path.join(CODE, "figures")))
        if f.endswith(".py")]

CSS = ("page.css", "css/theme.css", "css/widgets.css")
JS = ("sim.js", "plot.js", "ui.js", "i3.js", "i1.js", "i2.js")
NL = "\n"


def read(p):
    with open(p, encoding="utf-8") as fh:
        return fh.read()


def listing(relpath, label=None):
    src = read(os.path.join(CODE, relpath))
    lines = src.count(NL) + 1
    doc = ""
    if src.lstrip().startswith('"""'):
        body = src.split('"""')[1].strip().splitlines()
        doc = body[0].strip() if body else ""
    return (
        f'<details class="code"><summary><span>{html.escape(relpath)}</span>'
        f'<span class="meta">{lines} lines'
        + (f' &#183; {html.escape(doc)[:72]}' if doc else "")
        + '</span></summary>'
        f'<div style="padding:9px 15px 0;text-align:right">'
        f'<button class="copybtn" type="button">copy</button></div>'
        f'<pre><code>{html.escape(src)}</code></pre></details>'
    )


HEAD = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>HKU iGEM 2026 Modelling</title>
<meta name="description" content="Mathematical modelling of a dual-release wound-healing patch: release, tissue transport, receptor binding and wound closure.">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:ital,opsz,wght@0,9..144,500;0,9..144,600;1,9..144,400&family=Inter:wght@400;500;600;700&display=swap" rel="stylesheet">
<style>
%s
</style>
</head>
<body>
"""


def main():
    css = NL.join(read(os.path.join(HERE, p)) for p in CSS)
    js = NL.join(read(os.path.join(HERE, "js", f)) for f in JS)
    body = read(os.path.join(HERE, "body.html"))

    cards = NL.join(
        f'<div class="filecard"><div class="fn">{html.escape(n)}</div>'
        f'<div class="fd">{d}</div></div>' for n, d in FILES)
    cards += ('<div class="filecard"><div class="fn">figures/fig01 &#8230; fig18</div>'
              '<div class="fd">One script per figure, named in that figure\'s caption. '
              'Run any of them on its own.</div></div>')

    listings = NL.join(listing(n) for n, _ in FILES)
    listings += NL + '<h3 style="font-size:1rem">figures/</h3>' + NL
    listings += NL.join(listing(os.path.join("figures", f)) for f, _ in FIGS)

    body = (body.replace("__MODELMAP__", MODEL_MAP)
                .replace("__LEGEND__", LEGEND)
                .replace("__FILECARDS__", cards)
                .replace("__CODELISTINGS__", listings)
                .replace("__I1__", read(os.path.join(HERE, "widgets", "i1.html")))
                .replace("__I2__", read(os.path.join(HERE, "widgets", "i2.html")))
                .replace("__I3__", read(os.path.join(HERE, "widgets", "i3.html")))
                .replace("__BUILDDATE__",
                         datetime.date.today().strftime("%d %B %Y")))

    out = HEAD.replace("%s", css) + body + NL + "<script>" + NL + js + NL + "</script>" + NL \
        + "</body>" + NL + "</html>" + NL

    path = os.path.join(ROOT, "index.html")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(out)
    print(f"wrote {path}  ({len(out)/1024:.0f} KB)")


if __name__ == "__main__":
    main()
