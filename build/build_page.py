#!/usr/bin/env python3
"""
build_page.py -- assemble the Modelling wiki page.

Reads page.css + css/theme.css + body.html + js/ui.js and writes ../index.html.
"""
import os, datetime, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)


CSS = ("page.css", "css/theme.css")
JS = ("ui.js",)
NL = "\n"


def read(p):
    with open(p, encoding="utf-8") as fh:
        return fh.read()


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

    body = (body.replace("__BUILDDATE__",
                         datetime.date.today().strftime("%d %B %Y")))

    out = HEAD.replace("%s", css) + body + NL + "<script>" + NL + js + NL + "</script>" + NL \
        + "</body>" + NL + "</html>" + NL

    path = os.path.join(ROOT, "index.html")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(out)
    print(f"wrote {path}  ({len(out)/1024:.0f} KB)")


if __name__ == "__main__":
    main()
