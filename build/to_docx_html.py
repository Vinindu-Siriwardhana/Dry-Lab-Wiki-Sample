"""Turn build/body.html into a pandoc-ready HTML document.

The page's equations live as raw LaTeX inside <div class="eq">. Pandoc's HTML
reader recognises <span class="math display">\\[...\\]</span>, so rewriting them
that way makes pandoc emit real OMML equations in the .docx rather than text.
"""
import io, re, sys, html as htmllib, datetime

SRC, OUT = sys.argv[1], sys.argv[2]
s = io.open(SRC, encoding='utf-8').read()

s = s.replace('__BUILDDATE__', datetime.date.today().strftime('%d %B %Y'))

# ---- 1. equations -----------------------------------------------------------
def eq_block(m):
    # Leave the entities alone. The source escapes '<', '>' and the alignment
    # '&' inside equations, and pandoc's HTML reader decodes them back to LaTeX
    # before tex_math_dollars parses it. Decoding here instead would put a bare
    # '<' into the HTML and the parser would swallow the rest of the equation.
    latex = m.group(1).strip()
    return '<p>$$' + latex + '$$</p>'

s = re.sub(r'<div class="eq">(.*?)</div>', eq_block, s, flags=re.S)

# equation captions -> small italic line above the maths
s = re.sub(r'<div class="eqlabel">(.*?)</div>',
           r'<p><em>\1</em></p>', s, flags=re.S)
s = s.replace('<div class="eqwrap">', '<div>')

# ---- 2. drop chrome that means nothing on paper -----------------------------
s = re.sub(r'<nav class="toc">.*?</nav>', '', s, flags=re.S)
s = re.sub(r'<div class="num">.*?</div>', '', s, flags=re.S)
s = re.sub(r'<a class="dl".*?</a>', '', s, flags=re.S)

# ---- 3. section headings: fold the eyebrow into the h2 ----------------------
def head(m):
    eyebrow = re.sub(r'<[^>]+>', '', m.group(1)).strip()
    title = m.group(2).strip()
    if eyebrow and eyebrow.lower() not in title.lower():
        return '<h1>' + eyebrow + ' &#8212; ' + title + '</h1>'
    return '<h1>' + title + '</h1>'

s = re.sub(r'<p class="eyebrow">(.*?)</p>\s*<h2>(.*?)</h2>', head, s, flags=re.S)

# ---- 4. callouts -> blockquotes --------------------------------------------
s = re.sub(r'<div class="claim[^"]*">', '<blockquote>', s)
# the matching </div> is handled by the generic div strip below, so convert the
# close tags by walking the string instead
out, depth, i = [], None, 0
BQ = '<blockquote>'
while True:
    j = s.find(BQ, i)
    if j < 0:
        out.append(s[i:])
        break
    out.append(s[i:j])
    out.append(BQ)
    k = j + len(BQ)
    d = 1
    while d:
        nd = s.find('<div', k)
        ne = s.find('</div>', k)
        if ne < 0:
            break
        if 0 <= nd < ne:
            d += 1
            k = nd + 4
        else:
            d -= 1
            if d == 0:
                out.append(s[j + len(BQ):ne])
                out.append('</blockquote>')
                k = ne + len('</div>')
                break
            k = ne + len('</div>')
    i = k
s = ''.join(out)

# ---- 5. <details> -> heading + body ----------------------------------------
s = re.sub(r'<details class="[^"]*"><summary>(.*?)</summary>',
           r'<h3>\1</h3>', s, flags=re.S)
s = s.replace('</details>', '')
s = s.replace('<div class="inner">', '<div>')

# ---- 6. demote remaining headings so h2 sits under the section h1 ----------
s = re.sub(r'<h3( [^>]*)?>', '<h2>', s)
s = s.replace('</h3>', '</h2>')

# ---- 7. hero ---------------------------------------------------------------
s = re.sub(r'<p class="kicker">(.*?)</p>', r'<p><em>\1</em></p>', s, flags=re.S)
s = re.sub(r'<div class="keystrip">(.*?)</div>\s*</div></header>',
           r'<ul>\1</ul></header>', s, flags=re.S)
s = re.sub(r'<div><div class="v">(.*?)</div><div class="k">(.*?)</div></div>',
           r'<li><strong>\1</strong> &#8212; \2</li>', s, flags=re.S)

# ---- 8. strip the remaining layout divs & spans ----------------------------
s = re.sub(r'<div class="tablewrap">', '<div>', s)
s = re.sub(r'<span class="caveat">(.*?)</span>', r'<em>\1</em>', s, flags=re.S)
s = re.sub(r'<span class="(fid|src|tag|st)">(.*?)</span>', r'\2', s, flags=re.S)
s = re.sub(r'<div class="ref-banner">', '<blockquote>', s)

HEAD = ('<!doctype html><html><head><meta charset="utf-8">'
        '<title>HKU iGEM 2026 - Dry Lab Modelling</title></head><body>')
io.open(OUT, 'w', encoding='utf-8').write(HEAD + s + '</body></html>')

n_eq = s.count('$$') // 2
print('wrote %s  (%d display equations)' % (OUT, n_eq))
