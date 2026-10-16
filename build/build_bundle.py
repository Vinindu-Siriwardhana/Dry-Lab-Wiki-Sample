"""build_bundle.py -- rebuild assets/igem2026-hku-dry-lab-model-code.zip from the live sources.

Run after any change to the model package, a figure script or the wet-lab inputs:

    python build/build_bundle.py


The shipped zip had drifted: it carried the pre-correction sa_run.py (the one with
the spurious Lambda/delta factor), a params.py without the 145 uM target, a
verify.py that asserted the superseded GelMA gel, and only the retired fig01-18
scripts rather than the G-series that actually produces the figures on the page.

Layout written:

    igem2026-hku-model-v3-code/
        README.md  requirements.txt  constants_register.xlsx
        params.py  v3model.py  twodomain.py  sa_run.py  style.py  verify.py
        data/nitrite_dose_response_TEMPLATE.csv
        figures/common.py  figures/G*.py  figures/WL_*.py

The model modules sit ONCE, at the top level, so there is a single params.py in
the bundle.  figures/common.py is the project copy with its one path line
repointed at the parent directory.
"""
import io
import os
import re
import shutil
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))   # build -> _remote_check -> project folder
CODE = os.path.join(ROOT, '_remote_check', 'code')
NG = os.path.join(ROOT, 'New Graphs')
OUT = os.path.join(ROOT, '_remote_check', 'assets',
                   'igem2026-hku-dry-lab-model-code.zip')
TOP = 'igem2026-hku-model-v3-code'

MODEL = ['params.py', 'v3model.py', 'twodomain.py', 'sa_run.py', 'style.py']
FLAT = MODEL + ['verify.py', 'README.md', 'requirements.txt',
                'constants_register.xlsx']

written = []
with zipfile.ZipFile(OUT, 'w', zipfile.ZIP_DEFLATED) as z:
    for f in FLAT:
        src = os.path.join(CODE, f)
        if not os.path.exists(src):
            print('  !! missing', f)
            continue
        z.write(src, '%s/%s' % (TOP, f))
        written.append(f)

    tmpl = os.path.join(CODE, 'data', 'nitrite_dose_response_TEMPLATE.csv')
    if os.path.exists(tmpl):
        z.write(tmpl, '%s/data/nitrite_dose_response_TEMPLATE.csv' % TOP)
        written.append('data/nitrite_dose_response_TEMPLATE.csv')

    # common.py, repointed: model modules live one level up in the bundle
    common = io.open(os.path.join(NG, 'common.py'), encoding='utf-8').read()
    before = common
    common = common.replace(
        'sys.path.insert(0, os.path.join(HERE, "model"))',
        'sys.path.insert(0, ROOT)   # bundle layout: model modules sit one level up')
    common = common.replace(
        'The model code lives in ./model (params.py is the single source of truth;',
        'The model code lives one directory up (params.py is the single source of truth;')
    assert common != before, 'common.py path line not found -- layout changed?'
    z.writestr('%s/figures/common.py' % TOP, common)
    written.append('figures/common.py')

    for f in sorted(os.listdir(NG)):
        if re.match(r'^(G\d\d|WL)_.*\.py$', f):
            z.write(os.path.join(NG, f), '%s/figures/%s' % (TOP, f))
            written.append('figures/' + f)

    # raw wet-lab inputs, so the WL_* and G03 scripts run from the bundle as
    # shipped rather than failing on a missing path
    wl = os.path.join(ROOT, 'Wet Lab Data')
    for f in ('trials2_3_positive_gate_frequency_raw_data.csv',
              'trials2_3_student_raw_p_comparisons.csv',
              'Hydrogel BCA Sep 30.xlsx'):
        src = os.path.join(wl, f)
        if os.path.exists(src):
            z.write(src, '%s/Wet Lab Data/%s' % (TOP, f))
            written.append('Wet Lab Data/' + f)
        else:
            print('  !! missing wet-lab input', f)

print('wrote %s  (%d files, %.0f KB)'
      % (OUT, len(written), os.path.getsize(OUT) / 1024))
for f in written:
    print('   ', f)
