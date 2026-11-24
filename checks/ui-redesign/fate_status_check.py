"""Verify reviewed Fate hover statuses, including contexts absent from the screen."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'app'))
from data_store import ViewerData
from fate_view import FateCells, SOURCES, TIMES, STATUS_TEXT

data = ViewerData(ROOT / 'viewer_data')
fate = FateCells(ROOT / 'viewer_data')
values = data.expression(data.resolve_gene('MtrunA17_Chr5g0448621'))
# An absent fitted context must explain ineligibility on its descriptive curve.
missing = fate.figure(values, 'gene', 'detection', 'pooled', status={})
assert missing.data
assert all(STATUS_TEXT[None] in t.hovertemplate and t.line.dash == 'dash' for t in missing.data)
# Supplying a screen status still produces that status and its existing line style.
for status, dash in [('figure7', 'solid'), ('below_threshold', 'solid'), ('excluded_fit', 'dot')]:
    known = fate.figure(values, 'gene', 'detection', 'pooled',
                        status={(s, t): status for s in SOURCES for t in TIMES})
    assert all(STATUS_TEXT[status] in t.hovertemplate and t.line.dash == dash for t in known.data)
# Descriptive gene-set curves must not inherit a gene eligibility warning.
set_fig = fate.figure(values, 'set', 'detection', 'pooled')
assert all('not eligible' not in t.hovertemplate and t.line.dash == 'solid' for t in set_fig.data)
print('PASS: missing fitted contexts show approved eligibility explanation; known statuses and line styles preserved; gene-set curves remain descriptive.')
