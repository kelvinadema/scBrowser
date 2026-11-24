"""Check the actual Cells caption renderer against real sample-filtered expression."""
import ast
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'app'))
import app as a

# Execute the renderer body outside Shiny, supplying its reactive dependencies.
source = ast.parse((ROOT / 'app/app.py').read_text())
renderer = next(n for n in ast.walk(source) if isinstance(n, ast.FunctionDef) and n.name == 'plot_title')
renderer.decorator_list = []
module = ast.fix_missing_locations(ast.Module(body=[renderer], type_ignores=[]))
rows = [a.DATA.resolve_gene(g) for g in ('MtrunA17_Chr5g0448621', 'MtrunA17_Chr3g0134941')]
assert all(r is not None for r in rows)
score, controls = a.SCORER.score(rows)
state = {'mask': None, 'kind': 'gene', 'row': rows[0]}
env = dict(np=np, DATA=a.DATA, COPY=a.COPY, caption=lambda text: text,
           input=SimpleNamespace(colour_by=lambda: 'selection', sel_kind=lambda: state['kind']),
           cell_mask=lambda: state['mask'], gene_row=lambda: state['row'],
           geneset_score=lambda: (score, controls), geneset_parsed=lambda: (rows, []))
exec(compile(module, str(ROOT / 'app/app.py'), 'exec'), env)
render = env['plot_title']
# Independent sparse-expression reference, bypassing caption's dense expression API.
indptr, indices, _ = a.DATA._arrays()
detected_indices = indices[int(indptr[rows[0]]):int(indptr[rows[0] + 1])]
masks = [None, np.isin(a.SAMPLE_COLUMN, a.TIP_SAMPLES),
         a.SAMPLE_COLUMN == a.SAMPLE_IDS[0], np.zeros(a.DATA.n_cells, dtype=bool)]
for mask in masks:
    state['mask'] = mask
    n = a.DATA.n_cells if mask is None else int(mask.sum())
    for kind in ('gene', 'set'):
        state['kind'] = kind
        text = render()
        if not n:
            assert text == 'No cells in the selected samples.'
            continue
        scope = 'across the full dataset' if mask is None else 'in the selected samples'
        assert scope in text
        if kind == 'gene':
            detected = len(detected_indices) if mask is None else int(mask[detected_indices].sum())
            assert f'detected in {detected:,} cells ({detected / n * 100:.1f}%)' in text
        else:
            positive = np.count_nonzero((score > 0) if mask is None else ((score > 0) & mask))
            assert f'{positive / n * 100:.1f}%' in text
            assert f'{controls:,} expression-matched background genes' in text
state.update(mask=None, kind='gene', row=None)
assert render() == 'No gene(s) selected'
print('PASS: actual renderer; gene and set statistics for all samples, root tips, one sample and no samples; no-gene prompt; fixed control count.')
