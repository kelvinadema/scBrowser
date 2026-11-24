import sys
import time
from pathlib import Path

import numpy as np

APP = Path(__file__).resolve().parents[1] / 'app'
sys.path.insert(0, str(APP))
from data_store import ViewerData  # noqa: E402
from fate_assoc import FateAssociation  # noqa: E402
from gene_map import GeneMap  # noqa: E402
from gene_sets import GeneSetScorer, fate_heatmap  # noqa: E402

D = ViewerData(APP.parent / 'viewer_data')
S = GeneSetScorer(D, APP.parent / 'viewer_data')
F = FateAssociation(APP.parent / 'viewer_data')
G = GeneMap(APP.parent / 'viewer_data')

# small set: symbols + an id + a bogus token
rows, missing = S.resolve('MtNIN MtNF-YA1 MtENOD11, MtrunA17_Chr5g0448621; notagene')
print('resolved', [D.genes.label.iloc[r] for r in rows], 'missing', missing)

t = time.time(); v, nc = S.score(rows); t_small = time.time() - t
ctrl = S.controls(rows)
ref = (np.mean([D.expression(r) for r in sorted(rows)], 0)
       - np.mean([D.expression(r) for r in ctrl], 0))
print(f'small set: {len(rows)} genes, {nc} controls, {t_small:.2f}s, max |diff vs dense| {np.abs(v - ref).max():.2e}')
assert np.allclose(v, ref, atol=1e-5)
assert not set(ctrl) & set(rows)
assert np.array_equal(S.controls(rows), ctrl), 'controls not deterministic'
# controls drawn from the same expression bins as the set genes
assert set(S.bins[ctrl]) == set(S.bins[rows])
# a set gene that sits inside a control panel must be removed from it
inpanel = int(S.panel_rows[S.panel_rows >= 0][0])
rows_p = rows + [inpanel]
vp, ncp = S.score(rows_p); cp = S.controls(rows_p)
assert inpanel not in cp
refp = np.mean([D.expression(r) for r in rows_p], 0) - np.mean([D.expression(r) for r in cp], 0)
print('panel-overlap set: max |diff|', float(np.abs(vp - refp).max()))
assert np.allclose(vp, refp, atol=1e-4)

# score separates nodule-like cells for an NCR/nodule set? use module M07 (cortex-high) vs cortex
m07 = list(G.genes[G.genes.module_id == 'M07'].gene_id)
rows7, _ = S.resolve(' '.join(m07))
t = time.time(); v7, nc7 = S.score(rows7); t_big = time.time() - t
t = time.time(); S.score(rows7); t_cached = time.time() - t
tissue = D.cells.ch5_tissue.astype(str)
by = {k: float(v7[(tissue == k).to_numpy()].mean()) for k in tissue.unique()}
print(f'M07: {len(rows7)} genes, {nc7} controls, {t_big:.1f}s (cached {t_cached:.3f}s)')
print('  mean score by tissue:', {k: round(x, 3) for k, x in sorted(by.items(), key=lambda kv: -kv[1])})

genes = [(D.genes.gene_id.iloc[r], D.genes.label.iloc[r]) for r in rows] + [('MtrunA17_Chr3g0134941', 'x')]
fig = fate_heatmap(F, genes)
m01 = list(G.genes[G.genes.module_id == 'M01'].gene_id)
rows1, _ = S.resolve(' '.join(m01))
t = time.time(); S.score(rows1); print(f'M01: {len(rows1)} genes, {time.time() - t:.1f}s')
print('fate heatmap', fig.data[0].type, len(fig.data[0].y), 'genes x', len(fig.data[0].x), 'contexts')
print('DONE')
