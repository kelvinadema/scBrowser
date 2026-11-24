"""Real-data checks for omitting untested genes from dissociation heatmaps."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'app'))
from dissociation import Dissociation, NAMES

dissoc = Dissociation(ROOT / 'viewer_data')
viewer = pd.read_parquet(ROOT / 'viewer_data/genes.parquet')
missing = viewer.loc[~viewer.gene_id.isin(dissoc.table.index), 'gene_id'].head(15).tolist()
assert len(missing) == 15
available = dissoc.table.loc[dissoc.table.index.isin(viewer.gene_id)]
responsive = available.index[available.protoplast_call | available.nuclei_call][0]
nonresponsive = available.index[~(available.protoplast_call | available.nuclei_call)][0]
genes = [(missing[0], 'missing first'), (responsive, 'responsive')]
genes += [(gid, 'missing') for gid in missing[1:]]
genes += [(nonresponsive, 'nonresponsive')]
fig = dissoc.heatmap(genes)
assert list(fig.data[0].y) == ['responsive', 'nonresponsive']
expected = [[float(dissoc.row(g)[f'{k}_lfc']) for k in NAMES] for g in (responsive, nonresponsive)]
np.testing.assert_array_equal(np.array(fig.data[0].z), expected)
assert list(fig.data[0].text[1]) == ['', ''], 'Tested non-responsive genes must remain visible'
assert fig.layout.height == 220, 'Height should reflect only plotted rows'
assert dissoc.heatmap([]) is None
assert dissoc.heatmap([(g, g) for g in missing]) is None
print('PASS: missing genes omitted, tested non-responsive genes retained, values/order preserved, compact height, empty/all-untested plots absent.')
