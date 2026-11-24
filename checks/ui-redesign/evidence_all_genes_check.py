"""Evidence heatmaps retain full set membership beyond the former 100-gene cap."""
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'app'))
import app as a
ids=[g for g in a.DATA.genes.gene_id if a.DISSOC.row(g) is not None][:500]
genes=[(g,g) for g in ids]
d=a.DISSOC.heatmap(genes)
assert list(d.data[0].y)==ids
np.testing.assert_array_equal(np.array(d.data[0].z),[[a.DISSOC.row(g).protoplast_lfc,a.DISSOC.row(g).nuclei_lfc] for g in ids])
f=a.fate_heatmap(a.FATE,genes)
assert all(list(t.y)==ids for t in f.data)
assert f.layout.height<=3600 and d.layout.height<=3600
# The original 400-gene Expression heatmap cap remains independent.
assert a.HEAT_MAX==400
print('PASS: 500 genes retained in both evidence heatmaps, dissociation values unchanged, bounded figure height, independent Expression limit preserved.')
