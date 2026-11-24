#!/usr/bin/env python3
"""Per-gene mean log1p(CP10K) over all cells, from viewer_data/expr.h5.

Used to draw expression-matched control genes for gene-set (module) scores
(Seurat AddModuleScore / scanpy score_genes style). Row order = genes.parquet = expr.h5.
Writes viewer_data/gene_means.parquet.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent.parent / 'viewer_data'
BLOCK = 2000


def main() -> int:
    genes = pd.read_parquet(OUT / 'genes.parquet', columns=['gene_id', 'n_cells_detected'])
    with h5py.File(OUT / 'expr.h5', 'r') as h:
        n_cells = int(h.attrs['n_cells'])
        indptr = h['indptr'][:]
        data = h['data']
        assert len(indptr) == len(genes) + 1
        sums = np.zeros(len(genes))
        for g0 in range(0, len(genes), BLOCK):
            g1 = min(g0 + BLOCK, len(genes))
            lo, hi = int(indptr[g0]), int(indptr[g1])
            d = data[lo:hi].astype(np.float64)
            cs = np.concatenate([[0.0], np.cumsum(d)])
            sums[g0:g1] = cs[indptr[g0 + 1:g1 + 1] - lo] - cs[indptr[g0:g1] - lo]
            print(f'[{time.strftime("%H:%M:%S")}] {g1:,}/{len(genes):,}', flush=True)
        nnz = np.diff(indptr)
        assert (nnz == genes.n_cells_detected.to_numpy()).all(), 'nnz per gene != n_cells_detected'
        # spot check against a direct dense read
        rng = np.random.default_rng(0)
        for g in rng.choice(len(genes), 20, replace=False):
            v = np.zeros(n_cells, np.float64)
            v[h['indices'][indptr[g]:indptr[g + 1]]] = data[indptr[g]:indptr[g + 1]]
            assert np.isclose(v.mean(), sums[g] / n_cells), g

    out = pd.DataFrame({'gene_id': genes.gene_id, 'mean_log1p_cp10k': (sums / n_cells).astype('float32')})
    out.to_parquet(OUT / 'gene_means.parquet', index=False)
    print(out.mean_log1p_cp10k.describe())
    return 0


if __name__ == '__main__':
    sys.exit(main())
