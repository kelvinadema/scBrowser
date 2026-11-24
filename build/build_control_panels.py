#!/usr/bin/env python3
"""Fixed control-gene panels for gene-set (module) scores.

All genes are split into 25 bins by mean log1p(CP10K) (viewer_data/gene_means.parquet,
rank-based equal-width cut as in scanpy.tl.score_genes); from each bin 50 genes are drawn
once (seed 0). Their per-cell summed expression is stored, so a set's score only has to
read the set's own genes:
    score = mean(set genes) - mean(control genes of the bins the set occupies)
Set genes that happen to be in a panel are removed from it exactly at scoring time.
Writes viewer_data/control_panels.npz.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

APP = Path(__file__).resolve().parent.parent / 'app'
OUT = APP.parent / 'viewer_data'
sys.path.insert(0, str(APP))
from data_store import ViewerData  # noqa: E402

N_BINS, CTRL_SIZE, SEED = 25, 50, 0


def gene_bins(means: pd.Series) -> np.ndarray:
    return pd.cut(means.rank(method='min'), N_BINS, labels=False).to_numpy()


def main() -> int:
    data = ViewerData(OUT)
    means = pd.read_parquet(OUT / 'gene_means.parquet')
    assert (means.gene_id.to_numpy() == data.genes.gene_id.to_numpy()).all()
    bins = gene_bins(means.mean_log1p_cp10k)
    rng = np.random.default_rng(SEED)
    panel_rows = np.full((N_BINS, CTRL_SIZE), -1, dtype=np.int32)
    panel_sums = np.zeros((data.n_cells, N_BINS), dtype=np.float32)
    for b in range(N_BINS):
        pool = np.flatnonzero(bins == b)
        pick = np.sort(rng.choice(pool, min(CTRL_SIZE, pool.size), replace=False))
        panel_rows[b, :pick.size] = pick
        panel_sums[:, b] = data.expression_mean(list(pick)) * pick.size
        print(f'[{time.strftime("%H:%M:%S")}] bin {b}: {pool.size} genes, panel {pick.size}', flush=True)
    # spot check one bin against per-gene reads
    chk = panel_rows[N_BINS // 2]
    direct = np.sum([data.expression(int(r)) for r in chk if r >= 0], axis=0)
    assert np.allclose(direct, panel_sums[:, N_BINS // 2], atol=1e-3)
    np.savez(OUT / 'control_panels.npz', bins=bins.astype(np.int16), panel_rows=panel_rows,
             panel_sums=panel_sums, n_bins=N_BINS, ctrl_size=CTRL_SIZE, seed=SEED)
    print('written', OUT / 'control_panels.npz')
    return 0


if __name__ == '__main__':
    sys.exit(main())
