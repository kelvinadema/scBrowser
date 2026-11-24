#!/usr/bin/env python3
"""Precompute gene-set scores for the browser's precomputed sets: the ten fate-associated gene modules
and the top 200 genes of the 20 cNMF programmes (GEP01-GEP20).

GEP sets: notebooks/05-Mt-resource/sections/03_gene_programmes/data/program_top_genes.tsv,
rank <= 200 (the reference "top 200"; NOT the original_GEP_ranks table, whose top 200 is incomplete).

Uses the app's own scorer (app/gene_sets.py, control_panels.npz), so a preset gives exactly
the score a pasted copy of the same list would. Large modules (up to 3,482 genes) take up
to a minute to score from a cold read, hence the cache.
Writes viewer_data/module_scores.npz: for each module, its gene rows and per-cell score.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import json

import numpy as np
import pandas as pd

APP = Path(__file__).resolve().parent.parent / 'app'
OUT = APP.parent / 'viewer_data'
# Project root = the folder containing this scBrowser checkout (override: SCBROWSER_PROJECT_ROOT)
ROOT = Path(os.environ.get('SCBROWSER_PROJECT_ROOT') or Path(__file__).resolve().parents[2])
GEP_TABLE = ROOT / 'notebooks/05-Mt-resource/sections/03_gene_programmes/data/program_top_genes.tsv'
GEP_TOP = 200
sys.path.insert(0, str(APP))
from data_store import ViewerData  # noqa: E402
from gene_map import GeneMap  # noqa: E402
from gene_sets import GeneSetScorer  # noqa: E402


def main() -> int:
    data = ViewerData(OUT)
    scorer = GeneSetScorer(data, OUT)
    gmap = GeneMap(OUT)
    sets = {mid: (gmap.meta['modules'][mid]['name'], list(gmap.genes[gmap.genes.module_id == mid].gene_id))
            for mid in gmap.order}
    gep = pd.read_csv(GEP_TABLE, sep='\t')
    gep = gep[gep['rank'] <= GEP_TOP]
    for prog, g in gep.groupby('program', sort=True):
        assert len(g) == GEP_TOP and g.gene_id.is_unique, prog
        sets[prog] = (f'{prog} (top {GEP_TOP})', list(g.sort_values('rank').gene_id))
    assert sum(k.startswith('GEP') for k in sets) == 20
    arrays, labels = {}, {}
    for mid, (label, ids) in sets.items():
        labels[mid] = {'label': label, 'n_genes': len(ids),
                       'group': 'cnmf' if mid.startswith('GEP') else 'figure7'}
        rows, missing = scorer.resolve(' '.join(ids))
        assert not missing and len(rows) == len(ids), mid
        values, n_ctrl = scorer.score(rows)
        arrays[f'{mid}__rows'] = np.array(sorted(rows), dtype=np.int32)
        arrays[f'{mid}__score'] = values
        arrays[f'{mid}__n_ctrl'] = np.array(n_ctrl)
        print(mid, len(rows), 'genes,', n_ctrl, 'controls, mean score', float(values.mean()), flush=True)
    np.savez(OUT / 'module_scores.npz', **arrays)
    (OUT / 'module_scores_meta.json').write_text(json.dumps(
        {'sets': labels, 'gep_table': str(GEP_TABLE), 'gep_top': GEP_TOP}, indent=2))
    print('written', OUT / 'module_scores.npz')
    return 0


if __name__ == '__main__':
    sys.exit(main())
