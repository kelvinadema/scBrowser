#!/usr/bin/env python3
"""Gene UMAP tab data = fate-associated gene map.

Source: notebooks/05-Mt-resource/sections/09v2_b_mock_induced/tables/fig_gene_umap.tsv
(15,531 fate-associated genes changed vs matched mock; refitted 2D UMAP; modules M01-M10).
Module names and colours follow the final figure legend, not the
section table's `module_color` column (a tab10 draft palette).

Writes viewer_data/gene_umap.parquet and viewer_data/gene_umap_meta.json.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

import pandas as pd

# Project root = the folder containing this scBrowser checkout (override: SCBROWSER_PROJECT_ROOT)
ROOT = Path(os.environ.get('SCBROWSER_PROJECT_ROOT') or Path(__file__).resolve().parents[2])
SRC = ROOT / 'notebooks/05-Mt-resource/sections/09v2_b_mock_induced/tables/fig_gene_umap.tsv'
SRC3D = SRC.with_name('fig_gene_umap_3d.tsv')
OUT = ROOT / 'scBrowser/viewer_data'

# Module legend, top to bottom: name, colour, gene count
MODULES = {
    'M07': ('M07 cortex-high', '#008300', 840),
    'M05': ('M05 pericycle-high', '#4a3aa7', 1717),
    'M06': ('M06 procambium-high', '#eda100', 1094),
    'M08': ('M08 procambium-high', '#b07a00', 580),
    'M10': ('M10 procambium-high', '#f3c96a', 117),
    'M09': ('M09 meristem-high', '#e34948', 332),
    'M04': ('M04 mixed', '#52514e', 2353),
    'M03': ('M03 mixed', '#8a8984', 2298),
    'M02': ('M02 mixed', '#adaca6', 2718),
    'M01': ('M01 mixed', '#c4c3be', 3482),
}
MOCK_LABELS = {'induced': 'induced vs mock', 'repressed': 'repressed vs mock',
               'both_across_contexts': 'induced and repressed (different contexts)'}


def main() -> int:
    g = pd.read_csv(SRC, sep='\t')
    genes = pd.read_parquet(OUT / 'genes.parquet', columns=['gene_id', 'symbol'])
    assert len(g) == 15531 and g.gene_id.is_unique
    assert g.gene_id.isin(genes.gene_id).all(), 'gene missing from expr.h5'
    counts = g.module_id.value_counts()
    for mid, (_, _, n) in MODULES.items():
        assert counts[mid] == n, (mid, counts[mid], n)
    assert int(g.curated.sum()) == 96

    sym = g.symbol.where(g.symbol.notna() & (g.symbol != g.gene_id))
    viewer_sym = g.gene_id.map(genes.set_index('gene_id').symbol).replace('', pd.NA)
    out = pd.DataFrame({
        'gene_id': g.gene_id,
        'label': sym.fillna(viewer_sym).fillna(g.gene_id),
        'module_id': g.module_id,
        'module': g.module_id.map(lambda m: MODULES[m][0]),
        'color': g.module_id.map(lambda m: MODULES[m][1]),
        'mock_change': g.mock_change_class.map(MOCK_LABELS),
        'two_platform': g.two_platform.astype(bool),
        'curated': g.curated.astype(bool),
        'umap_1': g.umap_1.astype('float32'),
        'umap_2': g.umap_2.astype('float32'),
    })
    assert out.mock_change.notna().all()

    # optional 3D layout (same 64D input, n_components=3; README of 09v2_b, data/embedding_3d.json)
    g3 = pd.read_csv(SRC3D, sep='\t').set_index('gene_id').reindex(g.gene_id)
    assert len(g3) == len(g) and g3[['umap_1', 'umap_2', 'umap_3', 'module_id']].notna().all().all()
    for col in ('module_id', 'mock_change_class', 'curated'):
        assert (g3[col].to_numpy() == g[col].to_numpy()).all(), f'3D table disagrees on {col}'
    for i in (1, 2, 3):
        out[f'umap3d_{i}'] = g3[f'umap_{i}'].to_numpy('float32')
    out.to_parquet(OUT / 'gene_umap.parquet', index=False)
    (OUT / 'gene_umap_meta.json').write_text(json.dumps({
        'source': str(SRC),
        'source_sha256': hashlib.sha256(SRC.read_bytes()).hexdigest(),
        'source_3d': str(SRC3D),
        'source_3d_sha256': hashlib.sha256(SRC3D.read_bytes()).hexdigest(),
        'figure': 'fate-associated gene map',
        'module_order': list(MODULES),
        'modules': {k: {'name': v[0], 'color': v[1], 'n_genes': v[2]} for k, v in MODULES.items()},
        'n_genes': len(out), 'n_curated': int(out.curated.sum()),
        'built': time.strftime('%Y-%m-%d %H:%M:%S'),
        'script': str(Path(__file__).resolve()),
    }, indent=2))
    print(out.module.value_counts().to_dict(), 'curated', int(out.curated.sum()))
    return 0


if __name__ == '__main__':
    sys.exit(main())
