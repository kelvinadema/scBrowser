#!/usr/bin/env python3
"""Dissociation-response panel data = bulk RNA-seq isolation effect (Medicago).

Source: notebooks/03-sc-vs-sn/output/bulk/isolation-effects-deg-nuclear-only/
  deseq2_protoplast_vs_tissue.tsv, deseq2_nuclei_vs_tissue.tsv  (DESeq2 ~ method,
  T1-3 intact 1.5 cm root section vs P1-3 protoplasts / N6-8 nuclei; organellar genes dropped,
  genes with >=10 counts in >=3 samples tested; apeglm-shrunk log2FC)
  size_factors.tsv + input/bulk/medicago/gene_counts.tsv -> normalised counts per sample.

Call (as in the source analysis): DESeq2 lfcThreshold=1 (altHypothesis greaterAbs), padj < 0.01
(`is_DEG_lfcThreshold`).

Writes viewer_data/dissociation.parquet and viewer_data/dissociation_meta.json.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# Project root = the folder containing this scBrowser checkout (override: SCBROWSER_PROJECT_ROOT)
ROOT = Path(os.environ.get('SCBROWSER_PROJECT_ROOT') or Path(__file__).resolve().parents[2])
CH3 = ROOT / 'notebooks/03-sc-vs-sn'
DEG = CH3 / 'output/bulk/isolation-effects-deg-nuclear-only'
COUNTS = CH3 / 'input/bulk/medicago/gene_counts.tsv'
OUT = ROOT / 'scBrowser/viewer_data'
CONTRASTS = {'protoplast': 'deseq2_protoplast_vs_tissue.tsv', 'nuclei': 'deseq2_nuclei_vs_tissue.tsv'}
EXPECT_CALLS = {'protoplast': 3562, 'nuclei': 1184}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> int:
    out = None
    for name, fn in CONTRASTS.items():
        d = pd.read_csv(DEG / fn, sep='\t')
        assert d.gene.is_unique
        assert int(d.is_DEG_lfcThreshold.sum()) == EXPECT_CALLS[name], name
        d = d.rename(columns={
            'log2FoldChange': f'{name}_lfc', 'lfcSE': f'{name}_se',
            'padj_lfcThreshold': f'{name}_padj', 'is_DEG_lfcThreshold': f'{name}_call',
        })[['gene', 'baseMean', f'{name}_lfc', f'{name}_se', f'{name}_padj', f'{name}_call']]
        out = d if out is None else out.merge(d.drop(columns='baseMean'), on='gene', validate='one_to_one')
    out = out.rename(columns={'gene': 'gene_id'})
    for name in CONTRASTS:
        out[f'{name}_call'] = out[f'{name}_call'].astype(bool)

    sf = pd.read_csv(DEG / 'size_factors.tsv', sep='\t')
    counts = pd.read_csv(COUNTS, sep='\t', index_col=0)[sf['sample']]
    norm = counts.div(sf.set_index('sample').sizeFactor, axis=1).reindex(out.gene_id)
    assert norm.notna().all().all()
    # sanity: DESeq2 baseMean = mean of normalised (rounded) counts
    base = counts.round().div(sf.set_index('sample').sizeFactor, axis=1).reindex(out.gene_id).mean(axis=1)
    assert np.allclose(base.to_numpy(), out.baseMean.to_numpy(), rtol=1e-3), 'size factors do not reproduce baseMean'
    for s in sf['sample']:
        out[f'norm_{s}'] = norm[s].to_numpy(np.float32)

    viewer = pd.read_parquet(OUT / 'genes.parquet', columns=['gene_id'])
    out.to_parquet(OUT / 'dissociation.parquet', index=False)
    (OUT / 'dissociation_meta.json').write_text(json.dumps({
        'source': 'bulk RNA-seq isolation effect (Medicago, nuclear genes only)',
        'design': 'DESeq2 ~ method; tissue T1-T3 (intact 1.5 cm root section) vs protoplast P1-P3 / nuclei N6-N8',
        'call': 'lfcThreshold=1 greaterAbs, padj < 0.01 (is_DEG_lfcThreshold)',
        'lfc': 'apeglm-shrunk log2FoldChange vs intact tissue; error bar = lfcSE',
        'n_tested': len(out), 'n_calls': {k: int(out[f'{k}_call'].sum()) for k in CONTRASTS},
        'n_viewer_genes_tested': int(viewer.gene_id.isin(out.gene_id).sum()),
        'samples': {r.sample: r.method for r in sf.itertuples()},
        'sources_sha256': {str(p): sha(p) for p in [*(DEG / f for f in CONTRASTS.values()),
                                                    DEG / 'size_factors.tsv', COUNTS]},
        'built': time.strftime('%Y-%m-%d %H:%M:%S'),
        'script': str(Path(__file__).resolve()),
    }, indent=2))
    print(len(out), 'genes;', {k: int(out[f'{k}_call'].sum()) for k in CONTRASTS},
          '| viewer genes with bulk result', int(viewer.gene_id.isin(out.gene_id).sum()), '/', len(viewer))
    return 0


if __name__ == '__main__':
    sys.exit(main())
