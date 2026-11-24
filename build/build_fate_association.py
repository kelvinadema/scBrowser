#!/usr/bin/env python3
"""Fate-association panel data = per-gene fate-screen contexts.

Source: notebooks/05-Mt-resource/sections/09v2_b_mock_induced/data/association_mock_contexts.tsv.gz
(all 263,069 gene x 12-HPI source tissue x time rows of the 19,790-gene parent fate screen).

  effect      = 10th-to-90th weighted-fate-rank expression change (natural log in the source;
                stored here as log2, the plotted axis), pooled and per platform (BD / 10x protoplast)
  candidate   = original screen call (|effect| >= log 1.5, conditional noise ratio >= 3, valid fit)
  passes      = full criterion (`passes_1_5`): candidate AND same-direction >= 1.5-fold vs matched mock

Writes viewer_data/fate_association.parquet and viewer_data/fate_association_meta.json.
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
SRC = ROOT / 'notebooks/05-Mt-resource/sections/09v2_b_mock_induced/data/association_mock_contexts.tsv.gz'
PARENT = ROOT / 'notebooks/05-Mt-resource/sections/09v2_ot_gene_dynamics/data/fate_associations.tsv.gz'
GENE_MAP = ROOT / 'scBrowser/viewer_data/gene_umap.parquet'
OUT = ROOT / 'scBrowser/viewer_data'
COLS = ['gene_id', 'source', 'time_h', 'effect', 'effect_BD', 'effect_10xP', 'conditional_noise_score',
        'candidate', 'valid', 'passes_1_5', 'support_class', 'assessable_platforms']


def main() -> int:
    a = pd.read_csv(SRC, sep='\t', usecols=COLS, low_memory=False)
    assert len(a) == 263069 and a.gene_id.nunique() == 19790
    assert not a.duplicated(['gene_id', 'source', 'time_h']).any()
    assert int(a.passes_1_5.sum()) == 50589                      # README: 50,589 qualifying rows
    retained = set(a.loc[a.passes_1_5, 'gene_id'])
    assert len(retained) == 15531                                # fate-associated gene set
    assert retained == set(pd.read_parquet(GENE_MAP, columns=['gene_id']).gene_id)
    # effect is a natural log: every candidate clears log(1.5) but not necessarily log2(1.5)
    cand = a[a.candidate]
    assert (cand.effect.abs() >= np.log(1.5) - 1e-9).all()
    assert (cand.effect.abs() < np.log2(1.5)).any()

    # parent screen rows: detection support and fit safeguards (same 263,069 rows)
    par = pd.read_csv(PARENT, sep='\t', low_memory=False,
                      usecols=['gene_id', 'source', 'time_h', 'effect', 'valid', 'numerically_eligible',
                               'converged_BD', 'converged_10xP', 'detected_cells', 'detection'])
    a = a.merge(par, on=['gene_id', 'source', 'time_h'], how='left', validate='one_to_one',
                suffixes=('', '_parent'))
    assert a.detected_cells.notna().all()
    assert (a.valid == a.valid_parent).all()
    assert np.allclose(a.effect, a.effect_parent, equal_nan=True)
    # every fitted row clears the per-gene eligibility floor, so a missing context = below it
    assert a.detected_cells.min() >= 20 and a.detection.min() >= 0.01
    beta = a.effect.abs() / 0.8
    reason = np.select(
        [a.valid, beta > 10, (a.converged_BD == False) | (a.converged_10xP == False)],  # noqa: E712
        ['', 'unbounded fit (|fate coefficient| > 10)', 'fit did not converge'],
        default='failed or ill-conditioned fit')
    status = np.select([a.passes_1_5, a.candidate, ~a.valid],
                       ['figure7', 'fate_associated', 'excluded_fit'], default='below_threshold')

    out = pd.DataFrame({
        'gene_id': a.gene_id,
        'source': a.source,
        'time_h': a.time_h.astype('int16'),
        'log2_effect': (a.effect / np.log(2)).astype('float32'),
        'log2_effect_BD': (a.effect_BD / np.log(2)).astype('float32'),
        'log2_effect_10xP': (a.effect_10xP / np.log(2)).astype('float32'),
        'noise_ratio': a.conditional_noise_score.astype('float32'),
        'candidate': a.candidate.astype(bool),
        'passes': a.passes_1_5.astype(bool),
        'valid': a.valid.astype(bool),
        'status': status,
        'exclusion_reason': reason,
        'detected_cells': a.detected_cells.astype('int32'),
        'detection': a.detection.astype('float32'),
        'support_class': a.support_class,
        'platforms': a.assessable_platforms.fillna(''),
    })
    out.to_parquet(OUT / 'fate_association.parquet', index=False)
    (OUT / 'fate_association_meta.json').write_text(json.dumps({
        'source': str(SRC),
        'source_sha256': hashlib.sha256(SRC.read_bytes()).hexdigest(),
        'parent': str(PARENT),
        'parent_sha256': hashlib.sha256(PARENT.read_bytes()).hexdigest(),
        'eligibility': '>=100 retained cells, >=20 detecting cells, >=1% platform-balanced weighted detection; '
                       'contexts below this are absent (not fitted)',
        'excluded_fit': 'failed, nonfinite, ill-conditioned or unbounded (|fate coefficient| > 10) fits; '
                        'kept in the table, never candidates',
        'status_counts': {k: int(v) for k, v in pd.Series(status).value_counts().items()},
        'figure': 'per-gene fate-screen contexts',
        'effect': 'log2 of the 10th-to-90th weighted-fate-rank expression change (source stores ln)',
        'candidate': 'valid fit, |effect| >= ln 1.5, conditional noise ratio >= 3',
        'passes': 'candidate AND same-direction >= 1.5-fold vs matched source-tissue mock (passes_1_5)',
        'n_rows': len(out), 'n_genes': int(out.gene_id.nunique()),
        'n_passing_rows': int(out.passes.sum()), 'n_passing_genes': len(retained),
        'built': time.strftime('%Y-%m-%d %H:%M:%S'),
        'script': str(Path(__file__).resolve()),
    }, indent=2))
    print(len(out), 'rows,', out.gene_id.nunique(), 'genes,', int(out.passes.sum()), 'passing rows,',
          len(retained), 'genes with a fate-screen context')
    return 0


if __name__ == '__main__':
    sys.exit(main())
