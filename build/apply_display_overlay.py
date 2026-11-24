#!/usr/bin/env python3
"""Bring viewer_data/cells.parquet in line with the reference display and final calls.

Run after build_viewer_data.py (which takes coords/labels from the annotated-atlas h5ad,
i.e. the superseded `k200_soft` embedding and pre-2026-09-23 calls). This overlay replaces:

  umap_1..3        <- connect3d k200_paga_soft (reference display), fig_atlas_cells.tsv.gz umap3d_1..3
  ch5_tissue       <- calibration_identity (final identity; legacy nodule_cortex/pericycle returned)
  ch5_nodule_like  <- final NI258 nodule_identity_state (31,928 nodule-like)
  ch5_responder    <- final responder_state (56,830 responding)   [new column]

Row order of cells.parquet is preserved: expr.h5 indexes cells by row position.
expr.h5 / genes.parquet are untouched. Old cells.parquet + meta.json are moved to
viewer_data/_superseded_<date>/.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from build_viewer_data import harmonise_po, load_po_labels, relabel_po, tidy_categorical  # noqa: E402

# Project root = the folder containing this scBrowser checkout (override: SCBROWSER_PROJECT_ROOT)
ROOT = Path(os.environ.get('SCBROWSER_PROJECT_ROOT') or Path(__file__).resolve().parents[2])
CH5 = ROOT / 'notebooks/05-Mt-resource'
FIG = CH5 / 'sections/01_resource/tables/fig_atlas_cells.tsv.gz'
CALLS = CH5 / 'shared/cells/cell_response_calls.parquet'
NOD = CH5 / 'shared/cells/cell_nodule_identity_overlay.parquet'
OUT = ROOT / 'scBrowser/viewer_data'

EXPECT = {'n_cells': 351239, 'nodule_like': 31928, 'responding': 56830}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, 'rb') as fh:
        for b in iter(lambda: fh.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


def main() -> int:
    cells = pd.read_parquet(OUT / 'cells.parquet')
    meta = json.loads((OUT / 'meta.json').read_text())
    n = len(cells)
    assert n == EXPECT['n_cells'] == meta['n_cells'] and cells.cell_id.is_unique

    fig = pd.read_csv(FIG, sep='\t', low_memory=False).set_index('cell_id')
    calls = pd.read_parquet(CALLS, columns=['calibration_identity', 'response_state'])
    nod = pd.read_parquet(NOD, columns=['nodule_identity_state'])
    assert set(fig.index) == set(cells.cell_id) == set(calls.index) == set(nod.index), 'cell sets differ'

    f = fig.reindex(cells.cell_id)
    c = calls.reindex(cells.cell_id)
    d = nod.reindex(cells.cell_id)
    assert (f['sample'].astype(str).to_numpy() == cells['sample'].astype(str).to_numpy()).all()
    # the two final tables must agree with the figure table
    assert (f.nodule_like_state.to_numpy() == d.nodule_identity_state.to_numpy()).all()
    assert (f.responder_state.to_numpy() ==
            c.response_state.map({1: 'responding', 0: 'below_threshold', -1: 'unassessable'}).to_numpy()).all()
    assert (f.identity.to_numpy() != c.calibration_identity.to_numpy()).sum() == 1047  # the legacy nodule_* cells

    new = cells.copy()
    for i in (1, 2, 3):
        new[f'umap_{i}'] = f[f'umap3d_{i}'].to_numpy(np.float32)
    assert np.isfinite(new[['umap_1', 'umap_2', 'umap_3']].to_numpy()).all()

    po = load_po_labels(Path(__file__).with_name('po_terms.tsv'))
    tissue = harmonise_po(pd.Series(c.calibration_identity.astype(str).to_numpy()), po)
    new['ch5_tissue'] = relabel_po(tidy_categorical(tissue), po)

    new['ch5_nodule_like'] = tidy_categorical((f.nodule_like_state == 'nodule_like').map({True: 'True', False: 'False'}).to_numpy())
    resp = pd.Categorical(f.responder_state.to_numpy(), categories=['responding', 'below_threshold', 'unassessable'])
    new['ch5_responder'] = resp

    # sanity vs reference legend counts
    assert int((new.ch5_nodule_like == 'True').sum()) == EXPECT['nodule_like']
    assert int((new.ch5_responder == 'responding').sum()) == EXPECT['responding']
    assert not new.ch5_tissue.astype(str).str.startswith('nodule').any()
    assert (new.cell_id.to_numpy() == cells.cell_id.to_numpy()).all()

    stamp = time.strftime('%Y-%m-%d')
    keep = OUT / f'_superseded_{stamp}'
    keep.mkdir(exist_ok=True)
    for name in ('cells.parquet', 'meta.json'):
        if not (keep / name).exists():
            shutil.copy2(OUT / name, keep / name)

    cat = list(meta['categorical_columns'])
    if 'ch5_responder' not in cat:
        cat.insert(cat.index('ch5_nodule_like') + 1, 'ch5_responder')
    meta['categorical_columns'] = cat
    meta['column_labels'].update({
        'ch5_tissue': 'Cell type (Plant Ontology)',
        'ch5_nodule_like': 'Nodule-like state',
        'ch5_responder': 'Responding cells',
    })
    meta['umap_key'] = 'connect3d__k200_paga_soft (reference display)'
    meta['display_overlay'] = {
        'applied': time.strftime('%Y-%m-%d %H:%M:%S'),
        'script': str(Path(__file__).resolve()),
        'coords': f'{FIG} umap3d_1..3',
        'ch5_tissue': f'{CALLS} calibration_identity',
        'ch5_nodule_like': f'{NOD} nodule_identity_state == nodule_like',
        'ch5_responder': f'{FIG} responder_state',
        'sources_sha256': {str(p): sha256(p) for p in (FIG, CALLS, NOD)},
        'superseded_copy': str(keep),
        'counts': {'nodule_like': EXPECT['nodule_like'], 'responding': EXPECT['responding']},
    }

    tmp = OUT / 'cells.parquet.tmp'
    new.to_parquet(tmp, index=False)
    tmp.replace(OUT / 'cells.parquet')
    (OUT / 'meta.json').write_text(json.dumps(meta, indent=2))

    print('ch5_tissue:', new.ch5_tissue.value_counts().to_dict())
    print('ch5_nodule_like:', new.ch5_nodule_like.value_counts().to_dict())
    print('ch5_responder:', new.ch5_responder.value_counts().to_dict())
    print('written', OUT / 'cells.parquet', '| backup in', keep)
    return 0


if __name__ == '__main__':
    sys.exit(main())
