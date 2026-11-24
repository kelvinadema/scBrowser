#!/usr/bin/env python3
"""Per-cell design (condition, platform, time) and the sample list, for the browser's
Expression tab (dot plot) and the sample filter.

Source: notebooks/05-Mt-resource/sections/01_resource/tables/fig_atlas_cells.tsv.gz
        + fig_library_inventory.tsv (library inventory).
BD tip names are swapped at source (author-confirmed): BDTip10mm = 0-2 mm, BDTip2mm = 0-10 mm.
Writes viewer_data/cell_design.parquet (row order = cells.parquet) and viewer_data/samples.json.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pandas as pd

# Project root = the folder containing this scBrowser checkout (override: SCBROWSER_PROJECT_ROOT)
ROOT = Path(os.environ.get('SCBROWSER_PROJECT_ROOT') or Path(__file__).resolve().parents[2])
S01 = ROOT / 'notebooks/05-Mt-resource/sections/01_resource/tables'
OUT = ROOT / 'scBrowser/viewer_data'

CONDITION = {  # (treatment, time_h) -> display label, in display order
    ('untreated_tip', 0): 'Root tip', ('mock', 0): 'Mock', ('NF', 3): 'NF 3 h', ('BAP', 3): 'BAP 3 h',
    ('S_meliloti', 12): 'Rhizobia 12 h', ('S_meliloti', 24): 'Rhizobia 24 h',
    ('S_meliloti', 48): 'Rhizobia 48 h', ('S_meliloti', 72): 'Rhizobia 72 h',
    ('S_meliloti', 96): 'Rhizobia 96 h',
}
PLATFORM = {('BD', 'protoplast'): 'BD', ('10x', 'protoplast'): '10x protoplast',
            ('10xN', 'nuclei'): '10x nuclei'}
TIP_REGION = {'BDTip10mm': '0–2 mm', 'BDTip2mm': '0–10 mm'}


def main() -> int:
    cells = pd.read_parquet(OUT / 'cells.parquet', columns=['cell_id', 'sample'])
    f = pd.read_csv(S01 / 'fig_atlas_cells.tsv.gz', sep='\t', low_memory=False,
                    usecols=['cell_id', 'sample', 'scope', 'platform', 'modality', 'treatment', 'time_h'])
    f = f.set_index('cell_id').reindex(cells.cell_id)
    assert f['sample'].notna().all()
    assert (f['sample'].astype(str).to_numpy() == cells['sample'].astype(str).to_numpy()).all()
    cond = pd.Series([CONDITION[(t, int(h))] for t, h in zip(f.treatment, f.time_h)])
    plat = pd.Series([PLATFORM[(p, m)] for p, m in zip(f.platform, f.modality)])
    out = pd.DataFrame({
        'cell_id': cells.cell_id.to_numpy(),
        'condition': pd.Categorical(cond, categories=list(CONDITION.values())),
        'platform': pd.Categorical(plat, categories=list(PLATFORM.values())),
        'time_h': f.time_h.astype('int16').to_numpy(),
    })
    out.to_parquet(OUT / 'cell_design.parquet', index=False)

    inv = pd.read_csv(S01 / 'fig_library_inventory.tsv', sep='\t').set_index('sample')
    samples = []
    per = pd.DataFrame({'sample': cells['sample'].astype(str), 'condition': cond, 'platform': plat})
    for s, g in per.groupby('sample', sort=False):
        c, p = g.condition.iloc[0], g.platform.iloc[0]
        assert g.condition.nunique() == 1 and g.platform.nunique() == 1, s
        region = f" {TIP_REGION[s]}" if s in TIP_REGION else ''
        n = int(len(g))
        assert n == int(inv.loc[s, 'n_cells']), s
        samples.append({'sample': s, 'condition': c + region, 'platform': p, 'n_cells': n,
                        'label': f'{s} · {c}{region} · {p} ({n:,})'})
    order = {c: i for i, c in enumerate(CONDITION.values())}
    samples.sort(key=lambda d: (order[d['condition'].split(' 0')[0].strip()], d['platform'], d['sample']))
    (OUT / 'samples.json').write_text(json.dumps(samples, indent=2, ensure_ascii=False))
    print(out.condition.value_counts().to_dict()); print(len(samples), 'samples')
    return 0


if __name__ == '__main__':
    sys.exit(main())
