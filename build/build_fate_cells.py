#!/usr/bin/env python3
"""Per-cell nodule-fate data for the browser (Fate colouring + Fate tab decile curves).

Contexts = the fate screen (09v2; sideprojects/CH5screenV2 run_2026-10-28):
  per 12-HPI source tissue x time x platform; 12 HPI = ALL source-cohort cells (uniform weight,
  BD only); 24/48/72 HPI = forward-mass-weighted, not-yet-nodule-like descendants (the export's
  `w_from_<source>`, cell sets identical to the screen's context QC table).
  Fate decile = weighted mid-percentile rank of p_nodule96 within context (the screen's `wrank`),
  decile = clip(floor((rank + 0.5) * 10), 0, 9) + 1.
Expression scale for the curves = the screen's: CPM over common-gene counts (`screen_lib_size`).
The viewer stores log1p(CP10K over the full-depth total = X row sum, `corrected_counts_total_counts`), so
CPM_screen = expm1(v) * 100 * total_X / screen_lib_size  ->  stored here as `cpm_factor`.

Source: notebooks/05-Mt-resource/sections/_exports/viewer_ot_2026-10-09/cells_ot.parquet
0-HPI colouring probabilities are supplemented from the original arches_e05 endpoint export.
This does not add 0-HPI contexts to the expression curves or change their weights.
Writes viewer_data/fate_cells.parquet (row order = cells.parquet = expr.h5).
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

# Project root = the folder containing this scBrowser checkout (override: SCBROWSER_PROJECT_ROOT)
ROOT = Path(os.environ.get('SCBROWSER_PROJECT_ROOT') or Path(__file__).resolve().parents[2])
CH5 = ROOT / 'notebooks/05-Mt-resource'
EXPORT = CH5 / 'sections/_exports/viewer_ot_2026-10-09'
ENDPOINTS = CH5 / 'upstream/results/02_analyses/temporal_transport/meristem_signature_revision_2026-09-23/models/arches_e05/percell_endpoint_probabilities.parquet'
H5AD = CH5 / 'upstream/results/04_figures_and_exports/annotated_atlas/Mt_chapter5_annotated_3D.h5ad'
QC = ROOT / 'sideprojects/CH5screenV2/data/analyses/run_2026-10-28/tables/context_platform_QC.tsv'
APP = ROOT / 'scBrowser/app'
OUT = ROOT / 'scBrowser/viewer_data'
SOURCES = ['root_hair_cell', 'cortex', 'pericycle_cell', 'procambium']
sys.path.insert(0, str(CH5 / 'sections/_analyses/fate_window_2026-10-01'))
sys.path.insert(0, str(APP))
from lineage2_export import wrank  # noqa: E402  (the screen's own function)


def decile(rank: np.ndarray) -> np.ndarray:
    return np.clip(np.floor((rank + 0.5) * 10), 0, 9).astype(int) + 1


def main() -> int:
    cells = pd.read_parquet(OUT / 'cells.parquet', columns=['cell_id'])
    ot = pd.read_parquet(EXPORT / 'cells_ot.parquet').set_index('cell_id').reindex(cells.cell_id)
    assert ot.index.notna().all() and len(ot) == len(cells)
    with h5py.File(H5AD, 'r') as h:
        obs = h['obs']
        ids = obs[obs.attrs['_index']][:].astype(str)
        total = pd.Series(obs['corrected_counts_total_counts'][:], index=ids).reindex(cells.cell_id)
    assert total.notna().all()

    out = pd.DataFrame({'cell_id': cells.cell_id.to_numpy()})
    out['in_ot'] = ot.in_ot.fillna(False).astype(bool).to_numpy()
    out['ot_time_h'] = ot.time_h.astype('float').to_numpy().astype(np.float32)
    out['ot_platform'] = ot.platform.astype(object).where(ot.platform.notna(), '').to_numpy()
    out['is_ot_endpoint'] = ot.is_ot_endpoint.fillna(False).astype(bool).to_numpy()
    out['p_nodule96'] = ot.p_nodule96.to_numpy(np.float32)
    endpoints = pd.read_parquet(ENDPOINTS, columns=['cell', 'time_h', 'new_nodule'])
    assert endpoints.cell.is_unique
    endpoints = endpoints.set_index('cell')
    aligned = endpoints.reindex(out.cell_id)
    existing = out.p_nodule96.notna().to_numpy()
    assert np.allclose(out.p_nodule96.to_numpy()[existing],
                       aligned.new_nodule.to_numpy()[existing], rtol=0, atol=1e-6)
    mock = (out.in_ot & out.ot_time_h.eq(0)).to_numpy()
    assert int(mock.sum()) == int(endpoints.time_h.eq(0).sum())
    assert (aligned.time_h.to_numpy()[mock] == 0).all()
    mock_p = aligned.new_nodule.to_numpy()[mock]
    assert np.isfinite(mock_p).all() and ((mock_p >= 0) & (mock_p <= 1)).all()
    out.loc[mock, 'p_nodule96'] = mock_p.astype(np.float32)
    lib = ot.screen_lib_size.to_numpy()
    out['cpm_factor'] = np.where(lib > 0, 100.0 * total.to_numpy() / lib, np.nan).astype(np.float32)

    qc = pd.read_csv(QC, sep='\t')
    checks = {}
    for s in SOURCES:
        w = np.full(len(out), np.nan)
        dec = np.zeros(len(out), dtype=np.int8)
        fm, wf = ot[f'forward_mass_{s}'].to_numpy(), ot[f'w_from_{s}'].to_numpy()
        tt, p = ot.time_h.astype('float').to_numpy(), ot.p_nodule96.to_numpy()
        pl = ot.platform.astype(object).where(ot.platform.notna(), '').to_numpy()
        for t in (12, 24, 48, 72):
            for plat in ('BD', '10xP'):
                if t == 12:
                    m = (tt == 12) & (pl == plat) & np.isfinite(fm) & (fm > 0)
                    wt = np.ones(m.sum())
                else:
                    m = (tt == t) & (pl == plat) & np.isfinite(wf)
                    wt = wf[m]
                if not m.any():
                    continue
                want = qc[(qc.source == s) & (qc.time_h == t) & (qc.platform == plat)].retained_cells
                assert len(want) == 1 and int(want.iloc[0]) == int(m.sum()), (s, t, plat, int(m.sum()))
                assert np.isfinite(p[m]).all()
                dd = decile(wrank(p[m], wt))
                w[m], dec[m] = wt, dd
                if t > 12:   # must reproduce the export's stored deciles exactly
                    stored = ot[f'fate_decile_{s}'].to_numpy()[m].astype(int)
                    assert (stored == dd).all(), (s, t, plat)
                checks[f'{s}_{t}_{plat}'] = int(m.sum())
        out[f'w_{s}'] = w.astype(np.float32)
        out[f'decile_{s}'] = dec

    # known-answer: reproduce stored weighted-CPM decile profiles (24-72 HPI) from viewer data
    from data_store import ViewerData
    data = ViewerData(OUT)
    ka = pd.read_csv(EXPORT / 'known_answer_profiles.tsv', sep='\t')
    ka = ka[ka.time_h > 12]
    worst = 0.0
    for gid, g in ka.groupby('gene_id'):
        row = data.resolve_gene(gid)
        if row is None:
            continue
        cpm = np.expm1(data.expression(row).astype(np.float64)) * out.cpm_factor.to_numpy()
        for r in g.itertuples():
            m = ((out.ot_time_h == r.time_h) & (out.ot_platform == r.platform)
                 & (out[f'decile_{r.source}'] == r.fate_decile)).to_numpy()
            ww = out[f'w_{r.source}'].to_numpy()[m]
            val = float(np.sum(ww * cpm[m]) / np.sum(ww))
            ref = r.weighted_cpm_stored
            worst = max(worst, abs(val - ref) / max(ref, 1.0))
    print('known-answer max relative CPM error (float16 storage):', round(worst, 5))
    assert worst < 0.01, worst

    out.to_parquet(OUT / 'fate_cells.parquet', index=False)
    (OUT / 'fate_cells_meta.json').write_text(json.dumps({
        'source': str(EXPORT / 'cells_ot.parquet'), 'qc': str(QC), 'contexts': checks,
        'zero_hpi_probability_source': str(ENDPOINTS),
        'zero_hpi_probability_column': 'new_nodule',
        'zero_hpi_cells': int(mock.sum()),
        'known_answer_max_rel_error': worst,
        'p_nodule96': 'OT probability of reaching the 96-HPI nodule-like endpoint (NI258), arches_e05; modelled 0-96 HPI time-course cells; 96-HPI values are the endpoint call (0 or 1)',
        'built': time.strftime('%Y-%m-%d %H:%M:%S'), 'script': str(Path(__file__).resolve()),
    }, indent=2))
    print('written', OUT / 'fate_cells.parquet', '| p_nodule96 non-NA', int(out.p_nodule96.notna().sum()))
    return 0


if __name__ == '__main__':
    sys.exit(main())
