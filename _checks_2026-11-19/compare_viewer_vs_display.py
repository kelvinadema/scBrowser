"""Compare scBrowser viewer_data against the reference display sources.

Reference: notebooks/05-Mt-resource/sections/01_resource/tables/fig_atlas_cells.tsv.gz
(k200_paga_soft display coords + final 2026-09-23 calls) and shared/cells/*.parquet.
Also compares to the annotated-atlas h5ad obsm/X_umap_3d_ch5 (the viewer's declared source).
"""
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

# Project root = the folder containing this scBrowser checkout (override: SCBROWSER_PROJECT_ROOT)
ROOT = Path(os.environ.get('SCBROWSER_PROJECT_ROOT') or Path(__file__).resolve().parents[2])
CH5 = ROOT / 'notebooks/05-Mt-resource'
VIEW = ROOT / 'scBrowser/viewer_data'
H5AD = CH5 / 'upstream/results/04_figures_and_exports/annotated_atlas/Mt_chapter5_annotated_3D.h5ad'
pd.set_option('display.width', 250, 'display.max_rows', 200, 'display.max_columns', 30)


def hdr(t):
    print(f'\n{"=" * 100}\n{t}\n{"=" * 100}', flush=True)


def coord_report(name, a, b):
    d = np.abs(a - b)
    r = [np.corrcoef(a[:, i], b[:, i])[0, 1] for i in range(3)]
    # orthogonal Procrustes after centring: same shape up to rotation/reflection/translation?
    ac, bc = a - a.mean(0), b - b.mean(0)
    u, _, vt = np.linalg.svd(ac.T @ bc)
    resid = np.linalg.norm(ac @ (u @ vt) - bc) / np.linalg.norm(bc)
    print(f'{name}: max|d| {d.max():.4g}  mean|d| {d.mean():.4g}  exact(<1e-4) {(d.max(1) < 1e-4).mean():.4f}  '
          f'per-axis r {np.round(r, 4).tolist()}  procrustes rel.resid {resid:.4f}')


def knn_overlap(a, b, k=15, n=5000, seed=0):
    from sklearn.neighbors import NearestNeighbors
    rng = np.random.default_rng(seed)
    q = rng.choice(len(a), n, replace=False)
    ia = NearestNeighbors(n_neighbors=k + 1).fit(a).kneighbors(a[q], return_distance=False)[:, 1:]
    ib = NearestNeighbors(n_neighbors=k + 1).fit(b).kneighbors(b[q], return_distance=False)[:, 1:]
    return np.mean([len(set(x) & set(y)) / k for x, y in zip(ia, ib)])


hdr('Load')
v = pd.read_parquet(VIEW / 'cells.parquet')
f = pd.read_csv(CH5 / 'sections/01_resource/tables/fig_atlas_cells.tsv.gz', sep='\t', low_memory=False)
resp = pd.read_parquet(CH5 / 'shared/cells/cell_response_calls.parquet',
                       columns=['identity', 'calibration_identity', 'time_h', 'response_state', 'phase'])
nod = pd.read_parquet(CH5 / 'shared/cells/cell_nodule_identity_overlay.parquet', columns=['nodule_identity_state'])
print('viewer', v.shape, 'fig_atlas_cells', f.shape, 'resp', resp.shape)

hdr('Cell sets')
vs, fs = set(v.cell_id), set(f.cell_id)
print('viewer unique', v.cell_id.is_unique, '| shared', len(vs & fs), '| viewer-only', len(vs - fs), '| fig-only', len(fs - vs))
print('viewer-only examples', sorted(vs - fs)[:5]); print('fig-only examples', sorted(fs - vs)[:5])
# try the double-prefix normalisation used in 01_resource
ids = v.cell_id.astype(str); rest = ids.str.slice(3)
fix = ids.str.startswith('p1.') & (rest.str.split('.').str[0] == v['sample'].astype(str))
v['cid'] = ids.where(~fix, rest)
print('after p1. double-prefix normalisation: fixed', int(fix.sum()), '| shared', len(set(v.cid) & fs),
      '| viewer-only', len(set(v.cid) - fs))
m = v.merge(f, left_on='cid', right_on='cell_id', how='inner', suffixes=('_v', '_f'), validate='one_to_one')
print('joined rows', len(m), '| sample agreement', (m.sample_v.astype(str) == m.sample_f.astype(str)).mean())

hdr('Coordinates: viewer vs reference display (k200_paga_soft)')
V = m[['umap_1', 'umap_2', 'umap_3']].to_numpy(np.float64)
F = m[['umap3d_1', 'umap3d_2', 'umap3d_3']].to_numpy(np.float64)
coord_report('viewer vs fig k200_paga_soft', V, F)
print('15-NN overlap (5k cells):', round(knn_overlap(V, F), 4))

hdr('Coordinates: viewer vs annotated-atlas h5ad obsm/X_umap_3d_ch5 (declared viewer source)')
with h5py.File(H5AD, 'r') as h:
    obs = h['obs']; idx_key = obs.attrs.get('_index', '_index')
    hid = obs[idx_key][:].astype(str)
    X = h['obsm/X_umap_3d_ch5'][:].astype(np.float64)
    print('h5ad n', len(hid), 'obsm keys', list(h['obsm'].keys()))
hpos = pd.Series(np.arange(len(hid)), index=hid)
ok = v.cell_id.isin(hpos.index)
print('viewer ids found in h5ad', ok.mean())
H = X[hpos[v.cell_id[ok]].to_numpy()]
coord_report('viewer vs h5ad X_umap_3d_ch5', v.loc[ok, ['umap_1', 'umap_2', 'umap_3']].to_numpy(np.float64), H)

hdr('Annotations: viewer columns vs final calls')
r = resp.join(nod).reindex(m.cid)
m['identity'] = r.identity.to_numpy(); m['calibration_identity'] = r.calibration_identity.to_numpy()
m['nod_state'] = r.nodule_identity_state.to_numpy(); m['resp_state'] = r.response_state.to_numpy()
m['phase'] = r.phase.to_numpy()

print('\n-- ch5_tissue (viewer) vs calibration_identity (final calls)')
print(pd.crosstab(m.ch5_tissue.astype(str), m.calibration_identity.astype(str), margins=True))
print('\n-- ch5_tissue (viewer) vs identity (incl. legacy nodule_* labels)')
print(pd.crosstab(m.ch5_tissue.astype(str), m.identity.astype(str)).loc[:, lambda d: d.columns.str.startswith('nodule')])
print('\n-- ch5_tissue (viewer) vs lineage_display (display file)')
print(pd.crosstab(m.ch5_tissue.astype(str), m.lineage_display.astype(str)))

print('\n-- ch5_nodule_like (viewer) vs nodule_like_state (final NI258)')
print(pd.crosstab(m.ch5_nodule_like.astype(str), m.nod_state.astype(str), margins=True))
print('\n-- responder: viewer has a column?', [c for c in v.columns if 'respon' in c.lower()])
print(pd.Series(m.resp_state).value_counts())
print('\n-- ch5_time_h (viewer) vs time_h (fig)')
print(pd.crosstab(m.ch5_time_h.astype(str), m.time_h.astype(str)))
print('\n-- ch5_phase_calibrated (viewer) vs phase (final)')
print(pd.crosstab(m.ch5_phase_calibrated.astype(str), m.phase.astype(str)))
print('\n-- ch5_meristem_like / ch5_cortex_layer value counts')
print(v.ch5_meristem_like.value_counts(dropna=False)); print(v.ch5_cortex_layer.value_counts(dropna=False))

hdr('Reference legend numbers')
sz = f.scope.eq('susceptible_zone')
print('total', len(f), '(legend 351,239)')
print('nodule_like', int((f.nodule_like_state == 'nodule_like').sum()), '(legend 31,928)',
      '| viewer ch5_nodule_like True', int((v.ch5_nodule_like.astype(str) == 'True').sum()))
print('responding', int((f.responder_state == 'responding').sum()), '(legend 56,830)')
print('per time (treated, SZ):'); print(f[sz & f.treatment.ne('mock') & f.time_h.gt(0)].time_h.value_counts().sort_index())
print('legend: 12 34,917; 24 28,616; 48 39,997; 72 37,840; 96 14,041')
print('treatments:', f.treatment.value_counts().to_dict())
print('viewer ch5_time_h counts:', v.ch5_time_h.value_counts().to_dict())
print('\nDONE')
