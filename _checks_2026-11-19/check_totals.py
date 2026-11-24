import os
from pathlib import Path

import h5py, numpy as np, pandas as pd

# Project root = the folder containing this scBrowser checkout (override: SCBROWSER_PROJECT_ROOT)
ROOT = Path(os.environ.get('SCBROWSER_PROJECT_ROOT') or Path(__file__).resolve().parents[2])
h=h5py.File(ROOT / 'notebooks/05-Mt-resource/upstream/results/04_figures_and_exports/annotated_atlas/Mt_chapter5_annotated_3D.h5ad','r')
o=h['obs']; ids=o[o.attrs['_index']][:].astype(str)
cols={k:o[k][:] for k in ['foundation__total_counts_X','corrected_counts_total_counts','legacy_endpoint__total_counts_X']}
ot=pd.read_parquet(ROOT / 'notebooks/05-Mt-resource/sections/_exports/viewer_ot_2026-10-09/cells_ot.parquet').set_index('cell_id').reindex(ids)
screen=ot.screen_lib_size.notna().to_numpy()
for k,v in cols.items(): print(k,'NaN',int(np.isnan(v).sum()),'| NaN among screen cells', int(np.isnan(v[screen]).sum()))
ip=h['X/indptr'][:]; rng=np.random.default_rng(0); idx=np.sort(rng.choice(len(ids),2000,replace=False))
rs=np.array([h['X/data'][ip[i]:ip[i+1]].sum() for i in idx])
for k,v in cols.items():
    ok=~np.isnan(v[idx]); print(k,'== X row sum on', int(ok.sum()),'sampled:', bool(np.allclose(v[idx][ok], rs[ok], rtol=1e-5)))
