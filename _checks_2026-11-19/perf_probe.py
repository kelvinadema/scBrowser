import os, sys, time, json, resource
t0=time.time()
sys.path.insert(0,'app')
import app as A
t_start=time.time()-t0
import psutil
rss=lambda: psutil.Process().memory_info().rss/1e9
print(f'startup {t_start:.1f}s, RSS after startup {rss():.2f} GB')
import plotly.io as pio, numpy as np
from plotting import PlotSettings, build_figure
for n in (50000,150000,351239):
    for mode,kw in (('annotation',{}),('gene',{'gene_row':A.DATA.resolve_gene('MtrunA17_Chr5g0448621')})):
        t=time.time(); f=build_figure(A.DATA, PlotSettings(mode=mode, annotation='ch5_tissue', max_points=n, **kw), A.AXES); tb=time.time()-t
        t=time.time(); js=pio.to_json(f); tj=time.time()-t
        print(f'{mode:10s} {n:>7,} pts: build {tb:.2f}s, json {tj:.2f}s, payload {len(js)/1e6:.1f} MB')
# memory of big in-memory objects
def sz(o):
    import pandas as pd
    if isinstance(o, pd.DataFrame): return o.memory_usage(deep=True).sum()/1e9
    if isinstance(o, np.ndarray): return o.nbytes/1e9
    return float('nan')
print('cells df', round(sz(A.DATA.cells),3),'GB; panel_sums', round(sz(A.SCORER.panel_sums),3), 'GB')
before=rss(); 
import gc; gc.collect()
print('FATE by_gene dict entries', len(A.FATE.by_gene))
# exercise: heatmap GEP09 + fate profile
rows=sorted(A.SCORER.resolve(' '.join(A.MODULE_PRESETS['GEP09'][1]))[0])
t=time.time(); A.HEATMAP.build(rows,'ch5_tissue',True,None,'z',6); print(f'heatmap GEP09 split: {time.time()-t:.1f}s')
print(f'RSS after work {rss():.2f} GB; peak {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.2f} GB')
