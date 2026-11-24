import sys, time
from pathlib import Path
sys.path.insert(0, 'app')
T = time.perf_counter
t = T(); import numpy, pandas, plotly.io as pio, shiny; print(f'imports              {T()-t:5.2f}s')
D = Path('viewer_data')
from data_store import ViewerData
t = T(); DATA = ViewerData(D); print(f'ViewerData           {T()-t:5.2f}s')
from plotting import scene_axes
t = T(); scene_axes(DATA.coords); print(f'scene_axes           {T()-t:5.2f}s')
from gene_map import GeneMap
t = T(); GM = GeneMap(D); print(f'GeneMap              {T()-t:5.2f}s')
from dissociation import Dissociation
t = T(); Dissociation(D); print(f'Dissociation         {T()-t:5.2f}s')
from fate_assoc import FateAssociation
t = T(); FateAssociation(D); print(f'FateAssociation      {T()-t:5.2f}s')
from fate_view import FateCells
t = T(); FateCells(D); print(f'FateCells            {T()-t:5.2f}s')
from expression_view import DotPlot
t = T(); DP = DotPlot(DATA, D); print(f'DotPlot              {T()-t:5.2f}s')
from gene_sets import GeneSetScorer
t = T(); GeneSetScorer(DATA, D); print(f'GeneSetScorer        {T()-t:5.2f}s')
import json, numpy as np
t = T()
z = np.load(D / 'module_scores.npz'); meta = json.loads((D / 'module_scores_meta.json').read_text())
for k in meta['sets']: list(DATA.genes.gene_id.iloc[z[f'{k}__rows']])
print(f'module presets       {T()-t:5.2f}s')
t = T(); j2 = pio.to_json(GM.figure(dim=2)); print(f'gene map 2D json     {T()-t:5.2f}s ({len(j2)/1e6:.1f} MB)')
t = T(); j3 = pio.to_json(GM.figure(dim=3)); print(f'gene map 3D json     {T()-t:5.2f}s ({len(j3)/1e6:.1f} MB)')
