"""Regression checks with the shipped real data; no scientific fixture substitutes."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'app'))
import app as a
ids=['MtrunA17_Chr5g0448621','MtrunA17_Chr3g0134941']
rows=[a.DATA.resolve_gene(g) for g in ids]
assert all(r is not None for r in rows)
report=a.SCORER.resolve_report(' '.join(ids+ids+['not_a_gene']))
assert report['rows']==rows and len(report['duplicates'])==2 and report['missing']==['not_a_gene']
# New reporting must leave unique-ID module scoring unchanged.
score,n=a.SCORER.score(rows)
ctrl=a.SCORER.controls(rows)
reference=a.DATA.expression_mean(rows)-a.DATA.expression_mean(list(ctrl))
assert np.allclose(score,reference,atol=1e-5)
mask=np.isin(a.SAMPLE_COLUMN,a.TIP_SAMPLES)
full=a.HEATMAP.build(rows,'ch5_tissue',False,None,'z',0)
filtered=a.HEATMAP.build(rows,'ch5_tissue',False,mask,'raw',0)
assert full[0] is not None and filtered[0] is not None
preset=a.SCORER.resolve_report('\n'.join(a.MODULE_PRESETS['GEP12'][1]))['rows']
x=a.HEATMAP.build(preset,'ch5_tissue',False,None,'z',4)
y=a.HEATMAP.build(preset,'ch5_tissue',False,None,'raw',4)
assert x[1]==y[1], 'Clustering must not depend on display scale'
assert len(x[1])==4
fig=a.GENE_MAP.figure(dim=2)
assert len(fig.data)>0
values=np.expm1(a.DATA.expression(rows[0]).astype(float))*a.FATE_CELLS.cpm_factor
prof=a.FATE_CELLS.profiles(values,'detection','pooled')
assert all(np.nanmin(v[0])>=0 and np.nanmax(v[0])<=1 for v in prof.values())
assert a.FATE.figure(ids[0]) is not None and a.DISSOC.gene_figure(ids[0],'MtNIN') is not None
print(f'PASS: real-data parsing, score identity ({n} controls), sample-filtered heatmap, reproducible clusters ({len(x[1])}), gene map, fate detection units, both evidence figures')
# Edge cases that would alter scientific interpretation if the UI obscured them.
large=list(range(500))
big1=a.HEATMAP.build(large,'ch5_tissue',False,None,'z',0)
big2=a.HEATMAP.build(large,'ch5_tissue',False,None,'raw',0)
assert big1[0].layout.meta['genes']==big2[0].layout.meta['genes']
assert len(big1[0].layout.meta['genes'])==400 and big1[2] == '> 400 genes selected; showing a 400 gene subsample'
empty=a.HEATMAP.build(rows,'ch5_tissue',False,np.zeros(a.DATA.n_cells,dtype=bool),'z',0)
assert empty[0] is None and 'no cell group' in empty[2]
# Explicit fixture for ambiguous aliases (not an invented atlas finding).
from types import SimpleNamespace
import pandas as pd
probe=a.GeneSetScorer.__new__(a.GeneSetScorer)
probe.data=SimpleNamespace(genes=pd.DataFrame({'gene_id':['id1','id2'],'symbol':['alias','alias']}))
r=probe.resolve_report('alias id1 ID1 unknown')
assert r==dict(rows=[0],missing=['unknown'],duplicates=['ID1'],ambiguous=['alias'],omitted=0)
print('PASS: 400-gene deterministic sample, empty selection, and ambiguous-alias reporting')
