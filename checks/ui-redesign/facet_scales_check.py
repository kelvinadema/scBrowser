"""Shared facet geometry and colour scales, including an empty first fate facet."""
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'app'))
import app as a
from plotting import PlotSettings, build_figure
f=build_figure(a.DATA,PlotSettings(mode='fate',values=a.FATE_CELLS.p,facet='ch5_region',max_points=50000),a.AXES)
visible=[t for t in f.data if t.marker.showscale and len(t.x)]
assert len(visible)==1 and visible[0].marker.colorbar.x>1
assert f.layout.margin.r>=100
assert f.layout.scene.aspectratio.to_plotly_json()==f.layout.scene2.aspectratio.to_plotly_json()
# This atlas exposes categorical annotations only. Exercise the numeric branch
# with a test-only column using the actual UMAP coordinate values.
from copy import copy
numeric=copy(a.DATA)
numeric.cells=a.DATA.cells.assign(num__test=a.DATA.coords[:,0])
f=build_figure(numeric,PlotSettings(annotation='num__test',facet='ch5_region',max_points=50000),a.AXES)
coloured=[t for t in f.data if t.marker.color is not None and not isinstance(t.marker.color,str)]
assert len({(t.marker.cmin,t.marker.cmax) for t in coloured})==1
assert sum(bool(t.marker.showscale) for t in coloured)==1
v=numeric.cells['num__test'].to_numpy();v=v[np.isfinite(v)]
assert np.isclose(coloured[0].marker.cmin,v.min()) and np.isclose(coloured[0].marker.cmax,v.max())
print('PASS: fate scale survives unmodelled first facet, one external shared numeric colourbar with correct bounds, identical facet baseline geometry.')
