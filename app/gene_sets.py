"""Gene sets: parsing, module score and gene-set heatmaps for the sidebar.

Module score (Seurat AddModuleScore / scanpy score_genes style): per cell, the mean
log1p(CP10K) of the set's genes minus the mean of expression-matched control genes.
All genes are split into 25 bins by mean expression over all cells; each bin has a fixed
panel of 50 control genes (seed 0, build/build_control_panels.py). A set is compared with
the panels of every bin its genes fall in; set genes that sit in a panel are removed from it.
Positive = above background for genes of that expression level; centred, so comparable
between sets.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go

MAX_GENES = 5000

SOURCES = {"root_hair_cell": "Root hair", "cortex": "Cortex",
           "pericycle_cell": "Pericycle", "procambium": "Procambium"}
TIMES = [12, 24, 48, 72]
FATE_CAP = 4.0


def parse_tokens(text: str) -> list[str]:
    return [t for t in re.split(r"[\s,;]+", text or "") if t]


class GeneSetScorer:
    def __init__(self, data, directory: Path):
        self.data = data
        path = Path(directory) / "control_panels.npz"
        self.available = path.exists()
        self._cache: dict[tuple[int, ...], tuple[np.ndarray, int]] = {}
        if not self.available:
            return
        z = np.load(path)
        self.bins = z["bins"].astype(int)
        self.panel_rows = z["panel_rows"]
        self.panel_sums = z["panel_sums"]
        assert len(self.bins) == data.n_genes and self.panel_sums.shape[0] == data.n_cells, \
            "control_panels.npz is out of step with the viewer data; rerun build/build_control_panels.py"
        # precomputed scores of the fate-associated module presets (build/build_module_scores.py)
        pre = Path(directory) / "module_scores.npz"
        if pre.exists():
            z = np.load(pre)
            for key in z.files:
                if key.endswith("__rows"):
                    mid = key[:-6]
                    self._preset = getattr(self, "_preset", {})
                    self._preset[tuple(int(r) for r in z[key])] = (z[f"{mid}__score"],
                                                                   int(z[f"{mid}__n_ctrl"]))

    def resolve(self, text: str) -> tuple[list[int], list[str]]:
        rows, missing, seen = [], [], set()
        for token in parse_tokens(text):
            row = self.data.resolve_gene(token)
            if row is None:
                missing.append(token)
            elif row not in seen:
                seen.add(row)
                rows.append(row)
        return rows[:MAX_GENES], missing

    def resolve_report(self, text: str) -> dict:
        ids = {str(g).lower(): i for i, g in enumerate(self.data.genes.gene_id)}
        symbols = {}
        for i, symbol in enumerate(self.data.genes.symbol):
            if symbol:
                symbols.setdefault(str(symbol).lower(), []).append(i)
        rows, missing, duplicates, ambiguous, seen = [], [], [], [], set()
        for token in parse_tokens(text):
            key = token.lower()
            matches = [ids[key]] if key in ids else symbols.get(key, [])
            if len(matches) > 1:
                ambiguous.append(token)
            elif not matches:
                missing.append(token)
            elif matches[0] in seen:
                duplicates.append(token)
            else:
                seen.add(matches[0])
                rows.append(matches[0])
        return dict(rows=rows[:MAX_GENES], missing=missing, duplicates=duplicates,
                    ambiguous=ambiguous, omitted=max(0, len(rows) - MAX_GENES))

    def controls(self, rows: list[int]) -> np.ndarray:
        """Control genes actually used: panels of the occupied bins minus set genes."""
        used = np.unique(self.bins[rows])
        ctrl = self.panel_rows[used].ravel()
        ctrl = ctrl[ctrl >= 0]
        return np.setdiff1d(ctrl, rows)

    def score(self, rows: list[int]) -> tuple[np.ndarray, int]:
        key = tuple(sorted(rows))
        if key in getattr(self, "_preset", {}):
            return self._preset[key]
        if key in self._cache:
            return self._cache[key]
        used = np.unique(self.bins[list(key)])
        ctrl_sum = self.panel_sums[:, used].sum(axis=1, dtype=np.float64)
        panel = self.panel_rows[used].ravel()
        panel = panel[panel >= 0]
        overlap = np.intersect1d(panel, key)
        if overlap.size:
            ctrl_sum -= self.data.expression_mean(list(overlap)).astype(np.float64) * overlap.size
        n_ctrl = int(panel.size - overlap.size)
        values = (self.data.expression_mean(list(key)) - ctrl_sum / n_ctrl).astype(np.float32)
        if len(self._cache) > 16:
            self._cache.pop(next(iter(self._cache)))
        self._cache[key] = (values, n_ctrl)
        return values, n_ctrl


def fate_heatmap(fate, genes: list[tuple[str, str]]) -> go.Figure | None:
    """genes: [(gene_id, label)] -> genes x (source, time) heatmap of the fate association."""
    if not genes or fate is None or not fate.available:
        return None
    cols = [(s, t) for s in SOURCES for t in TIMES]
    z, text, hover = [], [], []
    for gid, label in genes:
        d = fate.rows(gid)
        d = d.set_index(["source", "time_h"]) if d is not None else None
        zr, tr, hr = [], [], []
        for s, t in cols:
            where = f"{label} · {SOURCES[s]} {t} HPI"
            if d is None or (s, t) not in d.index:
                zr.append(np.nan); tr.append("")
                hr.append(f"{where}<br>not eligible / not in the fate screen")
                continue
            r = d.loc[(s, t)]
            e = float(r.log2_effect)
            if r.status == "excluded_fit":
                zr.append(np.nan); tr.append("×")
                hr.append(f"{where}<br>excluded fit ({r.exclusion_reason})")
                continue
            zr.append(e)
            tr.append({"figure7": "●", "fate_associated": "○"}.get(r.status, ""))
            label_status = {"figure7": "Fate-associated, with ≥1.5-fold change vs matched mock in the same direction",
                            "fate_associated": "fate-associated, not changed vs mock",
                            "below_threshold": "below screen thresholds"}[r.status]
            hr.append(f"{where}<br>log2 change {e:+.2f}<br>{label_status}")
        z.append(zr); text.append(tr); hover.append(hr)
    # one framed block per source tissue, separated by white space
    from plotly.subplots import make_subplots
    z, text, hover = np.array(z, dtype=float), np.array(text, dtype=object), np.array(hover, dtype=object)
    labels = [g[1] for g in genes]
    fig = make_subplots(rows=1, cols=len(SOURCES), shared_yaxes=True, horizontal_spacing=0.035,
                        subplot_titles=list(SOURCES.values()))
    n_t = len(TIMES)
    for k in range(len(SOURCES)):
        sl = slice(k * n_t, (k + 1) * n_t)
        fig.add_trace(go.Heatmap(
            z=z[:, sl], x=[str(t) for t in TIMES], y=labels, text=text[:, sl], texttemplate="%{text}",
            textfont=dict(size=9, color="#111"), customdata=hover[:, sl],
            hovertemplate="%{customdata}<extra></extra>",
            colorscale="PRGn", zmin=-FATE_CAP, zmax=FATE_CAP, zmid=0, xgap=1, ygap=1,
            showscale=(k == len(SOURCES) - 1),
            colorbar=dict(title="log2", thickness=10, len=0.6),
        ), row=1, col=k + 1)
    frame = dict(showline=True, mirror=True, linecolor="#111", linewidth=1.2, showgrid=False,
                 zeroline=False, ticks="")
    fig.update_xaxes(type="category", tickfont=dict(size=9), tickangle=0, **frame)
    fig.update_yaxes(autorange="reversed", automargin=True, tickfont=dict(size=10), **frame)
    fig.update_xaxes(title_text="HPI", title_font=dict(size=10), title_standoff=2)
    for a in fig.layout.annotations:
        a.font = dict(size=10, color="#333")
    fig.update_layout(template="plotly_white", plot_bgcolor="white",
                      height=min(3600, max(240, 120 + 16 * len(genes))),
                      margin=dict(l=10, r=10, t=24, b=10), font=dict(size=11))
    return fig
