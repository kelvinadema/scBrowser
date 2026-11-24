"""Fate tab: expression along nodule fate (fate-screen contexts).

Per 12-HPI source tissue and time point, cells of the fate-screen context are split into
ten weighted fate deciles (1 = least, 10 = most nodule-directed; deciles hold equal forward
mass, not equal cell counts). Each curve is the forward-mass-weighted mean per decile:
  gene     -> CPM over common-gene counts (the screen's scale), or weighted detection
  gene set -> module score, or weighted share of cells with score > 0
Platforms: deciles are defined within platform; "pooled" averages the BD and 10x curves with
equal weight (the screen normalises platforms to equal total weight); 12 HPI is BD only.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

SOURCES = {"root_hair_cell": "Root hair", "cortex": "Cortex",
           "pericycle_cell": "Pericycle", "procambium": "Procambium"}
TIMES = [12, 24, 48, 72]
TIME_COLORS = {12: "#9ecae1", 24: "#6baed6", 48: "#3182bd", 72: "#08519c"}   # blues
PLATFORMS = {"BD": "BD", "10xP": "10x"}
STATUS_TEXT = {"figure7": "Fate-associated, with ≥1.5-fold change vs matched mock in the same direction", "fate_associated": "fate-associated, not changed vs mock",
               "below_threshold": "below screen thresholds (requires |log2 change| ≥ log2(1.5) and noise ratio ≥ 3)", "excluded_fit": "fit excluded by the screen (failed to converge, invalid coefficient or noise scale, ill-conditioned fit, or |fate coefficient| > 10)",
               None: "not eligible in the screen (not fitted: &lt;100 retained cells, detected in &lt;20 cells, or &lt;1% platform-balanced weighted detection)"}


class FateCells:
    def __init__(self, directory: Path):
        path = Path(directory) / "fate_cells.parquet"
        self.available = path.exists()
        if not self.available:
            return
        t = pd.read_parquet(path)
        self.p = t.p_nodule96.to_numpy(np.float32)
        self.cpm_factor = t.cpm_factor.to_numpy(np.float64)
        self.time = t.ot_time_h.to_numpy()
        self.platform = t.ot_platform.to_numpy()
        self.w = {s: t[f"w_{s}"].to_numpy(np.float64) for s in SOURCES}
        self.dec = {s: t[f"decile_{s}"].to_numpy() for s in SOURCES}
        self.meta = json.loads((Path(directory) / "fate_cells_meta.json").read_text())

    def profiles(self, values: np.ndarray, measure: str, platform: str) -> dict:
        """{(source, time): (deciles 1..10 means, cells per decile)} for the chosen platform view."""
        x = (values > 0).astype(np.float64) if measure == "detection" else values.astype(np.float64)
        out = {}
        plats = ["BD", "10xP"] if platform == "pooled" else [platform]
        for s in SOURCES:
            d, w = self.dec[s], self.w[s]
            for t in TIMES:
                curves, counts = [], np.zeros(10)
                for pl in plats:
                    m = (d > 0) & (self.time == t) & (self.platform == pl)
                    if not m.any():
                        continue
                    k = d[m] - 1
                    sw = np.bincount(k, weights=w[m], minlength=10)
                    swx = np.bincount(k, weights=w[m] * x[m], minlength=10)
                    with np.errstate(invalid="ignore", divide="ignore"):
                        curves.append(swx / sw)
                    counts += np.bincount(k, minlength=10)
                if curves:
                    out[(s, t)] = (np.nanmean(np.vstack(curves), axis=0), counts.astype(int), len(curves))
        return out

    def figure(self, values: np.ndarray, kind: str, measure: str, platform: str,
               status: dict | None = None) -> go.Figure:
        prof = self.profiles(values, measure, platform)
        fig = make_subplots(rows=2, cols=2, shared_xaxes=True, subplot_titles=list(SOURCES.values()),
                            horizontal_spacing=0.08, vertical_spacing=0.14)
        if kind == "set":
            ylab = "Score >0 (fraction of cells)" if measure == "detection" else "mean module score"
        else:
            ylab = "Expressing cells (weighted fraction)" if measure == "detection" else "Gene expression (weighted mean CPM)"
        deciles = list(range(1, 11))
        for i, (s, name) in enumerate(SOURCES.items()):
            r, c = divmod(i, 2)
            for t in TIMES:
                if (s, t) not in prof:
                    continue
                y, n, n_plat = prof[(s, t)]
                st = (status or {}).get((s, t)) if kind == "gene" else "n/a"
                dash = "solid"
                if kind == "gene":
                    dash = {"figure7": "solid", "fate_associated": "solid", "below_threshold": "solid",
                            "excluded_fit": "dot"}.get(st, "dash")
                note = STATUS_TEXT.get(st, "") if kind == "gene" else ""
                plat_note = {"pooled": "BD + 10x (equal weight)" if n_plat == 2 else "BD only",
                             "BD": "BD", "10xP": "10x only"}[platform]
                fig.add_trace(go.Scatter(
                    x=deciles, y=y, mode="lines+markers", name=f"{t} HPI", legendgroup=str(t),
                    showlegend=(i == 0),
                    line=dict(color=TIME_COLORS[t], width=2.2 if st == "figure7" else 1.6, dash=dash),
                    marker=dict(size=5, color=TIME_COLORS[t]),
                    customdata=np.stack([n], axis=1),
                    hovertemplate=(f"{name} · {t} HPI · decile %{{x}}<br>{ylab} %{{y:.3g}}"
                                   f"<br>%{{customdata[0]:,}} cells · {plat_note}"
                                   + (f"<br>{note}" if note else "") + "<extra></extra>"),
                ), row=r + 1, col=c + 1)
        fig.update_xaxes(tickmode="array", tickvals=deciles, range=[0.5, 10.5], fixedrange=True)
        fig.update_xaxes(title_text="fate decile (1 = least, 10 = most nodule-directed)", row=2, col=1)
        fig.update_xaxes(title_text="fate decile (1 = least, 10 = most nodule-directed)", row=2, col=2)
        fig.update_yaxes(title_text=ylab, col=1, rangemode="tozero" if kind == "gene" or measure == "detection" else None)
        fig.update_yaxes(rangemode="tozero" if kind == "gene" or measure == "detection" else None)
        fig.update_layout(template="plotly_white", margin=dict(l=60, r=20, t=40, b=50),
                          legend=dict(title="", orientation="h", x=1, xanchor="right", y=1.0, yref="container", yanchor="top"),
                          hovermode="closest", font=dict(size=12))
        for a in fig.layout.annotations:
            if a.text in SOURCES.values():
                a.font = dict(size=13, color="#333")
        return fig
