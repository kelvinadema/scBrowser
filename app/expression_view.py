"""Expression tab: dot plot of a gene (or gene-set score) per cell group x condition.

Dot size = % of cells expressing (gene: > 0 UMI; set: score > 0); dot colour = mean over all
cells of the group (gene: log1p CP10K; set: module score, diverging around 0).
Columns = condition (root tip, mock, NF/BAP 3 h, rhizobia 12-96 h), optionally split by platform
(BD / 10x protoplast / 10x nuclei). Groups with fewer than MIN_CELLS cells are left empty:
missing is not zero.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go

MIN_CELLS = 20
EXPRESSION_SCALE = [[0.00, "#dedede"], [0.25, "#f2c9b8"], [0.55, "#e8795c"], [1.00, "#c40000"]]
MUTED = {"not annotated", "Unassigned"}


class DotPlot:
    def __init__(self, data, directory: Path):
        path = Path(directory) / "cell_design.parquet"
        self.available = path.exists()
        if not self.available:
            return
        self.data = data
        d = pd.read_parquet(path)
        assert (d.cell_id.to_numpy() == data.cells.cell_id.to_numpy()).all()
        self.condition = d.condition
        self.platform = d.platform

    def figure(self, values: np.ndarray, kind: str, row_column: str, split_platform: bool,
               mask: np.ndarray | None, title: str) -> go.Figure:
        rows = self.data.cells[row_column]
        frame = pd.DataFrame({"row": rows.astype(str), "cond": self.condition.astype(str),
                              "plat": self.platform.astype(str), "v": values})
        if mask is not None:
            frame = frame[mask]
        frame["col"] = frame.cond + (" · " + frame.plat if split_platform else "")
        frame["pos"] = frame.v > 0
        g = frame.groupby(["row", "col"], observed=True).agg(cells=("v", "size"), frac=("pos", "mean"),
                                                             mean=("v", "mean")).reset_index()
        g = g[g.cells >= MIN_CELLS]
        if g.empty:
            return _empty("No group has at least 20 cells with the current sample selection.")
        # column order: condition order, then platform order
        cond_order = list(self.condition.cat.categories)
        plat_order = list(self.platform.cat.categories)
        cols = [c + (f" · {p}" if split_platform else "") for c in cond_order
                for p in (plat_order if split_platform else [None])]
        cols = [c for c in cols if c in set(g.col)]
        row_order = [r for r in rows.cat.categories if r in set(g.row)]
        row_order = [r for r in row_order if r not in MUTED] + [r for r in row_order if r in MUTED]

        size = 100 * g.frac.to_numpy()
        max_px = 22
        if kind == "set":
            lim = max(float(np.quantile(np.abs(g["mean"]), 0.98)), 1e-3)
            colour = dict(color=g["mean"], colorscale="RdBu_r", cmin=-lim, cmax=lim, cmid=0)
            ctitle, ptitle = "mean<br>score", "% score &gt; 0"
        else:
            colour = dict(color=g["mean"], colorscale=EXPRESSION_SCALE, cmin=0,
                          cmax=max(float(g["mean"].max()), 1e-3))
            ctitle, ptitle = "mean<br>log1p CP10K", "% expressing"
        fig = go.Figure(go.Scatter(
            x=g.col, y=g.row, mode="markers",
            marker=dict(size=size, sizemode="area", sizeref=2.0 * 100 / max_px ** 2, sizemin=1,
                        line=dict(color="#666", width=0.5), showscale=True,
                        colorbar=dict(title=ctitle, thickness=12, len=0.5, y=0.75), **colour),
            customdata=np.stack([g.cells, g.frac * 100, g["mean"]], axis=1),
            hovertemplate=("%{y} · %{x}<br>%{customdata[0]:,} cells<br>" + ptitle.replace("&gt;", ">")
                           + " %{customdata[1]:.1f}%<br>mean %{customdata[2]:.3f}<extra></extra>"),
            showlegend=False,
        ))
        for pct in (10, 50, 100):   # size key
            fig.add_trace(go.Scatter(
                x=[None], y=[None], mode="markers", name=f"{pct}%",
                marker=dict(size=pct, sizemode="area", sizeref=2.0 * 100 / max_px ** 2,
                            color="#bbb", line=dict(color="#666", width=0.5)),
                showlegend=True, hoverinfo="skip"))
        fig.update_xaxes(categoryorder="array", categoryarray=cols, tickangle=-40, showgrid=True,
                         gridcolor="#f0f0f0", automargin=True)
        fig.update_yaxes(categoryorder="array", categoryarray=row_order[::-1], showgrid=True,
                         gridcolor="#f0f0f0", automargin=True)
        fig.update_layout(template="plotly_white", height=max(360, 140 + 34 * len(row_order)),
                          margin=dict(l=10, r=10, t=40, b=10), title=dict(text=title, x=0, font=dict(size=13)),
                          legend=dict(title=ptitle, orientation="v", x=1.02, y=0.25, font=dict(size=11)))
        return fig


def _empty(text: str) -> go.Figure:
    fig = go.Figure()
    fig.update_layout(template="plotly_white", height=300, xaxis=dict(visible=False),
                      yaxis=dict(visible=False),
                      annotations=[dict(text=text, showarrow=False, font=dict(size=14, color="#666"))])
    return fig
