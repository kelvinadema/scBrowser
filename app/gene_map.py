"""Gene UMAP tab: the 15,531 fate-associated genes, modules M01-M10."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go

CURATED_COLOR = "#000000"
HIGHLIGHT_COLOR = "#f6cf33"


class GeneMap:
    def __init__(self, directory: Path):
        path = Path(directory) / "gene_umap.parquet"
        self.available = path.exists()
        if not self.available:
            return
        self.genes = pd.read_parquet(path)
        self.meta = json.loads((Path(directory) / "gene_umap_meta.json").read_text())
        self.order = self.meta["module_order"]  # legend order
        self.position = {g: i for i, g in enumerate(self.genes.gene_id)}

    @property
    def has_3d(self) -> bool:
        return self.available and "umap3d_1" in self.genes

    def coords(self, frame: pd.DataFrame, dim: int) -> dict:
        if dim == 3:
            return dict(x=frame.umap3d_1, y=frame.umap3d_2, z=frame.umap3d_3)
        return dict(x=frame.umap_1, y=frame.umap_2)

    def locate(self, gene_id: str | None, dim: int = 2) -> dict | None:
        if not self.available or gene_id is None or gene_id not in self.position:
            return None
        row = self.genes.iloc[[self.position[gene_id]]]
        xyz = {k: float(v.iloc[0]) for k, v in self.coords(row, dim).items()}
        return {**xyz, "z": xyz.get("z"), "label": str(row.label.iloc[0]),
                "module": str(row.module.iloc[0])}

    def figure(self, dim: int = 2, point_size: float | None = None) -> go.Figure:
        g = self.genes
        three = dim == 3
        Scatter = go.Scatter3d if three else go.Scattergl
        if point_size is None:
            point_size = 2.2 if three else 3.0
        fig = go.Figure()
        hover = ("<b>%{text}</b><br>%{customdata[0]}<br>%{customdata[1]}<br>"
                 "%{customdata[2]}<extra></extra>")

        # grey "mixed" modules underneath, tissue-high modules on top
        draw = sorted(self.order, key=lambda m: (not self.meta["modules"][m]["name"].endswith("mixed"),
                                                 -self.meta["modules"][m]["n_genes"]))
        legend_rank = {m: 100 + i for i, m in enumerate(self.order)}
        for mid in draw:
            q = g[g.module_id == mid]
            info = self.meta["modules"][mid]
            fig.add_trace(Scatter(
                **self.coords(q, dim), mode="markers",
                name=f"{info['name']} ({len(q):,})", legendrank=legend_rank[mid],
                marker=dict(size=point_size, color=info["color"], opacity=0.8,
                            **({} if three else {"line": dict(width=0)})),
                text=q.label,
                customdata=np.stack([q.gene_id, q.module, q.mock_change], axis=1),
                hovertemplate=hover,
            ))

        c = g[g.curated]
        fig.add_trace(Scatter(
            **self.coords(c, dim), mode="markers",
            name=f"known nodulation gene ({len(c)})", legendrank=50,
            marker=dict(size=point_size + (2 if three else 3), color=CURATED_COLOR),
            text=c.label,
            customdata=np.stack([c.gene_id, c.module, c.mock_change], axis=1),
            hovertemplate=hover,
        ))

        # selected-gene marker; filled in by the 'geneumap-highlight' message
        ring = (dict(size=9, color=HIGHLIGHT_COLOR, line=dict(color="#111", width=2)) if three else
                dict(size=16, color="rgba(0,0,0,0)", line=dict(color=HIGHLIGHT_COLOR, width=3)))
        fig.add_trace((go.Scatter3d if three else go.Scatter)(
            x=[], y=[], **({"z": []} if three else {}),
            mode="markers+text", name="selected gene", showlegend=False,
            marker=ring,
            text=[], textposition="top center",
            textfont=dict(size=13, color="#111"), hoverinfo="skip",
        ))

        labels = []
        for mid in self.order:
            q = g[g.module_id == mid]
            info = self.meta["modules"][mid]
            pos = {k: float(v.median()) for k, v in self.coords(q, dim).items()}
            labels.append(dict(**pos, text=f"<b>{mid}</b>", showarrow=False,
                               font=dict(size=13, color=info["color"]),
                               bgcolor="rgba(255,255,255,0.55)"))

        if three:
            axis = lambda t: dict(title=t, showgrid=True, gridcolor="#e8e8e8", zeroline=False,
                                  showticklabels=False, showbackground=False, showspikes=False)
            fig.update_layout(scene=dict(xaxis=axis("UMAP 1"), yaxis=axis("UMAP 2"),
                                         zaxis=axis("UMAP 3"), aspectmode="data",
                                         annotations=labels, uirevision="gene-map-3d"))
        else:
            for a in labels:
                fig.add_annotation(**a)
            fig.update_layout(
                xaxis=dict(title="UMAP 1", showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(title="UMAP 2", showgrid=False, zeroline=False, showticklabels=False,
                           scaleanchor="x", scaleratio=1),
                dragmode="pan")

        fig.update_layout(
            template="plotly_white",
            margin=dict(l=10, r=10, t=10, b=10),
            legend=dict(itemsizing="constant", font=dict(size=11), x=1.0, xanchor="right",
                        y=1.0, yanchor="top", bgcolor="rgba(255,255,255,0.8)",
                        bordercolor="rgba(0,0,0,0.12)", borderwidth=1),
            hovermode="closest",
            uirevision=f"gene-map-{dim}d",
        )
        return fig
