"""Dissociation-response panel: bulk RNA-seq isolation effect (Medicago).

Per gene: apeglm-shrunk log2FC vs intact tissue for protoplasting and nuclei extraction
(DESeq2), the source analysis' call (lfcThreshold = 1, padj < 0.01) and normalised counts per
sample. Several genes: heatmap of the two log2FCs, * = dissociation-responsive.
"""

from __future__ import annotations

import json
from html import escape
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Modality colours (col_modality in the figure code)
COLORS = {"tissue": "#8a8a8a", "protoplast": "#2a78d6", "nuclei": "#eb6834"}
NAMES = {"protoplast": "Protoplasting", "nuclei": "Nuclei extraction"}
LFC_CAP = 6.0


def _p(p: float) -> str:
    if p is None or not np.isfinite(p):
        return "padj n/a"
    return f"padj {p:.1e}" if p < 1e-3 else f"padj {p:.3f}"


class Dissociation:
    def __init__(self, directory: Path):
        path = Path(directory) / "dissociation.parquet"
        self.available = path.exists()
        if not self.available:
            return
        self.table = pd.read_parquet(path).set_index("gene_id")
        self.meta = json.loads((Path(directory) / "dissociation_meta.json").read_text())
        self.samples = self.meta["samples"]  # sample -> method, in design order

    def row(self, gene_id: str) -> pd.Series | None:
        if not self.available or gene_id not in self.table.index:
            return None
        return self.table.loc[gene_id]

    def verdict(self, gene_id: str) -> str:
        r = self.row(gene_id)
        if r is None:
            return "Not tested (organellar, or fewer than 3 samples with ≥10 counts)"
        parts = []
        for k in ("protoplast", "nuclei"):
            if r[f"{k}_call"]:
                parts.append(f"{'induced' if r[f'{k}_lfc'] > 0 else 'repressed'} by "
                             f"{NAMES[k].lower()}")
        return ("dissociation-responsive: " + "; ".join(parts)) if parts else \
            "not dissociation-responsive (neither contrast meets both |log2FC| > 1 and padj < 0.01)"

    def gene_figure(self, gene_id: str, label: str) -> go.Figure | None:
        r = self.row(gene_id)
        if r is None:
            return None
        fig = make_subplots(rows=2, cols=1, row_heights=[0.55, 0.45], vertical_spacing=0.2,
                            subplot_titles=("log2 fold change vs intact tissue",
                                            "Normalised counts per sample"))
        for k in ("protoplast", "nuclei"):
            lfc, se, call = float(r[f"{k}_lfc"]), float(r[f"{k}_se"]), bool(r[f"{k}_call"])
            fig.add_trace(go.Bar(
                x=[NAMES[k]], y=[lfc],
                error_y=dict(type="data", array=[se], color="#444", thickness=1.2, width=6),
                marker=dict(color=COLORS[k], opacity=1.0 if call else 0.35,
                            line=dict(color=COLORS[k], width=1.5)),
                text=[("responsive" if call else "criterion not met") + "<br>" + _p(r[f"{k}_padj"])],
                textposition="outside", textfont=dict(size=10), cliponaxis=False,
                hovertemplate=f"{NAMES[k]}<br>log2FC %{{y:.2f}} ± {se:.2f}<br>{_p(r[f'{k}_padj'])}"
                              "<extra></extra>",
                showlegend=False,
            ), row=1, col=1)
        lim = max(1.6, *(abs(float(r[f"{k}_lfc"])) + float(r[f"{k}_se"]) for k in NAMES)) * 1.45
        for y in (-1, 1):
            fig.add_hline(y=y, line=dict(color="#999", dash="dot", width=1), row=1, col=1)
        fig.add_hline(y=0, line=dict(color="#333", width=1), row=1, col=1)
        fig.update_yaxes(range=[-lim, lim], title_text="log2FC", row=1, col=1)

        groups = {"tissue": "Intact tissue", **NAMES}
        for method, name in groups.items():
            samples = [s for s, m in self.samples.items() if m == method]
            vals = [float(r[f"norm_{s}"]) for s in samples]
            fig.add_trace(go.Box(
                x=[name] * len(vals), y=vals, boxpoints="all", pointpos=0, jitter=0.25,
                fillcolor="rgba(0,0,0,0)", line=dict(color=COLORS[method], width=1.2),
                marker=dict(color=COLORS[method], size=7), text=samples,
                hovertemplate="%{text}: %{y:.1f}<extra></extra>", showlegend=False,
            ), row=2, col=1)
        fig.update_yaxes(title_text="counts / size factor", rangemode="tozero", row=2, col=1)
        fig.update_layout(template="plotly_white", height=510,
                          title=dict(text=escape(label or gene_id), x=0, xref="container",
                                     xanchor="left", y=1, yref="container", yanchor="top",
                                     pad=dict(t=8, l=10), font=dict(size=16)),
                          margin=dict(l=55, r=10, t=70, b=10), font=dict(size=11))
        for a in fig.layout.annotations:
            a.font = dict(size=11, color="#333")
        return fig

    def heatmap(self, genes: list[tuple[str, str]]) -> go.Figure | None:
        """genes: [(gene_id, label)]; omit genes without dissociation results."""
        if not genes:
            return None
        labels, z, star, hover = [], [], [], []
        for gid, label in genes:
            r = self.row(gid)
            if r is None:
                continue
            labels.append(label)
            z.append([float(r[f"{k}_lfc"]) for k in NAMES])
            star.append(["*" if r[f"{k}_call"] else "" for k in NAMES])
            hover.append([f"{NAMES[k]}: log2FC {float(r[f'{k}_lfc']):.2f}, {_p(r[f'{k}_padj'])}"
                          for k in NAMES])
        if not labels:
            return None
        fig = go.Figure(go.Heatmap(
            z=z, x=list(NAMES.values()), y=labels, text=star, texttemplate="%{text}",
            textfont=dict(size=16, color="#111"), customdata=hover,
            hovertemplate="%{y}<br>%{customdata}<extra></extra>",
            colorscale="RdBu_r", zmin=-LFC_CAP, zmax=LFC_CAP, zmid=0, xgap=2, ygap=1,
            colorbar=dict(title="log2FC", thickness=10, len=0.8),
        ))
        fig.update_yaxes(autorange="reversed", automargin=True)
        fig.update_layout(template="plotly_white", height=min(3600, max(220, 70 + 22 * len(labels))),
                          margin=dict(l=10, r=10, t=10, b=10), font=dict(size=11))
        return fig
