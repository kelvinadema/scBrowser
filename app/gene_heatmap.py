"""Expression tab, heatmap view: where each member of a gene set is expressed.

Rows = genes of the set (up to MAX_GENES; larger sets use a fixed, seeded random sample).
Columns = cell groups (the chosen annotation), optionally split by condition into framed blocks.
Value = mean log1p(CP10K) per column over the cells shown (sample filter applies); columns with
fewer than MIN_CELLS cells are empty (missing, not zero). Display = z-score per gene across the
columns (default) or the raw mean.
Clustering = k-means (scikit-learn KMeans, k-means++ init, 10 restarts, seed 0) on the per-gene z-scores;
genes are ordered by cluster, clusters and genes by the column where they peak.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

MAX_GENES = 400
NAMES_MAX = 250      # row labels up to this many genes (smaller font above 80)
MIN_CELLS = 20
SEED = 0
MUTED = {"not annotated", "Unassigned"}
EXPRESSION_SCALE = [[0.00, "#f2f2f2"], [0.25, "#f2c9b8"], [0.55, "#e8795c"], [1.00, "#c40000"]]
CLUSTER_COLORS = ["#4e79a7", "#f28e2b", "#59a14f", "#e15759", "#76b7b2", "#edc948",
                  "#b07aa1", "#ff9da7", "#9c755f", "#bab0ac", "#86bcb6", "#d37295"]
SHORT_COND = {"Root tip": "Tip", "Mock": "Mock", "NF 3 h": "NF 3h", "BAP 3 h": "BAP 3h",
              "Rhizobia 12 h": "Rh 12", "Rhizobia 24 h": "Rh 24", "Rhizobia 48 h": "Rh 48",
              "Rhizobia 72 h": "Rh 72", "Rhizobia 96 h": "Rh 96"}


def kmeans(x: np.ndarray, k: int, seed: int = SEED) -> np.ndarray:
    """scikit-learn KMeans (k-means++, 10 restarts, fixed seed): same set -> same clusters."""
    from sklearn.cluster import KMeans
    k = max(1, min(k, len(x)))
    return KMeans(n_clusters=k, n_init=10, random_state=seed).fit_predict(x)


class GeneHeatmap:
    def __init__(self, data, condition: pd.Series):
        self.data = data
        self.condition = condition

    def columns(self, group_col: str, by_condition: bool, mask: np.ndarray | None):
        """Column definitions: codes per cell, and [(group, condition|None, n_cells)] kept columns."""
        groups = self.data.cells[group_col]
        g_cats = [c for c in groups.cat.categories if c not in MUTED] + \
                 [c for c in groups.cat.categories if c in MUTED]
        g_code = pd.Categorical(groups.astype(str), categories=g_cats).codes.astype(np.int64)
        if by_condition:
            conds = list(self.condition.cat.categories)
            c_code = self.condition.cat.codes.to_numpy().astype(np.int64)
            codes = np.where(g_code >= 0, g_code * len(conds) + c_code, -1)
            labels = [(g, c) for g in g_cats for c in conds]
        else:
            codes = g_code
            labels = [(g, None) for g in g_cats]
        counted = np.where(mask, codes, -1) if mask is not None else codes
        n = np.bincount(counted[counted >= 0], minlength=len(labels))
        keep = [i for i in range(len(labels)) if n[i] >= MIN_CELLS]
        return codes, labels, n, keep

    def build(self, rows: list[int], group_col: str, by_condition: bool, mask, scale: str, k: int):
        """Returns (figure, {cluster name: [gene ids]}, note)."""
        note = ""
        if len(rows) > MAX_GENES:
            rng = np.random.default_rng(SEED)
            rows = sorted(rng.choice(rows, MAX_GENES, replace=False).tolist())
            note = f"> {MAX_GENES} genes selected; showing a {MAX_GENES} gene subsample"
        codes, labels, n, keep = self.columns(group_col, by_condition, mask)
        if not keep:
            return None, {}, "no cell group has at least 20 cells with the current sample selection"
        means = self.data.group_means(rows, codes, len(labels), keep=mask)[:, keep]
        cols = [labels[i] for i in keep]
        n_cols = n[keep]
        genes = self.data.genes.iloc[rows]
        ids, names = genes.gene_id.to_numpy(), genes.label.to_numpy()

        mu = np.nanmean(means, axis=1, keepdims=True)
        sd = np.nanstd(means, axis=1, keepdims=True)
        with np.errstate(invalid="ignore", divide="ignore"):
            z = np.where(sd > 0, (means - mu) / sd, 0.0)
        expressed = np.nanmax(means, axis=1) > 0
        peak = np.nanargmax(np.where(np.isfinite(z), z, -np.inf), axis=1)

        if k and k >= 2 and expressed.sum() > k:
            lab = np.full(len(rows), -1)
            lab[expressed] = kmeans(np.nan_to_num(z[expressed]), k)
            # order clusters by the mean peak position of their genes
            order_c = sorted(set(lab[lab >= 0]), key=lambda c: np.mean(peak[lab == c]))
            rename = {c: i for i, c in enumerate(order_c)}
            lab = np.array([rename.get(c, -1) for c in lab])
        else:
            lab = np.where(expressed, 0, -1)
        # undetected genes last; within a cluster, by peak column then by peak height
        key = np.lexsort((-np.nanmax(z, axis=1), peak, np.where(lab < 0, 999, lab)))
        z, means, ids, names, lab = z[key], means[key], ids[key], names[key], lab[key]
        clusters = {}
        if k and k >= 2:
            for c in range(lab.max() + 1):
                clusters[f"C{c + 1}"] = list(ids[lab == c])

        values = z if scale == "z" else means
        # spacer columns between groups when split by condition
        x_idx, x_lab, blocks, cur = [], [], [], None
        for j, (g, c) in enumerate(cols):
            if by_condition and g != cur:
                if cur is not None:
                    x_idx.append(None); x_lab.append(" " * (len(x_idx)))
                blocks.append([g, len(x_idx), len(x_idx)])
                cur = g
            x_idx.append(j)
            x_lab.append(f"{SHORT_COND.get(c, c)}·{len(x_idx)}" if by_condition else g)
            if by_condition:
                blocks[-1][2] = len(x_idx) - 1
        grid = np.full((len(ids), len(x_idx)), np.nan)
        hover = np.full(grid.shape, "", dtype=object)
        for xi, j in enumerate(x_idx):
            if j is None:
                continue
            grid[:, xi] = values[:, j]
            g, c = cols[j]
            where = g + (f" · {c}" if c else "")
            for r in range(len(ids)):
                hover[r, xi] = (f"{names[r]}<br>{where} ({n_cols[j]:,} cells)<br>mean {means[r, j]:.3f}"
                                f" · z {z[r, j]:+.2f}<br><i>click to select this gene</i>")

        fig = make_subplots(rows=1, cols=2, shared_yaxes=True, column_widths=[0.025, 0.975],
                            horizontal_spacing=0.004)
        if clusters:
            palette = [CLUSTER_COLORS[i % len(CLUSTER_COLORS)] for i in range(lab.max() + 1)]
            bar = np.where(lab >= 0, lab, np.nan)[:, None].astype(float)
            cs = [[i / max(1, len(palette) - 1), col] for i, col in enumerate(palette)] if len(palette) > 1 \
                else [[0, palette[0]], [1, palette[0]]]
            fig.add_trace(go.Heatmap(
                z=bar, x=["cluster"], y=list(range(len(ids))), colorscale=cs, zmin=0,
                zmax=max(1, len(palette) - 1), showscale=False,
                customdata=np.array([[f"C{l + 1}" if l >= 0 else "not detected"] for l in lab], dtype=object),
                hovertemplate="%{customdata}<extra></extra>"), row=1, col=1)
        if scale == "z":
            colour = dict(colorscale="RdBu_r", zmin=-2.5, zmax=2.5, zmid=0,
                          colorbar=dict(title="z-score<br>per gene", thickness=10, len=0.5))
        else:
            top = float(np.nanquantile(means, 0.99)) if np.isfinite(means).any() else 1.0
            colour = dict(colorscale=EXPRESSION_SCALE, zmin=0, zmax=max(top, 1e-3),
                          colorbar=dict(title="mean<br>log1p CP10K", thickness=10, len=0.5))
        fig.add_trace(go.Heatmap(
            z=grid, x=list(range(len(x_idx))), y=list(range(len(ids))), customdata=hover,
            hovertemplate="%{customdata}<extra></extra>", xgap=0, ygap=0, **colour), row=1, col=2)

        show_names = len(ids) <= NAMES_MAX
        row_px = 14 if len(ids) <= 80 else 10
        name_size = 9 if len(ids) <= 80 else 7
        fig.update_yaxes(autorange="reversed", tickmode="array", tickvals=list(range(len(ids))),
                         ticktext=list(names) if show_names else [""] * len(ids),
                         tickfont=dict(size=name_size), showgrid=False, zeroline=False, row=1, col=1)
        fig.update_yaxes(showgrid=False, zeroline=False, row=1, col=2)
        fig.update_xaxes(showticklabels=False, showgrid=False, row=1, col=1)
        fig.update_xaxes(tickmode="array", tickvals=[i for i, j in enumerate(x_idx) if j is not None],
                         ticktext=[x_lab[i].split("·")[0] for i, j in enumerate(x_idx) if j is not None],
                         tickangle=-90 if by_condition else -40, tickfont=dict(size=9),
                         showgrid=False, zeroline=False, row=1, col=2)
        frame = dict(type="rect", xref="x2", yref="y2", line=dict(color="#111", width=1.2),
                     fillcolor="rgba(0,0,0,0)")
        y0, y1 = -0.5, len(ids) - 0.5
        if by_condition:
            for g, a, b in blocks:
                fig.add_shape(x0=a - 0.5, x1=b + 0.5, y0=y0, y1=y1, **frame)
                fig.add_annotation(x=(a + b) / 2, y=1.0, xref="x2", yref="paper", yanchor="bottom",
                                   text=g, showarrow=False, textangle=-30 if len(blocks) > 5 else 0,
                                   font=dict(size=10, color="#333"))
        else:
            fig.add_shape(x0=-0.5, x1=len(x_idx) - 0.5, y0=y0, y1=y1, **frame)
        # cluster boundaries
        for c in range(1, (lab.max() + 1) if clusters else 0):
            edge = int(np.argmax(lab == c)) - 0.5
            fig.add_shape(type="line", xref="paper", yref="y2", x0=0, x1=1, y0=edge, y1=edge,
                          line=dict(color="#111", width=0.8, dash="dot"))
        height = int(min(3200, max(420, 170 + (row_px if show_names else max(2.0, 900 / len(ids))) * len(ids))))
        # row index -> gene id, so a click on a row can select that gene (app.py JS handler)
        fig.update_layout(meta={"genes": [str(g) for g in ids]})
        fig.update_layout(template="plotly_white", plot_bgcolor="white", height=height,
                          margin=dict(l=10, r=10, t=70 if by_condition else 20, b=10),
                          font=dict(size=11))
        return fig, clusters, note

    @staticmethod
    def highlight(fig, gene_id: str | None):
        """Copy of a heatmap with the row of gene_id outlined and labelled (no recompute)."""
        if fig is None or not gene_id:
            return fig
        genes = list((fig.layout.meta or {}).get("genes", []))
        if gene_id not in genes:
            return fig
        r = genes.index(gene_id)
        out = go.Figure(fig)
        out.add_shape(type="rect", xref="paper", yref="y2", x0=0, x1=1, y0=r - 0.5, y1=r + 0.5,
                      line=dict(color="#e0a800", width=2.5), fillcolor="rgba(0,0,0,0)")
        return out
