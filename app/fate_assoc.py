"""Fate-association panel: per-gene fate-screen contexts.

One small panel per 12-HPI source tissue (screen order), bars = log2 expression change
between the 10th and 90th weighted nodule-fate ranks at each time point (pooled platforms),
markers = per-platform estimates. Bar fill encodes the call:
  solid  = full criterion (fate-associated AND same-direction >= 1.5-fold vs matched mock)
  faded  = fate-associated, but not changed vs mock in the same direction
  grey   = below the fate-association thresholds
  hatched= fit excluded by the screen (unbounded |fate coefficient| > 10, failed/ill-conditioned)
  "–"    = context not eligible (gene detected in < 20 cells or < 1 % weighted detection)
"""

from __future__ import annotations

import json
from html import escape
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

SOURCES = {"root_hair_cell": "Root hair", "cortex": "Cortex",
           "pericycle_cell": "Pericycle", "procambium": "Procambium"}
TIMES = [12, 24, 48, 72]
HIGHER, LOWER, BELOW = "#1b7837", "#762a83", "#cfcfcc"   # PRGn
SCREEN_LOG2 = float(np.log2(1.5))
NOT_ELIGIBLE = ("not eligible: detected in &lt;20 cells or &lt;1% weighted detection<br>"
                "in this source/time (not fitted, not a zero effect)")


class FateAssociation:
    def __init__(self, directory: Path):
        path = Path(directory) / "fate_association.parquet"
        self.available = path.exists()
        if not self.available:
            return
        # one table sorted by gene + (start, stop) offsets: a cheap index instead of
        # ~20k small DataFrames (that split took ~2.5 s at startup)
        t = pd.read_parquet(path)
        t = t.sort_values("gene_id", kind="stable").reset_index(drop=True)
        ids = t.gene_id.to_numpy()
        starts = np.flatnonzero(np.r_[True, ids[1:] != ids[:-1]])
        stops = np.r_[starts[1:], len(ids)]
        self.table = t
        self.index = {ids[a]: (int(a), int(b)) for a, b in zip(starts, stops)}
        self.meta = json.loads((Path(directory) / "fate_association_meta.json").read_text())

    def rows(self, gene_id: str) -> pd.DataFrame | None:
        if not self.available:
            return None
        span = self.index.get(gene_id)
        return None if span is None else self.table.iloc[span[0]:span[1]]

    @property
    def by_gene(self):      # backwards compatible: membership tests / len()
        return self.index

    def verdict(self, gene_id: str, module: str | None) -> tuple[str, bool]:
        d = self.rows(gene_id)
        if d is None:
            return ("No fate-association result available", False)
        n_pass = int(d.passes.sum())
        if n_pass:
            best = d[d.passes].loc[lambda x: x.log2_effect.abs().idxmax()]
            way = "higher" if best.log2_effect > 0 else "lower"
            where = f"{SOURCES[best.source].lower()} {best.time_h} HPI"
            return (f"fate-associated in {n_pass} context{'s' if n_pass > 1 else ''}; "
                    f"strongest: {way} towards nodule fate in {where}", True)
        if d.candidate.any():
            return ("fate-associated, but never changed ≥1.5-fold vs matched mock in the same "
                    "direction", False)
        return "below the fate-association thresholds in every context", False

    def figure(self, gene_id: str, label: str | None = None) -> go.Figure | None:
        d = self.rows(gene_id)
        if d is None:
            return None
        d = d.set_index(["source", "time_h"])
        fig = make_subplots(rows=2, cols=2, shared_yaxes=True, shared_xaxes=True,
                            subplot_titles=list(SOURCES.values()),
                            horizontal_spacing=0.06, vertical_spacing=0.14)
        # scale to estimable fits; excluded (runaway) fits are clipped at the axis edge
        # pooled bars set the scale; single-platform markers beyond it are pinned to the edge
        vals = d.loc[d.status != "excluded_fit", ["log2_effect"]].abs().to_numpy()
        lim = max(1.2, float(np.nanmax(vals)) if np.isfinite(vals).any() else 1.0) * 1.2
        # fixed numeric slots so every panel shows 12/24/48/72 even when contexts are missing
        pos = {t: i for i, t in enumerate(TIMES)}
        x = [pos[t] for t in TIMES]

        for i, (src, name) in enumerate(SOURCES.items()):
            r, c = divmod(i, 2)
            r, c = r + 1, c + 1
            ys, colors, opac, hover, missing, hatch = [], [], [], [], [], []
            for t in TIMES:
                if (src, t) not in d.index:
                    ys.append(0); colors.append(BELOW); opac.append(0); missing.append(pos[t])
                    hatch.append(""); hover.append(f"{name} {t} HPI<br>{NOT_ELIGIBLE}")
                    continue
                row = d.loc[(src, t)]
                e = float(row.log2_effect)
                ys.append(e)
                hatch.append("")
                if row.status == "figure7":
                    colors.append(HIGHER if e > 0 else LOWER); opac.append(1.0); status = "Fate-associated, with ≥1.5-fold change vs matched mock in the same direction"
                elif row.status == "fate_associated":
                    colors.append(HIGHER if e > 0 else LOWER); opac.append(0.35)
                    status = "fate-associated, but not changed ≥1.5-fold vs mock in the same direction"
                elif row.status == "excluded_fit":
                    colors.append(HIGHER if e > 0 else LOWER); opac.append(0.8); hatch[-1] = "/"
                    if abs(e) > lim * 0.97:
                        ys[-1] = float(np.sign(e)) * lim * 0.97
                    status = f"excluded by the screen: {row.exclusion_reason}; effect not estimable"
                else:
                    colors.append(BELOW); opac.append(1.0)
                    fails = []
                    if abs(e) < SCREEN_LOG2:
                        fails.append("|change| &lt; log2 1.5")
                    if float(row.noise_ratio) < 3:
                        fails.append("noise ratio &lt; 3")
                    status = "below screen thresholds" + (f" ({', '.join(fails)})" if fails else "")
                plat = []
                for p, lab in (("BD", "BD"), ("10xP", "10x")):
                    v = row[f"log2_effect_{p}"]
                    if np.isfinite(v):
                        plat.append(f"{lab} {v:+.2f}")
                hover.append(f"{name} {t} HPI<br>log2 change {e:+.2f} ({', '.join(plat) or 'n/a'})"
                             f"<br>noise ratio {float(row.noise_ratio):.1f}; detected in "
                             f"{int(row.detected_cells):,} cells ({100 * float(row.detection):.1f}% weighted)"
                             f"<br>{status}"
                             f"<br>platforms assessable vs mock: {row.platforms or 'none'}")
            fig.add_trace(go.Bar(
                x=x, y=ys,
                marker=dict(color=[("white" if h else c) for c, h in zip(colors, hatch)],
                            opacity=opac,
                            line=dict(color=[(c if h else "rgba(0,0,0,0)") for c, h in zip(colors, hatch)],
                                      width=[1.5 if h else 0 for h in hatch]),
                            pattern=dict(shape=hatch, fgcolor=colors, solidity=0.35)),
                customdata=hover, hovertemplate="%{customdata}<extra></extra>", showlegend=False,
            ), row=r, col=c)
            edge = lim * 0.97
            for p, symbol in (("BD", "circle"), ("10xP", "square")):
                px, py, sym, txt = [], [], [], []
                for t in TIMES:
                    if (src, t) in d.index and d.loc[(src, t), "status"] != "excluded_fit":
                        v = float(d.loc[(src, t), f"log2_effect_{p}"])
                        if np.isfinite(v):
                            px.append(pos[t]); txt.append(f"{v:+.2f}")
                            if abs(v) > edge:
                                py.append(float(np.sign(v)) * edge)
                                sym.append("triangle-up" if v > 0 else "triangle-down")
                            else:
                                py.append(v); sym.append(symbol)
                fig.add_trace(go.Scatter(
                    x=px, y=py, mode="markers", name="BD" if p == "BD" else "10x protoplast",
                    marker=dict(symbol=sym, size=6, color="white", line=dict(color="#222", width=1.2)),
                    customdata=txt,
                    hovertemplate=f"{'BD' if p == 'BD' else '10x'} %{{customdata}}<extra></extra>",
                    showlegend=(i == 0), legendgroup=p,
                ), row=r, col=c)
            for y in (-SCREEN_LOG2, SCREEN_LOG2):
                fig.add_hline(y=y, line=dict(color="#999", dash="dot", width=1), row=r, col=c)
            fig.add_hline(y=0, line=dict(color="#333", width=1), row=r, col=c)
            for t in missing:
                fig.add_annotation(x=t, y=0, text="–", showarrow=False, yshift=8,
                                   hovertext=NOT_ELIGIBLE,
                                   font=dict(color="#999", size=12), row=r, col=c)

        fig.update_yaxes(range=[-lim, lim], zeroline=False)
        fig.update_yaxes(title_text="log2 change", row=1, col=1)
        fig.update_yaxes(title_text="log2 change", row=2, col=1)
        fig.update_xaxes(type="linear", range=[-0.6, len(TIMES) - 0.4], tickmode="array",
                         tickvals=x, ticktext=[str(t) for t in TIMES], showgrid=False,
                         fixedrange=True)
        fig.update_xaxes(title_text="HPI", row=2, col=1)
        fig.update_xaxes(title_text="HPI", row=2, col=2)
        fig.update_layout(template="plotly_white", height=470, bargap=0.25,
                          title=dict(text=escape(label or gene_id), x=0, xref="container",
                                     xanchor="left", y=1, yref="container", yanchor="top",
                                     pad=dict(t=8, l=10), font=dict(size=16)),
                          margin=dict(l=50, r=10, t=68, b=10), font=dict(size=11),
                          legend=dict(orientation="h", x=0, y=-0.16, font=dict(size=10)))
        for a in fig.layout.annotations:
            if a.text in SOURCES.values():
                a.font = dict(size=11, color="#333")
        return fig
