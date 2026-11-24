"""Figure builder
"""

from __future__ import annotations

import re
from dataclasses import dataclass, replace

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from data_store import ViewerData

# Positive expression: peach -> red; undetected cells use ZERO_COLOR.
EXPRESSION_SCALE = [
    [0.00, "#fbe4d8"],
    [0.25, "#f2c9b8"],
    [0.55, "#e8795c"],
    [1.00, "#c40000"],
]

# Colors for category variables
CATEGORY_COLORS = [
    "#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00",
    "#F0E442", "#7F3C8D", "#11A579", "#3969AC", "#F2B701", "#E73F74",
    "#80BA5A", "#E68310", "#008695", "#CF1C90", "#A5AA99", "#764E9F",
    "#ED645A", "#1D6996", "#73AF48", "#EDAD08", "#94346E", "#6F4070",
]

# Make not analyeses groups grey
MUTED_LABELS = {"not annotated", "Unassigned", "Unknown", "unassigned", "unknown"}
MISSING_LABEL = "not annotated"
MISSING_COLOR = "#cfcfcf"
CONTEXT_COLOR = "#e4e4e4"
ZERO_COLOR = "#dcdcdc"
# Opening view = the reference display (R project_view: az 120, el 45, view_mirror TRUE).
# The mirror is a reflection no camera can make, so data_store flips UMAP 1 for display
# (DISPLAY_MIRROR_UMAP1) and the camera looks from the matching direction:
# eye ∝ (sin az·cos el, −cos az·cos el, sin el), up = +UMAP 3, orthographic like the figure.
# Checked against the R projection on all cells: identical up to float error.
VIEW_AZ, VIEW_EL, VIEW_DISTANCE = 120.0, 45.0, 1.6
_az, _el = np.radians(VIEW_AZ), np.radians(VIEW_EL)
CAMERA = dict(
    eye=dict(x=float(VIEW_DISTANCE * np.sin(_az) * np.cos(_el)),
             y=float(-VIEW_DISTANCE * np.cos(_az) * np.cos(_el)),
             z=float(VIEW_DISTANCE * np.sin(_el))),
    up=dict(x=0, y=0, z=1),
    center=dict(x=0, y=0, z=0),
    projection=dict(type="orthographic"),
)


# Set legend over plot
COLORBAR = dict(
    title="log1p<br>CP10K",
    thickness=13,
    len=0.55,
    x=0.99,
    xanchor="right",
    y=0.5,
    bgcolor="rgba(255,255,255,0.78)",
)

OVERLAY_LEGEND = dict(
    itemsizing="constant",
    font=dict(size=11),
    x=0.99,
    xanchor="right",
    y=0.99,
    yanchor="top",
    bgcolor="rgba(255,255,255,0.78)",
    bordercolor="rgba(0,0,0,0.12)",
    borderwidth=1,
)


@dataclass
class PlotSettings:
    mode: str = "annotation"          
    annotation: str = ""
    gene_row: int | None = None
    facet: str = ""
    facet_context: bool = True     
    zero_display: str = "dim"      
    max_points: int = 150_000
    point_size: float = 1.6
    opacity: float = 0.75
    expr_quantile: float = 0.95
    # gene-set mode: precomputed per-cell score and its colour-bar title
    values: np.ndarray | None = None
    colorbar_title: str = "log1p<br>CP10K"
    # sample filter: boolean mask over all cells (None = all cells)
    cell_mask: np.ndarray | None = None


def scene_axes(coords: np.ndarray) -> dict:
    """Make facets share view"""
    lo, hi = coords.min(axis=0), coords.max(axis=0)
    pad = (hi - lo) * 0.03

    def axis(i: int, title: str) -> dict:
        return dict(
            title=title,
            range=[float(lo[i] - pad[i]), float(hi[i] + pad[i])],
            showspikes=False,
            # UMAP 1 is sign-flipped for display (reference-view mirror), so its tick values
            # would be negated; UMAP units are arbitrary, hide ticks on all three axes.
            showticklabels=False,
            gridcolor="#e8e8e8",
            zerolinecolor="#d0d0d0",
        )

    return dict(
        xaxis=axis(0, "UMAP 1"),
        yaxis=axis(1, "UMAP 2"),
        zaxis=axis(2, "UMAP 3"),
        # A shared explicit ratio also gives orthographic zoom a stable reset value.
        aspectmode="manual",
        aspectratio=dict(zip(("x", "y", "z"), (hi - lo) / np.cbrt(np.prod(hi - lo)))),
        camera=CAMERA,
        uirevision="keep-camera",
    )


# Identity palette (col_identity in the figure code). Keys are the
# label with any " (PO:...)" suffix removed and spaces -> underscores.
IDENTITY_COLORS = {
    "root_hair_cell": "#e87ba4",
    "cortex": "#008300",
    "pericycle_cell": "#4a3aa7",
    "procambium": "#eda100",
    "IC_like": "#008300",
    "OC_like": "#8fc98f",
    "root_meristem": "#e34948",
    "lateral_root_cap": "#9c8f7a",
    "endodermis": "#1baf7a",
    "xylem": "#8c510a",
    "companion_cell": "#2f8fbf",
    "sieve_element": "#b3246b",
    "nodule_cortex": "#0b4d1f",
    "nodule_pericycle_cell": "#1f1650",
    "Unassigned": "#cfcfcc",
    # viewer spelling of the cortex layer calls
    "cortex_IC": "#008300",
    "cortex_OC": "#8fc98f",
}
_PO_SUFFIX = re.compile(r"\s*\(PO:\d+\)\s*$")

# State layers: highlighted state in the figure colour, the rest is
# light-grey background, drawn first and fainter so the highlighted cells stay visible.
BACKGROUND_COLOR = "#d3d3d3"  # lightgrey
COLUMN_COLORS = {
    "ch5_nodule_like": {"True": "#0b4d1f", "False": BACKGROUND_COLOR},
    "ch5_responder": {
        "responding": "#762a83",
        "below_threshold": BACKGROUND_COLOR,
        "unassessable": BACKGROUND_COLOR,
    },
}
BACKGROUND_OPACITY = 0.25


def identity_key(category: str) -> str:
    return _PO_SUFFIX.sub("", category).strip().replace(" ", "_")


def is_background(column: str, category: str) -> bool:
    return COLUMN_COLORS.get(column, {}).get(category) == BACKGROUND_COLOR


def category_color(position: int, category: str, column: str = "") -> str:
    fixed = COLUMN_COLORS.get(column, {}).get(category) or IDENTITY_COLORS.get(identity_key(category))
    if fixed:
        return fixed
    if category in MUTED_LABELS:
        return MISSING_COLOR
    return CATEGORY_COLORS[position % len(CATEGORY_COLORS)]


def facet_axes(axes: dict) -> dict:
    stripped = dict(axes)
    for key in ("xaxis", "yaxis", "zaxis"):
        if key not in axes:
            continue
        stripped[key] = {
            **axes[key],
            "title": "",
            "showticklabels": False,
            "showbackground": False,
        }
    return stripped


def _context_trace(data, settings, n_facets=1):
    """All background cells in grey
    """
    budget = min(25_000, max(5_000, 120_000 // max(1, n_facets)))
    index = data.subsample(settings.cell_mask, min(settings.max_points, budget))
    coords = data.coords[index]
    return go.Scatter3d(
        x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
        mode="markers",
        marker=dict(size=max(0.5, settings.point_size * 0.6), color=CONTEXT_COLOR,
                    opacity=0.35),
        hoverinfo="skip",
        showlegend=False,
        name="all cells",
    )

# This is stupid
def facet_columns(n_facets: int) -> int:
    """Grid width"""
    if n_facets <= 1:
        return 1
    if n_facets <= 4:
        return 2
    if n_facets <= 9:
        return 3
    return 4


def figure_height(settings: PlotSettings, n_facets: int) -> str:
    """Make plotsize scale with screensize
    """
    if not settings.facet or n_facets <= 1:
        return "min(calc(100vh - 190px), 1100px)"
    rows = int(np.ceil(n_facets / facet_columns(n_facets)))
    return f"max({rows * 42}vh, {rows * 360 + 190}px)"

# Set limits on height to prevent extreme stretching
def plot_min_height(settings: PlotSettings, n_facets: int) -> str:
    if not settings.facet or n_facets <= 1:
        return "460px"
    rows = int(np.ceil(n_facets / facet_columns(n_facets)))
    return f"{rows * 360 + 190}px"


def _message(text: str) -> go.Figure:
    return go.Figure(
        layout=dict(
            annotations=[dict(text=text, showarrow=False, font=dict(size=16, color="#666"))],
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            template="plotly_white",
        )
    )


def _add_annotation_traces(
    figure, data, settings, index, row=None, col=None, legend=True, seen=None
) -> None:
    column = settings.annotation
    coords = data.coords[index]

    if column.startswith("num__"):
        values = data.cells[column].iloc[index].to_numpy(np.float32)
        scale_values = data.cells[column].to_numpy(np.float32)
        if settings.cell_mask is not None:
            scale_values = scale_values[settings.cell_mask]
        scale_values = scale_values[np.isfinite(scale_values)]
        lo, hi = (float(scale_values.min()), float(scale_values.max())) if scale_values.size else (0., 1.)
        figure.add_trace(
            go.Scatter3d(
                x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
                mode="markers",
                marker=dict(
                    size=settings.point_size,
                    color=values,
                    colorscale="Viridis",
                    cmin=lo, cmax=hi,
                    opacity=settings.opacity,
                    showscale=legend,
                    colorbar={**COLORBAR, "title": data.label_for(column[5:])},
                ),
                hoverinfo="skip",
                showlegend=False,
            ),
            row=row, col=col,
        )
        return

    series = data.cells[column].iloc[index]
    categories = list(enumerate(data.cells[column].cat.categories))
    # background states first so the highlighted state is drawn on top
    categories.sort(key=lambda pc: not is_background(column, str(pc[1])))
    for position, category in categories:
        selected = (series == category).to_numpy()
        if not selected.any():
            continue
        points = coords[selected]
        name = str(category)
        show = legend and (seen is None or name not in seen)
        if seen is not None:
            seen.add(name)
        background = is_background(column, name)
        figure.add_trace(
            go.Scatter3d(
                x=points[:, 0], y=points[:, 1], z=points[:, 2],
                mode="markers",
                name=name,
                legendgroup=name,
                showlegend=show,
                marker=dict(
                    size=settings.point_size * (0.8 if background else 1.0),
                    color=category_color(position, name, column),
                    opacity=min(settings.opacity, BACKGROUND_OPACITY) if background else settings.opacity,
                ),
                hovertemplate=f"{category}<extra></extra>",
            ),
            row=row, col=col,
        )


def _add_gene_trace(
    figure, data, settings, values, index, cmax, row=None, col=None, legend=True
) -> None:
    """Draw gene expression (or a gene-set score: cells <= 0 form the grey backdrop)
    """
    colorbar = {**COLORBAR, "title": settings.colorbar_title}
    coords = data.coords[index]
    subset = values[index]
    expressing = subset > 0

    # --- non-expressing backdrop ---------------------------------------
    if settings.zero_display != "hide" and (~expressing).any():
        flat = coords[~expressing]
        dim = settings.zero_display == "dim"
        figure.add_trace(
            go.Scatter3d(
                x=flat[:, 0], y=flat[:, 1], z=flat[:, 2],
                mode="markers",
                marker=dict(
                    size=settings.point_size * (0.75 if dim else 1.0),
                    color=ZERO_COLOR,
                    opacity=settings.opacity * (0.45 if dim else 1.0),
                ),
                hoverinfo="skip",
                showlegend=False,
                name="not detected",
            ),
            row=row, col=col,
        )

    # plot low to high to prevent hiding expression in point clouds
    if expressing.any():
        lit = coords[expressing]
        strength = subset[expressing]
        order = np.argsort(strength, kind="stable")
        lit, strength = lit[order], strength[order]

        figure.add_trace(
            go.Scatter3d(
                x=lit[:, 0], y=lit[:, 1], z=lit[:, 2],
                mode="markers",
                marker=dict(
                    size=settings.point_size * 1.7,
                    color=strength,
                    colorscale=EXPRESSION_SCALE,
                    cmin=0.0,
                    cmax=cmax,
                    opacity=settings.opacity,
                    showscale=legend,
                    colorbar=colorbar,
                ),
                hoverinfo="skip",
                showlegend=False,
                name="detected",
            ),
            row=row, col=col,
        )
    elif legend:
        figure.add_trace(
            go.Scatter3d(
                x=[], y=[], z=[], mode="markers",
                marker=dict(
                    color=[], colorscale=EXPRESSION_SCALE, cmin=0.0, cmax=cmax,
                    showscale=True,
                    colorbar=colorbar,
                ),
                showlegend=False, hoverinfo="skip",
            ),
            row=row, col=col,
        )


# Nodule-fate probability (OT): light grey-green -> PRGn green
FATE_SCALE = [[0.0, "#e6efe7"], [0.35, "#a6d1ad"], [0.7, "#4f9e5f"], [1.0, "#1b7837"]]
OUTSIDE_COLOR = "#e8e8e8"


def _add_fate_trace(figure, data, settings, values, index, cmax, row=None, col=None, legend=True):
    """Cells outside the fate model (NaN) form a faint backdrop; the rest are coloured by
    P(nodule fate), low to high so high-probability cells are drawn on top."""
    coords = data.coords[index]
    subset = values[index]
    inside = np.isfinite(subset)
    if (~inside).any() and settings.zero_display != "hide":
        flat = coords[~inside]
        figure.add_trace(go.Scatter3d(
            x=flat[:, 0], y=flat[:, 1], z=flat[:, 2], mode="markers",
            marker=dict(size=settings.point_size * 0.75, color=OUTSIDE_COLOR,
                        opacity=settings.opacity * 0.4),
            hoverinfo="skip", showlegend=False, name="not in the fate model",
        ), row=row, col=col)
    lit, strength = coords[inside], subset[inside]
    order = np.argsort(strength, kind="stable")
    figure.add_trace(go.Scatter3d(
        x=lit[order, 0], y=lit[order, 1], z=lit[order, 2], mode="markers",
        marker=dict(size=settings.point_size * 1.3, color=strength[order], colorscale=FATE_SCALE,
                    cmin=0.0, cmax=cmax, opacity=settings.opacity, showscale=legend,
                    colorbar={**COLORBAR, "title": settings.colorbar_title}),
        hoverinfo="skip", showlegend=False, name="P(nodule fate)",
    ), row=row, col=col)


def colour_max(values: np.ndarray, quantile: float) -> float:
    """Upper end of the colour scale from the expressing cells only
    """
    nonzero = values[values > 0]
    if nonzero.size == 0:
        return 1.0
    return max(float(np.quantile(nonzero, quantile)), 1e-3)


def build_figure(data: ViewerData, settings: PlotSettings, axes: dict | None = None):
    axes = axes if axes is not None else scene_axes(data.coords)

    values = cmax = None
    fate = settings.mode == "fate"
    if fate:
        values = settings.values
        finite = values[np.isfinite(values)]
        cmax = max(float(np.quantile(finite, settings.expr_quantile)), 1e-3) if finite.size else 1.0
        settings = replace(settings, mode="gene")
    elif settings.mode == "geneset":
        if settings.values is None:
            values = np.zeros(data.n_cells, dtype=np.float32)
            cmax = 1.0
            settings = replace(settings, zero_display="dim")
        else:
            values = settings.values
            cmax = colour_max(values, settings.expr_quantile)
        settings = replace(settings, mode="gene")   # same drawing as expression
    elif settings.mode == "gene":
        if settings.gene_row is None:
            values = np.zeros(data.n_cells, dtype=np.float32)
            cmax = 1.0
            settings = replace(settings, zero_display="dim")
        else:
            values = data.expression(settings.gene_row)
            cmax = colour_max(values, settings.expr_quantile)
    elif not settings.annotation:
        return _message("Choose an annotation level.")

    draw = _add_fate_trace if fate else _add_gene_trace
    if not settings.facet:
        index = data.subsample(settings.cell_mask, settings.max_points)
        figure = go.Figure()
        if settings.mode == "gene":
            draw(figure, data, settings, values, index, cmax)
        else:
            _add_annotation_traces(figure, data, settings, index)
        figure.update_layout(
            scene=axes,
            autosize=True,
            margin=dict(l=0, r=0, t=0, b=0),
            legend=OVERLAY_LEGEND,
            uirevision="keep-camera",
            template="plotly_white",
        )
        return figure

    # Facet parameters
    column = data.cells[settings.facet]
    present = (column[settings.cell_mask] if settings.cell_mask is not None else column).value_counts()
    facets = [c for c in column.cat.categories if present.get(c, 0) > 0]
    if not facets:
        return _message("Facet column has no populated levels.")

    n_cols = facet_columns(len(facets))
    n_rows = int(np.ceil(len(facets) / n_cols))

    figure = make_subplots(
        rows=n_rows,
        cols=n_cols,
        specs=[[{"type": "scene"} for _ in range(n_cols)] for _ in range(n_rows)],
        subplot_titles=[f"{c}  (n={present[c]:,})" for c in facets],
        horizontal_spacing=0.02,
        vertical_spacing=0.06,
    )

    # Split the point budget across panels so the total drawn stays constant
    per_panel = max(2_000, settings.max_points // len(facets))
    drawn_legend = False
    seen: set[str] = set()

    for position, facet in enumerate(facets):
        row, col = divmod(position, n_cols)
        mask = (column == facet).to_numpy()
        if settings.cell_mask is not None:
            mask = mask & settings.cell_mask
        panel_index = data.subsample(mask, per_panel)
        if panel_index.size == 0:
            continue
        if settings.facet_context and settings.mode != "gene":
            figure.add_trace(
                _context_trace(data, settings, len(facets)), row=row + 1, col=col + 1
            )
        legend = not drawn_legend
        if settings.mode == "gene":
            draw(
                figure, data, settings, values, panel_index, cmax,
                row + 1, col + 1, legend=legend,
            )
        else:
            _add_annotation_traces(
                figure, data, settings, panel_index, row + 1, col + 1,
                legend=legend if settings.annotation.startswith("num__") else True, seen=seen,
            )
        # A first facet with no detected/modelled cells has an empty colour
        # trace; Plotly does not render its colourbar. Use the next populated one.
        drawn_legend = any(t.marker.showscale and t.x is not None and len(t.x)
                           for t in figure.data)

    panel_axes = facet_axes(axes)
    for position in range(len(facets)):
        key = "scene" if position == 0 else f"scene{position + 1}"
        figure.update_layout(**{key: panel_axes})

    continuous = any(t.marker.showscale for t in figure.data)
    if continuous:
        for trace in figure.data:
            if trace.marker.showscale:
                trace.marker.colorbar.update(x=1.02, xanchor="left", y=.5, len=.75)
    figure.update_layout(
        autosize=True,
        margin=dict(l=0, r=110 if continuous else 8, t=30, b=30 if continuous else 190),
        legend=dict(orientation="h", x=0, xanchor="left", y=-.06, yanchor="top",
                    entrywidth=240, entrywidthmode="pixels", maxheight=150, itemsizing="constant", font=dict(size=11)),
        uirevision="keep-camera",
        template="plotly_white",
    )
    figure.update_annotations(font_size=12)
    return figure
