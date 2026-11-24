"""
shiny run app/app.py

Data directory: .../viewer_data

Layout: persistent gene search and active selection; Cells / Genes / Expression / Fate;
one contextual right panel; detailed evidence in non-destructive wide overlays.
Gene inspection is independent of the active gene set.
"""

from __future__ import annotations

import json
import os
import sys
import tempfile
from pathlib import Path

# `shiny run` puts this directory on sys.path; a bare ASGI server (uvicorn,
# gunicorn) does not. Doing it here means both invocations work.
sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
import pandas as pd
import plotly.io as pio
from plotly.offline import get_plotlyjs
from shiny import App, reactive, render, ui

from data_store import ViewerData
from dissociation import Dissociation
from expression_view import DotPlot
from fate_assoc import FateAssociation
from fate_view import FateCells
from gene_heatmap import MAX_GENES as HEAT_MAX, GeneHeatmap
from gene_map import GeneMap
from gene_sets import MAX_GENES, GeneSetScorer, fate_heatmap
from plotting import (CAMERA, PlotSettings, build_figure, figure_height, plot_min_height, scene_axes)

DATA_DIR = Path(os.environ.get("VIEWER_DATA", Path(__file__).parent.parent / "viewer_data"))
DATA = ViewerData(DATA_DIR)
AXES = scene_axes(DATA.coords)
GENE_MAP = GeneMap(DATA_DIR)
DISSOC = Dissociation(DATA_DIR)
FATE = FateAssociation(DATA_DIR)
FATE_CELLS = FateCells(DATA_DIR)
DOTPLOT = DotPlot(DATA, DATA_DIR)
HEATMAP = GeneHeatmap(DATA, DOTPLOT.condition) if DOTPLOT.available else None
SCORER = GeneSetScorer(DATA, DATA_DIR)

# Precomputed sets (build/build_module_scores.py): fate-associated gene modules + cNMF GEP top 200
MODULE_PRESETS: dict[str, tuple[str, list[str], str]] = {}
_meta_path = DATA_DIR / "module_scores_meta.json"
if _meta_path.exists():
    _npz = np.load(DATA_DIR / "module_scores.npz")
    for _key, _info in json.loads(_meta_path.read_text())["sets"].items():
        _ids = list(DATA.genes.gene_id.iloc[_npz[f"{_key}__rows"]])
        MODULE_PRESETS[_key] = (_info["label"], _ids, _info["group"])
PRESET_CHOICES = {"": "(choose a set)"}
for _group, _title in (("figure7", "gene modules"), ("cnmf", "cNMF programmes (GEP top 200)")):
    _opts = {k: f"{v[0]} ({len(v[1]):,} genes)" for k, v in MODULE_PRESETS.items() if v[2] == _group}
    if _opts:
        PRESET_CHOICES[_title] = _opts

# Samples (build/build_cell_design.py): grouped by condition for the sample filter
SAMPLES = json.loads((DATA_DIR / "samples.json").read_text()) if (DATA_DIR / "samples.json").exists() else []
SAMPLE_IDS = [s["sample"] for s in SAMPLES]
# one compact checkbox group per condition (tip region kept on the sample line)
SAMPLE_GROUPS: list[tuple[str, str, dict[str, str]]] = []
for _s in SAMPLES:
    _base = _s["condition"].split(" 0")[0].strip()
    _region = _s["condition"][len(_base):].strip()
    _short = f"{_s['sample']} · {_s['platform']}{' · ' + _region if _region else ''} ({_s['n_cells']:,})"
    if not SAMPLE_GROUPS or SAMPLE_GROUPS[-1][1] != _base:
        SAMPLE_GROUPS.append((f"samples_{len(SAMPLE_GROUPS)}", _base, {}))
    SAMPLE_GROUPS[-1][2][_s["sample"]] = _short
SAMPLE_COLUMN = DATA.cells["sample"].astype(str).to_numpy()
TIP_SAMPLES = [s["sample"] for s in SAMPLES if s["condition"].startswith("Root tip")]
SZ_SAMPLES = [s for s in SAMPLE_IDS if s not in TIP_SAMPLES]

ROW_CHOICES = {c: DATA.label_for(c) for c in DATA.categorical if c != "sample"}
GENE_MAP_JSON = {}
if GENE_MAP.available:
    GENE_MAP_JSON["2"] = json.loads(pio.to_json(GENE_MAP.figure(dim=2)))
    if GENE_MAP.has_3d:
        GENE_MAP_JSON["3"] = json.loads(pio.to_json(GENE_MAP.figure(dim=3)))


def _prepare_assets() -> Path:
    """Puts plotly.js somewhere"""
    candidates = []
    if os.environ.get("VIEWER_ASSETS"):
        candidates.append(Path(os.environ["VIEWER_ASSETS"]))
    candidates.append(Path(__file__).parent / "_assets")
    candidates.append(Path(tempfile.gettempdir()) / "mt_viewer_assets")

    for target in candidates:
        try:
            target.mkdir(parents=True, exist_ok=True)
            script = target / "plotly.min.js"
            if not script.exists() or script.stat().st_size == 0:
                script.write_text(get_plotlyjs(), encoding="utf-8")
            for name in ("browser.css", "browser.js"):
                source = Path(__file__).parent / "_assets" / name
                destination = target / name
                if source.resolve() != destination.resolve():
                    destination.write_bytes(source.read_bytes())
            return target
        except OSError:
            continue
    raise RuntimeError(
        "Couldn't write plotly.min.js")


ASSETS = _prepare_assets()

PLOT_CONFIG = {
    "responsive": True,
    "displaylogo": False,
    "toImageButtonOptions": {"format": "png", "scale": 2, "filename": "mt_umap3d"},
}

## Drawing all cells is sometimes a bit slow
POINT_BUDGETS = {
    "50000": "50,000 (fastest)",
    "100000": "100,000 (default)",
    "200000": "200,000",
    str(DATA.n_cells): f"All {DATA.n_cells:,} cells",
}

NOTE = "color:#777;font-size:0.78em;line-height:1.45;margin-top:6px"


def field(label: str, *controls, help: str | None = None):
    """One toolbar slot: small label on top, control below (fixed height so rows line up)."""
    return ui.div(
        ui.div(label, class_="tb-label"),
        ui.div(*controls, class_="tb-control"),
        ui.div(help, class_="tb-help") if help else None,
        class_="tb-field",
    )


def toolbar(*items, cls: str = "toolbar"):
    return ui.div(*items, class_=cls)


def caption(html: str) -> ui.HTML:
    return ui.HTML(f"<div style='padding:2px 2px 6px;line-height:1.35;color:#555;font-size:0.9em'>{html}</div>")


from ui_layout import build_ui, COPY

app_ui = build_ui(globals(),
    #// Fixing camera position was horrible, mostly copied claude/chatgpt suggestion fixes for this
    ui.head_content(
        ui.tags.script(src="assets/plotly.min.js"),
        ui.tags.script(
            r"""
            document.addEventListener('DOMContentLoaded', function () {
              // Orthographic zoom is stored in aspectratio, not camera.eye.
              var lastViews = {}, baseViews = {}, applyingView = false, activeScene = 'scene';
              function sceneKeys(obj) {
                return Object.keys(obj || {}).filter(k => /^scene[0-9]*$/.test(k));
              }
              function clone(value) { return JSON.parse(JSON.stringify(value)); }
              function currentView(el, key) {
                var layout = el._fullLayout && el._fullLayout[key];
                if (!layout) return null;
                var scene = layout._scene;
                var camera = scene && scene.getCamera ? scene.getCamera() : layout.camera;
                return clone({camera: camera, aspectratio: layout.aspectratio});
              }
              function rememberViews(el) {
                sceneKeys(el._fullLayout).forEach(k => { lastViews[k] = currentView(el, k); });
              }
              function syncEnabled() {
                var box = document.getElementById('sync_facets');
                return !!box && box.checked;
              }
              function movedScene(ev) {
                for (var key of Object.keys(ev || {})) {
                  var match = /^(scene[0-9]*)\.(camera|aspectratio|aspectmode)(\.|$)/.exec(key);
                  if (match) return match[1];
                }
                return null;
              }
              function viewUpdate(update, key, view) {
                update[key + '.camera'] = clone(view.camera);
                update[key + '.aspectratio'] = clone(view.aspectratio);
                update[key + '.aspectmode'] = 'manual';
              }
              function applyViews(el, update) {
                applyingView = true;
                return Plotly.relayout(el, update).finally(function () {
                  applyingView = false;
                  rememberViews(el);
                });
              }
              function alignPanels(el, source) {
                if (applyingView) return;
                var keys = sceneKeys(el._fullLayout);
                if (keys.length < 2) return;
                source = keys.includes(source) ? source : keys[0];
                var view = currentView(el, source), update = {};
                // Leave the scene being manipulated alone; update its peers.
                keys.filter(k => k !== source).forEach(k => viewUpdate(update, k, view));
                return applyViews(el, update);
              }
              Shiny.addCustomMessageHandler('umap-update', function (msg) {
                var wrap = document.getElementById('umap3d_wrap');
                var el = document.getElementById('umap3d');
                if (!el) return;
                if (msg.height) {
                  wrap.style.height = msg.height;
                  wrap.style.minHeight = msg.minHeight;
                }
                rememberViews(el);
                var fallback = lastViews[activeScene] || lastViews.scene;
                sceneKeys(msg.figure.layout).forEach(function (key) {
                  var layout = msg.figure.layout[key];
                  baseViews[key] = clone({camera: layout.camera, aspectratio: layout.aspectratio});
                  var saved = lastViews[key] || fallback;
                  if (saved) {
                    layout.camera = clone(saved.camera);
                    layout.aspectratio = clone(saved.aspectratio);
                    layout.aspectmode = 'manual';
                  }
                });
                applyingView = true;
                Plotly.react(el, msg.figure.data, msg.figure.layout, msg.config).then(function () {
                  applyingView = false;
                  if (!el._camListener) {
                    el._camListener = true;
                    el.on('plotly_relayout', function (ev) {
                      if (applyingView) return;
                      var source = movedScene(ev);
                      if (!source) return; // resizing is not a camera gesture
                      activeScene = source;
                      // A wheel event can arrive before _fullLayout has its new
                      // orthographic ratio. Commit the event's values, including
                      // the source, so redrawing a peer cannot undo that zoom.
                      var view = currentView(el, source), update = {};
                      if (ev[source + '.camera']) view.camera = clone(ev[source + '.camera']);
                      if (ev[source + '.aspectratio']) view.aspectratio = clone(ev[source + '.aspectratio']);
                      var keys = syncEnabled() ? sceneKeys(el._fullLayout) : [source];
                      keys.forEach(k => viewUpdate(update, k, view));
                      applyViews(el, update);
                    });
                  }
                  if (msg.align) alignPanels(el, activeScene);
                  Plotly.Plots.resize(el);
                }, function () { applyingView = false; });
              });

              // ---- Gene UMAP tab ----
              var pendingHighlight = null;
              function applyHighlight(msg) {
                var gel = document.getElementById('geneumap');
                if (!gel || !gel.data || gel._drawing) { pendingHighlight = msg; return; }
                var last = gel.data.length - 1;
                if (msg.xs) {
                  var many = {x: [msg.xs], y: [msg.ys], text: [msg.labels]};
                  if (gel.data[last].type === 'scatter3d') { many.z = [msg.zs]; }
                  Plotly.restyle(gel, many, [last]);
                  return;
                }
                var none = msg.x === null || msg.x === undefined;
                var upd = {
                  x: [none ? [] : [msg.x]],
                  y: [none ? [] : [msg.y]],
                  text: [msg.label && !none ? [msg.label] : []]
                };
                if (gel.data[last].type === 'scatter3d') {
                  upd.z = [none || msg.z === null ? [] : [msg.z]];
                }
                Plotly.restyle(gel, upd, [last]);
              }
              Shiny.addCustomMessageHandler('geneumap-init', function (msg) {
                var gel = document.getElementById('geneumap');
                if (!gel) { return; }
                var cfg = Object.assign({}, msg.config, {scrollZoom: true});
                gel._drawing = true;
                Plotly.newPlot(gel, msg.figure.data, msg.figure.layout, cfg).then(function () {
                  gel._drawing = false;
                  if (gel.removeAllListeners) { gel.removeAllListeners('plotly_click'); }
                  gel.on('plotly_click', function (ev) {
                    var p = ev && ev.points && ev.points[0];
                    if (p && p.customdata) {
                      Shiny.setInputValue('geneumap_click', p.customdata[0], {priority: 'event'});
                    }
                  });
                  if (pendingHighlight) { applyHighlight(pendingHighlight); pendingHighlight = null; }
                });
              });
              Shiny.addCustomMessageHandler('geneumap-highlight', applyHighlight);
              // hidden tabs / collapsing sidebars change plot sizes: resize on both
              function resizeAll() {
                ['geneumap', 'umap3d', 'fate_gene', 'dissoc_gene', 'dissoc_heat', 'dotplot', 'fateplot'].forEach(function (id) {
                  var p = document.getElementById(id);
                  if (p && p.data) { Plotly.Plots.resize(p); }
                });
              }
              document.addEventListener('shown.bs.tab', resizeAll);
              document.addEventListener('shown.bs.collapse', resizeAll);
              document.addEventListener('bslib.sidebar', function () { setTimeout(resizeAll, 350); });

              // ---- Dissociation sidebar (bulk isolation effect) ----
              Shiny.addCustomMessageHandler('dissoc-plot', function (msg) {
                var el = document.getElementById(msg.target);
                if (!el) { return; }
                if (!msg.figure) { Plotly.purge(el); el.style.height = '0px'; return; }
                el.style.height = (msg.figure.layout.height || 300) + 'px';
                Plotly.react(el, msg.figure.data, msg.figure.layout, msg.config).then(function () {
                  if (msg.target === 'dotplot') {
                    // gene heatmap: a click on a row selects that gene (layout.meta.genes).
                    // Re-bound on every draw: Plotly.purge (empty plot) drops listeners.
                    if (el.removeAllListeners) { el.removeAllListeners('plotly_click'); }
                    el.on('plotly_click', function (ev) {
                      var meta = el.layout && el.layout.meta;
                      var p = ev && ev.points && ev.points[0];
                      if (!meta || !meta.genes || !p || typeof p.y !== 'number') { return; }
                      var gene = meta.genes[p.y];
                      if (gene) { Shiny.setInputValue('heat_gene_click', gene, {priority: 'event'}); }
                    });
                  }
                });
              });

              // ---- custom gene sets: kept in this visitor's browser only ----
              var CUSTOM_KEY = 'mt_browser_custom_sets_v1';
              $(document).on('shiny:connected', function () {
                var stored = '{}';
                try { stored = localStorage.getItem(CUSTOM_KEY) || '{}'; } catch (e) { /* storage blocked */ }
                Shiny.setInputValue('stored_custom_sets', stored);
              });
              Shiny.addCustomMessageHandler('store-custom-sets', function (msg) {
                try { localStorage.setItem(CUSTOM_KEY, JSON.stringify(msg.sets || {})); } catch (e) { Shiny.setInputValue('storage_failed', Date.now(), {priority:'event'}); }
              });

              Shiny.addCustomMessageHandler('umap-reset', function (msg) {
                var el = document.getElementById('umap3d');
                if (!el || !el._fullLayout) { return; }
                var update = {};
                sceneKeys(el._fullLayout).forEach(function (key) {
                  var view = clone(baseViews[key]);
                  view.camera = clone(msg.camera);
                  viewUpdate(update, key, view);
                });
                applyViews(el, update);
              });
            });
            """
        ),
    ),
)


def server(input, output, session):
    DISSOC_CONFIG = {"displaylogo": False, "responsive": True,
                     "modeBarButtonsToRemove": ["select2d", "lasso2d", "autoScale2d"]}

    def _fig_json(fig):
        if fig is not None:
            fig.update_layout(font=dict(color="#202020"), paper_bgcolor="#FFFFFF")
        return json.loads(pio.to_json(fig)) if fig is not None else None

    async def push_plot(target: str, fig, filename: str):
        await session.send_custom_message(
            "dissoc-plot", {"target": target, "figure": _fig_json(fig),
                            "config": {**DISSOC_CONFIG, "toImageButtonOptions": {
                                "format": "png", "scale": 2, "filename": filename}}})

    # ================= selection (one gene or one gene set) =================
    @reactive.calc
    def gene_hits() -> pd.DataFrame:
        return DATA.gene_search(input.gene_query(), limit=100)

    @output(suspend_when_hidden=False)
    @render.ui
    def gene_picker():
        hits = gene_hits()
        if hits.empty:
            return ui.div(
                ui.tags.em(COPY["no_search_matches"]), style="color:#b00;padding:4px 0;"
            )
        choices = {"": "Choose a gene"} if not input.gene_query().strip() else {}
        for row in hits.itertuples():
            detail = row.description[:55] if row.description else f"{row.pct_cells_detected:.0f}% of cells"
            choices[row.gene_id] = f"{row.label} - {detail}"
        with reactive.isolate():
            previous = input.gene() if "gene" in input else None
        selected = previous if input.gene_query().strip() and previous in choices else next(iter(choices))
        return ui.input_select("gene", "Search results", choices, selected=selected, size=min(6, len(choices)))

    @reactive.calc
    def gene_row() -> int | None:
        try:
            token = input.gene()
        except Exception:
            return None
        return DATA.resolve_gene(token)

    @reactive.calc
    def selected_gene_id() -> str | None:
        row = gene_row()
        if row is None:
            return None
        return str(DATA.genes.gene_id.iloc[row])

    def activate_searched_gene():
        ui.update_radio_buttons("sel_kind", selected="gene")
        ui.update_radio_buttons("details_kind", selected="gene")
        ui.update_radio_buttons("fate_target", selected="active")
        ui.update_select("colour_by", selected="selection")

    @reactive.effect
    @reactive.event(input.gene_search_action)
    def search_activated():
        activate_searched_gene()

    @reactive.calc
    def geneset_parsed() -> tuple[list[int], list[str]]:
        if not SCORER.available:
            return [], []
        report = SCORER.resolve_report(input.geneset_text())
        return report["rows"], report["missing"]

    @reactive.calc
    def current_set() -> list[tuple[str, str]] | None:
        """[(gene_id, label)] of the active gene set, or None when one gene is selected."""
        if input.sel_kind() != "set":
            return None
        rows, _ = geneset_parsed()
        return [(str(DATA.genes.gene_id.iloc[r]), str(DATA.genes.label.iloc[r])) for r in rows]

    @reactive.calc
    def geneset_score():
        rows, _ = geneset_parsed()
        if input.sel_kind() != "set" or not rows or not SCORER.available:
            return None
        return SCORER.score(rows)

    @reactive.calc
    def selection():
        """(per-cell values, kind 'gene'|'set', label) of the selection, or None."""
        if input.sel_kind() == "set":
            scored = geneset_score()
            if scored is None:
                return None
            n = len(geneset_parsed()[0])
            return scored[0], "set", f"Gene set ({n:,} genes)"
        row = gene_row()
        if row is None:
            return None
        return DATA.expression(row), "gene", str(DATA.genes.label.iloc[row])

    # An explicit gene or selection change switches the colouring to the active selection.
    @reactive.effect
    @reactive.event(input.gene, input.sel_kind, ignore_init=True)
    def colour_by_selection():
        if input.colour_by() == "annotation" and (gene_row() is not None or input.sel_kind() == "set"):
            ui.update_select("colour_by", selected="selection")

    # ---- saved / precomputed sets ----
    custom_sets = reactive.value({})      # name -> [gene ids], persisted in the visitor's browser
    active_name = reactive.value("")      # name of the set currently in the text box ("" = unsaved)

    def preset_key() -> str:
        key = input.geneset_preset()
        if isinstance(key, (list, tuple)):
            key = key[0] if key else ""
        return str(key or "")

    def set_genes(key: str) -> list[str] | None:
        if key in MODULE_PRESETS:
            return MODULE_PRESETS[key][1]
        if key.startswith("custom:"):
            return custom_sets.get().get(key[7:])
        return None

    def set_label(key: str) -> str:
        if key in MODULE_PRESETS:
            return MODULE_PRESETS[key][0]
        return key[7:] if key.startswith("custom:") else ""

    def preset_choices() -> dict:
        choices = {"": "(choose a set)"}
        sets = custom_sets.get()
        if sets:
            choices["Custom sets"] = {f"custom:{n}": f"{n} ({len(g):,} genes)" for n, g in sets.items()}
        choices.update({k: v for k, v in PRESET_CHOICES.items() if k != ""})
        return choices

    def store_custom(sets: dict, select: str | None):
        custom_sets.set(sets)
        ui.update_select("geneset_preset", choices=preset_choices(),
                         selected=select if select is not None else preset_key())
        return session.send_custom_message("store-custom-sets", {"sets": sets})

    def unique_name(base: str) -> str:
        names, name, i = custom_sets.get(), base, 1
        while name in names:
            name, i = f"{base} ({i})", i + 1
        return name

    @reactive.effect
    @reactive.event(input.stored_custom_sets)
    def restore_custom_sets():
        try:
            raw = json.loads(input.stored_custom_sets() or "{}")
        except (TypeError, ValueError):
            return
        clean = {}
        for name, ids in (raw.items() if isinstance(raw, dict) else []):
            if isinstance(name, str) and isinstance(ids, list):
                ids = [str(g) for g in ids if DATA.resolve_gene(str(g)) is not None][:MAX_GENES]
                if ids:
                    clean[name[:80]] = ids
        if clean:
            custom_sets.set({**clean, **custom_sets.get()})
            ui.update_select("geneset_preset", choices=preset_choices(), selected=preset_key())

    @reactive.effect
    @reactive.event(input.geneset_preset)
    def load_preset():
        key = preset_key()
        genes = set_genes(key)
        if genes:
            ui.update_text_area("geneset_text", value="\n".join(genes))
            active_name.set(set_label(key))
            activate_set()

    def activate_set():
        # A newly activated set starts without a leftover inspected gene.
        ui.update_text("gene_query", value="")
        ui.update_select("gene", choices={"": "Choose a gene"}, selected="")
        ui.update_radio_buttons("sel_kind", selected="set")
        ui.update_radio_buttons("details_kind", selected="set")
        ui.update_radio_buttons("fate_target", selected="active")
        ui.update_select("colour_by", selected="selection")

    @reactive.effect
    @reactive.event(input.geneset_text, ignore_init=True)
    def text_edited():
        """Editing the box by hand detaches it from the selected set, so choosing that set
        again in the dropdown is a real change (and reloads it)."""
        key = preset_key()
        genes = set_genes(key) if key else None
        typed = set(t.lower() for t in __import__("re").split(r"[\s,;]+", input.geneset_text() or "") if t)
        if genes is not None and typed != set(g.lower() for g in genes):
            ui.update_select("geneset_preset", selected="")
            active_name.set("")

    @reactive.effect
    @reactive.event(input.custom_save)
    async def save_custom():
        rows, _ = geneset_parsed()
        if not rows:
            ui.notification_show(COPY["set_save_empty"], type="warning")
            return
        typed = (input.custom_name() or "").strip()[:80]
        # a typed name is used as given (saving under an existing name replaces that set);
        # otherwise assign the next available default name
        name = typed or unique_name("Custom set")
        ids = [str(DATA.genes.gene_id.iloc[r]) for r in rows]
        await store_custom({**custom_sets.get(), name: ids}, f"custom:{name}")
        active_name.set(name)
        ui.update_text("custom_name", value="")
        ui.notification_show(COPY["set_saved"].format(name=name, n=len(ids)), duration=4)

    @render.ui
    def custom_delete():
        key = preset_key()
        if not key.startswith("custom:"):
            return ui.HTML("")
        return ui.input_action_link("custom_remove", f"Delete '{key[7:]}' from custom sets.")

    @reactive.effect
    @reactive.event(input.custom_remove)
    async def remove_custom():
        key = preset_key()
        if key.startswith("custom:"):
            sets = {n: g for n, g in custom_sets.get().items() if n != key[7:]}
            await store_custom(sets, "")
            active_name.set("")

    @render.ui
    def geneset_note():
        if not SCORER.available:
            return ui.HTML("<p style='color:#b00'>control_panels.npz missing: run "
                           "build/build_control_panels.py.</p>")
        report = SCORER.resolve_report(input.geneset_text())
        rows = report["rows"]
        lines = [ui.p(COPY["parse_summary"].format(recognised=len(rows), duplicates=len(report["duplicates"]), missing=len(report["missing"]), ambiguous=len(report["ambiguous"])))]
        for key, label in (("missing", "Unmatched genes"), ("duplicates", "Duplicates"), ("ambiguous", "Ambiguous: Use locus tags (e.g. MtrunA17_Chr5g0448621)")):
            if report[key]:
                lines.append(ui.tags.details(ui.tags.summary(label), ui.p(", ".join(report[key]))))
        if report["omitted"]:
            lines.append(ui.p(COPY["parse_omitted"].format(omitted=report["omitted"], maximum=MAX_GENES)))
        lines.append(ui.p(COPY["score_help"], class_="help-copy"))
        return ui.div(*lines)

    @reactive.effect
    @reactive.event(input.use_set)
    def use_set():
        if geneset_parsed()[0]:
            activate_set()
        else:
            ui.notification_show(COPY["set_use_empty"], type="warning")

    @reactive.effect
    @reactive.event(input.use_gene)
    def use_gene():
        if gene_row() is not None:
            activate_searched_gene()

    @render.ui
    def selection_summary():
        row = gene_row()
        inspected = str(DATA.genes.label.iloc[row]) if row is not None else "No gene selected"
        if input.sel_kind() == "set":
            return ui.div(ui.div(COPY["selection_set"].format(name=active_name.get() or "Unsaved set", n=len(geneset_parsed()[0])), class_="selection-name"),
                          ui.div(COPY["selection_inspected"].format(gene=inspected), class_="selection-meta") if row is not None else None)
        return ui.div(ui.div(inspected, class_="selection-name"), ui.div("Active gene", class_="selection-meta"))

    @render.ui
    def gene_annotation():
        if input.sel_kind() == "set" and input.details_kind() == "set":
            return ui.HTML("")
        row = gene_row()
        if row is None:
            return ui.HTML("")
        return ui.p(str(DATA.genes.description.iloc[row]) or COPY["no_annotation"], class_="help-copy")

    @render.ui
    def expression_heading():
        return ui.h3("Gene expression heatmap" if input.expr_view() == "heat" else
                     "Gene-set score dot plot" if input.sel_kind() == "set" else "Gene expression")

    @reactive.effect
    @reactive.event(input.storage_failed)
    def storage_failed():
        ui.notification_show(COPY["storage_error"], type="warning", duration=None)

    # ================= sample filter =================
    @reactive.calc
    def chosen_samples() -> list[str]:
        out = []
        for gid, _, _ in SAMPLE_GROUPS:
            out += list(input[gid]() or ())
        return out

    def set_samples(keep: list[str]):
        for gid, _, choices in SAMPLE_GROUPS:
            ui.update_checkbox_group(gid, selected=[s for s in choices if s in keep])

    @reactive.calc
    def cell_mask() -> np.ndarray | None:
        chosen = set(chosen_samples())
        if not SAMPLE_IDS or chosen >= set(SAMPLE_IDS):
            return None
        return np.isin(SAMPLE_COLUMN, list(chosen))

    @reactive.effect
    @reactive.event(input.samples_all, input.samples_reset)
    def _all():
        set_samples(SAMPLE_IDS)

    @reactive.effect
    @reactive.event(input.samples_sz)
    def _sz():
        set_samples(SZ_SAMPLES)

    @reactive.effect
    @reactive.event(input.samples_tip)
    def _tip():
        set_samples(TIP_SAMPLES)

    @render.ui
    def sample_note():
        mask = cell_mask()
        n = DATA.n_cells if mask is None else int(mask.sum())
        chosen = len(chosen_samples())
        return ui.div(COPY["sample_summary"].format(selected=chosen, total=len(SAMPLE_IDS), cells=n),
                      class_="selection-meta", role="status")

    # ================= Cell map =================
    @reactive.calc
    def settings() -> PlotSettings:
        colour = input.colour_by()
        base = dict(
            annotation=input.annotation(), facet=input.facet(), max_points=int(input.max_points()),
            point_size=input.point_size(), opacity=input.opacity(), expr_quantile=input.expr_max(),
            zero_display=input.zero_display(), cell_mask=cell_mask(), facet_context=input.facet_context(),
        )
        if colour == "fate" and FATE_CELLS.available:
            return PlotSettings(mode="fate", values=FATE_CELLS.p, colorbar_title="P(nodule<br>fate)", **base)
        if colour == "selection":
            if input.sel_kind() == "set":
                scored = geneset_score()
                return PlotSettings(mode="geneset", values=scored[0] if scored else None,
                                    colorbar_title="module<br>score", **base)
            return PlotSettings(mode="gene", gene_row=gene_row(), **base)
        return PlotSettings(mode="annotation", **base)

    @reactive.calc
    def panel_count() -> int:
        column = input.facet()
        if not column:
            return 1
        col = DATA.cells[column]
        mask = cell_mask()
        counts = (col[mask] if mask is not None else col).value_counts()
        return sum(1 for c in col.cat.categories if counts.get(c, 0) > 0)

    rendered_counts = reactive.value((0, 0))

    @reactive.effect
    async def push_figure():
        # Pushed to the existing plot otherwise keeps resetting figure
        current = settings()
        panels = panel_count()
        sync = bool(input.sync_facets()) and bool(current.facet)
        figure = build_figure(DATA, current, AXES)
        figure.update_layout(font=dict(color="#202020"), paper_bgcolor="#FFFFFF")
        points = sum(len(t.x) if t.x is not None else 0 for t in figure.data)
        backdrop = sum(len(t.x) if t.x is not None else 0 for t in figure.data if t.name == "all cells")
        rendered_counts.set((points - backdrop, backdrop))
        await session.send_custom_message(
            "umap-update",
            {
                "figure": json.loads(pio.to_json(figure)),
                "config": PLOT_CONFIG,
                "height": figure_height(current, panels) if current.facet else "max(470px, calc(100vh - 280px))",
                "minHeight": plot_min_height(current, panels) if current.facet else "470px",
                "align": sync,
            },
        )

    @reactive.effect
    @reactive.event(input.reset_view)
    async def reset_view():
        await session.send_custom_message("umap-reset", {"camera": CAMERA})

    @render.ui
    def plot_title():
        colour = input.colour_by()
        if colour == "fate":
            return caption(COPY["plot_fate_caption"])
        if colour == "selection":
            mask = cell_mask()
            scope = "across the full dataset" if mask is None else "in the selected samples"
            if input.sel_kind() == "set":
                scored, (rows, _) = geneset_score(), geneset_parsed()
                if not rows:
                    return caption("Paste genes or load a precomputed set in the Gene sets panel.")
                values, n_ctrl = scored
                values = values if mask is None else values[mask]
                if not len(values):
                    return caption(COPY["plot_no_cells"])
                return caption(COPY["plot_set_caption"].format(
                    genes=len(rows), controls=n_ctrl, percent=(values > 0).mean() * 100, scope=scope))
            row = gene_row()
            if row is None:
                return caption(COPY["plot_gene_prompt"])
            gene = DATA.genes.iloc[row]
            values = DATA.expression(row)
            values = values if mask is None else values[mask]
            if not len(values):
                return caption(COPY["plot_no_cells"])
            detected = int(np.count_nonzero(values > 0))
            detail = f" - {gene.description}" if gene.description else ""
            return caption(
                f"<b style='color:#111'>{gene.label}</b> ({gene.gene_id}){detail}<br>detected in "
                f"{detected:,} cells ({detected / len(values) * 100:.1f}%) {scope}")
        column = input.annotation()
        return caption(f"<b style='color:#111'>{DATA.label_for(column.replace('num__', ''))}</b>")

    @render.ui
    def summary_note():
        _, backdrop = rendered_counts.get()
        if not backdrop or not COPY["render_background"]:
            return ui.HTML("")
        return ui.div(COPY["render_background"].format(points=backdrop), class_="help-copy")

    # ================= Expression tab (dot plot) =================
    @reactive.calc
    def heat_result():
        """(figure, clusters, note) of the gene heatmap, or None when it does not apply."""
        # The heatmap shows the set in the box even while one gene is selected (a row click
        # selects a gene; the set's heatmap stays, with that row highlighted).
        if input.expr_view() != "heat" or HEATMAP is None:
            return None
        rows, _ = geneset_parsed()
        if len(rows) < 2:
            return None
        k = int(input.heat_k() or 0) if input.heat_clustering() else 0
        return HEATMAP.build(rows, input.dot_rows(), bool(input.heat_by_condition()), cell_mask(),
                             input.heat_scale(), max(0, min(12, k)))

    sent: dict[str, object] = {}       # target -> key of what that hidden-tab plot last showed

    def tab_is(name: str) -> bool:
        return (input.main_tabs() or "cell_map") == name

    def selection_key():
        if input.sel_kind() == "set":
            return ("set", tuple(sorted(geneset_parsed()[0])))
        return ("gene", gene_row())

    @reactive.effect
    async def push_dotplot():
        if not tab_is("expression"):
            return
        if input.expr_view() == "heat":
            key = ("heat", tuple(sorted(geneset_parsed()[0])), selected_gene_id(), input.dot_rows(),
                   bool(input.heat_by_condition()), input.heat_scale(), input.heat_k(), input.heat_clustering(), input.heat_height(), input.heat_labels(),
                   tuple(sorted(chosen_samples())))
        else:
            key = (selection_key(), input.expr_view(), input.dot_rows(), bool(input.dot_split()),
                   tuple(sorted(chosen_samples())))
        if sent.get("dotplot") == key:
            return
        sent["dotplot"] = key
        if input.expr_view() == "heat":
            res = heat_result()
            fig = HEATMAP.highlight(res[0], selected_gene_id()) if res and res[0] is not None else None
            if fig is not None:
                import plotly.graph_objects as go
                fig = go.Figure(fig)
                n_rows = len((fig.layout.meta or {}).get("genes", []))
                base_height = max(fig.layout.height, 170 + 14 * n_rows) if input.heat_labels() and n_rows <= 250 else fig.layout.height
                fig.update_layout(height=max(320, int(base_height * input.heat_height() / 100)))
                fig.update_yaxes(tickfont=dict(size=11))
                fig.update_xaxes(tickfont=dict(size=11))
                if not input.heat_labels():
                    fig.update_yaxes(showticklabels=False)
            await push_plot("dotplot", fig, "gene_heatmap")
            return
        sel = selection()
        if sel is None or not DOTPLOT.available:
            await push_plot("dotplot", None, "dotplot")
            return
        values, kind, label = sel
        fig = DOTPLOT.figure(values, kind, input.dot_rows(), bool(input.dot_split()), cell_mask(), label)
        await push_plot("dotplot", fig, f"dotplot_{label}")

    @render.ui
    def heat_note():
        if len(geneset_parsed()[0]) < 2:
            return ui.p(COPY["empty_heatmap"], class_="help-copy")
        res = heat_result()
        if res is None:
            return ui.HTML("")
        _, clusters, note = res
        bits = [f"{len(clusters)} clusters" if clusters else "Clustering disabled"]
        if note:
            bits.append(note)
        if len(geneset_parsed()[0]) > 250:
            bits.append(COPY["heat_labels_omitted"])
        return ui.HTML(f"<span style='color:#666;font-size:0.9em'>{'; '.join(bits)}. "
                       f"{COPY['heat_inspection']}</span>")

    @render.ui
    def heat_cluster_picker():
        res = heat_result()
        if not res or not res[1]:
            return ui.HTML("<span style='color:#999;font-size:0.85em'>Clustering disabled</span>")
        choices = {"": "(choose)", **{c: f"{c} ({len(g)} genes)" for c, g in res[1].items()}}
        return ui.div(
            field("Select (sub)cluster", ui.input_select("heat_cluster", "Select (sub)cluster", choices)),
            ui.panel_conditional("input.heat_cluster && input.heat_clustering",
                ui.div(ui.input_text("heat_cluster_name", "Custom set name", placeholder="Name this set"),
                       ui.input_action_button("heat_use_cluster", "save as custom set"),
                       class_="cluster-save-fields"),
                class_="cluster-save"),
            class_="cluster-picker",
        )

    @reactive.effect
    @reactive.event(input.heat_use_cluster)
    async def use_cluster():
        """Load a heatmap cluster as the gene set and keep it under Custom sets,
        named after the set it came from (e.g. 'GEP09 (top 200) · C2')."""
        res = heat_result()
        key = input.heat_cluster()
        if not (res and key in res[1]):
            return
        parent = active_name.get() or "Custom set"
        name = (input.heat_cluster_name() or "").strip()[:80] or unique_name(f"{parent} · {key}")
        await store_custom({**custom_sets.get(), name: list(res[1][key])}, f"custom:{name}")
        activate_set()
        ui.notification_show(COPY["cluster_saved"].format(cluster=key, name=name),
                             duration=5)

    @reactive.effect
    @reactive.event(input.heat_gene_click)
    def heat_gene_selected():
        """A click on a heatmap row makes that gene the selection (the set stays in its box)."""
        gid = input.heat_gene_click()
        if not gid or DATA.resolve_gene(gid) is None:
            return
        ui.update_text("gene_query", value=gid)
        ui.update_radio_buttons("details_kind", selected="gene")

    # ================= Fate tab (expression along fate) =================
    @reactive.calc
    def fate_selection():
        if input.sel_kind() == "set" and input.fate_target() == "inspected" and gene_row() is not None:
            return DATA.expression(gene_row()), "gene", str(DATA.genes.label.iloc[gene_row()])
        return selection()

    @reactive.effect
    @reactive.event(input.sel_kind, input.fate_target, selected_gene_id)
    def update_fate_measure_labels():
        sel = fate_selection()
        kind = sel[1] if sel is not None else input.sel_kind()
        prefix = "fate_set" if kind == "set" else "fate"
        ui.update_radio_buttons(
            "fate_measure",
            choices={"cpm": COPY[f"{prefix}_expression_label"],
                     "detection": COPY[f"{prefix}_detection_label"]},
            selected=input.fate_measure(),
        )

    @render.ui
    def fate_help_text():
        sel = fate_selection()
        if sel is None:
            return ui.HTML("")
        kind = sel[1]
        measure = "detection" if input.fate_measure() == "detection" else "expression"
        key = f"fate_{kind}_{measure}_help"
        return ui.div(
            ui.p(COPY["fate_help"], class_="help-copy"),
            ui.p(COPY[key], class_="help-copy", **{"data-copy": key}),
        )

    @reactive.effect
    async def push_fateplot():
        if not tab_is("fate"):
            return
        key = (selection_key(), input.fate_measure(), input.fate_platform(), input.fate_target(), selected_gene_id())
        if sent.get("fateplot") == key:
            return
        sent["fateplot"] = key
        sel = fate_selection()
        if sel is None or not FATE_CELLS.available:
            await push_plot("fateplot", None, "fate")
            return
        values, kind, label = sel
        status = None
        if kind == "gene":
            values = np.expm1(values.astype(np.float64)) * FATE_CELLS.cpm_factor
            rows = FATE.rows(selected_gene_id()) if FATE.available else None
            status = {} if rows is None else {(r.source, int(r.time_h)): r.status for r in rows.itertuples()}
        fig = FATE_CELLS.figure(values, kind, input.fate_measure(), input.fate_platform(), status)
        fig.update_layout(height=640, margin=dict(t=80),
                          title=dict(text=f"<b>{label}</b>", x=0, xref="container", xanchor="left", y=1.0,
                                     yref="container", yanchor="top", pad=dict(t=8, l=10), font=dict(size=14)))
        await push_plot("fateplot", fig, f"fate_{label}")

    @reactive.effect
    @reactive.event(input.open_fate)
    def _open_fate():
        ui.update_radio_buttons("fate_target", selected="inspected" if input.details_kind() == "gene" else "active")
        ui.update_navs("main_tabs", selected="fate")

    @reactive.calc
    def evidence_genes():
        return current_set() if input.details_kind() == "set" or selected_gene_id() is None else None

    # ================= right sidebar: verdicts =================
    @reactive.effect
    async def push_dissoc_gene():
        gid = selected_gene_id()
        fig = None
        genes = evidence_genes()
        if genes and DISSOC.available:
            fig = DISSOC.heatmap(genes)
            gid = "gene_set"
        elif gid is not None and DISSOC.available:
            fig = DISSOC.gene_figure(gid, str(DATA.genes.label.iloc[gene_row()]))
        await push_plot("dissoc_gene", fig, f"dissociation_{gid}")

    @render.ui
    def gene_header():
        genes = evidence_genes()
        if genes is not None:
            if not genes:
                return ui.HTML("<p style='color:#666'>Paste genes or load a precomputed set in "
                               "the Gene sets panel.</p>")
            return ui.HTML(f"<div style='line-height:1.35;font-size:1.05em'><b>Gene set</b> "
                           f"<span style='color:#666;font-size:0.85em'>({len(genes):,} genes)"
                           f"</span></div>")
        gid = selected_gene_id()
        if gid is None:
            return ui.p(COPY["empty_gene"], class_="help-copy")
        label = DATA.genes.label.iloc[gene_row()]
        return ui.div(ui.strong(str(label), class_="gene-name" if str(label) != gid else "gene-id"),
                      ui.div(gid, class_="selection-meta") if str(label) != gid else None)


    @render.ui
    def dissoc_verdict():
        if not DISSOC.available:
            return ui.HTML("<p style='color:#b00'>dissociation.parquet missing: run "
                           "build/build_dissociation.py.</p>")
        genes = evidence_genes()
        if genes:
            rows = [DISSOC.row(g) for g, _ in genes]
            tested = [r for r in rows if r is not None]
            n_p = sum(bool(r.protoplast_call) for r in tested)
            n_n = sum(bool(r.nuclei_call) for r in tested)
            n_any = sum(bool(r.protoplast_call or r.nuclei_call) for r in tested)
            return ui.HTML(
                f"<div style='font-size:0.9em;line-height:1.35'>{n_any} of {len(tested)} tested "
                f"genes dissociation-responsive (protoplasting {n_p}, nuclei extraction {n_n})"
                f"{f'; {len(genes) - len(tested)} not tested in bulk' if len(tested) < len(genes) else ''}."
                f" <span style='color:#888'>* = responsive.</span></div>")
        gid = selected_gene_id()
        if gid is None:
            return ui.HTML("")
        verdict = DISSOC.verdict(gid)
        colour = "#b45309" if verdict.startswith("dissociation-responsive") else "#555"
        return ui.HTML(f"<div style='color:{colour};font-size:0.9em;line-height:1.35'>{verdict}</div>")

    @render.ui
    def fate_verdict():
        if not FATE.available:
            return ui.HTML("<p style='color:#b00'>fate_association.parquet missing: run "
                           "build/build_fate_association.py.</p>")
        genes = evidence_genes()
        if genes:
            qualifying = []
            for gid, _ in genes:
                rows = FATE.rows(gid)
                if rows is not None:
                    passed = rows.loc[rows.passes]
                    if not passed.empty:
                        qualifying.append(passed)
            n_negative = sum(bool((rows.log2_effect < 0).any()) for rows in qualifying)
            n_positive = sum(bool((rows.log2_effect > 0).any()) for rows in qualifying)
            return ui.div(
                ui.div(COPY["fate_set_summary"].format(
                    associated=len(qualifying), total=len(genes),
                    negative=n_negative, positive=n_positive)),
                ui.p(COPY["fate_set_criteria"], class_="help-copy"),
                ui.span(COPY["fate_set_legend"], class_="help-copy"),
                style="font-size:0.9em;line-height:1.35",
            )
        gid = selected_gene_id()
        if gid is None:
            return ui.HTML("")
        hit = GENE_MAP.locate(gid) if GENE_MAP.available else None
        text, strong = FATE.verdict(gid, hit["module"] if hit else None)
        colour = "#202020" if strong else "#666"
        return ui.HTML(f"<div style='color:{colour};font-size:0.9em;line-height:1.35'>{text}</div>")

    @reactive.effect
    async def push_fate_gene():
        gid = selected_gene_id()
        genes = evidence_genes()
        if genes:
            fig, gid = fate_heatmap(FATE, genes), "gene_set"
        else:
            fig = FATE.figure(gid, str(DATA.genes.label.iloc[gene_row()])) if (gid is not None and FATE.available) else None
        await push_plot("fate_gene", fig, f"fate_association_{gid}")

    # ================= Gene map tab =================
    @reactive.calc
    def gene_map_dim() -> int:
        dim = int(input.gene_map_dim() or 2)
        return dim if str(dim) in GENE_MAP_JSON else 2

    def highlight_message(dim: int) -> dict:
        """Marker(s) on the gene map: the whole active gene set, or the selected gene."""
        genes = evidence_genes()
        if genes and GENE_MAP.available:
            hits = [h for h in (GENE_MAP.locate(g, dim) for g, _ in genes) if h]
            return {"xs": [h["x"] for h in hits], "ys": [h["y"] for h in hits],
                    "zs": [h["z"] for h in hits],
                    "labels": [h["label"] for h in hits] if len(hits) <= 15 else [""] * len(hits)}
        hit = GENE_MAP.locate(selected_gene_id(), dim) if GENE_MAP.available else None
        return hit or {"x": None, "y": None, "z": None, "label": ""}

    gene_map_drawn = reactive.value(None)   # dimension currently drawn in the browser (None = not yet)

    @reactive.effect
    async def init_gene_map():
        if not GENE_MAP_JSON or not tab_is("gene_umap"):
            return
        dim = gene_map_dim()
        with reactive.isolate():
            if gene_map_drawn.get() == dim:
                return
        gene_map_drawn.set(dim)
        await session.send_custom_message(
            "geneumap-init",
            {"figure": GENE_MAP_JSON[str(dim)],
             "config": {**PLOT_CONFIG, "toImageButtonOptions": {
                 "format": "png", "scale": 2,
                 "filename": "gene_umap_2d" if dim == 2 else "gene_umap_3d"}}},
        )
        # re-place the ring on the freshly drawn layout
        with reactive.isolate():
            msg = highlight_message(dim)
        await session.send_custom_message("geneumap-highlight", msg)

    @reactive.effect
    async def highlight_gene_map():
        msg = highlight_message(gene_map_dim())
        if gene_map_drawn.get() is None:      # nothing drawn yet: the first draw places the ring
            return
        await session.send_custom_message("geneumap-highlight", msg)

    @reactive.effect
    @reactive.event(input.geneumap_click)
    def gene_map_clicked():
        gene_id = input.geneumap_click()
        if not gene_id or DATA.resolve_gene(gene_id) is None:
            return
        ui.update_text("gene_query", value=gene_id)
        ui.update_radio_buttons("details_kind", selected="gene")

    @render.ui
    def gene_map_note():
        if not GENE_MAP.available:
            return ui.HTML("<p style='color:#b00'>gene_umap.parquet missing: run "
                           "build/build_gene_umap.py.</p>")
        genes = evidence_genes()
        if genes:
            n_on = sum(GENE_MAP.locate(g) is not None for g, _ in genes)
            return ui.HTML(f"<span style='font-size:0.9em'><b>{n_on}</b> out of {len(genes)} "
                           "highlighted on map</span>")
        gid = selected_gene_id()
        if gid is None:
            return ui.HTML("")
        hit = GENE_MAP.locate(gid)
        label = DATA.genes.label.iloc[gene_row()]
        if hit is None:
            return ui.HTML(f"<span style='color:#666;font-size:0.9em'><b>{label}</b> is not among "
                           "the 15,531 fate-associated genes on this map.</span>")
        return ui.HTML(f"<span style='font-size:0.9em'>Selected: <b>{hit['label']}</b> "
                       f"<span style='color:#666'>({gid}) &middot; {hit['module']}</span></span>")


app = App(app_ui, server, static_assets={"/assets": ASSETS})
