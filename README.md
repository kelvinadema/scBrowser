# scBrowser

Interactive browser for a *Medicago truncatula* single-cell RNA-seq atlas, built with
[Shiny for Python](https://shiny.posit.co/py/) and Plotly.

Four views share one gene search and one active selection (a gene or a gene set):

- **Cells**: 3D cell map coloured by annotation, gene expression, gene-set score or nodule-fate probability, with sample filters and facets.
- **Genes**: 2D/3D map of genes by expression pattern.
- **Expression**: dot plot per cell group, or a gene heatmap with k-means clustering of the active set.
- **Fate association**: expression along predicted nodule fate, per source tissue.

The *Selection details* panel shows, for the inspected gene, its annotation, fate association and
dissociation response, with the underlying evidence plots.

Custom gene sets are stored in the visitor's browser (`localStorage`), not on the server.

## Running it

```bash
conda env create -f environment.yml
conda activate scBrowser
python -m shiny run --port 8000 app/app.py
```

Then open <http://127.0.0.1:8000>.

The app reads its data from `viewer_data/` next to `app/`, or from the folder in the `VIEWER_DATA`
environment variable. This data folder (about 5 GB) is **not** part of the repository. The scripts in
`build/` create it from the annotated AnnData object (`build/build_viewer_data.py` first, then the
other `build_*.py` scripts for the gene map, fate, dissociation and module-score panels).

## Data conventions

- Expression: log1p(counts / cell total × 10⁴), i.e. log1p(CP10K).
- Genes: kept if detected with more than 1 UMI in at least 10 cells.
- Cell types carry Plant Ontology (PO) terms.

## Repository layout

| Path | Contents |
|---|---|
| `app/app.py` | Shiny server: selection logic, plots, evidence panels |
| `app/ui_layout.py` | page layout, context panels and evidence dialogs |
| `app/ui_copy.json` | all explanatory text and wording templates (restart the app after editing) |
| `app/_assets/` | `browser.css`, `browser.js` and a local copy of `plotly.min.js` (no CDN needed) |
| `app/*.py` (others) | data access, plotting, gene sets, heatmap, fate and dissociation panels |
| `build/` | scripts that create `viewer_data/` |
| `checks/` | scripts that check the app against the data and drive a headless session |

## Checks

Run from the repository root with the app's environment and data available:

```bash
python checks/ui-redesign/science_check.py              # regression checks on the real data
python checks/ui-redesign/caption_filter_check.py       # cell-map captions with sample filters
python checks/ui-redesign/dissociation_missing_check.py # untested genes left out of dissociation heatmaps
python checks/ui-redesign/fate_status_check.py          # fate hover statuses
python checks/back_to_set_check.py 8000                 # "Back to set" flow; needs the app running on port 8000
```

`checks/ui-redesign/browser_*.py` and `set_activation_check.py` drive a real browser and need
Playwright with Chromium (`BROWSER_URL` sets the app's address).

## Offline Windows version

Releases tagged `win-…` are code-only updates for an offline Windows build of scBrowser. They
contain the app code only, no data, and are signed. Installed copies pick them up by themselves.
