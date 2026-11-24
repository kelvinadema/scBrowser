
```bash 
conda env create -f environment.yml
conda activate scBrowser
python -m shiny run --port 8000 app/app.py
```

open <http://127.0.0.1:8000>

# Data

- Expression counts defined as "log1p(counts / cell_total × 1e4)"
- Genes fitlered on > 1 UMI in >= 10 cells
- Cell types inherit PO terms

# Interface and copy review

The app remains Shiny for Python with Plotly and the existing `viewer_data` files.
The four views share a persistent gene search and active selection. Inspecting a
heatmap row or map gene retains the active set. Typing/choosing a search result or
**Plot this gene** activates that gene; choosing a preset or **Use set** activates the
set and its details, clearing the previous gene search and inspection. Samples affect Cells and Expression, while Fate and evidence retain
their original precomputed scope.

- `app/ui_layout.py`: reusable workspace, contextual panels and evidence dialogs.
- `app/_assets/browser.css` / `browser.js`: responsive styling and browser interaction.
- `app/ui_copy.json`: editable explanatory copy and dynamic wording templates.
- `checks/ui-redesign/control-map.md`: control mapping and scientific-scope decisions.
- `checks/ui-redesign/text-review.md` / `text-review.json`: completed editorial decisions.

The editorial pass is **complete**: all 197 listed items have approved decisions; none
are skipped or pending. Approved wording and removals are recorded in
`checks/ui-redesign/text-review.json`; retain computed values in wording templates.
After a copy change, restart Shiny to reload the JSON.

Validation scripts:

```bash
python checks/ui-redesign/science_check.py
python checks/ui-redesign/caption_filter_check.py
python checks/ui-redesign/dissociation_missing_check.py
python checks/ui-redesign/fate_status_check.py
# With Playwright/Chromium available and the app running on port 54032:
python checks/ui-redesign/browser_check.py
python checks/ui-redesign/browser_edges.py
```

Override the browser test URL with `BROWSER_URL`. Screenshots and browser results
contain real unpublished data and are ignored by Git. Browser dependencies used
for this review were installed temporarily under `/tmp`, not added to the app.

The post-review facet, selection and layout fixes are recorded in
`checks/ui-redesign/bugfixes.md`. Cluster selection and named saving are above the
Expression heatmap. Evidence heatmaps use the full set, with no 100-gene cap.
