# UI redesign validation

Implementation uses the existing Shiny/Plotly app and actual local viewer data:
351,239 cells and 40,532 genes. The initial redesign preserved data files and scientific fitting/scoring algorithms.
The subsequent approved 0-HPI export addition is documented below. `control-map.md` records the old/new controls and filter applicability.

## Passed

- Python compilation and `git diff --check`.
- Real-data module scores match direct means minus matched controls within 1e-5.
- Sample-filtered heatmaps; identical clustering across display scales; fixed seed,
  10 restarts and ordering retained.
- Deterministic 400-gene sampling, explicit sampling disclosure and empty-group handling.
- Gene-ID parsing, duplicate/unmatched reporting; a dedicated alias fixture verifies
  ambiguous symbols are reported instead of silently resolved.
- Fate weighted detection remains in [0, 1]; existing fate/dissociation figures build.
- Chromium workflow suite: 12 workflow groups passed with no browser JavaScript errors
  or Shiny output errors. Includes search, evidence/camera restoration, sample presets,
  search, partial/group checkboxes, reset, set creation/editing/save/reload, dot/heatmap,
  inspection retaining sets, clustering and saved cluster sets, gene map 2D/3D and clicks,
  Fate controls, drawing budgets, opacity, faceting, reset and fullscreen.
- Screens inspected at 1366 × 768, 1440 × 900 and 1920 × 1080. Panels scroll independently,
  plots resize when panels close, and the tested widths have no horizontal page overflow.
- Targeted checks passed for absent-map genes, unmatched search, set evidence, dialog focus
  restoration, evidence-to-Fate navigation retaining the active set, and empty sample selection.
- Active checkbox/radio colours checked against charcoal RGB(37,37,37).

- Editorial follow-up: fixed active-set evidence help being hidden when no individual
  gene had been selected. Verified visibility in a fresh Chromium session using a real
  two-gene set and no inspected gene; screenshot `review-set-evidence-help.png`.
  The explanation was subsequently removed at the user’s request; its empty paragraph
  is no longer rendered.

## Presentation corrections made during visual checks

- Fixed a Shiny custom-message callback signature that blocked initial output.
- Fixed sample searching for Shiny's actual checkbox wrapper structure.
- Removed Bootstrap blue from controls and slider handles.
- Made settings contextual when switching to Genes/Fate.
- Increased heatmap row-label size and height without changing gene ordering or values.
- Reported actual rendered cells and repeated facet backdrops separately from selected cells.
- Distinguished positive low expression with a peach scale; undetected cells remain grey.
- Removed duplicated IDs from details when a gene has no separate display symbol.

## Limits and editorial follow-up

- Mockup images were not present in the attachment; the implementation follows the written
  palette and interaction brief. Visual parity with the absent images cannot be certified.
- Browser QA used headless Chromium with software WebGL. No Safari/Firefox or physical
  touch-device validation was performed.
- Cells has the existing 3D embedding; there was no supported Cells 2D projection to preserve.
- Large heatmaps and multi-tissue plots require vertical scrolling at laptop sizes.
- The editorial pass is complete: all 197 listed items have approved decisions, with no
  skipped or pending items. Decisions are recorded in `text-review.json`. The colour-saturation label and
  shortened explanation are approved, with the original 0.80–1.00 range retained. The drawn/selected cell-count sentence
  was removed at the user’s request during editorial review. The approved gene-expression
  wording is shown only for that colouring mode; set and fate explanations remain separate.
  The approved search placeholder was verified to fit laptop and desktop widths.
- Screenshots contain real unpublished data and stay local; Git ignores image outputs.
- Validation used a separate local server on port 54032. Restart the usual app process to
  load the source changes. Browser dependencies are temporary under `/tmp`.

- Editorial naming follow-up: checked the actual save callback with isolated state.
  Repeated unnamed saves produce Custom set, Custom set (1), Custom set (2); occupied
  names are skipped, and explicitly entered existing names still replace that set.

- Fate help now follows the same selection function as the plot and switches by measure.
  Isolated execution with Shiny UI verified all four quantity explanations, inspecting a
  gene while a set remains active, missing-gene fallback, and empty selection.

- Fate Y-axis option labels now switch between gene CPM/detection and gene-set score/positive-score fraction.
  Isolated callback checks covered active genes, active sets, inspected genes, missing-gene
  fallback and empty selections, preserving both possible measure selections.

- Gene-set fate summary now reports qualifying genes and negative/positive counts using
  the existing full criterion. Isolated renderer checks covered repeated contexts,
  both directions in one gene, nonpassing contexts, missing genes and sets over 100 genes.

- 0-HPI export follow-up (2026-11-22): rebuilt `viewer_data/fate_cells.parquet` with
  39,808 original `arches_e05/new_nodule` probabilities matched by cell ID. All added
  values exactly match the original endpoint export after float32 conversion. The
  existing 155,411 probabilities and every other export column are exactly unchanged,
  including all curve weights and deciles. Non-modelled cells retain missing values;
  96-HPI endpoints retain binary calls. Total modelled probabilities: 195,219.
  The builder records the supplementary source and validates compatibility on rebuild.
  Context QC and known-answer CPM checks passed (maximum relative error 0.001258,
  unchanged). App loader verified; targeted Chromium check passed for the exact
  approved caption and rendered fate trace, with no Shiny or JavaScript errors.
  Local screenshot: `zero-hpi-fate-caption.png`. Compilation and diff checks passed.
  This updates the local data export (ignored by Git) as well as the tracked builder;
  restart an already-running app to load it.

- Sample-dependent Cells captions (2026-11-22): gene detection count/percentage and
  gene-set positive-score percentage now use all cells in the selected samples.
  The wording switches between “across the full dataset” and “in the selected samples”.
  Empty sample selections show an explicit message. Set membership, background genes
  and per-cell scores are unchanged; drawing limits do not enter these calculations.
  `caption_filter_check.py` checks the actual renderer with real data for all samples,
  root tips, one sample, no samples and no gene, using sparse detection counts as
  an independent reference. Targeted Chromium checks verified live updates for gene
  and set captions when selecting root tips and resetting to all samples; numerical
  values matched the real data and no Shiny or JavaScript errors occurred.
  Python compilation and `git diff --check` passed.

- Dissociation heatmaps now omit genes without results, including the comparison and
  active-set evidence views. Tested non-responsive genes remain; fold changes, response
  markers and input ordering are preserved. Height follows the plotted rows, and an
  entirely untested list returns no figure (the existing browser handler clears it).
  The comparison note reports “{n} genes not shown (no dissociation response available)”.
  `dissociation_missing_check.py` passed using real data for mixed, entirely untested,
  tested-only and empty inputs, including exact omission counts, row order, values,
  non-responsive retention and reduced height. Python compilation and diff checks passed.

## Final editorial validation — 2026-11-23

- All 197 listed items are approved (106 kept, 88 replaced, 3 removed); no skipped
  or pending items remain. The final dot plot hover template is unchanged.
- The shared full-criterion wording now also appears in the fate evidence bar and
  heatmap tooltips, removing the remaining user-facing figure-number criterion references.
- Fixed missing fitted contexts bypassing the approved eligibility tooltip. The
  existing descriptive curves, line styles and numerical values are unchanged.
  `fate_status_check.py` verifies absent contexts, full/below/excluded statuses,
  and descriptive gene-set curves.
- Final real-data science checks, caption-filter checks, dissociation omission checks,
  fate-status checks, Python compilation and diff checks passed.
- Chromium main suite: all 12 workflow groups passed, with no JavaScript or Shiny
  output errors. Updated test assertions for approved wording and Plotly's decoded
  point arrays; no drawing-budget behavior change was required.
- Edge checks passed: absent-map gene, unmatched search, active-set evidence, dialog
  focus restoration, evidence-to-Fate navigation retaining the set, and empty sample
  selection. The test now waits for a preset to load before clicking Use set.
- Captured and inspected the final 1366 × 768 and 1920 × 1080 Cells views, plus
  1440 × 900 heatmap and Fate views. Text wraps within panels; long heatmaps and
  settings scroll. No horizontal page overflow in the workflow checks.
  Screenshots: `final-cells-laptop.png`, `final-cells-desktop.png`,
  `heatmap-desktop.png`, `fate-desktop.png`.
- Validation used temporary local servers on ports 54032 and 54033, with Python 3.12
  Playwright under `/tmp/medicago_browser_tools` and Chromium under
  `/tmp/medicago_browser_bins`. App and scientific checks use the scBrowser environment.
  Browser coverage remains Chromium with software WebGL; Safari/Firefox were not tested.

## Browser bug follow-up — 2026-11-23

- Single-gene fate and dissociation evidence plots now show the gene name as a
  plot title (locus tag fallback), including in exported figures. Added top space
  preserves the subplot area. Chromium checks verified visible MtNIN titles above
  subplot headings and updates to MtTIP1e in both dialogs, with no JavaScript/Shiny
  errors. Screenshots: `fate-gene-title.png`, `dissoc-gene-title.png`.
  Python compilation and diff checks passed.

- Subsequent selection correction: activating a preset, Use set or a saved cluster
  clears the previous gene query and selected gene. The inspection label is absent
  until a gene is deliberately inspected again. `set_activation_check.py` passed in
  Chromium for MtNIN → GEP09, row inspection → Use set, row inspection → saved
  cluster, and a subsequent gene → set transition. Set colouring remains active;
  deliberate row inspection retains the set. No JavaScript/Shiny errors occurred.
  Python compilation and diff checks passed. Screenshot: `set-activation.png`.

- Implemented the follow-up requests recorded in `bugfixes.md`: facet zoom/reset
  and legends, explicit gene/set activation, aligned Expression and cluster controls,
  compact layout, removed comparison/full-grey controls and uncapped evidence rows.
- All 12 main Chromium workflow groups passed. Targeted checks passed for real
  wheel/drag gestures in either facet, independent/linked zoom, reset, redraw
  preservation, separate legends, named cluster saving and selection transitions.
- Both evidence plots were checked with 500 genes; all available rows and values
  were preserved. Untested dissociation genes are omitted. The separate Expression
  heatmap subsample remains unchanged.
- Real-data science, caption-filter, dissociation omission and facet-scale checks
  passed. Facet checks include an empty first fate panel and shared numeric bounds.
  Python/JavaScript syntax and diff checks passed.
- Laptop/desktop layouts fit without horizontal page overflow. Local screenshots:
  `bugfix-facets.png`, `bugfix-expression.png`, `bugfix-expression-1366.png` and
  `bugfix-expression-1920.png`. Browser tests used temporary port 54034 and Chromium
  with software WebGL; Safari/Firefox were not tested.
