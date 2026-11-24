# Browser follow-up fixes — 2026-11-23

This follow-up implements the user's bug/layout requests after the 197-item text review.
It supersedes the earlier 100-gene evidence limit and the removed comparison controls.

- Facet categorical legends now occupy reserved space below the panels. Continuous
  colour scales occupy a separate right margin. Legend marker sizes remain readable.
  Fate scales also render when the first facet has no modelled cells; numeric
  annotation facets share one scale and its bounds.
- Orthographic wheel zoom changes `scene.aspectratio`, not just `scene.camera`.
  Synchronization now commits both values from the originating relayout event to all
  linked scenes. Using the event values avoids copying stale layout state back over
  a new zoom. Resize events no longer initiate camera synchronization.
- With linking off, each scene preserves its own rotation and zoom. Plot updates
  retain both. Reset restores the original camera and shared aspect ratio in every
  scene. The checkbox is labelled “Rotate and zoom facets together”.
- Background cells offer only Dimmed grey and Hidden.
- The dissociation Compare genes input, secondary plot and server callbacks are removed.
- Expression controls use aligned grid rows. Clustering, cluster count and selection
  are above the heatmap; choosing a cluster reveals an optional name and save button.
  Blank names retain the existing unique parent/cluster naming behavior. Explicit
  names follow the existing custom-set save/replace convention.
- User typing or choosing a search result activates the individual gene and its
  details/UMAP colouring. Loading a preset or pressing Use set activates the set and
  its details/UMAP colouring and clears the previous gene search and inspection.
  Saving and activating a heatmap cluster also clears the previous inspected gene.
  Programmatic gene-map/heatmap inspection still retains the active set.
- Reduced page/header padding, removed the desktop width ceiling, and increased the
  available embedding height. Cluster naming spans a full grid row when revealed.
- Fate and dissociation evidence heatmaps use the full selected set. Untested
  dissociation genes are still omitted. Plot height is bounded at 3,600 pixels to avoid
  extremely tall browser canvases; rows are not sampled and ordinary plot zoom remains
  available. The separate Expression heatmap 400-gene subsample is unchanged.

## Verification

- `bugfix_browser_check.py`: real wheel/drag gestures on both facets, linked and
  independent zoom, reset and redraw preservation, categorical/continuous legend
  separation, selection transitions, preset activation, row inspection, inline named
  cluster saving, removed controls, >100-row evidence plots and responsive layouts.
- `evidence_all_genes_check.py`: 500 genes retained in both evidence heatmaps,
  dissociation values unchanged, bounded plot height and separate Expression limit.
- Existing real-data science, caption-filter and missing-dissociation checks.
- `facet_scales_check.py`: empty-first-facet colourbar, shared numeric scale/bounds
  and consistent facet geometry.
- Main browser workflow suite: all 12 groups passed. Targeted bug checks and all
  listed real-data checks passed, as did Python/JavaScript syntax and diff checks.

Screenshots stay local under `checks/ui-redesign/bugfix-*.png`. Browser coverage is
headless Chromium with software WebGL; Safari/Firefox are not covered by these checks.
