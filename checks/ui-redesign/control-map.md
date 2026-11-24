# UI control and scientific-scope map

| Existing control | Redesigned location | Scientific scope |
| --- | --- | --- |
| Gene search / one gene | Persistent search; Use this gene in details | Active gene drives Cells, Expression and Fate; inspected gene drives evidence |
| Set text / saved sets | Gene sets contextual panel | Same matched-control scoring and browser-local persistence |
| Colour by / annotation | Cells toolbar | Rendering only |
| Display accordion | Plot settings contextual panel | Drawing budget, size, opacity, facets, background, quantile |
| Sample checkbox groups | Samples & filters contextual panel | Masks Cells and Expression; does not refit precomputed analyses |
| Gene map dimension | Genes toolbar | Existing 2D and 3D coordinates and module colours |
| Dot / heatmap | Expression toolbar | Dot: selection summary; heatmap: individual genes |
| Heatmap scale, k, clusters | Plot settings | Existing means, z-scores, seed 0, 10 restarts, ordering and 400-gene cap |
| Fate measure/platform | Fate toolbar | All screen cells; forward-mass-weighted deciles, equally weighted platforms |
| Fate / dissociation sidebars | Concise details + wide evidence dialogs | Existing precomputed results and original plots |

## Resolved ambiguities

- Only Genes had a supported 2D view. Cells retains its actual 3D embedding and reset/fullscreen controls.
- Dot plots of sets represent module scores, not individual-gene means. The label now explicitly says “Gene-set score summary”.
- Inspection and active selection are separate: searching or clicking a heatmap/map gene does not discard an active set. Use this gene explicitly switches selection mode.
- The full-atlas detection and score-positive percentages in Cells captions are not recalculated by sample selection and are labelled accordingly.
- Colour-scale quantiles and expression-matched control panels retain their original full-atlas reference.
- Rendering budgets are not analysis cell counts. Count actual foreground trace points and report additional repeated facet backgrounds separately.
- Ambiguous symbols require a unique gene ID in the new set editor instead of silently choosing the first database match. Calculations on uniquely resolved IDs are unchanged.
- Fate eligibility metadata also requires >=100 retained cells. Methods now include that existing rule and use the absolute coefficient bound.
- “Criterion not met” replaces “n.s.” for dissociation bars failing the combined response call.
- No mockup images accompanied the attached brief; visual implementation follows its written colours and interaction requirements.
