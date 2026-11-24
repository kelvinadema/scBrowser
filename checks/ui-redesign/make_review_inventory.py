"""Create the initial review queue; never run over approvals without preserving them."""
import ast, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
existing=ROOT/'checks/ui-redesign/text-review.json'
if existing.exists() and any(i['status'] != 'awaiting review' for i in json.loads(existing.read_text())['items']):
    raise SystemExit('Review decisions exist: update the inventory without regenerating it.')
copy=json.loads((ROOT/'app/ui_copy.json').read_text())
locations={
 'intro':('Shared','Header subtitle','Explain the browser in one short line',75),
 'search_placeholder':('Shared','Gene search input','Guide searchable fields',50),
 'samples_scope':('Samples','Panel introduction','Explain filter applicability',180),
 'set_help':('Gene sets','Panel introduction','Explain gene input and review',160),
 'set_storage':('Gene sets','Below Save named set','Explain browser-local storage and replacement',160),
 'inspection_help':('Shared','Selection details','Distinguish inspection from replacing the set',140),
 'cells_help':('Cells','Plot footer and settings','Explain units and background cells',240),
 'genes_help':('Genes','Above the map','Explain precomputed map and markers',220),
 'dot_help':('Expression','Dot plot explanation','Explain summary calculation and small groups',320),
 'heat_help':('Expression','Heatmap settings','Explain scale, clustering, omissions and sampling',380),
 'fate_scope':('Fate','Above controls','Explain filter scope and independent axes',160),
 'fate_help':('Fate','Units, weighting & line styles','Explain scientific units and weighting',520),
 'fate_lines':('Fate','Units, weighting & line styles','Explain fit-status line styles',240),
 'fate_methods':('Evidence — fate','Methods disclosure','Retain thresholds, diagnostics and eligibility',None),
 'dissoc_methods':('Evidence — dissociation','Methods disclosure','Retain contrasts, units and response criteria',None),
 'evidence_scope':('Evidence — shared','Below dialog title','Identify precomputed independent evidence',110),
}
queue=[]
for key,value in copy.items():
    screen,location,purpose,limit=locations.get(key,('Shared' if key.startswith(('empty_gene','no_','selection_','storage')) else
        'Samples' if key.startswith('sample_') else 'Cells' if key.startswith('render_') else
        'Expression' if key.startswith(('heat_','cluster_','empty_heat')) else 'Gene sets',
        key.replace('_',' ').title(),'User-facing help, status or feedback',180))
    queue.append(dict(id=key,screen=screen,location=location,purpose=purpose,
        kind='dynamic wording template' if '{' in value else 'existing scientific text (source-checked, adapted)' if key.endswith('_methods') else 'draft placeholder',
        source=f'app/ui_copy.json#{key}',current=value,suggested_max_characters=limit,status='awaiting review'))
# Existing inline scientific descriptions remain at their established code locations.
for module,functions in {'app.py':['plot_title','heat_note','gene_header','fate_verdict','dissoc_verdict','dissoc_multi_note','gene_map_note','heat_cluster_picker'],
                         'fate_assoc.py':['verdict'], 'dissociation.py':['verdict'],
                         'fate_view.py':['figure'], 'gene_heatmap.py':['build'], 'expression_view.py':['figure']}.items():
    tree=ast.parse((ROOT/'app'/module).read_text())
    for n in ast.walk(tree):
        if isinstance(n,(ast.FunctionDef,ast.AsyncFunctionDef)) and n.name in functions:
            queue.append(dict(id=f'{module}:{n.name}',screen='Evidence' if module in ('fate_assoc.py','dissociation.py') or n.name.endswith('verdict') else 'Scientific plots',
                location=f'{module} / {n.name}',purpose='Review wording templates while preserving computed values, units, categories and thresholds',
                kind='existing dynamic scientific text',source=f'app/{module}:{n.lineno}',status='awaiting review'))
for screen,location in [('Shared','Navigation and search labels'),('Cells','Colouring and plot settings labels'),('Samples','Quick selection, search and group labels'),('Gene sets','Editor, saved-set and validation labels'),('Genes','Dimensions and selected-map-gene labels'),('Expression','View, grouping, scale and clustering labels'),('Fate','Measure, platform and curves-for labels'),('Evidence','Dialog titles, close and methods controls')]:
    queue.append(dict(id='labels:'+screen,screen=screen,location=location,purpose='Review each visible control label individually if wording changes are requested',kind='control-label inventory',source='app/ui_layout.py',status='awaiting review'))
# Shared items reviewed once; review proceeds by screen.
order=['Shared','Cells','Samples','Genes','Gene sets','Expression','Fate','Evidence — shared','Evidence — fate','Evidence — dissociation','Evidence','Scientific plots']
queue.sort(key=lambda x:order.index(x['screen']) if x['screen'] in order else 99)
(ROOT/'checks/ui-redesign/text-review.json').write_text(json.dumps(dict(phase='awaiting first response',current_item='intro',items=queue),indent=2,ensure_ascii=False)+'\n')
lines=['# Editorial review — pending','',
       'No wording has been approved by the user yet. Review one item at a time, starting with `intro`. Shared text is reviewed once. A skipped item remains unresolved; only an explicit “keep” or supplied replacement resolves it. Preserve supplied replacement wording verbatim. Keep computed values and scientific rules intact.','',
       'Central editable help, placeholders and wording templates: `app/ui_copy.json`. Existing scientific templates retain their established module locations, indexed below. Control-label inventory rows expand to individual labels during their screen’s review.','',
       '| Item | Screen / location | Text class | Status |','| --- | --- | --- | --- |']
for q in queue:lines.append(f"| `{q['id']}` | {q['screen']} / {q['location']} | {q['kind']} | Awaiting review |")
lines += ['', 'Data-derived gene identities, annotations, sample IDs/counts, module names/colours, plotted coordinates/values, statistical results and thresholds are not draft findings. Review only their surrounding labels/templates.','',
          'After each answer: edit the source, update this checklist and `text-review.json`, inspect the changed location, then ask for the next item. After all answers: capture final laptop and desktop screens and report any skipped/pending entries.']
(ROOT/'checks/ui-redesign/text-review.md').write_text('\n'.join(lines)+'\n')
print(len(queue),'review entries, all pending')
