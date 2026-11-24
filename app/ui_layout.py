"""Presentation only: one workspace, one contextual panel, persistent selections."""
import json
from pathlib import Path
from shiny import ui

COPY = json.loads(Path(__file__).with_name('ui_copy.json').read_text())


def copy(key):
    return ui.p(COPY[key], class_='help-copy', **{'data-copy': key})


def cell_help():
    """Explain grey cells according to the quantity currently plotted."""
    return ui.div(
        ui.panel_conditional("input.colour_by === 'selection' && input.sel_kind === 'gene'", copy('cells_help')),
        ui.panel_conditional("input.colour_by === 'selection' && input.sel_kind === 'set'", copy('cells_set_help')),
        ui.panel_conditional("input.colour_by === 'fate'", copy('cells_fate_help')),
    )


def button(label, **attrs):
    return ui.tags.button(label, type='button', class_='btn btn-outline-secondary', **attrs)


def build_ui(c, scripts):
    data = c['DATA']
    field, toolbar = c['field'], c['toolbar']

    def panel(key, title, *body):
        return ui.tags.section(ui.div(ui.h2(title), button('Close', **{'data-panel': 'none'}), class_='panel-heading'),
                          *body, id=f'panel-{key}', class_='context-panel', hidden=key != 'details')

    def evidence(key, title, *body):
        return ui.tags.dialog(
            ui.div(ui.h2(title), button('Close evidence', **{'data-close-evidence': key}), class_='panel-heading'),
            copy('evidence_scope'), *body, id=f'evidence-{key}', class_='evidence-dialog',
            **{'aria-label': title})

    samples = panel('samples', 'Samples & filters', copy('samples_scope'),
        ui.div(ui.input_action_button('samples_all', 'All samples'),
               ui.input_action_button('samples_sz', 'Susceptible zone'),
               ui.input_action_button('samples_tip', 'Root tip'), class_='quick-select'),
        ui.tags.label('Find samples or groups', **{'for': 'sample-search'}),
        ui.tags.input(id='sample-search', type='search', placeholder='Sample ID, condition or platform', class_='form-control'),
        *[ui.tags.details(ui.tags.summary(cond),
             ui.tags.label(ui.tags.input(type='checkbox', checked=True, **{'data-group': gid, 'aria-label': f'Select all {cond} samples'}), ' Select group', class_='group-toggle'),
             ui.input_checkbox_group(gid, 'Samples', choices, selected=list(choices)),
             class_='sample-group', open=True) for gid, cond, choices in c['SAMPLE_GROUPS']],
        ui.input_action_button('samples_reset', 'Reset to all samples'))
    sets = panel('sets', 'Gene sets', copy('set_help'),
        ui.input_select('geneset_preset', 'Saved & precomputed sets', c['PRESET_CHOICES']),
        ui.input_text_area('geneset_text', 'Gene IDs (Format: MtrunA17_ChrXgXXXXXXX) or gene names', rows=6, width='100%'),
        ui.output_ui('geneset_note'),
        ui.input_action_button('use_set', 'Use set', class_='btn-primary'),
        ui.input_text('custom_name', 'Set name', placeholder='Name this set', width='100%'),
        ui.input_action_button('custom_save', 'Save set'), ui.output_ui('custom_delete'), copy('set_storage'))
    settings = panel('settings', 'Plot settings',
        ui.panel_conditional("input.main_tabs === 'cell_map'",
            ui.input_select('facet', 'Facet by', data.facet_choices()),
            ui.input_checkbox('sync_facets', 'Rotate and zoom facets together', True),
            ui.input_checkbox('facet_context', 'Show background cells in facets', True),
            ui.input_select('max_points', COPY['cells_budget_label'], c['POINT_BUDGETS'], selected='100000'),
            ui.input_slider('point_size', 'Point size', .5, 5., 1.6, step=.1),
            ui.input_slider('opacity', 'Opacity', .1, 1., .75, step=.05),
            ui.input_slider('expr_max', COPY['cells_colour_limit_label'], .80, 1., .95, step=.01),
            copy('cells_colour_limit_help'),
            ui.input_select('zero_display', 'Background cells', {'dim':'Dimmed grey','hide':'Hidden'}),
            cell_help()),
        ui.panel_conditional("input.main_tabs === 'expression'",
            ui.input_radio_buttons('heat_scale', 'Color by', {'z':'Relative expression per gene (Z-score)','raw':'Mean expression (log1p(CP10K))'}, selected='z'),
            ui.input_slider('heat_height', 'Vertical heatmap size (%)', 60, 180, 100, step=10),
            ui.input_checkbox('heat_labels', 'Show gene labels when supported', True),
            copy('heat_help')))
    details = panel('details', 'Selection details',
        ui.panel_conditional("input.sel_kind === 'set'",ui.input_radio_buttons('details_kind','Details for',{'gene':'Inspected gene','set':'Active set'}, selected='gene',inline=True)),
        ui.output_ui('gene_header'), ui.output_ui('gene_annotation'),
        ui.panel_conditional("input.sel_kind === 'set'", copy('inspection_help')),
        ui.input_action_button('use_gene', 'Plot this gene'),
        ui.output_ui('back_to_set_ui'),
        ui.h3('Fate association'), ui.output_ui('fate_verdict'),
        button('View fate evidence', **{'data-evidence':'fate'}),
        ui.h3('Dissociation response'), ui.output_ui('dissoc_verdict'),
        button('View dissociation evidence', **{'data-evidence':'dissoc'}))
    tabs = ui.navset_tab(
        ui.nav_panel('Cells',
            toolbar(field('Colour cells by', ui.input_select('colour_by', 'Colour cells by',
                {'annotation':'Annotation','selection':'Selected gene / set','fate':'Nodule-fate probability'})),
                ui.panel_conditional("input.colour_by === 'annotation'", field('Annotation',ui.input_select('annotation','Annotation',data.annotation_choices()))),
                ui.div(ui.span('3D embedding',class_='selection-meta'),button('Plot settings', **{'data-panel':'settings'}),
                       ui.input_action_button('reset_view', 'Reset view'), button('Fullscreen', **{'data-fullscreen':'umap3d_wrap'}), class_='plot-actions')),
            ui.output_ui('plot_title'),ui.output_ui('summary_note'), ui.div(ui.div(id='umap3d'),id='umap3d_wrap',class_='embedding-wrap'),
            cell_help(), value='cell_map'),
        ui.nav_panel('Genes',
            toolbar(field(COPY['genes_dimensions_label'],ui.input_radio_buttons('gene_map_dim',COPY['genes_dimensions_label'],
                {'2':'2D', **({'3':'3D'} if c['GENE_MAP'].has_3d else {})}, selected='2',inline=True)),
                button('Fullscreen', **{'data-fullscreen':'geneumap_wrap'})),
            copy('genes_help'), ui.output_ui('gene_map_note'),
            ui.div(ui.div(id='geneumap'),id='geneumap_wrap',class_='embedding-wrap'),value='gene_umap'),
        ui.nav_panel('Expression',
            toolbar(field('View',ui.input_radio_buttons('expr_view','Expression view', {'dot':'Dot plot','heat':'Heatmap'},selected='dot',inline=True)),
                field('Group by',ui.input_select('dot_rows','Group by',c['ROW_CHOICES'],selected='ch5_tissue')),
                field('Split',
                    ui.panel_conditional("input.expr_view === 'dot'",ui.input_checkbox('dot_split','Split by platform (10x, BD)',False)),
                    ui.panel_conditional("input.expr_view === 'heat'",ui.input_checkbox('heat_by_condition','Split by condition',False))),
                ui.div(button('Plot settings', **{'data-panel':'settings'}), class_='toolbar-action'),
                cls='toolbar expression-toolbar'),
            ui.panel_conditional("input.expr_view === 'heat'",
                ui.div(field('Clustering',ui.input_checkbox('heat_clustering','Cluster genes',True)),
                       field('Number of clusters',ui.input_numeric('heat_k','Number of clusters',6,min=2,max=12)),
                       ui.output_ui('heat_cluster_picker'),class_='heat-controls')),
            ui.output_ui('expression_heading'),
            ui.panel_conditional("input.expr_view === 'dot' && input.sel_kind === 'gene'",copy('dot_help')),
            ui.panel_conditional("input.expr_view === 'dot' && input.sel_kind === 'set'",copy('dot_set_help')),
            ui.panel_conditional("input.expr_view === 'heat'",ui.output_ui('heat_note')),
            ui.div(id='dotplot'),value='expression'),
        ui.nav_panel(COPY['nav_fate'],copy('fate_scope'),
            toolbar(ui.panel_conditional("input.sel_kind === 'set'", field('Show curves for',ui.input_radio_buttons('fate_target','Show curves for',{'active':'Active set','inspected':'Inspected gene'},selected='active',inline=True))),
                    field(COPY['fate_measure_label'],ui.input_radio_buttons('fate_measure',COPY['fate_measure_label'],{'cpm':COPY['fate_expression_label'],'detection':COPY['fate_detection_label']},selected='cpm',inline=True)),
                    field('Platform',ui.input_radio_buttons('fate_platform','Platform',{'pooled':'Pooled (10x + BD)','BD':'BD','10xP':'10x'},selected='pooled',inline=True))),
            ui.div(id='fateplot'), ui.tags.details(ui.tags.summary('Units, weighting & line styles'),ui.output_ui('fate_help_text'),copy('fate_lines')),value='fate'),
        id='main_tabs', selected='cell_map')
    return ui.page_fluid(
        ui.head_content(ui.tags.link(rel='stylesheet',href='assets/browser.css'),ui.tags.script(src='assets/browser.js')),
        scripts,
        ui.tags.header(ui.div(ui.h1('Medicago ',ui.span('scRNA browser'))),
                       ui.div(f'{data.n_cells:,} cells',ui.span(' / '),f'{data.n_genes:,} genes',class_='atlas-count'), class_='app-header'),
        ui.div(ui.div(ui.input_text('gene_query',COPY['search_placeholder'],placeholder=COPY['search_placeholder']),
                      ui.div(ui.output_ui('gene_picker'),class_='search-results'),class_='gene-search'),
               button('Gene sets', **{'data-panel':'sets'}), ui.output_ui('selection_summary'),
               button('Selection details', **{'data-panel':'details'}), class_='selection-bar'),
        ui.div(ui.output_ui('sample_note'),button('Samples & filters', **{'data-panel':'samples'}),class_='sample-bar'),
        ui.div(ui.tags.main(tabs,class_='plot-card'), ui.tags.aside(details,samples,sets,settings,id='context-dock', **{'aria-label':'Contextual controls'}),class_='workspace'),
        ui.div(ui.input_radio_buttons('sel_kind','Active selection',{'gene':'Gene','set':'Set'},selected='gene'),hidden=True),
        evidence('fate','Fate association',ui.div(id='fate_gene'),ui.input_action_link('open_fate','View fate curves'),
                 ui.tags.details(ui.tags.summary('Fit statuses, thresholds & methods'),ui.HTML(COPY['fate_methods']),open=True)),
        evidence('dissoc','Dissociation response',ui.div(id='dissoc_gene'),
                 ui.tags.details(ui.tags.summary('Contrasts, thresholds & methods'),ui.HTML(COPY['dissoc_methods']),open=True)),
        title='Medicago scRNA browser')
