"""Activating a set clears the prior gene; deliberate row inspection still works."""
import os
from pathlib import Path
from playwright.sync_api import sync_playwright

with sync_playwright() as p:
    browser = p.chromium.launch(headless=True, args=['--no-sandbox', '--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader'])
    page = browser.new_page(viewport={'width': 1440, 'height': 1000})
    page.set_default_timeout(60000)
    errors = []
    page.on('pageerror', lambda e: errors.append(str(e)))
    page.goto(os.environ.get('BROWSER_URL', 'http://127.0.0.1:54034'), wait_until='networkidle')

    def cleared():
        page.wait_for_function("""(() => {
          const v=Shiny.shinyapp.$inputValues;
          return v.sel_kind==='set' && v.details_kind==='set' && v.gene_query==='' && v.gene==='' &&
            document.querySelector('#gene_query').value==='' &&
            !document.querySelector('#selection_summary').textContent.includes('Inspecting:');
        })()""")
        assert not errors, errors
        assert not page.locator('.shiny-output-error').count()

    page.locator('#gene_query').fill('MtNIN')
    page.wait_for_function("Shiny.shinyapp.$inputValues.sel_kind==='gene' && document.querySelector('#selection_summary').textContent.includes('MtNIN')")
    page.locator('[data-panel=sets]').click()
    page.locator('#geneset_preset').select_option('GEP09')
    cleared()
    page.wait_for_function("document.querySelector('#plot_title').textContent.includes('200 selected genes')")
    print('PASS MtNIN to preset clears search, gene selection and inspection label; set colouring remains active', flush=True)

    page.get_by_role('tab', name='Expression', exact=True).click()
    page.locator('input[name=expr_view][value=heat]').check()
    page.wait_for_function("document.getElementById('dotplot').layout.meta?.genes?.length===200")

    def inspect_row():
        gid = page.evaluate("document.getElementById('dotplot').layout.meta.genes[0]")
        page.evaluate("document.getElementById('dotplot').emit('plotly_click',{points:[{y:0}]})")
        page.wait_for_function("id=>Shiny.shinyapp.$inputValues.gene===id && document.querySelector('#selection_summary').textContent.includes('Inspecting:')", arg=gid)
        assert page.evaluate("Shiny.shinyapp.$inputValues.sel_kind") == 'set'

    inspect_row()
    page.locator('#use_set').click()
    cleared()
    print('PASS deliberate row inspection retains set; Use set clears inspection again', flush=True)
    inspect_row()
    page.locator('#heat_cluster').select_option('C1')
    page.locator('#heat_cluster_name').fill('Activation test cluster')
    page.locator('#heat_use_cluster').click()
    page.wait_for_function("document.querySelector('#selection_summary').textContent.includes('Activation test cluster')")
    cleared()
    assert page.evaluate("JSON.parse(localStorage.getItem('mt_browser_custom_sets_v1'))['Activation test cluster'].length") > 0
    print('PASS saving and activating a cluster clears the previously inspected gene', flush=True)

    page.locator('#gene_query').fill('MtNIN')
    page.wait_for_function("Shiny.shinyapp.$inputValues.sel_kind==='gene' && document.querySelector('#selection_summary').textContent.includes('MtNIN')")
    page.locator('#geneset_preset').select_option('GEP09')
    cleared()
    page.screenshot(path=str(Path(__file__).with_name('set-activation.png')))
    print('PASS gene search after set and subsequent set activation remain consistent', flush=True)
    browser.close()
