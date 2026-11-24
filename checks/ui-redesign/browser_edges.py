"""Targeted empty-state, set-evidence, accessibility and final screenshot checks."""
import json, os
from pathlib import Path
from playwright.sync_api import sync_playwright
OUT=Path(__file__).parent
with sync_playwright() as p:
    b=p.chromium.launch(headless=True,args=['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader'])
    page=b.new_page(viewport={'width':1366,'height':768});page.set_default_timeout(45000)
    errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
    page.goto(os.environ.get('BROWSER_URL','http://127.0.0.1:54032'),wait_until='networkidle')
    page.wait_for_function("document.getElementById('umap3d').data?.length>0")
    page.locator('.app-header').screenshot(path=str(OUT/'review-01-header.png'))
    page.locator('#gene_query').fill('MtrunA17_CPg0492181')
    page.get_by_role('tab',name='Genes',exact=True).click()
    page.wait_for_function("document.querySelector('#gene_map_note').textContent.includes('is not among')")
    assert page.locator('#gene_header').inner_text().count('MtrunA17_CPg0492181')==1
    page.locator('#gene_query').fill('no_gene_matches_this_query')
    page.wait_for_function("document.querySelector('#gene_picker').textContent.includes('Gene not found')")
    page.locator('#gene_query').fill('MtNIN')
    page.wait_for_function("document.querySelector('#selection_summary').textContent.includes('MtNIN')")
    page.locator('[data-panel=sets]').click()
    page.locator('#geneset_preset').select_option('GEP12')
    page.wait_for_function("document.querySelector('#geneset_note').textContent.includes('200 Genes matched')")
    page.locator('#use_set').click()
    page.wait_for_function("document.querySelector('#selection_summary').textContent.includes('Active set:')")
    page.locator('[data-panel=details]').click()
    page.locator('input[name=details_kind][value=set]').check()
    try:
        page.wait_for_function("document.querySelector('#gene_header').textContent.includes('200 genes')")
    except Exception:
        print(page.evaluate("({header:document.querySelector('#gene_header').textContent,selection:document.querySelector('#selection_summary').textContent,kind:Shiny.shinyapp.$inputValues.sel_kind,details:Shiny.shinyapp.$inputValues.details_kind,errors:[...document.querySelectorAll('.shiny-output-error')].map(e=>e.textContent)})"),flush=True)
        page.screenshot(path=str(OUT/'edge-failure.png'))
        raise
    page.locator('[data-evidence=fate]').click()
    page.wait_for_function("document.getElementById('fate_gene').data[0]?.type==='heatmap'")
    page.screenshot(path=str(OUT/'set-evidence-laptop.png'))
    page.keyboard.press('Escape')
    assert page.evaluate("document.activeElement.dataset.evidence")=='fate'
    page.locator('input[name=details_kind][value=gene]').check()
    page.locator('[data-evidence=fate]').click();page.locator('#open_fate').click()
    page.wait_for_function("document.getElementById('fateplot').layout?.yaxis?.title.text==='Gene expression (weighted mean CPM)'")
    assert page.evaluate("Shiny.shinyapp.$inputValues.sel_kind")=='set'
    page.locator('[data-panel=samples]').click()
    for control in page.locator('[data-group]').all():
        control.uncheck()
    page.wait_for_function("document.querySelector('#sample_note').textContent.includes('0 of 29')")
    page.get_by_role('tab',name='Expression',exact=True).click()
    page.wait_for_function("document.getElementById('dotplot').layout?.annotations?.[0]?.text.includes('No group')")
    page.locator('#samples_reset').click()
    page.get_by_role('tab',name='Cells',exact=True).click()
    page.locator('#panel-samples [data-panel=none]').click()
    page.wait_for_timeout(1200)
    page.screenshot(path=str(OUT/'final-cells-laptop.png'))
    page.set_viewport_size({'width':1920,'height':1080});page.wait_for_timeout(800)
    page.screenshot(path=str(OUT/'final-cells-desktop.png'))
    assert not errors,errors
    assert not page.locator('.shiny-output-error').count()
    print('PASS: absent map gene, unmatched search, set evidence, dialog focus restoration, evidence-to-Fate link retaining set, empty sample selection and final laptop/desktop rendering')
    b.close()
