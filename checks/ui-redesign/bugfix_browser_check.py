"""Regression checks for the user-reported facet, selection and layout glitches."""
import json, os
from pathlib import Path
from playwright.sync_api import sync_playwright
OUT=Path(__file__).parent
with sync_playwright() as p:
 b=p.chromium.launch(headless=True,args=['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader'])
 page=b.new_page(viewport={'width':1440,'height':1000});page.set_default_timeout(60000)
 errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
 page.goto(os.environ.get('BROWSER_URL','http://127.0.0.1:54034'),wait_until='networkidle')
 page.wait_for_function("document.getElementById('umap3d')._fullData?.length>0")
 def wait(js): page.wait_for_function(js)
 def check(name):
  assert not errors,errors
  assert not page.locator('.shiny-output-error').count(),page.locator('.shiny-output-error').all_text_contents()
  print('PASS',name,flush=True)
 def state(key):return page.evaluate("k=>Shiny.shinyapp.$inputValues[k]",key)
 def views():return page.evaluate("Object.fromEntries(Object.keys(document.getElementById('umap3d')._fullLayout).filter(k=>/^scene[0-9]*$/.test(k)).map(k=>{const l=document.getElementById('umap3d')._fullLayout[k];return [k,{camera:l._scene.getCamera(),ratio:l.aspectratio}]}))")
 def close(a,c):
  if isinstance(a,dict):return all(close(a[k],c[k]) for k in a)
  return a==c if isinstance(a,str) else abs(a-c)<1e-5
 def wheel(key):
  pos=page.evaluate("k=>{const e=document.getElementById('umap3d'),r=e.getBoundingClientRect(),f=e._fullLayout,d=f[k].domain;return {x:r.x+f._size.l+f._size.w*(d.x[0]+d.x[1])/2,y:r.y+f._size.t+f._size.h*(1-(d.y[0]+d.y[1])/2)}}",key)
  page.mouse.move(pos['x'],pos['y']);page.mouse.wheel(0,-300);page.wait_for_timeout(1500)
 page.get_by_role('button',name='Plot settings',exact=True).click();page.locator('#max_points').select_option('50000')
 assert page.locator('#zero_display option').all_text_contents()==['Dimmed grey','Hidden']
 page.locator('#facet').select_option('ch5_region');wait("!!document.getElementById('umap3d')._fullLayout?.scene2");page.wait_for_timeout(1200)
 baseline=views();wheel('scene2');v=views();assert not close(v['scene2']['ratio'],baseline['scene2']['ratio']);assert close(v['scene'],v['scene2']),v
 wheel('scene');v2=views();assert not close(v2['scene']['ratio'],v['scene']['ratio']);assert close(v2['scene'],v2['scene2'])
 # Real pointer drag changes rotation and synchronizes the other scene too.
 pos=page.evaluate("(()=>{const e=document.getElementById('umap3d'),r=e.getBoundingClientRect(),f=e._fullLayout,d=f.scene2.domain;return [r.x+f._size.l+f._size.w*(d.x[0]+d.x[1])/2,r.y+f._size.t+f._size.h*(1-(d.y[0]+d.y[1])/2)]})()")
 page.mouse.move(*pos);page.mouse.down();page.mouse.move(pos[0]+45,pos[1]+20,steps=8);page.mouse.up();page.wait_for_timeout(1000)
 assert close(views()['scene'],views()['scene2'])
 page.locator('#reset_view').click();page.wait_for_timeout(1200);assert all(close(views()[k],baseline[k]) for k in baseline)
 page.locator('#sync_facets').uncheck();page.wait_for_timeout(800);wheel('scene2');v=views();assert close(v['scene'],baseline['scene']);assert not close(v['scene2']['ratio'],baseline['scene2']['ratio'])
 # A redraw must retain independent zoom; reset must then restore both panels.
 page.locator('#facet_context').uncheck();page.wait_for_timeout(1200);assert close(views()['scene2'],v['scene2'])
 page.locator('#reset_view').click();page.wait_for_timeout(1000);assert all(close(views()[k],baseline[k]) for k in baseline)
 assert page.evaluate("(()=>{const e=document.getElementById('umap3d'),r=e.getBoundingClientRect(),f=e._fullLayout,l=e.querySelector('.legend').getBoundingClientRect();return l.top>=r.top+f._size.t+f._size.h})()")
 page.screenshot(path=str(OUT/'bugfix-facets.png'));check('wheel zoom on either facet, synchronized rotation/zoom, independent zoom, preserved zoom on redraw, reset and separate legend')
 page.locator('#colour_by').select_option('fate');wait("document.getElementById('umap3d')._fullData.some(t=>t.marker?.showscale)")
 page.wait_for_timeout(700)
 assert page.evaluate("(()=>{const e=document.getElementById('umap3d'),r=e.getBoundingClientRect(),f=e._fullLayout,c=e.querySelector('.colorbar').getBoundingClientRect();return c.left>=r.left+f._size.l+f._size.w})()")
 check('continuous facet colour scale also stays outside scene domains')
 page.locator('#facet').select_option('');page.locator('#gene_query').fill('MtNIN');wait("Shiny.shinyapp.$inputValues.sel_kind==='gene' && document.querySelector('#selection_summary').textContent.includes('MtNIN')")
 page.locator('[data-panel=sets]').click();page.locator('#geneset_preset').select_option('GEP12')
 wait("Shiny.shinyapp.$inputValues.sel_kind==='set' && Shiny.shinyapp.$inputValues.details_kind==='set' && document.querySelector('#selection_summary').textContent.includes('200 genes')")
 page.locator('[data-panel=details]').click();wait("document.querySelector('#gene_header').textContent.includes('200 genes')")
 page.locator('[data-evidence=fate]').click();wait("document.getElementById('fate_gene').data?.[0]?.y.length===200");page.keyboard.press('Escape')
 page.locator('[data-evidence=dissoc]').click();wait("document.getElementById('dissoc_gene').data?.[0]?.y.length>100")
 assert page.locator('#dissoc_genes').count()==0;assert 'Compare genes' not in page.locator('#evidence-dissoc').inner_text();page.keyboard.press('Escape')
 page.get_by_role('tab',name='Expression',exact=True).click();page.locator('input[name=expr_view][value=heat]').check()
 wait("document.getElementById('dotplot').layout.meta?.genes?.length===200")
 page.locator('#heat_k').fill('4');page.locator('#heat_k').press('Tab');wait("document.querySelector('#heat_note').textContent.includes('4 clusters')")
 assert page.locator('#heat_cluster').is_visible();assert not page.locator('#heat_cluster_name').is_visible()
 page.locator('#heat_cluster').select_option('C2');page.locator('#heat_cluster_name').fill('Named test cluster')
 assert page.locator('.cluster-save-fields').evaluate('(e)=>e.getBoundingClientRect().width>e.closest(".heat-controls").getBoundingClientRect().width*.9')
 page.screenshot(path=str(OUT/'bugfix-expression.png'))
 page.locator('#heat_use_cluster').click();wait("document.querySelector('#selection_summary').textContent.includes('Named test cluster')")
 assert page.evaluate("JSON.parse(localStorage.getItem('mt_browser_custom_sets_v1'))['Named test cluster'].length")>0
 check('preset activates set and details, both evidence heatmaps exceed 100 rows, comparison removed, visible cluster controls and named saving')
 page.get_by_role('tab',name='Cells',exact=True).click();page.locator('#gene_query').fill('MtNF-YA1')
 wait("Shiny.shinyapp.$inputValues.sel_kind==='gene' && Shiny.shinyapp.$inputValues.details_kind==='gene' && document.querySelector('#selection_summary').textContent.includes('MtNF-YA1')")
 wait("document.getElementById('plot_title').textContent.includes('detected in')")
 page.locator('[data-panel=sets]').click();page.locator('#geneset_preset').select_option('GEP12')
 wait("document.querySelector('#plot_title').textContent.includes('200 selected genes')")
 page.get_by_role('tab',name='Expression',exact=True).click();wait("document.getElementById('dotplot').layout.meta?.genes?.length===200")
 gid=page.evaluate("document.getElementById('dotplot').layout.meta.genes[0]")
 page.evaluate("document.getElementById('dotplot').emit('plotly_click',{points:[{y:0}]})")
 page.wait_for_function("id=>Shiny.shinyapp.$inputValues.gene===id",arg=gid)
 assert state('sel_kind')=='set'
 check('typing a gene replaces active set; selecting a set replaces active gene; heatmap inspection retains set')
 for width,height in [(1366,768),(1920,1080)]:
  page.set_viewport_size({'width':width,'height':height});page.wait_for_timeout(700)
  assert page.evaluate('document.documentElement.scrollWidth<=innerWidth')
  page.screenshot(path=str(OUT/f'bugfix-expression-{width}.png'))
 check('laptop/desktop layout fits without horizontal page overflow')
 b.close()
