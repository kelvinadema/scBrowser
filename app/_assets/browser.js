document.addEventListener('DOMContentLoaded', () => {
  let panelTrigger = null, evidenceTrigger = null;
  const resize = () => document.querySelectorAll('.js-plotly-plot').forEach(p => {
    if (p.offsetWidth && p.offsetHeight && window.Plotly) Plotly.Plots.resize(p);
  });
  function showPanel(key, trigger) {
    const dock = document.getElementById('context-dock');
    dock.hidden = key === 'none';
    document.querySelector('.workspace').classList.toggle('dock-closed', dock.hidden);
    document.querySelectorAll('.context-panel').forEach(p => p.hidden = p.id !== 'panel-' + key);
    document.querySelectorAll('[data-panel]').forEach(b => b.setAttribute('aria-expanded', String(b.dataset.panel === key && key !== 'none')));
    if (key === 'none') { if (panelTrigger) panelTrigger.focus(); }
    else { panelTrigger = trigger || panelTrigger; const heading = document.querySelector('#panel-' + key + ' h2'); heading.tabIndex = -1; heading.focus({preventScroll:true}); dock.scrollTop = 0; }
    setTimeout(resize, 50);
  }
  document.addEventListener('click', e => {
    const panel = e.target.closest('[data-panel]');
    if (panel) showPanel(panel.dataset.panel, panel);
    const ev = e.target.closest('[data-evidence]');
    if (ev) { evidenceTrigger = ev; document.getElementById('evidence-' + ev.dataset.evidence).showModal(); resize(); }
    const close = e.target.closest('[data-close-evidence]');
    if (close) document.getElementById('evidence-' + close.dataset.closeEvidence).close();
    const full = e.target.closest('[data-fullscreen]');
    if (full) { const el = document.getElementById(full.dataset.fullscreen); if (document.fullscreenElement) document.exitFullscreen(); else if (el.requestFullscreen) el.requestFullscreen(); }
    if (e.target.closest('#open_fate')) document.getElementById('evidence-fate').close();
  });
  // Only a user's search gesture replaces the active set. Server updates used
  // for map/heatmap inspection deliberately keep the set active.
  function searchAction(e) {
    if (!e.isTrusted || !window.Shiny) return;
    Shiny.setInputValue('gene_search_action', Date.now(), {priority: 'event'});
  }
  document.getElementById('gene_query').addEventListener('input', searchAction);
  document.getElementById('gene_picker').addEventListener('change', searchAction);
  document.querySelectorAll('dialog').forEach(d => d.addEventListener('close', () => { if (evidenceTrigger) evidenceTrigger.focus(); resize(); }));
  document.addEventListener('keydown', e => { if (e.key === 'Escape' && !document.querySelector('dialog[open]')) showPanel('none'); });
  document.addEventListener('fullscreenchange', resize);
  document.addEventListener('shown.bs.tab', e => {
    const target = e.target.getAttribute('data-value');
    if (['gene_umap','fate'].includes(target) && !document.getElementById('panel-settings').hidden) showPanel('details');
    resize();
  });
  const ro = new ResizeObserver(() => requestAnimationFrame(resize));
  ro.observe(document.querySelector('.plot-card'));
  document.querySelectorAll('[data-group]').forEach(group => {
    const root = document.getElementById(group.dataset.group);
    const boxes = () => [...root.querySelectorAll('input[type=checkbox]')];
    function state() { const b = boxes(), n = b.filter(x => x.checked).length; group.checked = n === b.length; group.indeterminate = n > 0 && n < b.length; group.setAttribute('aria-checked',group.indeterminate ? 'mixed' : String(group.checked)); }
    group.addEventListener('change', () => { boxes().forEach(b => b.checked = group.checked); $(root).trigger('change'); state(); });
    root.addEventListener('change', state);
    $(root).on('change', state);
    $(document).on('shiny:message', () => setTimeout(state, 0));
    state();
  });
  document.getElementById('sample-search').addEventListener('input', e => {
    const q = e.target.value.toLowerCase().trim();
    document.querySelectorAll('.sample-group').forEach(g => {
      const groupMatch = g.querySelector('summary').textContent.toLowerCase().includes(q);
      let any = false;
      g.querySelectorAll('.form-check, .checkbox').forEach(row => { const match = groupMatch || row.textContent.toLowerCase().includes(q); row.hidden = !match; any ||= match; });
      g.hidden = !any; if (q && any) g.open = true;
    });
  });
  $(document).on('shiny:connected', () => {
    Shiny.addCustomMessageHandler('show-details', (msg) => showPanel('details'));

  });
});
