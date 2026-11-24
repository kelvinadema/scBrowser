#!/bin/bash
cd "$(dirname "$0")/.."   # this scBrowser checkout
PY="${SCBROWSER_PYBIN:-${CONDA_PREFIX:?activate the scBrowser conda env or set SCBROWSER_PYBIN}/bin}"
$PY/shiny run app/app.py --host 127.0.0.1 --port 54112 > _checks_2026-11-19/latency_server.log 2>&1 &
P=$!
until curl -s -o /dev/null http://127.0.0.1:54112/ 2>/dev/null; do sleep 0.2; done
$PY/python -u _checks_2026-11-19/latency_probe.py 54112
echo "--- second session (warm caches)"
$PY/python -u _checks_2026-11-19/latency_probe.py 54112
kill $P
echo "tracebacks: $(grep -ic traceback _checks_2026-11-19/latency_server.log)"
