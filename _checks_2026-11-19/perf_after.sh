#!/bin/bash
cd "$(dirname "$0")/.."   # this scBrowser checkout
PY="${SCBROWSER_PYBIN:-${CONDA_PREFIX:?activate the scBrowser conda env or set SCBROWSER_PYBIN}/bin}"
T0=$(date +%s.%N)
$PY/shiny run app/app.py --host 127.0.0.1 --port 54113 > _checks_2026-11-19/perf_after_server.log 2>&1 &
P=$!
until curl -s -o /dev/null http://127.0.0.1:54113/ 2>/dev/null; do sleep 0.2; done
echo "launch -> first page $(echo "$(date +%s.%N) - $T0" | bc) s"
$PY/python -u _checks_2026-11-19/latency_probe.py 54113
echo "--- custom sets scenario"
$PY/python -u _checks_2026-11-19/custom_sets_test.py 54113 | grep -E "PASS|FAIL"
kill $P
echo "tracebacks: $(grep -ic traceback _checks_2026-11-19/perf_after_server.log)"
