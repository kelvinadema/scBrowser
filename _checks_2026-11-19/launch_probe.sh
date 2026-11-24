#!/bin/bash
cd "$(dirname "$0")/.."   # this scBrowser checkout
PY="${SCBROWSER_PYBIN:-${CONDA_PREFIX:?activate the scBrowser conda env or set SCBROWSER_PYBIN}/bin}"
for run in 1 2; do
  T0=$(date +%s.%N)
  $PY/shiny run app/app.py --host 127.0.0.1 --port 54111 > _checks_2026-11-19/launch_probe_$run.log 2>&1 &
  P=$!
  until curl -s -o /dev/null http://127.0.0.1:54111/ 2>/dev/null; do sleep 0.2; done
  T1=$(date +%s.%N)
  echo "run $run: launch -> first page $(echo "$T1 - $T0" | bc) s"
  kill $P; wait $P 2>/dev/null
done
$PY/python -X importtime -c "import sys; sys.path.insert(0,'app'); import app" 2> _checks_2026-11-19/importtime.txt
sort -t'|' -k2 -n -r _checks_2026-11-19/importtime.txt | head -12
