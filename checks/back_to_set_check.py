"""After "Plot this gene" the user can return to the gene set ("Back to set"), and the gene stays inspected.

Usage: python checks/back_to_set_check.py <port>     (app running: python -m shiny run --port <port> app/app.py)

Drives a headless Shiny session over the websocket protocol. Like a browser, it applies the
server's input updates (update_radio_buttons etc.) and sends them back as input changes.
Sequence (the reported bug): set active -> inspect a gene -> Plot this gene -> Cells tab ->
the Back to set button is there -> press it -> set active again, gene still inspected.
"""
import json
import re
import sys
import time

from websockets.sync.client import connect

PORT = int(sys.argv[1])
NIN = "MtrunA17_Chr5g0448621"
SET = ["MtrunA17_Chr5g0448621", "MtrunA17_Chr7g0253431", "MtrunA17_Chr3g0089541", "MtrunA17_Chr1g0182991"]

state = {
    "main_tabs": "expression", "sel_kind": "gene", "details_kind": "gene", "fate_target": "active",
    "gene_query": "", "gene": "", "geneset_text": "\n".join(SET), "geneset_preset": "",
    "colour_by": "annotation", "annotation": "ch5_tissue", "expr_max": 0.95, "zero_display": "dim",
    "facet": "", "facet_context": True, "sync_facets": True, "max_points": "20000", "point_size": 1.6,
    "opacity": 0.75, "gene_map_dim": "2", "dot_rows": "ch5_tissue", "dot_split": False,
    "fate_measure": "cpm", "fate_platform": "pooled", "expr_view": "dot", "heat_by_condition": False,
    "heat_scale": "z", "heat_k": 6, "heat_clustering": True, "heat_cluster": "", "stored_custom_sets": "{}",
    ".clientdata_url_protocol": "http:", ".clientdata_url_hostname": "127.0.0.1",
    ".clientdata_url_port": str(PORT), ".clientdata_url_pathname": "/", ".clientdata_url_search": "",
    ".clientdata_url_hash_initial": "", ".clientdata_url_hash": "", ".clientdata_pixelratio": 1,
    ".clientdata_singletons": "", ".clientdata_allowDataUriScheme": True,
    # Shiny only sends outputs the browser reports as visible
    ".clientdata_output_selection_summary_hidden": False, ".clientdata_output_back_to_set_ui_hidden": False,
}
outputs: dict = {}
errors: dict = {}
clicks: dict = {}


def pump(ws, seconds=60.0):
    """Receive until the server is idle; apply input updates as the browser would, and echo them back.

    Shiny wraps every update cycle in {"busy": "busy"} ... {"busy": "idle"}. An echoed input
    starts another cycle, so this waits for an idle with nothing echoed after it, then 1 s of quiet.
    """
    deadline = time.monotonic() + seconds
    start = time.monotonic()
    settled_at = None
    seen = False
    while time.monotonic() < deadline:
        if not seen and time.monotonic() - start > 5:
            return                       # no reaction at all (e.g. a button that does not exist)
        if settled_at is not None and time.monotonic() - settled_at > 1.0:
            return
        try:
            raw = ws.recv(timeout=0.2)
        except TimeoutError:
            continue
        if isinstance(raw, bytes):
            continue
        if raw[:1].isdigit() and "|" in raw[:12]:
            raw = raw.split("|", 1)[1]
        try:
            msg = json.loads(raw)
        except ValueError:
            continue
        seen = True
        outputs.update(msg.get("values") or {})
        for k, v in (msg.get("errors") or {}).items():
            errors[k] = v.get("message") if isinstance(v, dict) else str(v)
        echo = {}
        for im in msg.get("inputMessages") or []:
            value = im["message"].get("value")
            if isinstance(value, list) and len(value) == 1:      # update_select sends a list
                value = value[0]
            if "value" in im["message"] and state.get(im["id"]) != value:
                state[im["id"]] = echo[im["id"]] = value
        if echo:
            ws.send(json.dumps({"method": "update", "data": echo}))
            settled_at = None
        elif msg.get("busy") == "idle":
            settled_at = time.monotonic()
        elif msg.get("busy") == "busy":
            settled_at = None
    print("      (server did not become idle within the time limit)")


def send(ws, **values):
    state.update(values)
    ws.send(json.dumps({"method": "update", "data": values}))
    pump(ws)


def click(ws, button):
    clicks[button] = clicks.get(button, 0) + 1
    send(ws, **{button: clicks[button]})


def html(name):
    v = outputs.get(name)
    return (v.get("html") if isinstance(v, dict) else v) or ""


def check(cond, text):
    print(("PASS  " if cond else "FAIL  ") + text, flush=True)
    if not cond:
        print("      state:", {k: state[k] for k in ("main_tabs", "sel_kind", "details_kind", "colour_by", "gene", "gene_query")})
        print("      selection_summary:", re.sub("<[^>]+>", " ", html("selection_summary")).split())
        global failed
        failed = True


failed = False
with connect(f"ws://127.0.0.1:{PORT}/websocket/", max_size=2 ** 28, open_timeout=60) as ws:
    ws.send(json.dumps({"method": "init", "data": state}))
    pump(ws, 60)

    click(ws, "use_set")
    check(state["sel_kind"] == "set" and state["details_kind"] == "set", "Use set: the set is active")
    check("back_to_set" not in html("back_to_set_ui"), "no Back to set button while the set is active")

    send(ws, gene_query=NIN, gene=NIN, details_kind="gene")          # inspect a gene (heatmap click)
    check(state["sel_kind"] == "set" and "Inspecting:" in html("selection_summary"), "gene inspected, set still active")

    click(ws, "use_gene")                                             # "Plot this gene"
    check(state["sel_kind"] == "gene", "Plot this gene: the gene is the active selection")
    send(ws, main_tabs="cell_map")                                    # swap to the cell map
    button = html("back_to_set_ui")
    check('id="back_to_set"' in button and f"{len(SET)} genes" in button,
          f"Back to set button shown: {re.sub('<[^>]+>', ' ', button).split()[:12]}")

    click(ws, "back_to_set")
    check(state["sel_kind"] == "set" and state["details_kind"] == "set" and state["colour_by"] == "selection",
          "Back to set: the set is active again and colours the map")
    check(state["gene"] == NIN and "Inspecting:" in html("selection_summary"),
          "the gene stays inspected, so 'Details for: Inspected gene' works again")
    check("back_to_set" not in html("back_to_set_ui"), "button gone again while the set is active")

    send(ws, geneset_text="")                                         # no set loaded
    click(ws, "use_gene")
    check(state["sel_kind"] == "gene" and "back_to_set" not in html("back_to_set_ui"),
          "no Back to set button when the Gene sets box is empty")

check(not errors, f"no output errors {errors if errors else ''}")
print("RESULT", "FAILED" if failed else "PASSED")
sys.exit(1 if failed else 0)
