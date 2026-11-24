"""Browser-like session test for gene-set switching, custom sets and heatmap row clicks.

Echoes the server's input updates back (as a browser does), so chains like
'use cluster -> dropdown changes -> text box reloads' are exercised end to end.
Usage: python custom_sets_test.py <port>
"""
import asyncio
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import websockets

PORT = int(sys.argv[1])
VD = str(Path(__file__).resolve().parents[1] / 'viewer_data')   # this checkout's viewer_data
samples = json.load(open(f'{VD}/samples.json'))
groups = []
for s in samples:
    base = s["condition"].split(" 0")[0].strip()
    if not groups or groups[-1][1] != base:
        groups.append((f"samples_{len(groups)}", base, []))
    groups[-1][2].append(s["sample"])
z = np.load(f'{VD}/module_scores.npz')
genes = pd.read_parquet(f'{VD}/genes.parquet', columns=['gene_id']).gene_id
GEP09 = set(genes.iloc[z['GEP09__rows']])

state = {
    "sel_kind": "set", "gene_query": "", "geneset_text": "", "geneset_preset": "", "custom_name": "",
    "colour_by": "annotation", "annotation": "ch5_tissue", "expr_max": 0.95, "zero_display": "dim",
    "facet": "", "sync_facets": True, "max_points": "50000", "point_size": 1.6, "opacity": 0.75,
    "gene_map_dim": "2", "dot_rows": "ch5_tissue", "dot_split": False, "fate_measure": "cpm",
    "fate_platform": "pooled", "dissoc_genes": "", "expr_view": "heat", "heat_by_condition": False,
    "heat_scale": "z", "heat_k": 6, "heat_cluster": "", "main_tabs": "expression",
    **{gid: m for gid, _, m in groups},
    ".clientdata_url_protocol": "http:", ".clientdata_url_hostname": "127.0.0.1",
    ".clientdata_url_port": str(PORT), ".clientdata_url_pathname": "/", ".clientdata_url_search": "",
    ".clientdata_url_hash_initial": "", ".clientdata_url_hash": "", ".clientdata_pixelratio": 1,
    ".clientdata_singletons": "", ".clientdata_allowDataUriScheme": True,
}
log = {"stored": None, "notifications": [], "choices": None, "last_heat_genes": None, "last_heat_shapes": None}


async def pump(ws, seconds=15):
    """Read messages; echo input updates back like a browser would."""
    try:
        while True:
            raw = await asyncio.wait_for(ws.recv(), timeout=seconds)
            if isinstance(raw, bytes):
                continue
            raw = raw.split("|", 1)[-1] if raw[:1].isdigit() and "|" in raw[:12] else raw
            try:
                m = json.loads(raw)
            except json.JSONDecodeError:
                continue
            echo = {}
            for im in m.get("inputMessages") or []:
                iid, msg = im.get("id"), im.get("message", {})
                if "options" in msg:
                    log["choices" if iid == "geneset_preset" else iid + "_options"] = msg["options"]
                if "value" in msg and msg["value"] is not None:
                    v = msg["value"]
                    if isinstance(v, list) and iid in ("geneset_preset", "heat_cluster", "sel_kind"):
                        v = v[0] if v else ""          # single selects report a string, as in a browser
                    state[iid] = v; echo[iid] = v
            for k, v in (m.get("custom") or {}).items():
                if k == "store-custom-sets":
                    log["stored"] = v["sets"]
                if k == "dissoc-plot" and v.get("target") == "dotplot" and v.get("figure"):
                    meta = v["figure"]["layout"].get("meta") or {}
                    log["last_heat_genes"] = meta.get("genes")
                    log["last_heat_shapes"] = v["figure"]["layout"].get("shapes") or []
            if m.get("notification"):
                log["notifications"].append(str(m["notification"])[:120])
            if echo:
                await ws.send(json.dumps({"method": "update", "data": echo}))
    except asyncio.TimeoutError:
        pass


async def send(ws, upd, seconds=15):
    state.update(upd)
    await ws.send(json.dumps({"method": "update", "data": upd}))
    await pump(ws, seconds)


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""))


async def main():
    async with websockets.connect(f"ws://127.0.0.1:{PORT}/websocket/", max_size=2 ** 28) as ws:
        await ws.send(json.dumps({"method": "init", "data": {**state, "stored_custom_sets": "{}"}}))
        await pump(ws, 25)

        await send(ws, {"geneset_preset": "GEP09"}, 40)
        typed = set(state["geneset_text"].split())
        check("GEP09 preset loads its 200 genes", typed == GEP09, f"{len(typed)} genes")
        check("heatmap carries row gene ids", bool(log["last_heat_genes"]) and set(log["last_heat_genes"]) <= GEP09,
              f"{len(log['last_heat_genes'] or [])} rows")

        await send(ws, {"heat_cluster": "C2"}, 5)
        await send(ws, {"heat_use_cluster": 1}, 30)
        stored = log["stored"] or {}
        cname = next((n for n in stored if n.startswith("GEP09") and n.endswith("C2")), None)
        check("cluster saved as a custom set named after its parent", cname is not None, str(list(stored)))
        check("dropdown now points at the saved cluster", state["geneset_preset"] == f"custom:{cname}",
              state["geneset_preset"])
        typed = set(state["geneset_text"].split())
        check("text box holds the cluster genes", cname and typed == set(stored[cname]), f"{len(typed)} genes")

        await send(ws, {"geneset_preset": "GEP09"}, 40)
        typed = set(state["geneset_text"].split())
        check("re-selecting GEP09 reloads GEP09 (the reported bug)", typed == GEP09, f"{len(typed)} genes")

        await send(ws, {"geneset_text": state["geneset_text"] + "\nMtNIN"}, 10)
        check("editing the box by hand detaches it from the preset", state["geneset_preset"] == "",
              repr(state["geneset_preset"]))
        await send(ws, {"geneset_preset": "GEP09"}, 30)
        check("...after which GEP09 can be chosen again", set(state["geneset_text"].split()) == GEP09)

        gene = (log["last_heat_genes"] or [None])[0]
        await send(ws, {"heat_gene_click": gene}, 20)
        check("heatmap row click selects that gene", state["sel_kind"] == "gene" and state["gene_query"] == gene,
              f"{state['sel_kind']} {state['gene_query']}")
        check("the set is still in its box", set(state["geneset_text"].split()) == GEP09)
        await send(ws, {"gene": gene}, 15)   # a browser's gene picker reports its (re-rendered) choice
        check("heatmap stays, selected row outlined",
              any(sh.get("line", {}).get("color") == "#e0a800" for sh in (log["last_heat_shapes"] or [])),
              f"{len(log['last_heat_genes'] or [])} rows")

        await send(ws, {"sel_kind": "set"}, 10)
        await send(ws, {"custom_name": "my GEP09 copy"}, 2)
        await send(ws, {"custom_save": 1}, 15)
        check("named save stored", "my GEP09 copy" in (log["stored"] or {}), str(list(log["stored"] or {})))
        groups_ = [o for o in (log["choices"] or []) if isinstance(o, dict)]
        check("custom sets appear in the dropdown", "Your custom sets" in json.dumps(log["choices"]))
        saved = log["stored"]

    # reconnect with what the browser stored: the sets must come back
    log["choices"] = None
    async with websockets.connect(f"ws://127.0.0.1:{PORT}/websocket/", max_size=2 ** 28) as ws:
        await ws.send(json.dumps({"method": "init", "data": {**state, "stored_custom_sets": json.dumps(saved)}}))
        await pump(ws, 25)
        check("custom sets restored after reload", "my GEP09 copy" in json.dumps(log["choices"]))
    print("notifications:", log["notifications"][:4])


asyncio.run(main())
