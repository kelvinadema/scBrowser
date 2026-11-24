"""Drive a headless Shiny session over the websocket protocol and report what the server sends.

Usage: python session_test.py <port>
Sends an init with browser-like inputs, then a sequence of input updates (one gene, gene set,
fate colouring, sample filter, tab controls); collects custom messages and output errors.
"""
import asyncio
import json
from pathlib import Path
import sys

import websockets

PORT = int(sys.argv[1])
NIN = "MtrunA17_Chr5g0448621"
VD = str(Path(__file__).resolve().parents[1] / 'viewer_data')   # this checkout's viewer_data
ALL_SAMPLES = json.load(open(f'{VD}/samples.json'))
ALL = [s['sample'] for s in ALL_SAMPLES]
import numpy as _np
_z = _np.load(f'{VD}/module_scores.npz')
import pandas as _pd
_g = _pd.read_parquet(f'{VD}/genes.parquet', columns=['gene_id'])
GEP12 = "\n".join(_g.gene_id.iloc[_z['GEP12__rows']])
# sample checkbox groups, same grouping as the app (one group per base condition)
GROUPS = []
for _s in ALL_SAMPLES:
    _base = _s["condition"].split(" 0")[0].strip()
    if not GROUPS or GROUPS[-1][1] != _base:
        GROUPS.append((f"samples_{len(GROUPS)}", _base, []))
    GROUPS[-1][2].append(_s["sample"])


def samples(keep):
    return {gid: [x for x in members if x in keep] for gid, _, members in GROUPS}

BASE = {
    "sel_kind": "gene", "gene_query": "", "geneset_text": "", "geneset_preset": "",
    "colour_by": "annotation", "annotation": "ch5_tissue", "expr_max": 0.95, "zero_display": "dim",
    **samples(ALL), "facet": "", "sync_facets": True, "max_points": "50000", "point_size": 1.6,
    "opacity": 0.75, "gene_map_dim": "2", "dot_rows": "ch5_tissue", "dot_split": False,
    "fate_measure": "cpm", "fate_platform": "pooled", "dissoc_genes": "",
    "expr_view": "dot", "heat_by_condition": False, "heat_scale": "z", "heat_k": 6, "heat_cluster": "",
    ".clientdata_url_protocol": "http:", ".clientdata_url_hostname": "127.0.0.1",
    ".clientdata_url_port": str(PORT), ".clientdata_url_pathname": "/", ".clientdata_url_search": "",
    ".clientdata_url_hash_initial": "", ".clientdata_url_hash": "", ".clientdata_pixelratio": 1,
    ".clientdata_singletons": "", ".clientdata_allowDataUriScheme": True,
}

STEPS = [
    ("pick NIN (one gene)", {"gene": NIN, "gene_query": "MtNIN", "colour_by": "selection"}),
    ("Expression split by platform", {"dot_split": True}),
    ("Fate detection, BD only", {"fate_measure": "detection", "fate_platform": "BD"}),
    ("colour by nodule fate", {"colour_by": "fate"}),
    ("samples: rhizobia BD only", samples(["SMu007", "SMu009", "SMu010", "SMu011", "SMu012"])),
    ("facet by time", {"facet": "ch5_time_h"}),
    ("facet by treatment", {"facet": "ch5_treatment_group"}),
    ("facet by region, colour by treatment", {"facet": "ch5_region", "colour_by": "annotation",
                                              "annotation": "ch5_treatment_group"}),
    ("gene set: GEP12 preset", {"sel_kind": "set", "geneset_preset": "GEP12", "colour_by": "selection",
                                "facet": "", **samples(ALL)}),
    ("gene set: pasted 3 genes", {"geneset_text": "MtNIN MtNF-YA1 MtrunA17_Chr3g0134941"}),
    ("gene map 3D", {"gene_map_dim": "3"}),
    ("heatmap view, single gene (should explain)", {"sel_kind": "gene", "expr_view": "heat"}),
    ("heatmap: GEP12 pasted", {"sel_kind": "set", "geneset_text": GEP12}),
    ("heatmap: split by condition", {"heat_by_condition": True}),
    ("heatmap: k = 0", {"heat_k": 0}),
    ("heatmap: k = 4, raw", {"heat_k": 4, "heat_scale": "raw"}),
    ("heatmap: use cluster C2", {"heat_cluster": "C2", "heat_use_cluster": 1}),
]


def summarise(msgs):
    out = {}
    for m in msgs:
        if "custom" in m:
            for k, v in m["custom"].items():
                tgt = v.get("target") if isinstance(v, dict) else None
                key = f"{k}:{tgt}" if tgt else k
                fig = v.get("figure") if isinstance(v, dict) else None
                n = len(fig["data"]) if isinstance(fig, dict) else ("null" if tgt else "")
                out[key] = n
        if "errors" in m and m["errors"]:
            for k, v in m["errors"].items():
                out[f"ERROR {k}"] = (v.get("message") if isinstance(v, dict) else str(v))[:200]
        if "inputMessages" in m and m["inputMessages"]:
            for im in m["inputMessages"]:
                out[f"input->{im.get('id')}"] = str(im.get("message"))[:80]
        if "values" in m and m["values"]:
            for k in ("heat_note", "heat_cluster_picker"):
                if k in m["values"]:
                    v = m["values"][k]
                    out[f"output {k}"] = (v.get("html") if isinstance(v, dict) else str(v))[:160]
    return out


async def collect(ws, seconds):
    msgs = []
    try:
        while True:
            raw = await asyncio.wait_for(ws.recv(), timeout=seconds)
            if isinstance(raw, bytes):
                continue
            raw = raw.split("|", 1)[-1] if raw[:1].isdigit() and "|" in raw[:12] else raw
            try:
                msgs.append(json.loads(raw))
            except json.JSONDecodeError:
                pass
    except asyncio.TimeoutError:
        pass
    return msgs


async def main():
    async with websockets.connect(f"ws://127.0.0.1:{PORT}/websocket/", max_size=2 ** 28) as ws:
        await ws.send(json.dumps({"method": "init", "data": BASE}))
        print("== init");  [print("  ", k, v) for k, v in summarise(await collect(ws, 25)).items()]
        for name, upd in STEPS:
            await ws.send(json.dumps({"method": "update", "data": upd}))
            print("==", name); [print("  ", k, v) for k, v in summarise(await collect(ws, 20)).items()]


asyncio.run(main())
