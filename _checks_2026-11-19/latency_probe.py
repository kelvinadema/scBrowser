"""Interaction latency: time from sending an input change until the server goes quiet.

Usage: python latency_probe.py <port>
Reports, per step, time to the last server message and bytes received (browser-like echo of
input updates). 'quiet' = no message for QUIET seconds.
"""
import asyncio
import json
from pathlib import Path
import sys
import time

import websockets

PORT = int(sys.argv[1])
QUIET = 1.5
VD = str(Path(__file__).resolve().parents[1] / 'viewer_data')   # this checkout's viewer_data
samples = json.load(open(f'{VD}/samples.json'))
groups = []
for s in samples:
    base = s["condition"].split(" 0")[0].strip()
    if not groups or groups[-1][1] != base:
        groups.append((f"samples_{len(groups)}", base, []))
    groups[-1][2].append(s["sample"])

state = {
    "sel_kind": "gene", "gene_query": "", "geneset_text": "", "geneset_preset": "", "custom_name": "",
    "colour_by": "annotation", "annotation": "ch5_tissue", "expr_max": 0.95, "zero_display": "dim",
    "facet": "", "sync_facets": True, "max_points": "100000", "point_size": 1.6, "opacity": 0.75,
    "gene_map_dim": "2", "dot_rows": "ch5_tissue", "dot_split": False, "fate_measure": "cpm",
    "fate_platform": "pooled", "dissoc_genes": "", "expr_view": "dot", "heat_by_condition": False,
    "heat_scale": "z", "heat_k": 6, "heat_cluster": "", "main_tabs": "cell_map",
    "stored_custom_sets": "{}",
    **{gid: m for gid, _, m in groups},
    ".clientdata_url_protocol": "http:", ".clientdata_url_hostname": "127.0.0.1",
    ".clientdata_url_port": str(PORT), ".clientdata_url_pathname": "/", ".clientdata_url_search": "",
    ".clientdata_url_hash_initial": "", ".clientdata_url_hash": "", ".clientdata_pixelratio": 1,
    ".clientdata_singletons": "", ".clientdata_allowDataUriScheme": True,
}
STEPS = [
    ("select a gene (MtNIN)", {"gene": "MtrunA17_Chr5g0448621", "gene_query": "MtNIN", "colour_by": "selection"}),
    ("select another gene (MtLBD16 query)", {"gene_query": "MtLBD16"}),
    ("colour by annotation", {"colour_by": "annotation"}),
    ("facet by time", {"facet": "ch5_time_h"}),
    ("no facet", {"facet": ""}),
    ("open Expression tab", {"main_tabs": "expression"}),
    ("open Fate tab", {"main_tabs": "fate"}),
    ("open Gene map tab", {"main_tabs": "gene_umap"}),
    ("back to Cell map", {"main_tabs": "cell_map"}),
    ("load GEP09 set", {"sel_kind": "set", "geneset_preset": "GEP09"}),
]


async def pump(ws):
    t0, last, nbytes, n = time.perf_counter(), None, 0, 0
    try:
        while True:
            raw = await asyncio.wait_for(ws.recv(), timeout=QUIET)
            last = time.perf_counter()
            nbytes += len(raw); n += 1
            if isinstance(raw, bytes):
                continue
            body = raw.split("|", 1)[-1] if raw[:1].isdigit() and "|" in raw[:12] else raw
            try:
                m = json.loads(body)
            except json.JSONDecodeError:
                continue
            echo = {}
            for im in m.get("inputMessages") or []:
                iid, msg = im.get("id"), im.get("message", {})
                if "value" in msg and msg["value"] is not None:
                    v = msg["value"]
                    if isinstance(v, list) and iid in ("geneset_preset", "heat_cluster", "sel_kind", "colour_by"):
                        v = v[0] if v else ""
                    state[iid] = v; echo[iid] = v
            if echo:
                await ws.send(json.dumps({"method": "update", "data": echo}))
    except asyncio.TimeoutError:
        pass
    return (last - t0) if last else 0.0, nbytes, n


async def main():
    async with websockets.connect(f"ws://127.0.0.1:{PORT}/websocket/", max_size=2 ** 28) as ws:
        t = time.perf_counter()
        await ws.send(json.dumps({"method": "init", "data": state}))
        dt, nb, n = await pump(ws)
        print(f"{'session start':38s} {dt:6.2f}s  {nb / 1e6:6.2f} MB  {n:3d} msgs")
        for name, upd in STEPS:
            state.update(upd)
            await ws.send(json.dumps({"method": "update", "data": upd}))
            dt, nb, n = await pump(ws)
            print(f"{name:38s} {dt:6.2f}s  {nb / 1e6:6.2f} MB  {n:3d} msgs")


asyncio.run(main())
