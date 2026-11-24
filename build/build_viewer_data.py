#!/usr/bin/env python3
"""
Take the full h5ad object and write required files for browser, should only be run once and creates view_data
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


# Select columns to keep
CATEGORICAL_COLUMNS = [
    "ch5_tissue",
    "ch5_cortex_layer",
    "ch5_phase_calibrated",
    "ch5_nodule_like",
    "ch5_meristem_like",
    "leiden_scArches",
    "sample",
]

# Select columns for facet
FACET_COLUMNS = ["ch5_time_h"]

CONTINUOUS_COLUMNS: list[str] = []

# Harmonise PO terms and actual readable names
PO_HARMONISE = {"ch5_tissue", "annotation_combined"}

COLUMN_LABELS = {
    "ch5_tissue": "Cell type (Plant Ontology)",
    "ch5_cortex_layer": "Cortex layer (inner/outer)",
    "ch5_phase_calibrated": "Cell cycle phase",
    "ch5_nodule_like": "Nodule-like",
    "ch5_meristem_like": "Meristem-like",
    "leiden_scArches": "Leiden cluster (scArches)",
    "sample": "Sample",
    "ch5_time_h": "Timepoint (HPI)",
}

UMAP_KEY = "X_umap_3d_ch5"


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _decode(values) -> np.ndarray:
    """h5py gives back bytes for string, normalise to str."""
    out = np.asarray(values)
    if out.dtype.kind in ("S", "O"):
        return np.array(
            [v.decode() if isinstance(v, bytes) else str(v) for v in out], dtype=object
        )
    return out


def read_obs_column(obs: h5py.Group, key: str, limit: int | None = None):
    """Return an obs column as a pandas series"""
    node = obs[key]
    stop = slice(None) if limit is None else slice(0, limit)
    if isinstance(node, h5py.Group) and "categories" in node and "codes" in node:
        categories = _decode(node["categories"][:])
        codes = node["codes"][stop]
        return pd.Categorical.from_codes(codes, categories=list(categories))

    values = node[stop]
    if values.dtype.kind in ("S", "O"):
        return pd.Series(_decode(values))
    if values.dtype == bool:
        return pd.Series(values)
    return pd.Series(values)


def load_po_labels(path: Path) -> dict[str, str]:
    if not path.exists():
        log(f"WARNING: {path} not found")
        return {}
    table = pd.read_csv(path, sep="\t", dtype=str).dropna()
    return dict(zip(table["po_id"], table["label"]))


MISSING_LABEL = "not annotated"


def tidy_categorical(values) -> pd.Categorical:
    """Normalise obs columns"""
    if isinstance(values, pd.Categorical):
        series = pd.Series(values.astype(object))
    else:
        series = pd.Series(values)

    if series.dtype.kind in ("i", "u", "f"):
        formatted = series.map(
            lambda v: "" if pd.isna(v) else (f"{v:g}" if float(v) % 1 else f"{int(v)}")
        )
    else:
        formatted = series.astype(object).map(lambda v: "" if pd.isna(v) else str(v))

    formatted = formatted.replace({"": MISSING_LABEL, "nan": MISSING_LABEL,
                                   "NA": MISSING_LABEL, "None": MISSING_LABEL})

    categories = sorted(set(formatted) - {MISSING_LABEL})
    # Numeric should sort numerically (timepoints, leiden ids)
    try:
        categories.sort(key=float)
    except ValueError:
        pass
    if (formatted == MISSING_LABEL).any():
        categories.append(MISSING_LABEL)
    return pd.Categorical(formatted, categories=categories)


def harmonise_po(values: pd.Series, po_labels: dict[str, str]) -> pd.Series:
    """Make PO terms readable and map to annotation
    """
    name_to_po = {label: po_id for po_id, label in po_labels.items()}
    bare_id = re.compile(r"^PO:\d+$")
    already_labelled = re.compile(r"\(PO:\d+\)\s*$")

    def convert(value: str) -> str:
        value = value.strip()
        if already_labelled.search(value):
            return value
        if bare_id.match(value):
            label = po_labels.get(value)
            return f"{label} ({value})" if label else value
        plain = value.replace("_", " ").strip()
        if plain in name_to_po:
            return f"{plain} ({name_to_po[plain]})"
        return plain

    return values.map(convert)


def relabel_po(series: pd.Categorical, po_labels: dict[str, str]) -> pd.Categorical:
    """Rewrite 'PO:0020131' to 'lateral root cap (PO:0020131)' in categories."""
    if not po_labels:
        return series
    renamed = [
        f"{po_labels[c]} ({c})" if c in po_labels else c for c in series.categories
    ]
    if len(set(renamed)) != len(renamed):
        return series
    return series.rename_categories(renamed)


def build_cells(h5: h5py.File, n_cells: int, po_labels: dict[str, str]) -> pd.DataFrame:
    obs = h5["obs"]
    frame = pd.DataFrame(index=np.arange(n_cells))

    index_key = obs.attrs.get("_index", "_index")
    if isinstance(index_key, bytes):
        index_key = index_key.decode()
    if index_key in obs:
        frame["cell_id"] = _decode(obs[index_key][:n_cells])

    used_categorical, used_continuous, used_facet = [], [], []

    for key in CATEGORICAL_COLUMNS + FACET_COLUMNS:
        if key not in obs:
            log(f"  skip (absent): {key}")
            continue

        raw = read_obs_column(obs, key, n_cells)
        if key in PO_HARMONISE:
            as_text = pd.Series(np.asarray(raw).astype(str))
            before = as_text.nunique()
            as_text = harmonise_po(as_text, po_labels)
            log(f"  {key}: harmonised {before} -> {as_text.nunique()} levels")
            raw = as_text
        values = relabel_po(tidy_categorical(raw), po_labels)
        frame[key] = values
        (used_facet if key in FACET_COLUMNS else used_categorical).append(key)

    for key in CONTINUOUS_COLUMNS:
        if key not in obs:
            log(f"  skip (absent): {key}")
            continue
        values = read_obs_column(obs, key, n_cells)
        if isinstance(values, pd.Categorical):
            values = pd.to_numeric(values.astype(str), errors="coerce")
        frame[f"num__{key}"] = np.asarray(values, dtype=np.float32)
        used_continuous.append(key)

    coords = h5["obsm"][UMAP_KEY][:n_cells]
    frame["umap_1"] = coords[:, 0].astype(np.float32)
    frame["umap_2"] = coords[:, 1].astype(np.float32)
    frame["umap_3"] = coords[:, 2].astype(np.float32)

    frame.attrs["categorical"] = used_categorical
    frame.attrs["continuous"] = used_continuous
    frame.attrs["facet"] = used_facet
    return frame


def build_genes(h5: h5py.File) -> pd.DataFrame:
    var = h5["var"]
    index_key = var.attrs.get("_index", "_index")
    if isinstance(index_key, bytes):
        index_key = index_key.decode()

    genes = pd.DataFrame({"gene_id": _decode(var[index_key][:])})

    for source, target in [("ch5_symbol", "symbol"), ("ch5_description", "description")]:
        if source in var:
            values = read_obs_column(var, source)
            values = pd.Series(np.asarray(values).astype(str))
            genes[target] = values.replace({"NA": "", "nan": "", "None": ""}).values
        else:
            genes[target] = ""

    # What the search box shows and matches against.
    genes["label"] = np.where(genes["symbol"] != "", genes["symbol"], genes["gene_id"])
    return genes

def transpose_matrix(
    h5: h5py.File,
    out_path: Path,
    n_cells: int,
    n_genes: int,
    chunk_cells: int,
    compression: str | None,
    min_cells: int,
) -> tuple[np.ndarray, np.ndarray]:
    data_ds = h5["X/data"]
    indices_ds = h5["X/indices"]
    indptr = h5["X/indptr"][: n_cells + 1].astype(np.int64)
    nnz = int(indptr[-1])

    log(f"matrix: {n_cells:,} cells x {n_genes:,} genes, {nnz:,} nonzeros")

    # How many cells express gene
    log("pass 1/2: counting nonzeros per gene")
    gene_counts = np.zeros(n_genes, dtype=np.int64)
    step = 50_000_000
    for start in range(0, nnz, step):
        stop = min(start + step, nnz)
        block = indices_ds[start:stop]
        gene_counts += np.bincount(block, minlength=n_genes).astype(np.int64)
        log(f"  {stop / nnz:6.1%}")

    assert gene_counts.sum() == nnz, "nonzero count mismatch between passes"

    # Filter genes
    keep = gene_counts >= max(0, min_cells)
    n_kept = int(keep.sum())
    kept_nnz = int(gene_counts[keep].sum())
    log(
        f"gene filter (>=1 UMI in >={min_cells} cells): keeping {n_kept:,} of "
        f"{n_genes:,} genes, {kept_nnz:,} of {nnz:,} nonzeros "
        f"({kept_nnz / nnz:.4%})"
    )

    remap = np.full(n_genes, -1, dtype=np.int32)
    remap[keep] = np.arange(n_kept, dtype=np.int32)

    out_indptr = np.zeros(n_kept + 1, dtype=np.int64)
    np.cumsum(gene_counts[keep], out=out_indptr[1:])

    log(f"pass 2/2: transposing (allocating ~{kept_nnz * 6 / 1e9:.1f} GB)")
    out_indices = np.empty(kept_nnz, dtype=np.int32)
    out_data = np.empty(kept_nnz, dtype=np.float16)
    cursor = out_indptr[:-1].copy()

    for cell_start in range(0, n_cells, chunk_cells):
        cell_stop = min(cell_start + chunk_cells, n_cells)
        lo, hi = int(indptr[cell_start]), int(indptr[cell_stop])
        if hi == lo:
            continue

        counts = data_ds[lo:hi].astype(np.float32)
        genes_of_nz = indices_ds[lo:hi]

        # Per-cell library size then CP10K + log1p.
        row_lengths = np.diff(indptr[cell_start : cell_stop + 1])
        cell_of_nz = np.repeat(
            np.arange(cell_start, cell_stop, dtype=np.int32), row_lengths
        )
        totals = np.bincount(
            cell_of_nz - cell_start, weights=counts, minlength=cell_stop - cell_start
        ).astype(np.float32)
        totals[totals == 0] = 1.0
        # Library size uses every gene including to be dropped
        values = np.log1p(counts / totals[cell_of_nz - cell_start] * 1e4)

        mapped = remap[genes_of_nz]
        surviving = mapped >= 0
        if not surviving.all():
            mapped = mapped[surviving]
            cell_of_nz = cell_of_nz[surviving]
            values = values[surviving]
        if mapped.size == 0:
            continue

        order = np.argsort(mapped, kind="stable")
        genes_sorted = mapped[order]
        chunk_counts = np.bincount(genes_sorted, minlength=n_kept).astype(np.int64)
        chunk_starts = np.zeros(n_kept, dtype=np.int64)
        np.cumsum(chunk_counts[:-1], out=chunk_starts[1:])
        rank = np.arange(genes_sorted.size, dtype=np.int64) - chunk_starts[genes_sorted]

        destination = cursor[genes_sorted] + rank
        out_indices[destination] = cell_of_nz[order]
        out_data[destination] = values[order].astype(np.float16)
        cursor += chunk_counts

        log(f"  {cell_stop / n_cells:6.1%}")

    assert np.array_equal(cursor, out_indptr[1:]), "scatter did not fill every slot"

    log(f"writing {out_path}")
    kwargs = {"compression": compression} if compression else {}
    with h5py.File(out_path, "w") as out:
        out.create_dataset("indptr", data=out_indptr, dtype=np.int64)
        out.create_dataset("indices", data=out_indices, dtype=np.int32, **kwargs)
        out.create_dataset("data", data=out_data, dtype=np.float16, **kwargs)
        out.attrs["n_cells"] = n_cells
        out.attrs["n_genes"] = n_kept
        out.attrs["layout"] = "csc_gene_major"
        out.attrs["transform"] = "log1p(counts / cell_total * 1e4)"

    return gene_counts, keep


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path, help="output directory")
    parser.add_argument(
        "--po-terms",
        type=Path,
        default=Path(__file__).with_name("po_terms.tsv"),
        help="Plant Ontology id -> label mapping",
    )
    parser.add_argument("--chunk-cells", type=int, default=20_000)
    parser.add_argument(
        "--min-cells",
        type=int,
        default=10,
        help="drop genes with >=1 UMI in fewer than this many cells",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=0,
        help="Use for quick testing of small subsets",
    )
    parser.add_argument(
        "--compression",
        default="none",
        choices=["none", "lzf", "gzip"],
        help="Compress to make smaller, but slower",
    )
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    compression = None if args.compression == "none" else args.compression
    po_labels = load_po_labels(args.po_terms)

    started = time.time()
    with h5py.File(args.h5ad, "r") as h5:
        n_cells = int(h5["X/indptr"].shape[0]) - 1
        if args.max_cells:
            n_cells = min(n_cells, args.max_cells)
            log(f"building from the first {n_cells:,} cells only")
        genes = build_genes(h5)
        n_genes = len(genes)

        log("building cell metadata")
        cells = build_cells(h5, n_cells, po_labels)
        categorical = cells.attrs["categorical"]
        continuous = cells.attrs["continuous"]
        facet = cells.attrs["facet"]

        gene_counts, keep = transpose_matrix(
            h5, args.out / "expr.h5", n_cells, n_genes, args.chunk_cells,
            compression, args.min_cells,
        )

    genes["n_cells_detected"] = gene_counts
    genes["pct_cells_detected"] = (gene_counts / n_cells * 100).round(2)
    # Row order here must match expr.h5's gene axis exactly.
    genes = genes[keep].reset_index(drop=True)
    n_genes = len(genes)

    cells.to_parquet(args.out / "cells.parquet", index=False)
    genes.to_parquet(args.out / "genes.parquet", index=False)

    (args.out / "meta.json").write_text(
        json.dumps(
            {
                "source_h5ad": str(args.h5ad),
                "subset_build": bool(args.max_cells),
                "n_cells": n_cells,
                "n_genes": n_genes,
                "umap_key": UMAP_KEY,
                "expression_transform": "log1p(counts / cell_total * 1e4)",
                "expression_layer": "corrected_counts (= X)",
                "categorical_columns": categorical,
                "continuous_columns": continuous,
                "facet_columns": facet,
                "min_cells_detected": args.min_cells,
                "column_labels": COLUMN_LABELS,
                "built": time.strftime("%Y-%m-%d %H:%M:%S"),
            },
            indent=2,
        )
    )

    log(f"done in {(time.time() - started) / 60:.1f} min -> {args.out}")
    for path in sorted(args.out.iterdir()):
        log(f"  {path.name:16s} {path.stat().st_size / 1e9:8.2f} GB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
