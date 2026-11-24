from __future__ import annotations

import functools
import json
import threading
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


# The reference view is drawn mirrored (view_mirror = TRUE); see plotting.CAMERA.
DISPLAY_MIRROR_UMAP1 = True


class ViewerData:
    def __init__(self, directory: Path):
        self.dir = Path(directory)
        if not (self.dir / "meta.json").exists():
            raise FileNotFoundError(
                f"No viewer data in {self.dir}. Run build/build_viewer_data.py first."
            )

        self.meta = json.loads((self.dir / "meta.json").read_text())
        self.cells = pd.read_parquet(self.dir / "cells.parquet")
        self.genes = pd.read_parquet(self.dir / "genes.parquet")

        self.n_cells = int(self.meta["n_cells"])
        self.n_genes = int(self.meta["n_genes"])
        self.labels = self.meta.get("column_labels", {})

        self.categorical = [c for c in self.meta["categorical_columns"] if c in self.cells]
        self.continuous = [
            c for c in self.meta["continuous_columns"] if f"num__{c}" in self.cells
        ]
        self.facet_columns = [
            c for c in self.meta.get("facet_columns", []) if c in self.cells
        ]
        self._add_design_columns()

        # Display coordinates. UMAP 1 is negated so a (non-mirroring) camera reproduces the
        # mirrored reference view; see plotting.CAMERA. cells.parquet is unchanged.
        self.coords = self.cells[["umap_1", "umap_2", "umap_3"]].to_numpy(np.float32, copy=True)
        if DISPLAY_MIRROR_UMAP1:
            self.coords[:, 0] *= -1

        # Gene lookup by id or symbol
        self._gene_row = {}
        for row, (gene_id, symbol) in enumerate(
            zip(self.genes["gene_id"], self.genes["symbol"])
        ):
            self._gene_row.setdefault(gene_id.lower(), row)
            if symbol:
                self._gene_row.setdefault(symbol.lower(), row)

        # Keep subsample stable
        self._shuffle = np.random.default_rng(0).permutation(self.n_cells)

        self._expr_path = self.dir / "expr.h5"
        self._expr: h5py.File | None = None
        self._expr_lock = threading.Lock()


    # condition (cell_design.parquet) -> treatment / sampled region
    TREATMENT = {"Root tip": "Untreated", "Mock": "Untreated", "NF 3 h": "NF", "BAP 3 h": "BAP"}
    TREATMENT_ORDER = ["Rhizobia", "NF", "BAP", "Untreated"]
    REGION_ORDER = ["Root tip", "Susceptible zone"]

    def _add_design_columns(self) -> None:
        """Treatment and sampled region, derived from the per-cell design table
        (build/build_cell_design.py), offered as annotations, facets and dot-plot groups."""
        path = self.dir / "cell_design.parquet"
        if not path.exists():
            return
        d = pd.read_parquet(path, columns=["cell_id", "condition"])
        assert (d.cell_id.to_numpy() == self.cells.cell_id.to_numpy()).all(), "cell_design out of step"
        cond = d.condition.astype(str)
        treat = cond.map(lambda c: "Rhizobia" if c.startswith("Rhizobia") else self.TREATMENT[c])
        region = cond.map(lambda c: "Root tip" if c == "Root tip" else "Susceptible zone")
        self.cells["ch5_treatment_group"] = pd.Categorical(treat.to_numpy(), categories=self.TREATMENT_ORDER)
        self.cells["ch5_region"] = pd.Categorical(region.to_numpy(), categories=self.REGION_ORDER)
        self.labels = {**self.labels, "ch5_treatment_group": "Treatment",
                       "ch5_region": "Region (root tip / susceptible zone)"}
        for col in ("ch5_treatment_group", "ch5_region"):
            self.categorical.append(col)
            self.facet_columns.append(col)

    def label_for(self, column: str) -> str:
        return self.labels.get(column, column)

    def annotation_choices(self) -> dict[str, str]:
        choices = {c: self.label_for(c) for c in self.categorical}
        choices.update({f"num__{c}": self.label_for(c) for c in self.continuous})
        return choices

    def facet_choices(self) -> dict[str, str]:
        usable = {
            c: f"{self.label_for(c)} ({self.cells[c].nunique()})"
            for c in self.facet_columns
            if self.cells[c].nunique() >= 2
        }
        return {"": "No faceting", **usable}


    def gene_search(self, query: str, limit: int = 100) -> pd.DataFrame:
        query = (query or "").strip().lower()
        genes = self.genes
        if not query:
            return genes.nlargest(limit, "n_cells_detected")

        label = genes["label"].str.lower()
        gene_id = genes["gene_id"].str.lower()
        description = genes["description"].str.lower()

        exact = label.eq(query) | gene_id.eq(query)
        prefix = label.str.startswith(query) | gene_id.str.startswith(query)
        anywhere = label.str.contains(query, regex=False) | gene_id.str.contains(
            query, regex=False
        )
        in_text = description.str.contains(query, regex=False)

        rank = np.select([exact, prefix, anywhere, in_text], [0, 1, 2, 3], default=9)
        hits = genes.assign(_rank=rank).query("_rank < 9")
        return hits.sort_values(["_rank", "n_cells_detected"], ascending=[True, False]).head(
            limit
        )

    def resolve_gene(self, token: str) -> int | None:
        if not token:
            return None
        return self._gene_row.get(str(token).strip().lower())

    @functools.lru_cache(maxsize=64)
    def expression(self, gene_row: int) -> np.ndarray:
        """Dense log1p-CP10K vector over all cells for one gene."""
        with self._expr_lock:
            indptr, indices, data = self._arrays()
            start, stop = int(indptr[gene_row]), int(indptr[gene_row + 1])

            values = np.zeros(self.n_cells, dtype=np.float32)
            if stop > start:
                values[indices[start:stop]] = data[start:stop].astype(np.float32)
        return values

    def _arrays(self):
        """(indptr, indices, data) of expr.h5. The builder writes contiguous, uncompressed
        datasets, so they are memory-mapped (fast random gene access, OS page cache);
        otherwise fall back to h5py datasets. Call with _expr_lock held."""
        if getattr(self, "_mapped", None) is None:
            with h5py.File(self._expr_path, "r") as h:
                specs = {k: (h[k].id.get_offset(), h[k].dtype, h[k].shape, h[k].chunks, h[k].compression)
                         for k in ("indptr", "indices", "data")}
                indptr = h["indptr"][:]
            if all(v[0] is not None and v[3] is None and v[4] is None for v in specs.values()):
                mm = {k: np.memmap(self._expr_path, mode="r", dtype=v[1], offset=v[0], shape=v[2])
                      for k, v in specs.items() if k != "indptr"}
                self._mapped = (indptr, mm["indices"], mm["data"])
            else:
                if self._expr is None:
                    self._expr = h5py.File(self._expr_path, "r")
                self._mapped = (indptr, self._expr["indices"], self._expr["data"])
        return self._mapped

    def group_means(self, gene_rows: list[int], codes: np.ndarray, n_groups: int,
                    keep: np.ndarray | None = None) -> np.ndarray:
        """Mean log1p-CP10K per cell group for each gene: (len(gene_rows), n_groups).

        codes: group index per cell (-1 = excluded); keep: optional cell mask (sample filter).
        Uses only the non-zero entries of each gene (zeros contribute nothing to the sum)."""
        codes = np.where(keep, codes, -1) if keep is not None else codes
        sizes = np.bincount(codes[codes >= 0], minlength=n_groups).astype(np.float64)
        out = np.zeros((len(gene_rows), n_groups), dtype=np.float64)
        with self._expr_lock:
            indptr, indices, data = self._arrays()
            for i, row in enumerate(gene_rows):
                start, stop = int(indptr[row]), int(indptr[row + 1])
                if stop == start:
                    continue
                c = codes[indices[start:stop]]
                ok = c >= 0
                out[i] = np.bincount(c[ok], weights=data[start:stop][ok].astype(np.float64),
                                     minlength=n_groups)
        with np.errstate(invalid="ignore", divide="ignore"):
            return out / sizes

    def expression_mean(self, gene_rows: list[int]) -> np.ndarray:
        """Mean log1p-CP10K over several genes for every cell (gene-set scores).

        Reads the genes straight from expr.h5 (bypasses the per-gene LRU cache, which
        would be flushed by hundreds of control genes)."""
        total = np.zeros(self.n_cells, dtype=np.float64)
        if not gene_rows:
            return total.astype(np.float32)
        with self._expr_lock:
            indptr, indices, data = self._arrays()
            for row in sorted(gene_rows):
                start, stop = int(indptr[row]), int(indptr[row + 1])
                if stop > start:
                    total[indices[start:stop]] += data[start:stop].astype(np.float64)
        return (total / len(gene_rows)).astype(np.float32)

    def subsample(self, mask: np.ndarray | None, max_points: int) -> np.ndarray:
        order = self._shuffle
        if mask is not None:
            order = order[mask[order]]
        if 0 < max_points < order.size:
            order = order[:max_points]
        return np.sort(order)
