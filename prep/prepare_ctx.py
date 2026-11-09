"""Build ctx_umap.parquet: joint context UMAP (500-8,000) with AIT33 labels and colours.

Run on the analysis host (paths default to the house layout there):

    /opt/dlami/nvme/envs/bloom-figures/bin/python prep/prepare_ctx.py

Inputs
- bloom-figures data/formatted/figure2_context_umap_coords.csv: one UMAP fitted on the stack
  of all five context embeddings (28_figure2_context.py; n_neighbors=10, min_dist=0.6,
  spread=1.0, euclidean, seed 0). Rows are keyed by (context, row_in_dataset).
- hyacinth-v1 eval_balanced_50/allen_multiome metadata.h5ad: cell_id (obs index) and the
  per-cell AIT33_class / AIT33_subclass. row_in_dataset is the positional row in this file.
- bloom-umaps data/ait33_subclass_taxonomy.csv: class_color / subclass_color.

Colours follow bloom-umaps: a cell is coloured by its OWN AIT33_class (metadata assigns class
and subclass independently, so a few cells carry a class that is not their subclass's parent).
`color` is the default display colour, i.e. the class colour.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

HOME = Path("/home/ubuntu/code")
DEF_COORDS = HOME / "bloom-figures/data/formatted/figure2_context_umap_coords.csv"
DEF_META = Path(
    "/mnt/bloom-output/datasets/hyacinth-v1/eval_balanced_50/allen_multiome/metadata.h5ad"
)
DEF_TAX = HOME / "bloom-umaps/data/ait33_subclass_taxonomy.csv"
CONTEXTS = [500, 1000, 2000, 4000, 8000]
COLUMNS = [
    "cell_id", "row_in_dataset", "ctxlen", "umap_1", "umap_2", "color",
    "AIT33_class", "AIT33_subclass", "class_color", "subclass_color",
]


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--coords", type=Path, default=DEF_COORDS)
    ap.add_argument("--metadata", type=Path, default=DEF_META)
    ap.add_argument("--taxonomy", type=Path, default=DEF_TAX)
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parents[1] / "ctx_umap.parquet")
    a = ap.parse_args()

    xy = pd.read_csv(a.coords)
    obs = ad.read_h5ad(a.metadata, backed="r").obs[["AIT33_class", "AIT33_subclass"]].copy()
    obs = obs.astype(str)
    obs["cell_id"] = obs.index.astype(str)
    obs = obs.reset_index(drop=True)
    assert obs.cell_id.is_unique, "duplicate cell_id in metadata"

    assert sorted(xy.context.unique()) == CONTEXTS, sorted(xy.context.unique())
    for n, g in xy.groupby("context"):
        assert len(g) == len(obs) and g.row_in_dataset.is_unique, f"context {n}: rows"
        assert set(g.row_in_dataset) == set(range(len(obs))), f"context {n}: row set"
    assert np.isfinite(xy[["umap_1", "umap_2"]].to_numpy()).all()

    df = xy.merge(obs, left_on="row_in_dataset", right_index=True, how="left",
                  suffixes=("_coords", ""), validate="many_to_one")
    bad = (df.AIT33_subclass_coords != df.AIT33_subclass).sum()
    assert bad == 0, f"{bad} rows: subclass in coords csv != metadata"

    tax = pd.read_csv(a.taxonomy)
    ccol = dict(zip(tax.class_id_label, tax.class_color))
    scol = dict(zip(tax.subclass_id_label, tax.subclass_color))
    miss_c = sorted(set(df.AIT33_class) - set(ccol))
    miss_s = sorted(set(df.AIT33_subclass) - set(scol))
    assert not miss_c and not miss_s, (miss_c, miss_s)
    df["class_color"] = df.AIT33_class.map(ccol).str.upper()
    df["subclass_color"] = df.AIT33_subclass.map(scol).str.upper()
    df["color"] = df.class_color
    df = df.rename(columns={"context": "ctxlen"})
    df["ctxlen"] = df.ctxlen.astype("int16")
    df["row_in_dataset"] = df.row_in_dataset.astype("int32")
    df[["umap_1", "umap_2"]] = df[["umap_1", "umap_2"]].astype("float32")
    for col in ["color", "AIT33_class", "AIT33_subclass", "class_color", "subclass_color"]:
        df[col] = df[col].astype("category")
    df = df[COLUMNS].sort_values(["ctxlen", "row_in_dataset"]).reset_index(drop=True)

    # identity / labels must not change across contexts
    piv = df.pivot(index="cell_id", columns="ctxlen", values="AIT33_subclass")
    assert piv.notna().all().all() and (piv.nunique(axis=1) == 1).all()

    off = (
        df[df.ctxlen == CONTEXTS[0]]
        .assign(parent=lambda d: d.AIT33_subclass.map(dict(zip(tax.subclass_id_label, tax.class_id_label))))
        .query("parent != AIT33_class")
    )
    df.to_parquet(a.out, index=False)
    prov = dict(
        output=a.out.name, output_sha256=sha(a.out), n_rows=len(df), columns=COLUMNS,
        n_cells=int(df.cell_id.nunique()), contexts=CONTEXTS,
        n_classes=int(df.AIT33_class.nunique()), n_subclasses=int(df.AIT33_subclass.nunique()),
        n_cells_class_not_subclass_parent=int(len(off)),
        umap="single fit on the stack of all five contexts; umap-learn n_neighbors=10, "
             "min_dist=0.6, spread=1.0, euclidean, random_state=0 (bloom-figures "
             "scripts/prepare/28_figure2_context.py)",
        color="color == class_color (cell's own AIT33_class); colours from ait33_subclass_taxonomy.csv",
        sources=[dict(source=str(p), sha256=sha(p)) for p in (a.coords, a.taxonomy)]
        + [dict(source=str(a.metadata), size_bytes=a.metadata.stat().st_size)],
    )
    a.out.with_suffix(".prov.json").write_text(json.dumps(prov, indent=2) + "\n")
    print(json.dumps({k: v for k, v in prov.items() if k != "sources"}, indent=1))


if __name__ == "__main__":
    main()
