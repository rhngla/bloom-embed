"""Build mod_umap.parquet: Figure 4i RNA-fit UMAP with projected ATAC, AIT33 labels and colours.

Run on the analysis host (paths default to the house layout there):

    /opt/dlami/nvme/envs/bloom-figures/bin/python prep/prepare_mod.py

Inputs
- bloom-figures data/formatted/figure4_modality_umap_coords.csv, variant "project"
  (41_build_modality_umaps.py; n_neighbors=10, min_dist=0.6, spread=1.0, euclidean, seed 0).
  The UMAP is fit on Bloom's RNA embedding alone (top4096, genes only) and Bloom's ATAC
  embedding of the same cells (resamp_atac_8192_d16, peaks only) is placed with
  umap.transform(). This is Figure 4i. Rows are keyed by (source, cell_row).
- hyacinth-v1 eval_balanced_50/allen_multiome metadata.h5ad: cell_id (obs index) and the
  per-cell AIT33_class / AIT33_subclass. cell_row is the positional row in this file.
- bloom-umaps data/ait33_subclass_taxonomy.csv: class_color / subclass_color.

Colours follow prepare_ctx.py: a cell is coloured by its OWN AIT33_class.
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
DEF_COORDS = HOME / "bloom-figures/data/formatted/figure4_modality_umap_coords.csv"
DEF_META = Path(
    "/mnt/bloom-output/datasets/hyacinth-v1/eval_balanced_50/allen_multiome/metadata.h5ad"
)
DEF_TAX = HOME / "bloom-umaps/data/ait33_subclass_taxonomy.csv"
VARIANT = "project"
MODALITIES = ["RNA", "ATAC"]  # display order: fitted modality first
ROLE = {"RNA": "fit", "ATAC": "projected"}
COLUMNS = [
    "cell_id", "row_in_dataset", "modality", "role", "umap_1", "umap_2",
    "AIT33_class", "AIT33_subclass", "class_color", "subclass_color",
]


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--coords", type=Path, default=DEF_COORDS)
    ap.add_argument("--metadata", type=Path, default=DEF_META)
    ap.add_argument("--taxonomy", type=Path, default=DEF_TAX)
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parents[1] / "mod_umap.parquet")
    a = ap.parse_args()

    xy = pd.read_csv(a.coords)
    xy = xy[xy.variant == VARIANT].copy()
    obs = ad.read_h5ad(a.metadata, backed="r").obs[["AIT33_class", "AIT33_subclass"]].copy()
    obs = obs.astype(str)
    obs["cell_id"] = obs.index.astype(str)
    obs = obs.reset_index(drop=True)
    assert obs.cell_id.is_unique, "duplicate cell_id in metadata"

    assert sorted(xy.source.unique()) == sorted(MODALITIES), sorted(xy.source.unique())
    for m, g in xy.groupby("source"):
        assert len(g) == len(obs) and g.cell_row.is_unique, f"{m}: rows"
        assert set(g.cell_row) == set(range(len(obs))), f"{m}: row set"
    assert np.isfinite(xy[["umap_1", "umap_2"]].to_numpy()).all()

    df = xy.merge(obs, left_on="cell_row", right_index=True, how="left", validate="many_to_one")
    bad = (df.subclass != df.AIT33_subclass).sum()
    assert bad == 0, f"{bad} rows: subclass in coords csv != metadata"

    tax = pd.read_csv(a.taxonomy)
    ccol = dict(zip(tax.class_id_label, tax.class_color))
    scol = dict(zip(tax.subclass_id_label, tax.subclass_color))
    miss_c = sorted(set(df.AIT33_class) - set(ccol))
    miss_s = sorted(set(df.AIT33_subclass) - set(scol))
    assert not miss_c and not miss_s, (miss_c, miss_s)
    df["class_color"] = df.AIT33_class.map(ccol).str.upper()
    df["subclass_color"] = df.AIT33_subclass.map(scol).str.upper()
    csv_col = xy.set_index(["source", "cell_row"]).subclass_color.str.upper()
    key = df.set_index(["source", "cell_row"])
    bad = (key.subclass_color != csv_col.reindex(key.index)).sum()
    assert bad == 0, f"{bad} rows: subclass colour in coords csv != taxonomy"

    df = df.rename(columns={"source": "modality", "cell_row": "row_in_dataset"})
    df["role"] = df.modality.map(ROLE)
    df["row_in_dataset"] = df.row_in_dataset.astype("int32")
    df[["umap_1", "umap_2"]] = df[["umap_1", "umap_2"]].astype("float32")
    df["_m"] = df.modality.map({m: k for k, m in enumerate(MODALITIES)})
    df = df.sort_values(["_m", "row_in_dataset"])[COLUMNS].reset_index(drop=True)
    for col in ["modality", "role", "AIT33_class", "AIT33_subclass", "class_color", "subclass_color"]:
        df[col] = df[col].astype("category")

    piv = df.pivot(index="cell_id", columns="modality", values="AIT33_subclass")
    assert piv.notna().all().all() and (piv.nunique(axis=1) == 1).all()

    off = (
        df[df.modality == MODALITIES[0]]
        .assign(parent=lambda d: d.AIT33_subclass.map(dict(zip(tax.subclass_id_label, tax.class_id_label))))
        .query("parent != AIT33_class")
    )
    df.to_parquet(a.out, index=False)
    prov = dict(
        output=a.out.name, output_sha256=sha(a.out), n_rows=len(df), columns=COLUMNS,
        n_cells=int(df.cell_id.nunique()), modalities=MODALITIES,
        n_classes=int(df.AIT33_class.nunique()), n_subclasses=int(df.AIT33_subclass.nunique()),
        n_cells_class_not_subclass_parent=int(len(off)),
        umap=f"variant={VARIANT}: umap-learn fit on Bloom's RNA embedding alone (top4096, genes "
             "only); Bloom's ATAC embedding of the same cells (resamp_atac_8192_d16, peaks only) "
             "placed with umap.transform(). n_neighbors=10, min_dist=0.6, spread=1.0, euclidean, "
             "random_state=0 (bloom-figures scripts/prepare/41_build_modality_umaps.py). Figure 4i.",
        color="class_color = cell's own AIT33_class; colours from ait33_subclass_taxonomy.csv",
        sources=[dict(source=str(p), sha256=sha(p)) for p in (a.coords, a.taxonomy)]
        + [dict(source=str(a.metadata), size_bytes=a.metadata.stat().st_size)],
    )
    a.out.with_suffix(".prov.json").write_text(json.dumps(prov, indent=2) + "\n")
    print(json.dumps({k: v for k, v in prov.items() if k != "sources"}, indent=1))


if __name__ == "__main__":
    main()
