"""Build mod.html: self-contained, iframe-embeddable RNA -> ATAC UMAP page (Figure 4i).

Input: mod_umap.parquet (prep/prepare_mod.py). Template: src/mod_template.html.
"""
import base64
import collections
import json
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parent
rows = pq.read_table(ROOT / "mod_umap.parquet").to_pylist()
frames = ["RNA", "ATAC"]
assert sorted({r["modality"] for r in rows}) == sorted(frames)
assert all(r["role"] == ("fit" if r["modality"] == "RNA" else "projected") for r in rows)
by = {f: {r["cell_id"]: r for r in rows if r["modality"] == f} for f in frames}
assert all(len(by[f]) == sum(r["modality"] == f for r in rows) for f in frames), "duplicate cell_id"
ids = sorted(by[frames[0]])
assert all(set(by[f]) == set(ids) for f in frames), "cells differ across modalities"
base = by[frames[0]]
keys = ("row_in_dataset", "AIT33_class", "AIT33_subclass", "class_color", "subclass_color")
for f in frames:
    assert all(tuple(by[f][i][k] for k in keys) == tuple(base[i][k] for k in keys) for i in ids), f

classes = sorted({base[i]["AIT33_class"] for i in ids})
subs = sorted({base[i]["AIT33_subclass"] for i in ids})
ccol = {base[i]["AIT33_class"]: base[i]["class_color"] for i in ids}
scol = {base[i]["AIT33_subclass"]: base[i]["subclass_color"] for i in ids}
assert len({(base[i]["AIT33_class"], base[i]["class_color"]) for i in ids}) == len(classes)
assert len({(base[i]["AIT33_subclass"], base[i]["subclass_color"]) for i in ids}) == len(subs)
# legend grouping: each subclass under its most frequent class
par = collections.defaultdict(collections.Counter)
for i in ids:
    par[base[i]["AIT33_subclass"]][base[i]["AIT33_class"]] += 1
ci = {c: k for k, c in enumerate(classes)}
si = {s: k for k, s in enumerate(subs)}
coords = np.array(
    [[[by[f][i]["umap_1"], by[f][i]["umap_2"]] for i in ids] for f in frames], dtype="<f4"
)
assert np.isfinite(coords).all()
data = {
    "frames": frames,
    "classes": [{"name": c, "color": ccol[c]} for c in classes],
    "subclasses": [
        {"name": s, "color": scol[s], "parent": ci[par[s].most_common(1)[0][0]]} for s in subs
    ],
    "cells": [[i, ci[base[i]["AIT33_class"]], si[base[i]["AIT33_subclass"]]] for i in ids],
    "coords": base64.b64encode(coords.tobytes()).decode(),
}
tpl = (ROOT / "src/mod_template.html").read_text()
assert tpl.count("__DATA__") == 1
out = ROOT / "mod.html"
out.write_text(tpl.replace("__DATA__", json.dumps(data, separators=(",", ":")).replace("</", "<\\/")))
print(f"Built {out.name}: {out.stat().st_size / 1e6:.2f} MB; {len(ids):,} cells x {len(frames)} modalities; "
      f"{len(classes)} classes, {len(subs)} subclasses")
