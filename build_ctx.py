"""Build ctx.html: self-contained, iframe-embeddable context-length UMAP page.

Input: ctx_umap.parquet (prep/prepare_ctx.py). Template: src/ctx_template.html.
"""
import base64
import collections
import json
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parent
rows = pq.read_table(ROOT / "ctx_umap.parquet").to_pylist()
frames = sorted({r["ctxlen"] for r in rows})
assert frames == [500, 1000, 2000, 4000, 8000], frames
by = {f: {r["cell_id"]: r for r in rows if r["ctxlen"] == f} for f in frames}
assert all(len(by[f]) == sum(r["ctxlen"] == f for r in rows) for f in frames), "duplicate cell_id"
ids = sorted(by[frames[0]])
assert all(set(by[f]) == set(ids) for f in frames), "cells differ across contexts"
base = by[frames[0]]
keys = ("AIT33_class", "AIT33_subclass", "class_color", "subclass_color", "color")
for f in frames:
    assert all(tuple(by[f][i][k] for k in keys) == tuple(base[i][k] for k in keys) for i in ids), f
assert all(base[i]["color"] == base[i]["class_color"] for i in ids)

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
tpl = (ROOT / "src/ctx_template.html").read_text()
assert tpl.count("__DATA__") == 1
out = ROOT / "ctx.html"
out.write_text(tpl.replace("__DATA__", json.dumps(data, separators=(",", ":")).replace("</", "<\\/")))
print(f"Built {out.name}: {out.stat().st_size / 1e6:.2f} MB; {len(ids):,} cells x {len(frames)} contexts; "
      f"{len(classes)} classes, {len(subs)} subclasses")
