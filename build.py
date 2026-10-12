"""Build a portable, dependency-free HTML visualization from the source parquet."""
import base64, json, collections
from pathlib import Path
import numpy as np
import pyarrow.parquet as pq
rows=pq.read_table('figure2_umap_ctxlen_metadata.parquet').to_pylist()
frames=sorted(set(r['ctxlen'] for r in rows))
groups={f:{r['cell_id']:r for r in rows if r['ctxlen']==f} for f in frames}
assert all(len(g)==15750 for g in groups.values()), 'Duplicate or missing cells'
ids=sorted(groups[frames[0]])
assert all(set(g)==set(ids) for g in groups.values())
classes=sorted(set(r['AIT33_class'] for r in rows))
subs=sorted(set(r['AIT33_subclass'] for r in rows))
base=groups[frames[0]]
for g in groups.values():
    assert all((g[i]['AIT33_class'],g[i]['AIT33_subclass'])==(base[i]['AIT33_class'],base[i]['AIT33_subclass']) for i in ids)
assert all(base[i]['rna_umi_counts'] > 0 for i in ids), 'Log coloring requires positive counts'
assert all(g[i]['rna_umi_counts']==base[i]['rna_umi_counts'] for g in groups.values() for i in ids), 'RNA counts vary across frames'
coords=np.array([[[groups[f][i]['umap_1'],groups[f][i]['umap_2']] for i in ids] for f in frames],dtype='<f4')
assert np.isfinite(coords).all()
data={'frames':frames,'classes':classes,'subclasses':subs,'parents':[classes.index(next(r['AIT33_class'] for r in rows if r['AIT33_subclass']==s)) for s in subs], 'cells':[[i,classes.index(base[i]['AIT33_class']),subs.index(base[i]['AIT33_subclass']),base[i]['roi'],base[i]['donor_name'],base[i]['rna_umi_counts']] for i in ids], 'coords':base64.b64encode(coords.tobytes()).decode(), 'bounds':[float(coords[:,:,0].min()),float(coords[:,:,0].max()),float(coords[:,:,1].min()),float(coords[:,:,1].max())]}
template=Path('src/template.html').read_text()
Path('index.html').write_text(template.replace('__DATA__',json.dumps(data,separators=(',',':')).replace('</','<\\/')))
print(f'Built index.html: {Path("index.html").stat().st_size/1e6:.2f} MB; {len(ids):,} cells × {len(frames)} frames')
