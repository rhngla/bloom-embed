"""Build Figure 4h by reusing the original page's layout, palette and interactions.

Never writes index.html or src/template.html. Guarded substitutions fail if the
shared template changes unexpectedly, rather than silently producing a bad page.
"""
from collections import Counter, defaultdict
import base64
import json
from pathlib import Path
import numpy as np
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parent
rows = pq.read_table(ROOT / 'figure4h_modality_umap_metadata.parquet').to_pylist()
modalities = ['RNA (genes)', 'ATAC (peaks)']
assert set(r['modality'] for r in rows) == set(modalities)
groups = []
for modality in modalities:
    subset = [r for r in rows if r['modality'] == modality]
    group = {r['cell_id']: r for r in subset}
    assert len(group) == len(subset), f'Duplicate cell IDs in {modality}'
    groups.append(group)
assert groups[0].keys() == groups[1].keys(), 'Unpaired cells'
ids = sorted(groups[0])
base = groups[0]
classes = sorted({r['AIT33_class'] for r in rows})
subclasses = sorted({r['AIT33_subclass'] for r in rows})
parent_counts = defaultdict(Counter)
for r in base.values():
    parent_counts[r['AIT33_subclass']][r['AIT33_class']] += 1
# Some subclass labels span multiple recorded classes. Use their modal parent
# for palette organization only; preserve each cell's original class metadata.
parents = {s: sorted(c, key=lambda k: (-c[k], k))[0] for s, c in parent_counts.items()}
for i in ids:
    assert all(base[i][k] == groups[1][i][k] for k in ['AIT33_class', 'AIT33_subclass'])
coords = np.array([[[g[i]['umap_1'], g[i]['umap_2']] for i in ids] for g in groups], dtype='<f4')
counts = [[g[i][key] for i in ids] for g, key in zip(groups, ['rna_umi_counts', 'atac_umi_counts'])]
assert np.isfinite(coords).all()
assert all(v is not None and v > 0 for frame in counts for v in frame), 'Counts must be positive for logarithmic colors'
data = {
    'frames': ['RNA', 'ATAC'], 'classes': classes, 'subclasses': subclasses,
    'parents': [classes.index(parents[s]) for s in subclasses],
    'cells': [[i, classes.index(base[i]['AIT33_class']), subclasses.index(base[i]['AIT33_subclass']), base[i]['roi'], base[i]['donor_name'], counts[0][j], counts[1][j]] for j, i in enumerate(ids)],
    'coords': base64.b64encode(coords.tobytes()).decode(),
    'bounds': [float(coords[:, :, 0].min()), float(coords[:, :, 0].max()), float(coords[:, :, 1].min()), float(coords[:, :, 1].max())],
}
html = (ROOT / 'src/template.html').read_text()
def replace(old, new, count=1):
    global html
    assert html.count(old) == count, f'Shared template changed: {old[:80]}'
    html = html.replace(old, new)

def replace_function(name, replacement):
    global html
    start = html.index(f'function {name}(')
    end = html.index('\nfunction ', start + 1)
    html = html[:start] + replacement.rstrip() + html[end:]

replace('Subclass shades progress in numerical order within their parent class.', 'Subclass shades follow numerical order within the most common parent class. Recorded per-cell class labels are preserved.')
replace('Context in motion', 'Modalities in motion', 2)
replace('One cell. Five perspectives.', 'One cell. Two modalities.')
replace('<strong>5</strong><span>Context lengths</span>', '<strong>2</strong><span>Modalities</span>')
replace('500 ctxlen · observed embedding', 'RNA · observed embedding')
replace('<span id="ctx">500</span><small>context length</small>', '<span id="ctx">RNA</span><small>modality blend</small>')
replace('max="4" step="0.001"', 'max="1" step="0.001"')
replace('aria-label="Context length"', 'aria-label="RNA to ATAC interpolation"')
a=html.index('<div class="ticks">'); b=html.index('</div>', a)
html=html[:a]+'<div class="ticks"><button data-frame="0">RNA · genes</button><button data-frame="0.5">50 / 50</button><button data-frame="1">ATAC · peaks</button>'+html[b:]
replace('Equal spacing per context doubling. Intermediate positions are visual interpolations, not measured embeddings.', 'RNA → ATAC, matched cell by cell. Intermediate positions and count colors are visual blends, not measured values.')
replace('<option value="ctxlen">ctxlen</option><option value="rna">log(rna_umi_counts)</option>', '<option value="modality">modality</option><option value="umi">UMI counts</option>')
replace('</select><p class="note" id="color-note">', '</select><div id="umi-scale-wrap" hidden style="margin-top:14px"><label for="umi-scale" class="label" style="display:block">Count scale</label><select id="umi-scale" style="width:100%"><option value="linear">Raw counts</option><option value="log">Log counts (ln)</option></select></div><p class="note" id="color-note">')
replace('across all five embeddings.', 'across both modalities.', 2)
start=html.index('const rnaLogs='); end=html.index('const countsClass=',start)
html=html[:start]+(ROOT/'src/modality-colors.js').read_text()+'\n'+html[end:]
replace("function ctxColor(v){return `hsl(${180-v*32} 74% ${34+v*2}%)`}", "function modalityColor(v){const u=ease(v);return `rgb(${[39,120,142].map((c,j)=>Math.round(c+([170,87,150][j]-c)*u)).join(',')})`}")
replace_function('position', "function position(i,v=t){const u=ease(v),a=i*2,b=a+N*2;return [coords[a]+(coords[b]-coords[a])*u,coords[a+1]+(coords[b+1]-coords[a+1])*u]}")
replace("mode==='rna'?rnaColors[i]:ctxColor(t)", "mode==='umi'?umiColor(i):modalityColor(t)")
replace_function('setTime', """function setTime(v){t=Math.max(0,Math.min(1,v));$('#slider').value=t;const pct=Math.round(t*100),label=t===0?'RNA':t===1?'ATAC':`${pct}% ATAC`;$('#ctx').textContent=label;$('#slider').setAttribute('aria-valuetext',t===0?'RNA':t===1?'ATAC':`${100-pct}% RNA, ${pct}% ATAC transition`);$('#frame-label').textContent=t===0?'RNA · observed embedding':t===1?'ATAC · observed embedding':`RNA → ATAC · ${pct}% transition · interpolated`;$('#tooltip').style.display='none';dirty=true}""")
# Retain the original ordered categorical legend and filtering implementation.
replace("mode==='rna'?'Natural logarithm of RNA UMI counts per cell. One fixed color scale across all five context settings.':'All cells share the current context length. Their color changes together as the context grows.'", "mode==='umi'?'RNA endpoint: rna_umi_counts. ATAC endpoint: atac_umi_counts. Colors blend on one shared scale during the transition.':'Teal for RNA, orchid for ATAC. The shared color blends as cells move between modalities.'")
replace("$('#search').hidden=mode==='ctxlen'||mode==='rna';$('#clear-filter').hidden=mode==='ctxlen'||mode==='rna';$('#legend-label').textContent=mode==='ctxlen'?'Context scale':mode==='rna'?'ln(RNA UMI counts)':'Click a label to spotlight';", "$('#umi-scale-wrap').hidden=mode!=='umi';$('#search').hidden=mode==='modality'||mode==='umi';$('#clear-filter').hidden=mode==='modality'||mode==='umi';$('#legend-label').textContent=mode==='modality'?'Modality':mode==='umi'?(umiScale==='log'?'ln(UMI counts)':'UMI counts'):'Click a label to spotlight';")
start=html.index("if(mode==='rna'){html="); end=html.index('else D.classes.forEach',start)
html=html[:start]+"""if(mode==='umi'){const [lo,hi]=countExtent();html=`<div class="gradient" style="background:linear-gradient(to right,${umiStops.map((_,i)=>paletteColor(i/4)).join(',')})"></div><div class="scale-labels"><span>${formatCount(lo)}</span><span>${formatCount((lo+hi)/2)}</span><span>${formatCount(hi)}</span></div><p class="note">Low → high counts · fixed across modalities<br>${umiScale==='log'?'Natural logarithm (ln)':'Raw UMI counts'}<br>Intermediate colors are illustrative blends.</p>`}else if(mode==='modality'){html=`<div class="gradient" style="background:linear-gradient(to right,${modalityColor(0)},${modalityColor(1)})"></div><div class="scale-labels"><span>RNA · genes</span><span>ATAC · peaks</span></div>`}"""+html[end:]
replace('RNA UMIs: ${c[5].toLocaleString()} · ln: ${Math.log(c[5]).toFixed(2)}', 'RNA UMIs: ${c[5].toLocaleString()} · ln: ${Math.log(c[5]).toFixed(2)}<br>ATAC UMIs: ${c[6].toLocaleString()} · ln: ${Math.log(c[6]).toFixed(2)}', 2)
replace("$('#play').onclick=", "$('#umi-scale').onchange=e=>{umiScale=e.target.value;legend();dirty=true};\n$('#play').onclick=")
replace('/6000*Number', '/12000*Number')
replace('if(next>=4){next=4;', 'if(next>=1){next=1;')
replace('.ctxvalue{font-size:29px', '.ctxvalue{font-size:25px')
replace('.controls{gap:5px}.controls select', '.ctxvalue{font-size:20px}.ctxvalue small{display:block;margin:3px 0 0}.controls{gap:5px}.controls select')
output=ROOT/'figure4h.html'
output.write_text(html.replace('__DATA__',json.dumps(data,separators=(',',':')).replace('</','<\\/')))
print(f'Built {output.name}: {output.stat().st_size/1e6:.2f} MB; {len(ids):,} paired cells')
