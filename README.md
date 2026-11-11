# Bloom — interactive figures

Self-contained, embeddable pages for the Bloom manuscript figures. Each page embeds its data, styles, scripts and wordmark font and makes no network requests.

| Page | Figure | Template | Builder | Data |
|---|---|---|---|---|
| `ctx.html` | 2e, embeddings across context length | `src/ctx_template.html` | `build_ctx.py` | `ctx_umap.parquet` (`prep/prepare_ctx.py`) |
| `mod.html` | 4i, RNA-fit UMAP with ATAC projected | `src/mod_template.html` | `build_mod.py` | `mod_umap.parquet` (`prep/prepare_mod.py`) |

Setup: `python3 -m venv .venv && .venv/bin/pip install -r requirements.txt`.

## ctx.html — context length with the new UMAP settings

`ctx.html` is a self-contained page (1.5 MB, no network requests) built from `ctx_umap.parquet`.

- Data: one UMAP fitted jointly on the five context embeddings (500, 1,000, 2,000, 4,000, 8,000) of the same 15,750 cells, with the figure settings n_neighbors=10, min_dist=0.6, spread=1.0, euclidean, seed 0 (bloom-figures `scripts/prepare/28_figure2_context.py` → `data/formatted/figure2_context_umap_coords.csv`).
- `ctx_umap.parquet` (78,750 rows, one per cell × context): `cell_id`, `row_in_dataset`, `ctxlen`, `umap_1`, `umap_2`, `color`, `AIT33_class`, `AIT33_subclass`, `class_color`, `subclass_color`. Colours come from the AIT33 taxonomy table used in bloom-umaps. As in bloom-umaps, each cell takes the colour of its own `AIT33_class`, and `color` equals `class_color`. 128 cells have an annotated class that differs from their subclass's usual parent class. Input hashes are in `ctx_umap.prov.json`.
- Autoplay runs once from 500 to 8,000 and stops (6 s per doubling at 1×; 1× or 2×). Pressing play at 8,000 restarts from 500. The arrow keys step between measured contexts.
- The bottom bar shows the context readout (a measured value, or `500 → 1,000` while interpolating), a short timeline, and the figure caption. In frames narrower than 760 px (phones) the caption sits below the controls, clamped to two lines with a More toggle. `caption=0` removes it.
- Touch: drag to pan, pinch to zoom, tap a cell to see its labels and trace it.
- The header wordmark follows the bloom-preview site (orange `b`, `#ff6e00`). It uses Figtree ExtraBold (SIL Open Font License 1.1), embedded as a 1.9 KB subset of the glyphs b, l, o, m, so the page still makes no network requests.
- The legend panel is closed by default (**Legend** button or `panel=1`). It holds the colour mode (Taxonomy: AIT33 class or subclass, or context), label filtering and highlighting, and the selected cell. Click a cell to trace its five positions.

### Embedding

Inside an iframe, or with `?embed=1`, the page hides its header, fills the frame, and shows a **Full page ↗** link. When the frame is narrower than 760 px, the legend opens as a drawer over the plot. The standalone page's **Embed** button shows this snippet:

```html
<iframe src="https://rhngla.github.io/bloom-embed/ctx.html?embed=1" title="Bloom: cell embeddings across context length"
        width="100%" height="640" style="border:0;border-radius:8px" loading="lazy"></iframe>
```

URL parameters: `color=class|subclass|ctxlen`, `ctx=500|1000|2000|4000|8000` (starting frame), `autoplay=1`, `panel=1` (open the legend), `caption=0` (omit the description text, e.g. when the host page supplies its own caption), `embed=0|1`.

### Rebuild

The prep step runs on the analysis host, where its inputs are:

```sh
/opt/dlami/nvme/envs/bloom-figures/bin/python prep/prepare_ctx.py   # -> ctx_umap.parquet + ctx_umap.prov.json
.venv/bin/python build_ctx.py                                       # -> ctx.html
```

## mod.html — RNA → ATAC (Figure 4i)

`mod.html` is a self-contained page (1.0 MB, no network requests) built from `mod_umap.parquet` with `src/mod_template.html`. It uses the same layout, styling and interactions as `ctx.html` and is embeddable in the same way.

- Data: Figure 4i, i.e. bloom-figures `scripts/prepare/41_build_modality_umaps.py` → `data/formatted/figure4_modality_umap_coords.csv`, variant `project`. A UMAP is fit on Bloom's RNA embedding alone (`top4096`, genes only) of the 15,750 allen_multiome eval_balanced_50 cells, and Bloom's ATAC embedding of the same cells (`resamp_atac_8192_d16`, peaks only) is placed with `umap.transform()`. Settings: n_neighbors=10, min_dist=0.6, spread=1.0, euclidean, seed 0.
- `mod_umap.parquet` (31,500 rows, one per cell × modality): `cell_id`, `row_in_dataset`, `modality` (`RNA`, `ATAC`), `role` (`fit`, `projected`), `umap_1`, `umap_2`, `AIT33_class`, `AIT33_subclass`, `class_color`, `subclass_color`. Labels come from metadata.h5ad and colours from the AIT33 taxonomy table, as for `ctx.html`. Input hashes are in `mod_umap.prov.json`.
- Two frames, RNA (fit) and ATAC (projected), matched by `cell_id`; positions in between are smoothstep-interpolated for display. Autoplay runs once from RNA to ATAC (6 s at 1×).
- Colour modes: AIT33 class, AIT33 subclass, or modality (green RNA `#5DC863` → indigo ATAC `#3B528B`, bloom-figures `palette.MODALITY_COLORS`, as in Figure 4i).
- Click a cell to trace its RNA and ATAC positions.

URL parameters: `color=class|subclass|modality`, `mod=rna|atac` (starting frame), `autoplay=1`, `panel=1`, `caption=0`, `embed=0|1`. Embed snippet:

```html
<iframe src="https://rhngla.github.io/bloom-embed/mod.html?embed=1" title="Bloom: cell embeddings across modalities"
        width="100%" height="640" style="border:0;border-radius:8px" loading="lazy"></iframe>
```

### Rebuild

```sh
/opt/dlami/nvme/envs/bloom-figures/bin/python prep/prepare_mod.py   # -> mod_umap.parquet + mod_umap.prov.json
.venv/bin/python build_mod.py                                       # -> mod.html
```

## GitHub Pages

The `.github/workflows/pages.yml` workflow publishes the self-contained HTML pages on every push to `main`. It uploads only the generated pages and the landing page, not the parquet inputs or build environment.

- Landing page: https://rhngla.github.io/bloom-embed/ (source `src/landing.html`)
- Context length (Figure 2e): https://rhngla.github.io/bloom-embed/ctx.html
- Modalities (Figure 4i): https://rhngla.github.io/bloom-embed/mod.html

After editing a template, rerun its builder and commit the generated HTML.
