# Bloom — Context in motion

Open `figure2e.html` directly in a modern browser. The 2 MB page contains all data, styles, and scripts, with no external dependencies or network requests.

- Autoplay moves back and forth between 500 and 8,000 context length (24 seconds each way at 1×).
- Drag the slider or click a measured context setting. Choose 0.5×, 1×, or 2× playback.
- Color by AIT33_subclass, AIT33_class, ctxlen, or log(rna_umi_counts). Search or select annotation labels to spotlight cells.
- Scroll to zoom, drag to pan, and click a cell to see its metadata and five-position trajectory. Escape clears the trajectory.

## Data and interpolation

All 15,750 unique `cell_id` values are present at 500, 1,000, 2,000, 4,000, and 8,000 context length. Class and subclass annotations are checked for consistency across frames. Positions are joined by `cell_id`, never by row order. Coordinates are stored as float32 and displayed at a fixed aspect-preserving scale across frames; no alignment, rotation, or normalization is applied.

Each adjacent pair is interpolated using `u²(3−2u)` along its straight-line segment. This passes through each recorded coordinate and eases velocity to zero at each measured setting, without overshoot. Intermediate embeddings are illustrative, not model outputs. Each doubling occupies an equal slider interval; the displayed context value changes linearly within that interval.

The 41 classes receive separated hues. The 315 subclasses are sorted by their zero-padded numeric labels within each parent class, with gradual changes in hue, saturation, and lightness. Context coloring assigns every cell the same color for the current shared context length.

## Rebuild

The source is `src/template.html`; `build.py` validates and embeds `figure2_umap_ctxlen_metadata.parquet`.

```sh
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python build.py
```

Browser checks cover all 78,750 cell/frame positions, slider, playback, four color modes, search, filtering, trajectory selection, and zoom/reset. Desktop and narrow layouts were visually inspected. The browser test tool blocks file URLs, so direct-file opening was not tested; the page uses only embedded resources and was tested through a local HTTP server.

RNA coloring uses the natural logarithm of positive RNA UMI counts (1,452–150,445), with a fixed purple–teal–yellow scale. Counts are verified identical for each cell across all five frames. The page uses a light ivory background with darker annotation colors for contrast.

## Figure 4h — RNA ↔ ATAC

Open `figure4h.html` for the separate modality visualization. It is a self-contained HTML file; the context-length page remains `figure2e.html`.

Rebuild with:

```sh
.venv/bin/python build_figure4h.py
```

The builder reuses `src/template.html` for the light layout, categorical palettes, canvas rendering, filtering, cell inspection, and playback controls. Figure-specific count logic lives in `src/modality-colors.js`. Guarded template substitutions detect incompatible changes. The Figure 4h builder never writes the original template or `figure2e.html`.

- 15,750 unique cells are paired by `cell_id` between `RNA (genes)` and `ATAC (peaks)`. Duplicate IDs, missing pairs, annotation changes, invalid coordinates, and nonpositive counts fail the build.
- The original coordinates are retained on one fixed scale, with smoothstep interpolation from RNA to ATAC and back (12 seconds each way at 1×). Intermediate positions are visual blends, not measured embeddings.
- Color modes: `AIT33_class`, `AIT33_subclass`, modality, and UMI counts. UMI coloring offers raw counts (default) and natural-log counts, with a fixed scale spanning both modalities (664–203,786 raw UMIs).
- RNA endpoint counts come from the RNA rows' `rna_umi_counts`; ATAC endpoint counts come from the ATAC rows' `atac_umi_counts`. Between endpoints, values interpolate on the selected raw/log scale with the same easing as coordinates. These intermediate values are for coloring only; cell details always report both original endpoint counts.
- Some subclass labels have cells assigned to different classes in the supplied metadata. Subclass palettes and legend groups use the most frequent parent class, while class coloring and cell details preserve every recorded per-cell class.

Validated in-browser: all 31,500 endpoint positions and raw/log endpoint counts match the embedded data; all four color modes, slider endpoints, playback/pause, annotation search/filtering, cell trajectories, and zoom/reset work. Desktop and mobile layouts were inspected with no console errors. The original context-length visualization is retained as `figure2e.html`.

## ctx.html — context length with the new UMAP settings

`ctx.html` is a self-contained page (1.5 MB, no network requests) built from `ctx_umap.parquet`. `figure2e.html` is unchanged.

- Data: one UMAP fitted jointly on the five context embeddings (500, 1,000, 2,000, 4,000, 8,000) of the same 15,750 cells, with the figure settings n_neighbors=10, min_dist=0.6, spread=1.0, euclidean, seed 0 (bloom-figures `scripts/prepare/28_figure2_context.py` → `data/formatted/figure2_context_umap_coords.csv`).
- `ctx_umap.parquet` (78,750 rows, one per cell × context): `cell_id`, `row_in_dataset`, `ctxlen`, `umap_1`, `umap_2`, `color`, `AIT33_class`, `AIT33_subclass`, `class_color`, `subclass_color`. Colours come from the AIT33 taxonomy table used in bloom-umaps. As in bloom-umaps, each cell takes the colour of its own `AIT33_class`, and `color` equals `class_color`. 128 cells have an annotated class that differs from their subclass's usual parent class. Input hashes are in `ctx_umap.prov.json`.
- Autoplay runs once from 500 to 8,000 and stops (6 s per doubling at 1×; 1× or 2×). Pressing play at 8,000 restarts from 500. The arrow keys step between measured contexts.
- The bottom bar shows the context readout (a measured value, or `500 → 1,000` while interpolating), a short timeline, and the figure caption. In frames narrower than 760 px (phones) the caption sits below the controls, clamped to two lines with a More toggle. `caption=0` removes it.
- Touch: drag to pan, pinch to zoom, tap a cell to see its labels and trace it.
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

## GitHub Pages

The `.github/workflows/pages.yml` workflow publishes the three self-contained HTML pages on every push to `main`. It uploads only the generated pages and a small landing page, not the parquet inputs or build environment.

- Landing page: https://rhngla.github.io/bloom-embed/
- Figure 2e: https://rhngla.github.io/bloom-embed/figure2e.html
- Figure 4h: https://rhngla.github.io/bloom-embed/figure4h.html
- Context length (new UMAP settings): https://rhngla.github.io/bloom-embed/ctx.html

After editing the shared template, run both builders and commit their generated HTML outputs. The landing page source is `src/landing.html`.
