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

## GitHub Pages

The `.github/workflows/pages.yml` workflow publishes the two self-contained HTML pages on every push to `main`. It uploads only the generated pages and a small landing page, not the parquet inputs or build environment.

- Landing page: https://rhngla.github.io/bloom-embed/
- Figure 2e: https://rhngla.github.io/bloom-embed/figure2e.html
- Figure 4h: https://rhngla.github.io/bloom-embed/figure4h.html

After editing the shared template, run both builders and commit their generated HTML outputs. The landing page source is `src/landing.html`.
