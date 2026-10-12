# Bloom — Context in motion

Open `index.html` directly in a modern browser. The 2 MB page contains all data, styles, and scripts, with no external dependencies or network requests.

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
