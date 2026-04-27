# Architecture: Spatial

High‑level layers (bottom‑up):

1. **IO sub‑package** (`spatial.io`) — platform‑specific loaders build AnnData objects
   with consistent schema (X, obs, var, obsm['spatial'], uns['platform']).
2. **Analysis sub‑packages** — independent algorithms:
   - `spatial.analysis.autocorrelation`  Moran, LISA, Geary, Getis‑Ord
   - `spatial.analysis.clustering`      spatially‑constrained Leiden / k‑means
   - `spatial.analysis.neighborhood`    Ripley's K, interaction index, NES
   - `spatial.analysis.deconvolution`   Stereoscope, Tangram
   - `spatial.analysis.communication`   ligand‑receptor scoring
3. **Visualisation** (`spatial.visualization`) — Matplotlib/Plotly wrappers for
   spot scatter, expression overlays, domain boundaries.
4. **Integration** (`spatial.integration`) — batch correction, label transfer,
   coordinate registration, cross‑platform downsampling.
5. **Configuration** (`spatial.config`) — YAML‑driven knobs: graph (k, radius),
   clustering (resolution, min_size), deconvolution (epochs, device).

## Data flow invariants

```
Loader → AnnData
  ↓
(Optional: QC, normalisation)
  ↓
neighbors() → spatial graph in adata.obsp['connectivities']
  ↓
Any analysis module (reads graph + coordinates)
  ↓
Results stored in adata.obs['<method>'] or adata.uns['<result_dict>']
  ↓
visualization.* reads adata + result key
```

All analysis functions are **pure** (no side‑effects on global state), accept and
return AnnData (or add `.obs` columns in‑place for convenience), and respect
`spatial.graph.*` config keys set via `config.set()`.

## Subpackage boundaries

```
spatial/
  io/                  loaders only
  analysis/
    __init__.py        public façade (exposes top‑level functions)
    autocorrelation.py
    clustering.py
    …
  visualization/
    __init__.py
    plots.py           matplotlib helpers
    interactive.py     plotly / bokeh
  integration/
    __init__.py
    batch.py           batch_integrate()
    registration.py    register_coordinates()
    label_transfer.py  transfer_labels()
  config.py            schema + defaults
```

No analysis module imports `visualization`; no `integration` depends on `io` except
for loader type hints. Dependencies flow upward only.
