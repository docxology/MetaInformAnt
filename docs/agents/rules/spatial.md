# Spatial Transcriptomics Module Rules

## Purpose
Spatial transcriptomics analysis supporting Visium, MERFISH, and Xenium platforms with spatial clustering, deconvolution, autocorrelation statistics, neighborhood analysis, and scRNA-seq integration.

## Source Structure
```
src/metainformant/spatial/
├── io/
│   ├── visium.py          # 10x Visium data loading, tissue positions, spatial images
│   ├── merfish.py         # MERFISH transcript spots, cell metadata, aggregation
│   └── xenium.py          # 10x Xenium transcripts, cell features, boundaries
├── analysis/
│   ├── clustering.py      # Spatial clustering (Leiden/Louvain), graph construction
│   ├── deconvolution.py   # Cell type deconvolution (NNLS, NMF), enrichment scoring
│   ├── autocorrelation.py # Moran's I, Geary's C, Getis-Ord G, variograms
│   └── neighborhood.py    # Neighborhood enrichment, niche detection, Ripley's K
├── integration/
│   └── scrna_mapping.py   # scRNA-seq mapping, anchor-based transfer, gene imputation
└── visualization/
    └── plots.py           # Spatial scatter, tissue overlay, gene maps, interaction heatmaps
```

## Dependencies
- **Required**: numpy, pandas, scipy, anndata
- **Optional**: scanpy, squidpy, libpysal (spatial weights)

## Import Patterns
```python
from metainformant.spatial.io import visium, merfish, xenium
from metainformant.spatial.analysis import clustering, deconvolution, autocorrelation
from metainformant.spatial.integration import scrna_mapping
```

## Configuration
- Environment prefix: `SPATIAL_` (e.g., `SPATIAL_THREADS`, `SPATIAL_WORK_DIR`)
- Output path: `output/spatial/<analysis_type>/` (e.g., `clustering/`, `deconvolution/`, `integration/`)

## Integration
- **Spatial → Single-Cell**: scRNA-seq reference for deconvolution and mapping
- **Spatial → Visualization**: Spatial plots, tissue overlays
- **Spatial → Networks**: Spatial interaction networks, ligand-receptor

## Testing
- Generate synthetic spatial coordinate and expression data programmatically
- Use `@pytest.mark.slow` for large spatial dataset tests
- All test outputs to `tmp_path`
