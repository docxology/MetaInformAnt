# Spatial Transcriptomics Module

Spatial transcriptomics analysis supporting Visium, MERFISH, and Xenium platforms.

## Submodules

### I/O (`io/`)
- **visium.py** - 10x Visium data loading, tissue positions, spatial images
- **merfish.py** - MERFISH transcript spots, cell metadata, aggregation
- **xenium.py** - 10x Xenium transcripts, cell features, boundaries

### Analysis (`analysis/`)
- **clustering.py** - Spatial clustering (Leiden/Louvain), spatial graph construction, domain detection
- **deconvolution.py** - Cell type deconvolution (NNLS, NMF), reference profiles, enrichment scoring
- **autocorrelation.py** - Moran's I, Geary's C, Getis-Ord G, variograms, spatial weights
- **neighborhood.py** - Neighborhood enrichment, niche detection, Ripley's K, ligand-receptor signaling

### Integration (`integration/`)
- **scrna_mapping.py** - scRNA-seq mapping, anchor-based transfer, gene imputation, correlation

### Visualization (`visualization/`)
- **plots.py** - Spatial scatter, tissue overlay, gene/cell-type maps, deconvolution pies, interaction heatmaps

## Usage

```python
from metainformant.spatial.io import visium, merfish, xenium
from metainformant.spatial.analysis import clustering, deconvolution, autocorrelation
from metainformant.spatial.integration import scrna_mapping
```
