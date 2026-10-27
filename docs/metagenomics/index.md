# Metagenomics Module

Metagenomic analysis tools for amplicon profiling, shotgun metagenomics, and functional annotation.

## Submodules

### Amplicon (`amplicon/`)
- **otu_clustering.py** - OTU clustering, identity calculation, chimera filtering
- **asv_denoising.py** - ASV denoising, error rate estimation, paired-read merging
- **taxonomy.py** - Taxonomic classification, tree building, confidence scoring

### Shotgun (`shotgun/`)
- **assembly.py** - Metagenome assembly, scaffolding, assembly statistics
- **binning.py** - Contig binning, tetranucleotide frequency, marker-gene bin quality (completeness, contamination) and refinement
- **profiling.py** - Community profiling, relative abundance, k-mer indexing

### Diversity (`diversity/`)
- **metrics.py** - Alpha and beta diversity, rarefaction, PERMANOVA, ordination (PCoA, NMDS)

### Functional (`functional/`)
- **annotation.py** - Gene annotation, ORF prediction, gene family classification
- **pathways.py** - Pathway reconstruction, completeness calculation, profile comparison

### Comparative (`comparative/`)
- **differential_abundance.py** - Differential abundance, CLR transform, indicator species, LEfSe-style effect sizes, biomarker discovery

### Visualization (`visualization/`)
- **plots.py** - Krona charts, stacked bars, rarefaction curves, ordination (PCoA or NMDS), alpha diversity plots, abundance heatmaps

## Usage

```python
from metainformant.metagenomics.amplicon import otu_clustering, taxonomy
from metainformant.metagenomics.shotgun import assembly, binning, profiling
from metainformant.metagenomics.functional import annotation, pathways
from metainformant.metagenomics.diversity import metrics
from metainformant.metagenomics.comparative import differential_abundance
```
