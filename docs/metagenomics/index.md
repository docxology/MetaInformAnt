# Metagenomics Module

Metagenomic analysis tools for amplicon profiling, shotgun metagenomics, and functional annotation.

## Submodules

### Amplicon (`amplicon/`)
- **otu_clustering.py** - OTU clustering, identity calculation, chimera filtering
- **asv_denoising.py** - ASV denoising, error rate estimation, paired-read merging
- **taxonomy.py** - Taxonomic classification, tree building, confidence scoring

### Shotgun (`shotgun/`)
- **assembly.py** - Metagenome assembly, scaffolding, assembly statistics
- **binning.py** - Contig binning, tetranucleotide frequency, bin refinement/QC
- **profiling.py** - Community profiling, relative abundance, k-mer indexing

### Functional (`functional/`)
- **annotation.py** - Gene annotation, ORF prediction, gene family classification
- **pathways.py** - Pathway reconstruction, completeness calculation, profile comparison

### Visualization (`visualization/`)
- **plots.py** - Krona charts, stacked bars, rarefaction curves, ordination, alpha diversity heatmaps

## Usage

```python
from metainformant.metagenomics.amplicon import otu_clustering, taxonomy
from metainformant.metagenomics.shotgun import assembly, binning, profiling
from metainformant.metagenomics.functional import annotation, pathways
```
