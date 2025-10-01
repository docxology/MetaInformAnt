# Epigenome Module

The `epigenome` module provides tools for epigenetic modification analysis, including DNA methylation, chromatin structure, and histone modifications.

## Overview

This module handles epigenetic data processing, from raw sequencing data to biological interpretation of epigenetic modifications.

## Submodules

### Methylation Analysis (`methylation.py`)
DNA methylation pattern detection and analysis.

**Key Features:**
- Methylation calling from sequencing data
- Differential methylation analysis
- Methylation quantitative trait locus (mQTL) mapping
- Methylation age prediction

**Usage:**
```python
from metainformant.epigenome import methylation

# Methylation analysis
methylation_data = methylation.load_bismark_output("methylation_calls.txt")
differential = methylation.find_differential_methylation(methylation_data, groups)
```

### Chromatin Tracks (`tracks.py`)
Chromatin state and accessibility analysis.

**Key Features:**
- ATAC-seq and DNase-seq processing
- Chromatin accessibility quantification
- Peak calling and annotation
- Chromatin state segmentation

**Usage:**
```python
from metainformant.epigenome import tracks

# Chromatin analysis
atac_data = tracks.load_atac_peaks("peaks.bed")
accessibility = tracks.calculate_accessibility(atac_data, genome_annotation)
```

## Integration with Other Modules

### With DNA Module
```python
from metainformant.dna import sequences
from metainformant.epigenome import methylation

# Epigenetic regulation of DNA sequences
promoter_sequences = sequences.get_promoter_sequences(genes)
methylation_levels = methylation.analyze_promoter_methylation(promoter_sequences)
```

### With Visualization Module
```python
from metainformant.epigenome import methylation
from metainformant.visualization import plots

# Visualize methylation patterns
methylation_profile = methylation.get_methylation_profile(genomic_region)
plots.lineplot(methylation_profile, title="DNA Methylation Profile")
```

## Data Formats

- **Bisulfite sequencing**: Bismark, BS-Seeker output formats
- **ATAC-seq**: BED, narrowPeak, broadPeak formats
- **ChIP-seq**: BED, SAM, BAM formats
- **WGBS**: Methylation call formats

## Performance Features

- Memory-efficient processing of large epigenomic datasets
- Parallel peak calling and analysis
- Streaming processing for genome-scale data

## Testing

Comprehensive tests cover:
- Methylation calling accuracy
- Peak detection algorithms
- Integration with genomic annotations

## Dependencies

- PyBedTools for genomic interval operations
- Optional: specialized epigenomic analysis packages

This module provides essential tools for epigenetic research and chromatin biology.
