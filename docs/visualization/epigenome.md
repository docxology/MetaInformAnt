# Epigenome Visualization

This document provides comprehensive documentation for epigenetic data visualization capabilities in METAINFORMANT, including DNA methylation, histone modifications, and chromatin accessibility analysis.

## Overview

Epigenome visualization includes specialized plots for genomic epigenetic marks, including methylation profiles, ChIP-seq peaks, ATAC-seq signals, and chromatin state analysis. These tools create publication-quality figures for epigenetic research.

## Module Functions

### DNA Methylation

#### Methylation Profiles
```python
from metainformant.epigenome import visualization as epi_viz
import numpy as np

# Plot methylation levels across CpG sites
methylation_data = np.random.rand(1000)  # 0-1 methylation levels
positions = np.arange(1000)  # Genomic positions

ax = epi_viz.plot_methylation_profile(methylation_data, positions, figsize=(12, 6))
```

#### Differential Methylation
```python
# Volcano plot for differential methylation analysis
dm_results = [
    {'delta': 0.3, 'pvalue': 1e-10, 'position': 1000},
    {'delta': -0.2, 'pvalue': 1e-5, 'position': 2000},
    # ... more results
]

ax = epi_viz.plot_differential_methylation(dm_results)
```

#### Methylation Clustering
```python
# Cluster methylation data across samples
methylation_matrix = np.random.rand(50, 1000)  # Samples x CpG sites
cluster_labels = np.random.randint(0, 3, 50)  # Cluster assignments

ax = epi_viz.plot_dna_methylation_clusters(methylation_matrix, cluster_labels)
```

### Histone Modifications

#### ChIP-seq Peaks
```python
# Visualize ChIP-seq peaks along chromosome
peaks = [
    {'start': 1000, 'end': 1500, 'score': 25.0},
    {'start': 3000, 'end': 3500, 'score': 18.5},
    # ... more peaks
]

ax = epi_viz.plot_chipseq_peaks(peaks, chromosome="chr1")
```

#### Histone Modification Heatmaps
```python
# Plot multiple histone modifications
histone_data = {
    'H3K4me3': np.random.rand(1000),
    'H3K27me3': np.random.rand(1000),
    'H3K36me3': np.random.rand(1000)
}

ax = epi_viz.plot_histone_modification_heatmap(histone_data)
```

### Chromatin Accessibility

#### ATAC-seq Signals
```python
# Plot chromatin accessibility signal
signal_data = np.random.rand(2000) + 0.5  # Accessibility scores
positions = np.arange(2000)

ax = epi_viz.plot_atacseq_signal(signal_data, positions)
```

### Chromatin States

#### State Segmentation
```python
# Visualize chromatin state annotations
states = np.random.randint(1, 6, 1000)  # 1-5 state assignments
positions = np.arange(1000)
state_labels = {1: 'Active', 2: 'Repressed', 3: 'Heterochromatin', 4: 'Enhancer', 5: 'Promoter'}

ax = epi_viz.plot_chromatin_states(states, positions, state_labels)
```

### Multi-omic Integration

#### Genome Browser Tracks
```python
# Create genome browser-style multi-track view
tracks = {
    'methylation': {
        'data': np.random.rand(1000),
        'color': 'blue',
        'label': 'Methylation'
    },
    'H3K4me3': {
        'data': np.random.rand(1000),
        'color': 'red',
        'label': 'H3K4me3'
    },
    'accessibility': {
        'data': np.random.rand(1000),
        'color': 'green',
        'label': 'ATAC-seq'
    }
}

ax = epi_viz.plot_genome_browser_tracks(tracks, region_start=0, region_end=1000)
```

#### Epigenetic Correlation Analysis
```python
# Analyze correlations between epigenetic marks
epigenetic_data = {
    'DNA_methylation': np.random.rand(100),
    'H3K4me3': np.random.rand(100),
    'H3K27me3': np.random.rand(100),
    'ATAC_signal': np.random.rand(100)
}

ax = epi_viz.plot_epigenetic_correlation_heatmap(epigenetic_data)
```

### Advanced Features

#### Interactive Genome Browser
```python
# Create interactive epigenetic data browser
epigenetic_tracks = {
    'methylation': {'positions': list(range(1000)), 'values': np.random.rand(1000)},
    'chipseq': {'positions': list(range(1000)), 'values': np.random.rand(1000)},
    'atacseq': {'positions': list(range(1000)), 'values': np.random.rand(1000)}
}

fig = epi_viz.create_interactive_epigenome_browser(epigenetic_tracks)
```

## Integration with Epigenome Module

### With Methylation Analysis
```python
from metainformant.epigenome import methylation, visualization as epi_viz

# Load methylation data
methylation_sites = methylation.load_methylation_bedgraph("methylation.bedgraph")

# Convert to visualization format
positions = [site['position'] for site in methylation_sites.values()]
levels = [site['methylation_level'] for site in methylation_sites.values()]

# Visualize
ax = epi_viz.plot_methylation_profile(np.array(levels), np.array(positions))
```

### With ChIP-seq Analysis
```python
from metainformant.epigenome import chipseq, visualization as epi_viz

# Load ChIP-seq peaks
peaks = chipseq.load_chip_peaks("peaks.narrowPeak")

# Visualize peak distribution
ax = epi_viz.plot_chipseq_peaks(peaks, chromosome="chr1")
```

### With ATAC-seq Analysis
```python
from metainformant.epigenome import atacseq, visualization as epi_viz

# Load ATAC-seq peaks
peaks = atacseq.load_atac_peaks("atac_peaks.narrowPeak")

# Create accessibility profile
# (Note: would need signal data for full profile)
signal_data = np.random.rand(1000)  # Placeholder
ax = epi_viz.plot_atacseq_signal(signal_data, np.arange(1000))
```

## Output Options

All visualization functions support:
```python
# Save to file
ax = epi_viz.plot_methylation_profile(methylation_data, output_path="methylation.png")

# Interactive browser
fig = epi_viz.create_interactive_epigenome_browser(tracks, output_path="browser.html")
```

## Data Format Support

- **BED Files**: Standard genomic interval format
- **BedGraph Files**: Continuous genomic signal data
- **NarrowPeak Files**: ChIP-seq and ATAC-seq peak calls
- **Methylation Data**: Bismark or custom methylation formats
- **Chromatin States**: Segway, ChromHMM, or custom state annotations

## Performance Considerations

- **Large Genomic Regions**: For whole chromosomes, consider tiling or subsampling
- **Multiple Samples**: Heatmaps work best with 10-100 samples
- **Interactive Plots**: Require Plotly for web-based genome browsers
- **Memory Usage**: High-resolution tracks can be memory-intensive

## Dependencies

- **Required**: matplotlib, numpy, pandas
- **Optional**: seaborn (enhanced styling), plotly (interactive plots)
- **Integration**: pybedtools (genomic intervals), pyBigWig (signal files)

## Examples

### Complete Epigenetic Analysis Workflow
```python
from metainformant.epigenome import methylation, chipseq, visualization as epi_viz
import numpy as np

# Load epigenetic data
methylation_data = methylation.load_methylation_bedgraph("sample.methylation.bedgraph")
chip_peaks = chipseq.load_chip_peaks("H3K4me3_peaks.narrowPeak")

# Create comprehensive visualization
fig, axes = plt.subplots(2, 2, figsize=(15, 12))

# Methylation profile
positions = [site['position'] for site in methylation_data.values()]
levels = [site['methylation_level'] for site in methylation_data.values()]
epi_viz.plot_methylation_profile(np.array(levels), np.array(positions), ax=axes[0,0])

# ChIP-seq peaks
epi_viz.plot_chipseq_peaks(chip_peaks, chromosome="chr1", ax=axes[0,1])

# Multi-track browser
tracks = {
    'methylation': {'data': np.array(levels), 'color': 'blue', 'label': 'DNAme'},
    'H3K4me3': {'data': np.random.rand(len(levels)), 'color': 'red', 'label': 'H3K4me3'}
}
epi_viz.plot_genome_browser_tracks(tracks, region_start=min(positions), region_end=max(positions), ax=axes[1,0])

# Correlation analysis
epi_data = {
    'DNA_methylation': np.array(levels),
    'H3K4me3_signal': np.random.rand(len(levels)),
    'accessibility': np.random.rand(len(levels))
}
epi_viz.plot_epigenetic_correlation_heatmap(epi_data, ax=axes[1,1])

plt.tight_layout()
plt.savefig("epigenetic_analysis.png", dpi=300, bbox_inches='tight')
```

## Color Schemes

Recommended color schemes for epigenetic data:
```python
# DNA methylation
methylation_cmap = 'RdYlBu_r'  # Blue for unmethylated, red for methylated

# Histone modifications
histone_colors = {
    'H3K4me3': '#FF6B6B',    # Red for active promoters
    'H3K27me3': '#4ECDC4',   # Teal for repressed regions
    'H3K36me3': '#45B7D1',   # Blue for gene bodies
    'H3K9me3': '#96CEB4',    # Green for heterochromatin
}

# Chromatin states
state_colors = {
    1: '#FF6B6B',   # Active
    2: '#4ECDC4',   # Repressed
    3: '#45B7D1',   # Heterochromatin
    4: '#96CEB4',   # Enhancer
    5: '#FFEAA7'    # Promoter
}
```

## Troubleshooting

### Common Issues

1. **Empty Plots**: Check data ranges and genomic coordinates
2. **Memory Errors**: Subsample large genomic regions for visualization
3. **File Format Errors**: Ensure BED/BedGraph files have correct column formats
4. **Color Scaling**: Use appropriate colormaps for different data types

### Data Requirements

- **Genomic Coordinates**: Must be sorted and non-overlapping for proper display
- **Signal Values**: Should be normalized for comparison across samples
- **Peak Scores**: Higher scores typically indicate stronger enrichment
- **Methylation Levels**: Values should be in 0-1 range for proper color scaling

## Related Documentation

- **[Epigenome Analysis](../epigenome/)**: Core epigenetic analysis functions
- **[Methylation Analysis](../epigenome/methylation.md)**: DNA methylation workflows
- **[ChIP-seq Analysis](../epigenome/chipseq.md)**: Histone modification analysis
- **[Visualization Integration](integration.md)**: Cross-module visualization patterns

This module provides comprehensive epigenetic visualization capabilities integrated with METAINFORMANT's epigenomic analysis workflows.
