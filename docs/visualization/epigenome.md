# Epigenome Visualization

This document provides comprehensive documentation for epigenetic data visualization capabilities in METAINFORMANT, including DNA methylation, histone modifications, and chromatin accessibility analysis.

## Overview

Epigenome visualization includes specialized plots for genomic epigenetic marks, including methylation profiles, ChIP-seq peaks, ATAC-seq signals, and chromatin state analysis. These tools create publication-quality figures for epigenetic research.

All plotting functions live in `metainformant.epigenome.visualization.visualization`; the examples import that module as `epi_viz`. They take in-memory data (NumPy arrays, lists of dicts, dicts of arrays) rather than file paths, and return the `matplotlib.axes.Axes` (the interactive browser returns a Plotly figure). Each plotting function accepts `ax=` to draw on existing axes, `output_path=` to save the figure, and `figsize=` to size a figure it creates. Colors are fixed by the functions (for example `RdYlBu_r` for methylation and `Set3` for chromatin states), and `**kwargs` is not used except where noted below.

## Module Functions

### DNA Methylation

#### Methylation Profiles

`plot_methylation_profile(methylation_data, positions=None, *, ax=None, output_path=None, figsize=(12, 6), **kwargs)` scatters methylation levels (a NumPy array of values between 0 and 1) along genomic `positions` (a NumPy array, default the site indices), colored by level with a colorbar. A smoothed trend line is added when there are more than 10 sites.

```python
from metainformant.epigenome.visualization import visualization as epi_viz
import numpy as np

# Plot methylation levels across CpG sites
methylation_data = np.random.default_rng(0).random(1000)  # 0-1 methylation levels
positions = np.arange(1000)  # Genomic positions

ax = epi_viz.plot_methylation_profile(methylation_data, positions, figsize=(12, 6))
```

#### Differential Methylation

`plot_differential_methylation(dm_results, *, ax=None, output_path=None, figsize=(10, 6), **kwargs)` draws a volcano plot from a list of dicts that each have `delta` (methylation difference) and `pvalue` (raw p-value, converted with -log10). Points with a p-value below 0.05 are red and the rest gray; pass `sig_threshold=` to change that cutoff. Other keys (such as `position`) are ignored.

```python
from metainformant.epigenome.visualization import visualization as epi_viz

# Volcano plot for differential methylation analysis
dm_results = [
    {'delta': 0.3, 'pvalue': 1e-10, 'position': 1000},
    {'delta': -0.2, 'pvalue': 1e-5, 'position': 2000},
    {'delta': 0.05, 'pvalue': 0.4, 'position': 3000},
]

ax = epi_viz.plot_differential_methylation(dm_results)
```

#### Methylation Clustering

`plot_dna_methylation_clusters(methylation_matrix, cluster_labels=None, *, ax=None, output_path=None, figsize=(12, 8), **kwargs)` draws a samples x CpG sites heatmap (values between 0 and 1). With `cluster_labels` (a NumPy array, one label per sample) the rows are sorted by cluster and black lines mark the cluster boundaries.

```python
from metainformant.epigenome.visualization import visualization as epi_viz
import numpy as np

# Cluster methylation data across samples
rng = np.random.default_rng(0)
methylation_matrix = rng.random((50, 1000))  # Samples x CpG sites
cluster_labels = rng.integers(0, 3, 50)  # Cluster assignments

ax = epi_viz.plot_dna_methylation_clusters(methylation_matrix, cluster_labels)
```

### Histone Modifications

#### ChIP-seq Peaks

`plot_chipseq_peaks(peaks, chromosome, *, ax=None, output_path=None, figsize=(12, 4), **kwargs)` draws each peak as a rectangle from `start` to `end` whose height is its `score` (default 1.0), colored by score. `peaks` is a list of dicts with `start` and `end` keys. `chromosome` is only used for the labels, so pass the peaks of a single chromosome.

```python
from metainformant.epigenome.visualization import visualization as epi_viz

# Visualize ChIP-seq peaks along chromosome
peaks = [
    {'start': 1000, 'end': 1500, 'score': 25.0},
    {'start': 3000, 'end': 3500, 'score': 18.5},
    {'start': 5200, 'end': 5900, 'score': 31.2},
]

ax = epi_viz.plot_chipseq_peaks(peaks, chromosome="chr1")
```

#### Histone Modification Heatmaps

`plot_histone_modification_heatmap(histone_data, positions=None, *, ax=None, output_path=None, figsize=(12, 8), **kwargs)` draws a marks x positions heatmap from a dict that maps each modification name to a signal array (all the same length). Columns are plotted by index; `positions` only switches the x-axis label text.

```python
from metainformant.epigenome.visualization import visualization as epi_viz
import numpy as np

# Plot multiple histone modifications
rng = np.random.default_rng(0)
histone_data = {
    'H3K4me3': rng.random(1000),
    'H3K27me3': rng.random(1000),
    'H3K36me3': rng.random(1000),
}

ax = epi_viz.plot_histone_modification_heatmap(histone_data)
```

### Chromatin Accessibility

#### ATAC-seq Signals

`plot_atacseq_signal(signal_data, positions, *, ax=None, output_path=None, figsize=(12, 6), **kwargs)` draws a filled line of accessibility values; both arguments are required NumPy arrays of the same length.

```python
from metainformant.epigenome.visualization import visualization as epi_viz
import numpy as np

# Plot chromatin accessibility signal
signal_data = np.random.default_rng(0).random(2000) + 0.5  # Accessibility scores
positions = np.arange(2000)

ax = epi_viz.plot_atacseq_signal(signal_data, positions)
```

### Chromatin States

#### State Segmentation

`plot_chromatin_states(states, positions, state_labels=None, *, ax=None, output_path=None, figsize=(12, 4), **kwargs)` draws a colored segment for each run of identical state ids along `positions` (both NumPy arrays of the same length, with `positions` sorted). `state_labels` is an optional dict from state id to name that is written on the segments and used in the legend.

```python
from metainformant.epigenome.visualization import visualization as epi_viz
import numpy as np

# Visualize chromatin state annotations
states = np.repeat([1, 3, 2, 5, 4, 1], [150, 200, 100, 250, 200, 100])  # state ids 1-5 per position
positions = np.arange(len(states))
state_labels = {1: 'Active', 2: 'Repressed', 3: 'Heterochromatin', 4: 'Enhancer', 5: 'Promoter'}

ax = epi_viz.plot_chromatin_states(states, positions, state_labels)
```

### Multi-omic Integration

#### Genome Browser Tracks

`plot_genome_browser_tracks(tracks, region_start, region_end, *, ax=None, output_path=None, figsize=(12, 8), **kwargs)` stacks one track per entry of `tracks` inside the axes, spacing each track's `data` values evenly over `region_start`-`region_end`. Each track is a dict with `data` (a list or array) and optional `color` (default blue) and `label` (default the track name); the first track is drawn at the bottom.

```python
from metainformant.epigenome.visualization import visualization as epi_viz
import numpy as np

# Create genome browser-style multi-track view
rng = np.random.default_rng(0)
tracks = {
    'methylation': {
        'data': rng.random(1000),
        'color': 'blue',
        'label': 'Methylation'
    },
    'H3K4me3': {
        'data': rng.random(1000),
        'color': 'red',
        'label': 'H3K4me3'
    },
    'accessibility': {
        'data': rng.random(1000),
        'color': 'green',
        'label': 'ATAC-seq'
    }
}

ax = epi_viz.plot_genome_browser_tracks(tracks, region_start=0, region_end=1000)
```

#### Epigenetic Correlation Analysis

`plot_epigenetic_correlation_heatmap(epigenetic_data, *, ax=None, output_path=None, figsize=(10, 8), **kwargs)` draws the Pearson correlation matrix between marks. `epigenetic_data` maps each mark name to a NumPy array; arrays are cut to the shortest length. Extra keyword arguments are passed to `seaborn.heatmap` (which is used when seaborn is available), whose `annot`, `fmt` and `cmap` are fixed by the function.

```python
from metainformant.epigenome.visualization import visualization as epi_viz
import numpy as np

# Analyze correlations between epigenetic marks
rng = np.random.default_rng(0)
epigenetic_data = {
    'DNA_methylation': rng.random(100),
    'H3K4me3': rng.random(100),
    'H3K27me3': rng.random(100),
    'ATAC_signal': rng.random(100)
}

ax = epi_viz.plot_epigenetic_correlation_heatmap(epigenetic_data)
```

### Advanced Features

#### Interactive Genome Browser

`create_interactive_epigenome_browser(epigenetic_tracks, *, output_path=None, **kwargs)` builds a Plotly figure and returns it; it requires `plotly` and raises `ImportError` without it. Each entry of `epigenetic_tracks` is a dict with `positions`, `values` and an optional `color`. A track is drawn as a line at its own row on the y axis, and its `values` appear in the hover text. Extra keyword arguments go to `fig.update_layout`, and `output_path` writes an HTML file (the suffix is forced to `.html`).

```python
from metainformant.epigenome.visualization import visualization as epi_viz
import numpy as np

# Create interactive epigenetic data browser (requires plotly)
rng = np.random.default_rng(0)
epigenetic_tracks = {
    'methylation': {'positions': list(range(1000)), 'values': rng.random(1000)},
    'chipseq': {'positions': list(range(1000)), 'values': rng.random(1000)},
    'atacseq': {'positions': list(range(1000)), 'values': rng.random(1000)}
}

fig = epi_viz.create_interactive_epigenome_browser(epigenetic_tracks)
```

## Integration with Epigenome Module

### With Methylation Analysis
```python
import numpy as np
from metainformant.epigenome.assays import methylation
from metainformant.epigenome.visualization import visualization as epi_viz

# Load methylation data: {chromosome: [MethylationSite, ...]}
methylation_sites = methylation.load_methylation_bedgraph("methylation.bedgraph")

# Convert one chromosome to visualization format
sites = methylation_sites["chr1"]
positions = np.array([site.position for site in sites])
levels = np.array([site.methylation_level for site in sites])

# Visualize
ax = epi_viz.plot_methylation_profile(levels, positions)
```

### With ChIP-seq Analysis
```python
from metainformant.epigenome.assays import chipseq
from metainformant.epigenome.visualization import visualization as epi_viz

# Load ChIP-seq peaks: a list of ChIPPeak objects
peaks = chipseq.load_chip_peaks("peaks.narrowPeak")

# Visualize the peak distribution on one chromosome (to_dict gives start, end and score)
chr1_peaks = [peak.to_dict() for peak in peaks if peak.chromosome == "chr1"]
ax = epi_viz.plot_chipseq_peaks(chr1_peaks, chromosome="chr1")
```

### With ATAC-seq Analysis
```python
import numpy as np
from metainformant.epigenome.assays import atacseq
from metainformant.epigenome.visualization import visualization as epi_viz

# Load ATAC-seq peaks: a list of ATACPeak objects
peaks = [p for p in atacseq.load_atac_peaks("atac_peaks.narrowPeak") if p.chromosome == "chr1"]

# Peaks are intervals, not a signal track, so build a binned signal from their signal values
bin_size = 1000
positions = np.arange(0, max(p.end for p in peaks) + bin_size, bin_size)
signal = np.zeros(len(positions))
for p in peaks:
    first, last = p.start // bin_size, p.end // bin_size
    signal[first : last + 1] = np.maximum(signal[first : last + 1], p.signal_value)

# Create accessibility profile
ax = epi_viz.plot_atacseq_signal(signal, positions)
```

## Output Options

`output_path` saves the current pyplot figure (it calls `plt.savefig`) at 300 dpi and creates missing parent directories. That is the figure just drawn unless you created other figures in between; when you draw on your own `ax=`, save with `fig.savefig(...)` instead.

```python
from metainformant.epigenome.visualization import visualization as epi_viz
import numpy as np

methylation_data = np.random.default_rng(0).random(500)

# Save to file
ax = epi_viz.plot_methylation_profile(methylation_data, output_path="methylation.png")
```

The interactive browser writes an HTML file instead (it requires plotly, and the file suffix is forced to `.html`):

```python
from metainformant.epigenome.visualization import visualization as epi_viz

epigenetic_tracks = {'methylation': {'positions': [1, 2, 3], 'values': [0.1, 0.5, 0.9]}}
fig = epi_viz.create_interactive_epigenome_browser(epigenetic_tracks, output_path="browser.html")
```

## Data Sources

The plotting functions do not read files; load data with the epigenome module first.

- **BedGraph Files**: `metainformant.epigenome.assays.methylation.load_methylation_bedgraph` returns `{chromosome: [MethylationSite, ...]}`; each site has `position` and `methylation_level`
- **Peak Files**: `chipseq.load_chip_peaks` and `atacseq.load_atac_peaks` read `narrowpeak` (default), `broadpeak` or `bed` files and return `ChIPPeak` / `ATACPeak` objects; call `.to_dict()` on a `ChIPPeak` to get the `start`, `end` and `score` dict that `plot_chipseq_peaks` expects
- **Signal Tracks**: `metainformant.epigenome.analysis.tracks.load_bedgraph_track` reads BedGraph files and `load_bigwig_track` reads bigWig files (which needs `pyBigWig`)
- **Chromatin States**: integer state ids per position, plus an optional `{id: name}` dict

## Performance Considerations

- **Large Genomic Regions**: For whole chromosomes, consider tiling or subsampling
- **Multiple Samples**: Heatmaps work best with 10-100 samples
- **Interactive Plots**: Require Plotly for web-based genome browsers
- **Memory Usage**: High-resolution tracks can be memory-intensive

## Dependencies

- **Required**: matplotlib, numpy, scipy (smoothing in the methylation profile)
- **Optional**: seaborn (the correlation heatmap uses its `heatmap`; a plain matplotlib heatmap is used without it), plotly (interactive browser), pyBigWig (only for `load_bigwig_track`)

## Examples

### Complete Epigenetic Analysis Workflow
```python
import matplotlib.pyplot as plt
import numpy as np
from metainformant.epigenome.assays import chipseq, methylation
from metainformant.epigenome.visualization import visualization as epi_viz

# Load epigenetic data
methylation_data = methylation.load_methylation_bedgraph("sample.methylation.bedgraph")
chip_peaks = chipseq.load_chip_peaks("H3K4me3_peaks.narrowPeak")

sites = methylation_data["chr1"]
positions = np.array([site.position for site in sites])
levels = np.array([site.methylation_level for site in sites])
chr1_peaks = [peak.to_dict() for peak in chip_peaks if peak.chromosome == "chr1"]

# Create comprehensive visualization
fig, axes = plt.subplots(2, 2, figsize=(15, 12))

# Methylation profile
epi_viz.plot_methylation_profile(levels, positions, ax=axes[0, 0])

# ChIP-seq peaks
epi_viz.plot_chipseq_peaks(chr1_peaks, chromosome="chr1", ax=axes[0, 1])

# Multi-track browser (the H3K4me3 values here are random stand-ins; use real signal tracks)
rng = np.random.default_rng(0)
tracks = {
    'methylation': {'data': levels, 'color': 'blue', 'label': 'DNAme'},
    'H3K4me3': {'data': rng.random(len(levels)), 'color': 'red', 'label': 'H3K4me3'}
}
epi_viz.plot_genome_browser_tracks(
    tracks, region_start=int(positions.min()), region_end=int(positions.max()), ax=axes[1, 0]
)

# Correlation analysis (random stand-ins for the other marks)
epi_data = {
    'DNA_methylation': levels,
    'H3K4me3_signal': rng.random(len(levels)),
    'accessibility': rng.random(len(levels))
}
epi_viz.plot_epigenetic_correlation_heatmap(epi_data, ax=axes[1, 1])

plt.tight_layout()
plt.savefig("epigenetic_analysis.png", dpi=300, bbox_inches='tight')
```

## Color Schemes

Recommended color schemes for your own epigenetic figures (the plotting functions above use fixed colormaps and take no color arguments):
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
4. **Color Scaling**: The methylation plots fix the color range at 0-1, so rescale values outside that range first

### Data Requirements

- **Genomic Coordinates**: `positions` for chromatin states must be sorted; plot one chromosome at a time
- **Signal Values**: Should be normalized for comparison across samples
- **Peak Scores**: Higher scores typically indicate stronger enrichment
- **Methylation Levels**: Values should be in 0-1 range for proper color scaling

## Related Documentation

- **[Epigenome Analysis](../epigenome/)**: Core epigenetic analysis functions
- **[Methylation Analysis](../epigenome/methylation.md)**: DNA methylation workflows
- **[ChIP-seq Analysis](../epigenome/chip_seq.md)**: Histone modification analysis
- **[Visualization Integration](integration.md)**: Cross-module visualization patterns

This module provides comprehensive epigenetic visualization capabilities integrated with METAINFORMANT's epigenomic analysis workflows.
