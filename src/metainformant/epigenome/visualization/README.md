# Epigenome Visualization

Plotting functions for epigenomic data: methylation profiles, ChIP-seq peak landscapes, ATAC-seq signals, histone modification heatmaps, and interactive genome browser views.

## Contents

| File | Purpose |
|------|---------|
| `visualization.py` | All epigenome visualization functions (matplotlib, seaborn, plotly) |

## Key Functions

| Function | Description |
|----------|-------------|
| `plot_methylation_profile()` | Per-base methylation level across a genomic region |
| `plot_chipseq_peaks()` | ChIP-seq peak signal landscape with annotations |
| `plot_atacseq_signal()` | ATAC-seq accessibility signal across regions |
| `plot_histone_modification_heatmap()` | Multi-mark histone modification comparison |
| `plot_differential_methylation()` | Volcano/MA plot of differentially methylated regions |
| `plot_chromatin_states()` | Color-coded chromatin state segmentation diagram |
| `plot_epigenetic_correlation_heatmap()` | Correlation matrix across epigenetic marks |
| `plot_genome_browser_tracks()` | Multi-track genome browser style figure |
| `plot_dna_methylation_clusters()` | Clustered methylation patterns across samples |
| `create_interactive_epigenome_browser()` | Plotly-based interactive genome browser |

## Usage

```python
from metainformant.epigenome.visualization.visualization import (
    plot_methylation_profile, plot_chipseq_peaks
)

fig = plot_methylation_profile(meth_data, region="chr1:1M-2M", output_path="output/meth.png")
fig = plot_chipseq_peaks(peaks, output_path="output/chip_peaks.png")
```
