# Visualization

Publication-quality plotting for long-read sequencing data, including read length distributions, quality scatter plots, sequence dotplots, alignment views, methylation tracks, and phasing block diagrams.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports plots submodule |
| `plots.py` | All visualization functions using matplotlib/seaborn |

## Key Functions

| Function | Description |
|----------|-------------|
| `plots.plot_read_length_histogram()` | Plot read length distribution as histogram |
| `plots.plot_quality_vs_length()` | Scatter plot of quality score versus read length |
| `plots.plot_dotplot()` | Sequence self-dotplot or pairwise dotplot |
| `plots.plot_alignment_view()` | Visualize read alignments against a reference |
| `plots.plot_methylation_track()` | Plot methylation levels along a genomic region |
| `plots.plot_phasing_blocks()` | Visualize phase blocks with haplotype assignments |

## Usage

```python
from metainformant.longread.visualization import plots

plots.plot_read_length_histogram(reads, output_path="output/length_dist.png")
plots.plot_quality_vs_length(reads, output_path="output/qv_length.png")
plots.plot_methylation_track(sites, region="chr1:1000-5000", output_path="output/meth.png")
```
