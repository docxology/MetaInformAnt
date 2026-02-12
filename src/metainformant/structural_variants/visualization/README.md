# Visualization

Publication-quality plots for structural variant analysis, including Circos-style genome-wide views, coverage tracks with SV overlays, size distributions, breakpoint detail views, and CNV profiles.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports plots submodule |
| `plots.py` | All SV visualization functions using matplotlib/seaborn |

## Key Functions

| Function | Description |
|----------|-------------|
| `plots.plot_circos()` | Circos-style genome-wide SV visualization |
| `plots.plot_coverage_track()` | Coverage track with SV overlay annotations |
| `plots.plot_sv_size_distribution()` | Histogram of SV sizes by type |
| `plots.plot_sv_type_summary()` | Bar chart summary of SV types and counts |
| `plots.plot_breakpoint_detail()` | Detailed breakpoint view with read support |
| `plots.plot_cnv_profile()` | Genome-wide CNV log2 ratio profile |

## Usage

```python
from metainformant.structural_variants.visualization import plots

plots.plot_circos(variants, chrom_sizes, output_path="output/circos.png")
plots.plot_coverage_track(coverage, variants, region="chr1:1000-50000",
                          output_path="output/coverage.png")
plots.plot_sv_size_distribution(variants, output_path="output/sv_sizes.png")
plots.plot_cnv_profile(segments, output_path="output/cnv_profile.png")
```
