# Specification: visualization

## Scope
Comprehensive documentation for the visualization domain producing 70+ plot types (Manhattan, heatmaps, network graphs, phylogenetic trees, Circos, dashboards, animations) with Matplotlib, Seaborn, Plotly, and Graph-tool backends. Publication-quality vector and raster outputs, accessible color palettes, and reproducible figure generation.

## Architecture
- **Dependency Level**: Domain (depends on core.io, core.utils.config; used by all analysis modules)
- **Component Type**: Presentation Layer
- **Location**: `docs/visualization/`
- **Source Root**: `src/metainformant/visualization/`

## Data Structures

### Documentation Files
- `README.md`: Module overview, installation, style guide
- `index.md`: Entry point with quick-start guide and plot taxonomy
- `AGENTS.md`: AI attribution, maintenance guidelines, cross-module dependencies
- `SPEC.md`: Technical specification (this file)
- `plots.md`: Function reference by category (genomic, analysis, network, phylogeny, dashboard)
- `gallery.md`: Sample figures (images + full runnable code)
- `export.md`: Output format guidelines (PDF/PNG/SVG/PPTX), DPI requirements, journal specs
- `themes.md`: Color palettes (colorbrewer, viridis, custom), style sheets, font guidelines
- `animation.md`: Animation generation (MP4, GIF) via Matplotlib/Plotly
- `dashboards.md`: Interactive widget construction (Plotly/Streamlit)

### Source Code Structure
```
src/metainformant/visualization/
├── plots/                  # Core plotting functions
│   ├── basic.py            # lineplot, scatter_plot, heatmap, bar_plot, pie_chart, area_plot, step_plot
│   ├── general.py          # volcano_plot, manhattan_plot, qq_plot, pca_plot, expression/correlation heatmaps
│   ├── specialized.py      # Venn, Sankey, chord, alluvial, circular bar, UpSet
│   ├── multidim.py         # Pairwise, parallel coordinates, radar, 3D scatter
│   └── animations.py       # Time-series, evolution, clustering, network, trajectory animations
├── genomics/               # genomics.py (Manhattan, volcano, regional, ideogram), expression.py, networks.py, trees.py
├── analysis/               # dimred, statistical, timeseries, information and quality-control plots
├── dashboards/             # composite.py (multi-panel figures), interactive.py (Plotly figures)
├── interactive_dashboards/ # dashboards.py (JSON-serializable interactive plot data, HTML export)
├── config/                 # themes.py, palettes.py, conventions.py (Okabe-Ito colors, deterministic export)
└── __init__.py             # Exposes the subpackages plus the WONG and OKABE_ITO palettes
```

## Integration
- **Source**: `src/metainformant/visualization/{plots,genomics,analysis,dashboards,interactive_dashboards,config}/`
- **Tests**: `tests/visualization/test_visualization_plots.py`, `tests/visualization/test_visualization_dashboards.py`, `tests/visualization/test_visualization_animation.py`
- **Visual regression**: `tests/visualization/baseline/` (PNG references compared via `pytest-mpl`)
- **Dependencies**:
  - **Core**: `metainformant.core.io` — Data loading (TSV/CSV/Parquet)
  - **Core**: `metainformant.core.utils.config` — Theme configuration (colors, fonts, DPI)
  - **DNA**: `metainformant.dna.phylogeny` — Tree objects (Newick)
  - **GWAS**: `metainformant.gwas.visualization` — Manhattan/QQ (wrappers)
  - **RNA**: `metainformant.rna.analysis` — PCA/heatmap integration

## Key Functions

### Genomic Plots
| Function | Module | Purpose | Returns |
|----------|--------|---------|---------|
| `manhattan_plot()` | genomics.genomics | GWAS Manhattan plot with genome-wide significance line | `matplotlib.axes.Axes` |
| `qq_plot()` | analysis.statistical | Q-Q plot against a normal or uniform distribution | `Axes` |
| `regional_plot()` | genomics.genomics | Regional association plot for a chromosome interval | `Axes` |
| `chromosome_ideogram()` | genomics.genomics | Chromosome length ideogram | `Axes` |

### Expression Analysis
| Function | Module | Purpose |
|----------|--------|---------|
| `pca_plot()` | plots.general | Scatter of pre-computed principal components colored by metadata (tissue, condition) |
| `expression_heatmap()` | plots.general | Gene expression heatmap (Seaborn, matplotlib fallback) |
| `volcano_plot()` | genomics.genomics | Volcano plot (logFC vs -log10 p) with significance coloring |

### Network
| Function | Module | Purpose |
|----------|--------|---------|
| `plot_network_force_directed()` | genomics.networks | Force-directed Fruchterman-Reingold layout |
| `plot_network_circular()` | genomics.networks | Circular layout for small networks |

### Phylogeny
| Function | Module | Purpose |
|----------|--------|---------|
| `plot_phylo_tree()` | genomics.trees | Tree (NetworkX DiGraph, dict-based or Bio.Phylo) in a hierarchical layout |
| `circular_tree_plot()` | genomics.trees | Circular phylogram |
| `tree_comparison_plot()` | genomics.trees | Side-by-side comparison of two trees |

### Animation
| Function | Module | Purpose |
|----------|--------|---------|
| `animate_trajectory()` | plots.animations | Animated point trajectories with trails |
| `animate_evolution()` | plots.animations | Sequence evolution across generations |

### Dashboards
| Function | Module | Purpose |
|----------|--------|---------|
| `multi_panel()` | dashboards.composite | Arrange several plots into one matplotlib figure |
| `interactive_scatter()` | dashboards.interactive | Interactive Plotly scatter plot with hover data |

## Testing Policy
- **Real Implementation**: All tests use real data; small synthetic figures allowed for edge/error cases only.
- **Visual regression**: Baseline images committed to `tests/visualization/baseline/`; updated via `pytest --mpl-update-baseline`.
- **Performance**: Large-plot benchmarks in `tests/benchmark/test_plot_speed.py`; ensure <5s for 10k×10k heatmap.

## Public API
```python-snippet
# Top-level imports: subpackages and the WONG / OKABE_ITO palettes
from metainformant.visualization import analysis, config, dashboards, genomics, interactive_dashboards, plots

# Plot functions live in their defining submodules
from metainformant.visualization.genomics.genomics import (
    # Genomic plots
    manhattan_plot, regional_plot, chromosome_ideogram, volcano_plot,
)
from metainformant.visualization.analysis.statistical import qq_plot
from metainformant.visualization.plots.general import pca_plot, expression_heatmap
from metainformant.visualization.genomics.networks import (
    # Networks
    plot_network_force_directed, plot_network_circular,
)
from metainformant.visualization.genomics.trees import (
    # Phylogeny
    plot_phylo_tree, circular_tree_plot, tree_comparison_plot,
)
from metainformant.visualization.plots.animations import (
    # Animation
    animate_trajectory, animate_evolution,
)
from metainformant.visualization.dashboards.composite import multi_panel
from metainformant.visualization.config.themes import apply_theme
from metainformant.visualization.config.palettes import categorical
from metainformant.visualization.config.conventions import save_figure_deterministic
```

## Examples
```python-snippet
from metainformant.visualization.genomics.genomics import manhattan_plot
from metainformant.visualization.plots.general import pca_plot, expression_heatmap
import pandas as pd

# 1. Manhattan (columns CHR, BP, P by default; saved at 300 DPI)
gwas = pd.read_csv("gwas_results.tsv", sep='\t')
ax = manhattan_plot(gwas, output_path="manhattan.pdf")

# 2. PCA (pre-computed components with their variance explained, in percent)
pcs = pd.read_csv("pca_scores.tsv", sep='\t')  # PC1, PC2, PC1_variance, PC2_variance, tissue
ax = pca_plot(pcs, hue="tissue", output_path="pca.png")

# 3. Heatmap
deseq = pd.read_csv("deseq2_results.tsv", sep='\t')
top100 = deseq.nsmallest(100, "padj")["gene_id"]
counts = pd.read_csv("norm_counts.tsv", sep='\t', index_col=0).loc[top100]
ax = expression_heatmap(counts, output_path="heatmap.pdf")
```

## Related Specifications
- **Core I/O**: [../../core/SPEC.md](../core/SPEC.md) — Data loading and DataFrame I/O
- **DNA phylogeny**: [../dna/SPEC.md](../dna/SPEC.md) — Tree structures used by phylogeny plots
- **GWAS**: [../gwas/SPEC.md](../gwas/SPEC.md) — Manhattan/QQ generation interface
- **Task Reference**: [../../tasks/visualize_results.md](../tasks/visualize_results.md) — Quick examples for common plots
