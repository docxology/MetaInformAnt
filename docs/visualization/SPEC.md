# Specification: visualization

## Scope
Comprehensive documentation for the visualization domain producing 70+ plot types (Manhattan, heatmaps, network graphs, phylogenetic trees, Circos, dashboards, animations) with Matplotlib, Seaborn, Plotly, and Graph-tool backends. Publication-quality vector and raster outputs, accessible color palettes, and reproducible figure generation.

## Architecture
- **Dependency Level**: Domain (depends on core.io, core.config; used by all analysis modules)
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
├── plots/            # Core plotting functions (manhattan, pca, heatmap, volcano...)
│   ├── genomic/      # Manhattan, QQ, regional, ideogram
│   ├── analysis/     # ROC, MA, correlation, enrichment
│   ├── network/      # Spring layout, circular, hive plots
│   ├── phylogeny/    # Tree plotting with bootstrap, rectangular/circular
│   └── composite/    # Multi-panel figure assembly
├── dashboards/       # Interactive Plotly/Streamlit apps
├── animation/        # Trajectory/evolution animations
├── export/           # Format conversion (PDF/PNG/SVG), figure assembly
├── themes.py         # Color palette definitions, style setters
└── __init__.py       # Public API exports
```

## Integration
- **Source**: `src/metainformant/visualization/{plots,dashboards,animation,export}/`
- **Tests**: `tests/test_visualization_plots.py`, `tests/test_visualization_dashboards.py`, `tests/test_visualization_animation.py`
- **Visual regression**: `tests/visualization/baseline/` (PNG references compared via `pytest-mpl`)
- **Dependencies**:
  - **Core**: `metainformant.core.io` — Data loading (TSV/CSV/Parquet)
  - **Core**: `metainformant.core.config` — Theme configuration (colors, fonts, DPI)
  - **DNA**: `metainformant.dna.phylogeny` — Tree objects (Newick)
  - **GWAS**: `metainformant.gwas.visualization` — Manhattan/QQ (wrappers)
  - **RNA**: `metainformant.rna.analysis` — PCA/heatmap integration

## Key Functions

### Genomic Plots
| Function | Module | Purpose | Returns |
|----------|--------|---------|---------|
| `manhattan()` | plots.genomic | GWAS Manhattan plot with significance threshold | `matplotlib.figure.Figure` |
| `qq()` | plots.genomic | QQ plot with genomic inflation λ annotation | `Figure` |
| `regional()` | plots.genomic | Regional association with LD heatmap and gene tracks | `Figure` |
| `ideogram()` | plots.genomic | Chromosome ideogram with cytobands | `Figure` |

### Expression Analysis
| Function | Module | Purpose |
|----------|--------|---------|
| `pca()` | plots.analysis | PCA scatter colored by metadata (tissue, condition) |
| `heatmap()` | plots.analysis | Clustered heatmap with dendrogram (Seaborn clustermap) |
| `volcano()` | plots.analysis | Volcano plot (logFC vs -log10 p) with DE gene labeling |
| `ma()` | plots.analysis | MA plot (mean expression vs log fold-change) |

### Network
| Function | Module | Purpose |
|----------|--------|---------|
| `spring_layout()` | plots.network | Force-directed Fruchterman-Reingold layout |
| `circular()` | plots.network | Circular layout for small networks |
| `hive()` | plots.network | Hive plot for bipartite networks |

### Phylogeny
| Function | Module | Purpose |
|----------|--------|---------|
| `plot_tree()` | plots.phylogeny | Newick tree with branch support, rectangular layout |
| `plot_circular()` | plots.phylogeny | Circular phylogram (Radial) |
| `draw_tanglegram()` | plots.phylogeny | Compare two trees with tanglegram |

### Animation
| Function | Module | Purpose |
|----------|--------|---------|
| `trajectory()` | animation | Pseudotime flow in UMAP/t-SNE embedding |
| `evolution()` | animation | Tree evolution over time (BEAST-style) |

### Dashboards
| Class | Module | Purpose |
|-------|--------|---------|
| `Dashboard` | dashboards | Grid of linked Plotly widgets with callbacks |
| `Explorer` | dashboards | Single-page app with dropdowns/sliders |

## Testing Policy
- **Zero Mock**: All tests use real data; small synthetic figures allowed for edge/error cases only.
- **Visual regression**: Baseline images committed to `tests/visualization/baseline/`; updated via `pytest --mpl-update-baseline`.
- **Performance**: Large-plot benchmarks in `tests/benchmark/test_plot_speed.py`; ensure <5s for 10k×10k heatmap.

## Public API
```python
# Top-level imports
from metainformant.visualization import (
    # Genomic plots
    manhattan, qq, regional, ideogram,
    # Analysis plots
    pca, heatmap, volcano, ma,
    # Networks
    spring_layout, circular, hive,
    # Phylogeny
    plot_tree, plot_circular, tanglegram,
    # Animation
    trajectory, evolution,
    # Dashboards
    Dashboard, Explorer,
    # Utils
    set_style, set_palette, save_figure, to_caption
)
```

## Examples
```python
from metainformant.visualization import manhattan, pca, heatmap
import pandas as pd

# 1. Manhattan
gwas = pd.read_csv("gwas_results.tsv", sep='\t')
fig = manhattan(gwas, genome="amellifera", highlight_genes=["Amel_004123"])
fig.save("manhattan.pdf", dpi=300)

# 2. PCA
expr = pd.read_csv("counts.tsv", sep='\t', index_col=0)
fig = pca(expr, color_by="tissue", shape_by="individual")
fig.save("pca.png", dpi=150)

# 3. Heatmap
deseq = pd.read_csv("deseq2_results.tsv", sep='\t')
top100 = deseq.nsmallest(100, "padj")["gene_id"]
counts = pd.read_csv("norm_counts.tsv", sep='\t', index_col=0).loc[top100]
fig = heatmap(counts, z_score="row", cluster=True)
fig.save("heatmap.pdf")
```

## Related Specifications
- **Core I/O**: [../../core/SPEC.md](../core/SPEC.md) — Data loading and DataFrame I/O
- **DNA phylogeny**: [../dna/SPEC.md](../dna/SPEC.md) — Tree structures used by phylogeny plots
- **GWAS**: [../gwas/SPEC.md](../gwas/SPEC.md) — Manhattan/QQ generation interface
- **Task Reference**: [../../tasks/visualize_results.md](../tasks/visualize_results.md) — Quick examples for common plots
