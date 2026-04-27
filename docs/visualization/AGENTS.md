# Agent Directives: docs/visualization

## Role

Documentation for the visualization module generating 70+ plot types including Manhattan plots, heatmaps, network graphs, phylogenetic trees, multi-omics Circos plots, interactive dashboards, animations, and publication-quality figures across all biological domains.

## Module Scope

End-to-end plotting workflows:
- **Genomic visualizations**: Manhattan (GWAS), QQ plots, regional association, LocusZoom-style, circos (multi-omics)
- **Expression analysis**: PCA, MDS, heatmaps, volcano plots, MA plots, time-series
- **Network graphics**: Protein-protein interaction (STRING), regulatory networks, pathway diagrams, force-directed layouts
- **Phylogenetics**: Tree plotting with bootstrap support, clade coloring, rectangular/circular layouts
- **Statistical plots**: Violin, box, scatter, correlation matrix, ROC/AUC, precision-recall
- **Dashboarding**: Interactive Plotly/Streamlit apps, linked brushing, widget grids
- **Animations**: Pseudotime trajectories, evolutionary simulations, 3D rotations

## Key Source Files

| Path | Description |
|------|-------------|
| `src/metainformant/visualization/plots/` | Core Matplotlib/Seaborn-based plotting functions (manhattan, pca, heatmap, volcano) |
| `src/metainformant/visualization/genomics/` | Genomic-specific plots: Manhattan, LocusZoom, ideograms, coverage tracks |
| `src/metainformant/visualization/analysis/` | Statistical plots: ROC, MA, correlation, enrichment, survival curves |
| `src/metainformant/visualization/dashboards/` | Interactive widgets: Plotly-based exploration tools |
| `src/metainformant/visualization/animation/` | Matplotlib/Plotly animation generators (trajectories, evolutionary sequences) |
| `src/metainformant/visualization/export/` | Output formatters: PDF (vector), PNG (raster), SVG (web), PowerPoint (presentation) |

## Cross-Module Dependencies

### Upstream
- **Core I/O**: `metainformant.core.io` — Load TSV/CSV/Parquet for plotting
- **Core config**: `metainformant.core.config` — Theme/style/config (DARK/LIGHT, COLORMAPS)
- **DNA**: `metainformant.dna.phylogeny` — Tree data structures

### Downstream (consumers of visualizations)
- **RNA**: `metainformant.rna.analysis` — PCA/heatmap plots for expression
- **GWAS**: `metainformant.gwas.visualization` — Manhattan/QQ outputs
- **Quality**: `metainformant.quality.reporting` — QC report figures
- **MCP**: Planning to expose plot-generating tools to LLMs

## Rule Constraints
1. **Vector first** — Default output format is PDF (vector) for publications; PNG for web.
2. **Accessible colors** — Use colorblind-safe palettes (`colorbrewer` Set2, viridis); avoid red-green alone.
3. **Reproducible style** — Theme via `set_style("meta")` for consistent branding.
4. **Caption-ready** — Every figure must support `fig.to_caption()` generating LaTeX/GMarkdown caption snippets.

## Related Tasks & Guides
- **Task reference**: [../../tasks/visualize_results.md](../tasks/visualize_results.md) — Quick examples: manhattan, heatmap, network, tree
- **API reference**: `docs/visualization/plots.md` — Function signatures for all plot types
- **Gallery**: `docs/visualization/gallery.md` — Image thumbnails and full sample code
- **Export**: `docs/visualization/export.md` — DPI, size, format guidelines for journals
- **Community**: QA for plot issues in `docs/quality/plots.md` (if any)

## Maintenance Notes
- **Matplotlib API**: Use object-oriented interface (not pyplot) for composability.
- **Testing**: Visual regression tests via `pytest-mpl`; baseline images in `tests/visualization/baseline/`.
- **Performance**: Large heatmaps (>10k genes) require `rasterized=True`; see `performance_tuning.md`.
- **Dependencies**: Matplotlib 3.8+, Seaborn 0.13+, Plotly 5.18+ (optional).

## AI Assistant Guidance
When adding new plot type:
1. Create function in `src/metainformant/visualization/plots/<type>.py`
2. Add to `__init__.py` registry
3. Write docstring with `Parameters`, `Returns`, `Example` sections
4. Add sample figure to `docs/visualization/gallery.md`
5. Link from task page `docs/tasks/visualize_results.md`
6. Add pytest test with `@pytest.mark.mpl_image_compare`
