# Personal AI Infrastructure (PAI) - visualization

## Context & Intent
- **Path**: `docs/visualization`
- **Purpose**: Comprehensive documentation for producing publication-ready figures, interactive dashboards, animations, and multi-panel composite figures covering all METAINFORMANT analysis outputs (genomics, transcriptomics, proteomics, networks).
- **Domain**: docs

## Virtual Hierarchy
- **Type**: Documentation
- **Parent**: `docs`

## Maintenance Notes
- **System**: Part of the METAINFORMANT Domain layer (visualization and reporting).
- **Style**: Consistent Matplotlib/Seaborn themes, accessible color palettes, reproducible figures.
- **Stability**: API stable; new plot types added as submodules (e.g., `visualization/circos/`).

## AI Workflows
- **Modification**: Run visual regression tests (`pytest tests/visualization/`) before committing figure changes.
- **Documentation**: Each plot type needs a gallery entry (image + code) in `docs/visualization/gallery.md`.

## Cross-References
- **Main docs**: [index.md](../visualization/index.md) — Module overview, installation, style guide
- **API reference**: [plots.md](../visualization/plots.md) — Function signatures by plot category (genomic, analysis, network, phylogeny, dashboard)
- **Gallery**: [gallery.md](../visualization/gallery.md) — Visual examples with full source code
- **Export**: export.md — Format/size/DPI guidelines for journals and posters
- **Related modules**: [DNA](../dna/index.md) for phylogeny trees; [GWAS](../gwas/index.md) for Manhattan/QQ; [RNA](../rna/index.md) for heatmaps/PCA; [Quality](../quality/index.md) for QC plots
