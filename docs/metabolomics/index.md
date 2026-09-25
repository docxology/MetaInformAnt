# Metabolomics

## Overview

Metabolomics analysis module for METAINFORMANT. Provides metabolite identification, MGF/CSV mass-spectrometry I/O, metabolite set enrichment, and metabolomics plots; cross-module metabolite-gene integration runs through `multiomics`.

## Sub-packages

| Sub-package | Description |
|-------------|-------------|
| `analysis` | Metabolite identification, normalization, fold change, differential abundance |
| `io` | Mass spectrometry I/O: MGF spectra and CSV tables (no mzML/mzXML readers) |
| `pathways` | Metabolite set enrichment (hypergeometric + Benjamini-Hochberg) |
| `visualization` | Metabolomics-specific plots and figures |

## Output

- Generated outputs belong under `output/` (program-generated results are ephemeral; no `config/metabolomics/` templates exist in this checkout).

## Integration

- **Multi-Omics**: metabolite layers are accepted by `multiomics.analysis.integration.integrate_omics_data`
- **Visualization**: metabolomics plots live in `visualization` (`volcano_plot_data`, `pca_metabolomics`, `intensity_heatmap_data`)

## See Also

- [Source Code](../../src/metainformant/metabolomics/)
- [Agent Rules](../agents/rules/metabolomics.md)
