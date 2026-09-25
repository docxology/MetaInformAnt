# Agent Directives: docs/metabolomics

## Role
Documentation for the metabolomics analysis module.

## Module Scope
- Mass spectrometry file I/O: MGF spectra and CSV metabolite × sample tables (no mzML/mzXML readers)
- Metabolite identification, normalization, fold change, and differential abundance
- Metabolite set enrichment analysis (hypergeometric + Benjamini-Hochberg)
- Cross-module metabolite-gene integration via `multiomics` (`integrate_omics_data` accepts metabolomics layers)
- Visualization (volcano plots, PCA ordination, concentration heatmaps)

## Key Source Files
- `src/metainformant/metabolomics/io/` - MGF and CSV readers/writers, spectra filtering, chromatogram extraction
- `src/metainformant/metabolomics/analysis/` - Identification, normalization, fold change, differential abundance
- `src/metainformant/metabolomics/pathways/` - Metabolite set enrichment
- `src/metainformant/metabolomics/visualization/` - Plotting utilities
