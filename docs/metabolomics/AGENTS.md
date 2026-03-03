# Agent Directives: docs/metabolomics

## Role
Documentation for the metabolomics analysis module.

## Module Scope
- Mass spectrometry data processing (mzML, mzXML, CSV)
- Metabolite identification and peak quantification
- Data normalization and scaling
- KEGG/Reactome pathway mapping
- Metabolite set enrichment analysis
- Metabolite-gene integration
- Visualization (volcano plots, PCA ordination, concentration heatmaps)

## Key Source Files
- `src/metainformant/metabolomics/io/` - MS file reading and format conversion
- `src/metainformant/metabolomics/analysis/` - Identification, quantification, normalization
- `src/metainformant/metabolomics/pathways/` - Pathway mapping and enrichment
- `src/metainformant/metabolomics/visualization/` - Plotting utilities
