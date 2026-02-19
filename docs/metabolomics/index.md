# Metabolomics

## Overview

Metabolomics analysis module for METAINFORMANT. Provides tools for metabolite identification, mass spectrometry data processing, pathway mapping, and metabolite-gene integration analysis.

## Sub-packages

| Sub-package | Description |
|-------------|-------------|
| `analysis` | Metabolite identification and quantification |
| `io` | Mass spectrometry data I/O (mzML, mzXML, CSV) |
| `pathways` | Metabolic pathway mapping and enrichment |
| `visualization` | Metabolomics-specific plots and figures |

## Configuration

- **Config Prefix**: `METAB_` (e.g., `METAB_THREADS`, `METAB_WORK_DIR`)
- **Config Path**: `config/metabolomics/default.yaml`
- **Output Path**: `output/metabolomics/<analysis_type>/`

## Integration

- **Multi-Omics**: Metabolite-gene integration via `multiomics`
- **Networks**: Metabolic pathway networks via `networks`
- **Visualization**: Metabolomics plots via `visualization`
- **Quality**: QC metrics for mass spectrometry data via `quality`

## See Also

- [Source Code](../../src/metainformant/metabolomics/)
- [Agent Rules](../agents/rules/metabolomics.md)
