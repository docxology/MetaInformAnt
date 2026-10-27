# Specification: eqtl

## Scope
Documentation for the eqtl domain in MetaInformAnt.

## Architecture
- **Dependency Level**: Documentation
- **Component Type**: Integration
- **Location**: `docs/eqtl/`

## Data Structures
- **Files**:
  - `README.md`: Overview
  - `AGENTS.md`: AI Attribution
  - `SPEC.md`: This file
  - `*.md`: Topic-specific guides

## Integration
- **Source**: `src/metainformant/gwas/finemapping/eqtl.py` (primary), `src/metainformant/gwas/finemapping/colocalization.py`, `src/metainformant/gwas/analysis/eqtl.py`, `src/metainformant/eqtl/` (input builders, transcriptome SNP calling)
- **Scripts**: `scripts/eqtl/` (pipeline orchestrators)
- **Tests**: `tests/gwas/test_gwas_analysis_eqtl.py`, `tests/eqtl/` (`test_eqtl_failed_fits.py`, `test_eqtl_unique_variant_ids.py`, `test_eqtl_real_script_labels.py`), `tests/integration/test_eqtl_integration.py`

## Contracts

- `scripts/eqtl/run_eqtl_real.py` scans real expression against synthetic genotypes and synthetic gene positions. The chromosome accession is parsed from the transcript ID and the TSS positions are invented (no GFF/GTF). The gene positions are written as `gene_positions_synthetic.tsv`, and the scan, the summaries and the plots are named `*_synthetic_genotypes.*`; the derived tables (`variant_positions.tsv`, `cis_eqtl_annotated.tsv`, `top_100_eqtls.tsv`, `top_effect_sizes.tsv`) keep plain names. The summaries carry `expression_data`, `genotype_data` and `gene_position_data`.
- Synthetic variant IDs are `var_<full gene id>_<i>`; `cis_eqtl_scan()` raises `ValueError` on a repeated ID in the gene positions, the variant positions or the genotype index.
- A failed fit has NaN statistics, `status` and `error` (`cis_eqtl_scan`, `conditional_eqtl`, `eqtl_effect_sizes`); `eqtl_summary_stats()` reports `n_tests` (tests with a p-value) and `n_tests_failed`; `write_eqtl_results()` and the scripts write `NA` for missing values.

## Testing Policy
- **Real Implementation**: All tests must use real implementations. Mocks are strictly prohibited.
