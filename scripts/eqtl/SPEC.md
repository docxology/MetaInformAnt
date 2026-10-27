# Specification: eqtl

## 🎯 Scope
Functionality for eqtl.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Orchestration Script

## 💾 Data Structures
- **Modules**: 3 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `rna_snp_pipeline.py`
- `run_eqtl_demo.py`
- `run_eqtl_real.py`

## Provenance contract

`run_eqtl_real.py` scans real expression against synthetic genotypes and synthetic gene positions (the chromosome accession is parsed from the transcript ID, the TSS positions are invented, no GFF/GTF is read). It writes `gene_positions_synthetic.tsv` and the `*_synthetic_genotypes.*` outputs (`cis_eqtl_results_synthetic_genotypes.tsv`, `summary_stats_synthetic_genotypes.json`, `run_summary_synthetic_genotypes.json` and three plots), plus `real_expression_matrix.tsv`, `synthetic_genotype_matrix.tsv`, `variant_positions.tsv`, `gene_expression_stats.tsv`, `cis_eqtl_annotated.tsv`, `top_100_eqtls.tsv` and `top_effect_sizes.tsv`. Both JSON summaries carry `expression_data`, `genotype_data` and `gene_position_data`; the summary also carries `data_source`, and the run summary `analysis_type` and `note`. Variant IDs are `var_<full gene id>_<i>`. The scripts write `NA` for missing values in their result tables.
