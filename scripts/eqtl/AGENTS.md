# Agent Directives: scripts/eqtl

## Role
Thin orchestrator scripts for eQTL integration pipelines.

## Contents
- `run_eqtl_real.py` — Scan real Amalgkit RNA-seq expression against synthetic genotypes at synthetic gene positions (the TSS positions are invented; no GFF/GTF is read)
- `run_eqtl_demo.py` — Demonstrate eQTL logic with fully synthetic data
- `rna_snp_pipeline.py` — Call SNP variants from transcriptome RNA-seq data

## Rules
- Scripts are thin wrappers that call `metainformant.gwas.finemapping.eqtl` library code
- Label provenance truthfully: only the expression of `run_eqtl_real.py` is real. The synthetic gene positions, the scan, the summaries and the plots carry `synthetic` in their names (`gene_positions_synthetic.tsv`, `*_synthetic_genotypes.*`), the tables derived from them (`variant_positions.tsv`, `cis_eqtl_annotated.tsv`, `top_100_eqtls.tsv`, `top_effect_sizes.tsv`) do not, and the summaries record `expression_data`, `genotype_data` and `gene_position_data`
- Write tables with `na_rep="NA"` so a failed fit (NaN statistics) is never written as a number
- Follow REAL IMPLEMENTATION policy
- Use `uv` for Python dependencies
