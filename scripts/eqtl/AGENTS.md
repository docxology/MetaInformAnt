# Agent Directives: scripts/eqtl

## Role
Thin orchestrator scripts for eQTL integration pipelines.

## Contents
- `run_eqtl_real.py` — Run eQTL analysis with real Amalgkit RNA-seq quantification data
- `run_eqtl_demo.py` — Demonstrate eQTL logic with synthetic data
- `rna_snp_pipeline.py` — Call SNP variants from transcriptome RNA-seq data

## Rules
- Scripts are thin wrappers that call `metainformant.gwas.finemapping.eqtl` library code
- Follow REAL IMPLEMENTATION policy
- Use `uv` for Python dependencies
