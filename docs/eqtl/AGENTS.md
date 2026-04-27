# Agent Directives: docs/eqtl

## Role

Documentation for the eQTL integration pipeline.

## Module Scope

Expression quantitative trait locus analysis bridging GWAS variants with Amalgkit RNA-seq expression data. Cross-cutting module with logic in gwas/finemapping/ and multiomics/.

## Key Source Files

| Path | Description |
|------|-------------|
| `src/metainformant/gwas/finemapping/eqtl.py` | cis/trans eQTL scanning |
| `src/metainformant/gwas/finemapping/colocalization.py` | eQTL-GWAS colocalization |
| `scripts/eqtl/run_eqtl_real.py` | Real-data pipeline orchestrator |
