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
| `src/metainformant/eqtl/workflow/synthetic.py` | Synthetic genotype and gene-position builders, real expression loader |
| `scripts/eqtl/run_eqtl_real.py` | Real expression with synthetic genotypes and gene positions (outputs named `*_synthetic_genotypes.*` and `gene_positions_synthetic.tsv`) |

## Rules

- Never describe `run_eqtl_real.py` output as an analysis of real genotypes or real gene positions: only the expression is real. The summaries carry `expression_data`, `genotype_data` and `gene_position_data` to say so.
- IDs passed to `cis_eqtl_scan()` must be unique (`ValueError` otherwise). A failed fit has NaN statistics with `status` and `error`; tables write missing values as `NA`.
