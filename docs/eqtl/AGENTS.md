# Agent Directives: docs/eqtl

## Role
Documentation for eQTL (expression Quantitative Trait Loci) integration pipeline.

## Contents
- `README.md` - eQTL pipeline overview and architecture
- `pipeline_guide.md` - Step-by-step transcriptome SNP calling walkthrough
- `configuration.md` - YAML configuration reference

## Key Topics
- cis-eQTL scanning and association testing
- Transcriptome SNP calling from RNA-seq data (HISAT2, bcftools)
- Expression matrix loading and normalization
- Multiple testing correction
- Volcano plots and effect-size boxplots
- Integration with Amalgkit (RNA) and GWAS (DNA) pipelines
- Synthetic genotype generation for method validation

## Key Source Files
- `src/metainformant/gwas/finemapping/eqtl.py` - Core statistical scanning
- `src/metainformant/gwas/visualization/eqtl_visualization.py` - Plotting utilities
- `scripts/eqtl/` - Workflow orchestration scripts
