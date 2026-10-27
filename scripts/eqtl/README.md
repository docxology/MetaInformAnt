# eQTL Integration Scripts

Thin orchestrator scripts for the eQTL integration pipeline, bridging GWAS variants with Amalgkit RNA-seq expression data.

## Scripts

| Script | Purpose |
|--------|---------|
| `run_eqtl_real.py` | cis-eQTL scan of real Amalgkit expression against synthetic genotypes at synthetic gene positions |
| `run_eqtl_demo.py` | Demo with fully synthetic data |
| `rna_snp_pipeline.py` | SNP calling from transcriptome RNA-seq |

## Usage

```bash
# Scan real Amalgkit expression against synthetic genotypes and gene positions
uv run python scripts/eqtl/run_eqtl_real.py

# Explore with synthetic data
uv run python scripts/eqtl/run_eqtl_demo.py

# Call SNPs from RNA-seq
uv run python scripts/eqtl/rna_snp_pipeline.py --species amellifera --n-samples 3
```

## What is real and what is synthetic

`run_eqtl_real.py` reads real expression (Amalgkit/kallisto TPM). The genotypes are synthetic random dosages, because no genotype data are matched to these RNA-seq samples. The gene positions are synthetic too: only the chromosome accession is parsed from each transcript ID (`NC_...`), and the TSS positions are invented (evenly spaced 10 kb apart in the order of the IDs; no GFF/GTF is read). The associations demonstrate the workflow and are not biological findings.

The synthetic gene positions, the scan, the two summaries and the three plots say so in their names: `gene_positions_synthetic.tsv`, `cis_eqtl_results_synthetic_genotypes.tsv`, `summary_stats_synthetic_genotypes.json`, `run_summary_synthetic_genotypes.json`, `volcano_plot_synthetic_genotypes.png`, `summary_plot_synthetic_genotypes.png` and `top_eqtl_boxplot_synthetic_genotypes.png`. `variant_positions.tsv`, `cis_eqtl_annotated.tsv`, `top_100_eqtls.tsv` and `top_effect_sizes.tsv` are derived from the same synthetic inputs and keep their plain names. The two JSON summaries carry the provenance keys `expression_data` (starts with `real`), `genotype_data` and `gene_position_data` (both start with `synthetic`); the summary also has `data_source`, and the run summary `analysis_type` and `note`.

Both analysis scripts write their result tables (the scan, the annotated and top-hit tables, the effect sizes) with `NA` for missing values, so a fit that could not be made (NaN statistics, `status` `"error"`) is never written as a number. Variant IDs are `var_<full gene id>_<i>`, unique whenever the gene IDs are.

## Related

- [eQTL Docs](../../docs/eqtl/README.md) — Integration pipeline documentation
- [GWAS Fine-Mapping](../../src/metainformant/gwas/finemapping/) — Library code
