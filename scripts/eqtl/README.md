# eQTL Integration Scripts

Thin orchestrator scripts for the eQTL integration pipeline, bridging GWAS variants with Amalgkit RNA-seq expression data.

## Scripts

| Script | Purpose |
|--------|---------|
| `run_eqtl_real.py` | Full eQTL analysis with real quantification data |
| `run_eqtl_demo.py` | Demo with synthetic data |
| `rna_snp_pipeline.py` | SNP calling from transcriptome RNA-seq |

## Usage

```bash
# Run eQTL with real Amalgkit quantification data
uv run python scripts/eqtl/run_eqtl_real.py

# Explore with synthetic data
uv run python scripts/eqtl/run_eqtl_demo.py

# Call SNPs from RNA-seq
uv run python scripts/eqtl/rna_snp_pipeline.py --species amellifera --n-samples 3
```

## Related

- [eQTL Docs](../../docs/eqtl/README.md) — Integration pipeline documentation
- [GWAS Fine-Mapping](../../src/metainformant/gwas/finemapping/) — Library code
