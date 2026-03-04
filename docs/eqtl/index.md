# eQTL Integration Pipeline

> **Note**: eQTL is a **cross-cutting integration pipeline** — it does not have its own `src/metainformant/eqtl/` module. Core logic lives in `metainformant.gwas.finemapping.eqtl`, `metainformant.gwas.visualization.eqtl_visualization`, and `metainformant.multiomics.analysis.integration`. Scripts are in `scripts/eqtl/`.

## Documentation

| Doc | Description |
|-----|-------------|
| [README.md](README.md) | Overview, architecture, integration scripts |
| [Pipeline Guide](pipeline_guide.md) | Step-by-step walkthrough of transcriptome SNP calling |
| [Configuration Reference](configuration.md) | All YAML config options |

## Quick Start

```bash
# Run real-data eQTL analysis (uses Amalgkit expression data)
uv run python scripts/eqtl/run_eqtl_real.py

# Run synthetic demo (for testing/validation)
uv run python scripts/eqtl/run_eqtl_demo.py

# Run transcriptome SNP calling
uv run python scripts/eqtl/rna_snp_pipeline.py --species amellifera --n-samples 3
```

## Related

- [GWAS Pipeline](../gwas/index.md)
- [RNA-seq Pipeline](../rna/index.md)
- [Multi-Omics Integration](../multiomics/index.md)
