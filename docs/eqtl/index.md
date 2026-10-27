# eQTL Integration Pipeline

> **Note**: eQTL is a **cross-cutting integration pipeline**. Scan statistics live in `metainformant.gwas.finemapping.eqtl`, plots in `metainformant.gwas.visualization.eqtl_visualization` and multi-omics helpers in `metainformant.multiomics.analysis.integration`; the `metainformant.eqtl` package holds the synthetic and real input builders and the transcriptome SNP-calling pipeline. Scripts are in `scripts/eqtl/`.

## Documentation

| Doc | Description |
|-----|-------------|
| [README.md](README.md) | Overview, architecture, integration scripts |
| [Pipeline Guide](pipeline_guide.md) | Step-by-step walkthrough of transcriptome SNP calling |
| [Configuration Reference](configuration.md) | All YAML config options |

## Quick Start

```bash
# Real Amalgkit expression scanned against SYNTHETIC genotypes at SYNTHETIC gene positions
# (a workflow demonstration: the associations are not biological findings)
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
