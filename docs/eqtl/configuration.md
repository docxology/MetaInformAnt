# eQTL Pipeline: Configuration Reference

This document describes all configuration options for the eQTL Transcriptome SNP Calling Pipeline.

## Configuration File

Config files live in `config/eqtl/` and are YAML format. Example: [`eqtl_amellifera.yaml`](../../config/eqtl/eqtl_amellifera.yaml).

## Options

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `species` | string | *required* | Species name matching the amalgkit output directory |
| `samples.mode` | string | `auto` | `auto` to discover from completed quant, or `explicit` |
| `samples.max_samples` | int | `5` | Max samples when mode=auto |
| `samples.explicit_ids` | list | — | SRR/ERR/DRR IDs when mode=explicit |
| `alignment.tool` | string | `hisat2` | Aligner (only hisat2 supported) |
| `alignment.threads` | int | `4` | CPU threads for alignment |
| `variant_calling.min_base_qual` | int | `20` | Minimum base quality for pileup |
| `variant_calling.min_map_qual` | int | `20` | Minimum mapping quality |
| `variant_calling.ploidy` | int | `2` | Ploidy for genotype calling |
| `filtering.min_qual` | int | `30` | Minimum variant QUAL score |
| `filtering.min_depth` | int | `10` | Minimum read depth |
| `output.base_dir` | string | `output/eqtl/{species}` | Output directory |
| `output.cleanup_fastq` | bool | `true` | Delete FASTQs after alignment |
| `eqtl_scan.enabled` | bool | `false` | Auto-run eQTL scan after calling |
| `eqtl_scan.cis_window` | int | `500000` | cis-eQTL window size (bp) |

## Example Usage

```bash
# Run with config file
uv run python scripts/eqtl/rna_snp_pipeline.py --config config/eqtl/eqtl_amellifera.yaml

# Override species and sample count from CLI
uv run python scripts/eqtl/rna_snp_pipeline.py --species amellifera --n-samples 3
```
