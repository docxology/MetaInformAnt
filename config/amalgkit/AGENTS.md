# Agent Directives: config/amalgkit

## Role

Production-ready amalgkit RNA-seq workflow configurations for automated transcript quantification pipelines.

## Contents

| File | Description |
| :--- | :--- |
| `amalgkit_template.yaml` | **Reference**: 400+ line template with all options documented |
| `amalgkit_test.yaml` | Minimal test configuration for validation |
| `amalgkit_pbarbatus.yaml` | **Production**: Full P. barbatus dataset (95/110 quantified) |
| `amalgkit_apis_mellifera_all.yaml` | **Production**: Full A. mellifera dataset (~7,270 samples) |
| `amalgkit_cross_species.yaml` | Cross-species TMM normalization config |
| `tissue_mapping.yaml` | Canonical tissue name synonyms for normalization |
| `tissue_patches.yaml` | Per-bioproject/sample tissue overrides |

## Configuration Structure

```yaml
# Core paths (relative to repo root)
work_dir: output/amalgkit/{species}/work
log_dir: output/amalgkit/{species}/logs
threads: 12

# Species identification
species_list:
  - Pogonomyrmex_barbatus
taxon_id: 144034

# Reference genome
genome:
  accession: GCF_000187915.1
  dest_dir: output/amalgkit/shared/genome/Pogonomyrmex_barbatus

# Step-specific parameters
steps:
  getfastq:
    redo: no           # Skip already-downloaded
    keep_fastq: no     # Delete after quant
  quant:
    redo: no           # Skip already-quantified
    index_dir: ...     # Reuse kallisto index
```

## Critical Patterns

### Stream-and-Clean (Disk Management)

For large datasets with limited disk space:

```yaml
steps:
  getfastq:
    redo: no           # Resume capability
  quant:
    keep_fastq: no     # Immediate cleanup
    redo: no           # Idempotent
```

### Shared Resources

Reuse genome/index across configs:

```yaml
genome:
  dest_dir: output/amalgkit/shared/genome/Pogonomyrmex_barbatus
steps:
  quant:
    index_dir: output/amalgkit/shared/genome/Pogonomyrmex_barbatus/index
```

### Metadata Filtering

Filter to RNA-Seq + Illumina to prevent genomic samples leaking in:

```yaml
steps:
  metadata:
    search_string: '"Species"[Organism] AND "RNA-Seq"[Strategy] AND "Illumina"[Platform]'
```

## Adding New Species

1. Copy `amalgkit_template.yaml` → `amalgkit_{species}.yaml`
2. Update `species_list`, `taxon_id`, and `genome.accession`
3. Adjust paths: `work_dir`, `log_dir`, `genome.dest_dir`
4. **Validation**: Run `python3 scripts/rna/validate_configs.py` to ensure schema compliance.
5. Test with small sample subset first (use `max_sample: 5`)
6. Scale to full dataset after validation

## Validation and Testing

- **Config Validation**: usage of `scripts/rna/validate_configs.py` is mandatory for all new configurations.
- **Zero-Mock Policy**: All Amalgkit tests strictly adhere to the Zero-Mock policy, ensuring real functional verification of the CLI and environment.

## Environment Overrides

Prefix with `AK_`:

```bash
export AK_THREADS=16
export AK_WORK_DIR=/fast/storage/amalgkit
export NCBI_EMAIL=your@email.com
```
