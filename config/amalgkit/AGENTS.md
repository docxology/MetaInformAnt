# Agent Directives: config/amalgkit

## Role

Production-ready amalgkit RNA-seq workflow configurations for automated transcript quantification pipelines.

## Contents

| File | Description |
|------|-------------|
| `amalgkit_template.yaml` | **Reference**: 400+ line template with all options documented |
| `amalgkit_test.yaml` | Minimal test configuration for validation |
| `amalgkit_pbarbatus_5sample.yaml` | 5-sample quick test |
| `amalgkit_pbarbatus_25sample.yaml` | 25-sample robustness validation |
| `amalgkit_pbarbatus_all.yaml` | **Production**: Full 110-sample P. barbatus dataset |
| `amalgkit_pogonomyrmex_barbatus.yaml` | Species-specific reference config |

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

## Adding New Species

1. Copy `amalgkit_template.yaml` → `amalgkit_{species}.yaml`
2. Update `species_list`, `taxon_id`, and `genome.accession`
3. Adjust paths: `work_dir`, `log_dir`, `genome.dest_dir`
4. Test with small sample subset first (use `max_sample: 5`)
5. Scale to full dataset after validation

## Environment Overrides

Prefix with `AK_`:

```bash
export AK_THREADS=16
export AK_WORK_DIR=/fast/storage/amalgkit
export NCBI_EMAIL=your@email.com
```
