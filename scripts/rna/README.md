# RNA workflow scripts

These scripts operate the current MetaInformAnt/Amalgkit workflow and keep
large data outside the repository. The producer uses ENA-first acquisition
with the configured NCBI/SRA fallback and records resumable evidence.

## Execution

| Script | Purpose |
|---|---|
| `run_all_species.py` | Configuration-derived multi-species launcher |
| `process_species.py` | Current single-species wrapper using the same orchestrator as the all-species run |
| `check_pipeline_status.py` | Inspect SQLite state and downstream evidence |
| `clean_external_artifacts.py` | Manifest-backed archival of superseded artifacts and explicitly selected legacy quant directories |

## Validation and evidence

| Script | Purpose |
|---|---|
| `validate_all_species_workflow.py` | Validate configuration and plan order |
| `validate_configs.py` | Validate YAML schema and required inputs |
| `check_environment.py` | Check Python, Amalgkit, and external tools |
| `check_pipeline_status.py` | Inspect run-state and downstream evidence |
| `verify_rna.py` | Run RNA module checks |
| `sync_quant_results.py` | Copy compact review outputs to a review directory |

## Preparation and analysis

| Script | Purpose |
|---|---|
| `setup_genome.py` | Download and prepare transcriptome/index assets |
| `batch_genome_index.py` | Prepare indexes for a declared configuration set |
| `normalize_tissue_metadata.py` | Apply declared tissue mappings |
| `generate_orthologs.py` | Build an ortholog mapping with provenance |
| `verify_tissue_coverage.py` | Report mapped and unmapped tissue labels |

## Quick start

```bash
export AMALGKIT_DATA_ROOT=/Volumes/blue/data/amalgkit
uv run python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run
uv run python scripts/rna/process_species.py \
  --species apis_mellifera \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run
uv run python scripts/rna/validate_all_species_workflow.py
bash projects/hymenoptera_amalgkit/scripts/verify_setup.sh \
  --data-root "$AMALGKIT_DATA_ROOT" --require-data --report
```

Use `--dry-run` and inspect the species inventory before starting downloads.
For the full scientific analysis, continue with the executable workflow and
evidence commands in `projects/hymenoptera_amalgkit`.
