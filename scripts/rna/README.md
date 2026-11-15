# RNA workflow scripts

These scripts operate the current MetaInformAnt/Amalgkit workflow and keep
large data outside the repository. The producer uses ENA-first acquisition
with the configured NCBI/SRA fallback and records resumable evidence.

## Execution

| Script | Purpose |
|---|---|
| `run_all_species.py` | Multi-species streaming producer (the campaign launcher is `projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh`) |
| `process_species.py` | Current single-species wrapper using the same orchestrator as the all-species run |
| `check_pipeline_status.py` | Inspect SQLite state and downstream evidence |
| `clean_external_artifacts.py` | Manifest-backed archival of superseded artifacts and, explicitly, quant directories without a provenance sidecar |

## Cohort and data-root maintenance

| Script | Purpose |
|---|---|
| `prepare_cohort_host.py` | Prepare a fresh data root for frozen-cohort quantification on a new host |
| `import_quantifications.py` | Import verified quantifications from other data roots into a campaign data root |
| `export_quant_snapshot.py` | Export verified, reusable quantification outputs to a portable snapshot |
| `reclaim_quantified_raw.py` | Reclaim raw inputs only for verified, reusable quantification outputs |
| `reconcile_quantification_compatibility.py` | Classify and safely restore existing quantification outputs |
| `campaign_watcher.py` | Detached, read-only campaign watcher |
| `benchmark_campaign_scale.py` | Campaign-scale performance benchmark on disposable fixtures |
| `new_clade_scaffold.py` | Scaffold a sibling clade project from the Hymenoptera template pattern |

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
| `batch_genome_index.py` | Prepare indexes for a declared configuration set |
| `normalize_tissue_metadata.py` | Apply declared tissue mappings |
| `generate_orthologs.py` | Build an ortholog mapping with provenance |
| `verify_tissue_coverage.py` | Report mapped and unmapped tissue labels |

## Quick start

```bash
export AMALGKIT_DATA_ROOT=/Volumes/external_drive/Data/amalgkit
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
