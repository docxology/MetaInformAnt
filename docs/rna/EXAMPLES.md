# Current RNA examples

All examples use the canonical project configuration directory and an
explicit external data root.

## Inspect and validate

```bash
export AMALGKIT_DATA_ROOT=/Volumes/external_drive/Data/amalgkit
uv run python scripts/rna/validate_configs.py
uv run python scripts/rna/validate_all_species_workflow.py
uv run python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run
```

## Run one species

```bash
env -u VIRTUAL_ENV uv run --frozen --no-sync python scripts/rna/process_species.py \
  --species apis_mellifera \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" \
  --cohort "$AMALGKIT_DATA_ROOT/cohort/hymenoptera27-20261005" \
  --workers 4 --threads 8
```

## Run the complete cohort

```bash
AMALGKIT_DATA_ROOT="$AMALGKIT_DATA_ROOT" \
bash projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh
```

The producer derives its worker and thread budgets from the host; override
them with `PIPELINE_WORKERS` and `PIPELINE_THREADS` only after measuring a
bottleneck.

## Observe progress

```bash
uv run python projects/hymenoptera_amalgkit/scripts/report_campaign_status.py \
  --data-root "$AMALGKIT_DATA_ROOT"
uv run python projects/hymenoptera_amalgkit/scripts/generate_pipeline_report.py \
  --data-root "$AMALGKIT_DATA_ROOT"
```

## Cross-species readiness

```bash
uv run python projects/hymenoptera_amalgkit/scripts/prepare_cross_species_inputs.py \
  --data-root "$AMALGKIT_DATA_ROOT" --check
uv run python projects/hymenoptera_amalgkit/scripts/run_cross_species_analysis.py \
  --data-root "$AMALGKIT_DATA_ROOT" --check
```

Use the generated evidence and manuscript status files as the source of
analysis claims. The examples do not imply that every configured species is
complete.
