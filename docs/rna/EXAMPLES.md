# Current RNA examples

All examples use the canonical project configuration directory and an
explicit external data root.

## Inspect and validate

```bash
export AMALGKIT_DATA_ROOT=/Volumes/blue/data/amalgkit
uv run python scripts/rna/validate_configs.py
uv run python scripts/rna/validate_all_species_workflow.py
uv run python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run
```

## Run one species

```bash
uv run python scripts/rna/process_species.py \
  --species apis_mellifera \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" \
  --workers 4 --threads 8 --quant-slots 2
```

## Run the complete cohort

```bash
AMALGKIT_DATA_ROOT="$AMALGKIT_DATA_ROOT" \
AMALGKIT_PIPELINE_WORKERS=8 \
AMALGKIT_PIPELINE_THREADS=8 \
bash projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh
```

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
