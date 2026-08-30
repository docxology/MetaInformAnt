# Getting started with current RNA production

## 1. Verify the environment

```bash
uv run python scripts/rna/check_environment.py
uv run python scripts/rna/validate_configs.py
uv run python scripts/rna/validate_all_species_workflow.py
```

The pinned project release is Amalgkit 0.16.60. Confirm `amalgkit --help`,
`amalgkit merge --help`, `amalgkit wsfilter --help`, and
`amalgkit finalize --help` before a real run.

## 2. Select and inspect the data root

```bash
export AMALGKIT_DATA_ROOT=/Volumes/external_drive/Data/amalgkit
mkdir -p "$AMALGKIT_DATA_ROOT"
df -h "$AMALGKIT_DATA_ROOT"
uv run python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run
```

The data root contains raw reads, indexes, abundance tables, logs, SQLite
state, and evidence. Do not put these artifacts in Git.

## 3. Run the project campaign

```bash
bash projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh
```

The script is the only full-campaign entrypoint. It coordinates one producer,
one lock, and the downstream checkpoint sequence for the canonical 27-species
configuration inventory.

For a single-species bounded run:

```bash
uv run python scripts/rna/process_species.py \
  --species pogonomyrmex_barbatus \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT"
```

## 4. Check progress

```bash
uv run python projects/hymenoptera_amalgkit/scripts/report_campaign_status.py \
  --data-root "$AMALGKIT_DATA_ROOT"
```

The report distinguishes configured, discovered, selected, quantified,
finalized, and analyzed species. It preserves failures and missing evidence.

## 5. Prepare cross-species analysis

Only finalized species that pass the current matrix and provenance gates may
enter the cross-species lane:

```bash
uv run python projects/hymenoptera_amalgkit/scripts/prepare_cross_species_inputs.py \
  --data-root "$AMALGKIT_DATA_ROOT" --check
uv run python projects/hymenoptera_amalgkit/scripts/run_cross_species_analysis.py \
  --data-root "$AMALGKIT_DATA_ROOT" --check
```

Use the project manuscript checklist before making biological claims.
