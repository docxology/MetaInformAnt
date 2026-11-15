# Amalgkit integration

The RNA package integrates the pinned Amalgkit 0.16.108 command line through a
small command registry and the shared streaming producer.

## Current command chain

```text
metadata → select → getfastq → integrate → quant
         → merge → wsfilter → finalize → sanity
```

The producer owns acquisition and quantification. In frozen-cohort mode
(`--cohort`, the launcher's only mode) `metadata` and `select` are fixed by the
cohort's frozen selection table and never run; without `--cohort` they run
only with `--refresh-metadata`. The project checkpoint runner owns merge,
filtering, finalization, and sanity after the producer has stopped.

## Inspect the installed release

```bash
amalgkit --version
amalgkit --help
uv run python scripts/rna/validate_configs.py
```

## Inspect the configured cohort

```bash
export AMALGKIT_DATA_ROOT=/Volumes/external_drive/Data/amalgkit
uv run python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run
```

## Execute

```bash
bash projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh
```

For a single diagnostic species:

```bash
env -u VIRTUAL_ENV uv run --frozen --no-sync python scripts/rna/process_species.py \
  --species apis_mellifera \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" \
  --cohort "$AMALGKIT_DATA_ROOT/cohort/hymenoptera27-20261005"
```

## Evidence and resume

The producer records SQLite sample state, current metadata and quantification
hashes, and resumable download progress. The downstream runner verifies those
receipts before reusing any matrix. See the [command reference](commands.md),
[workflow stages](steps/README.md), and [project storage contract](../../../projects/hymenoptera_amalgkit/doc/01_infrastructure/02_storage_contract.md).
