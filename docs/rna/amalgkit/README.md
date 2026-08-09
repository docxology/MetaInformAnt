# Amalgkit integration

The RNA package integrates the pinned Amalgkit 0.16.38 command line through a
small command registry and the shared streaming producer.

## Current command chain

```text
metadata → select → getfastq → integrate → quant
         → merge → wsfilter → finalize → sanity
```

The producer owns acquisition and quantification. The project checkpoint
runner owns merge, filtering, finalization, and sanity after the producer lock
is released.

## Inspect the installed release

```bash
amalgkit --version
amalgkit --help
uv run python scripts/rna/validate_configs.py
```

## Inspect the configured cohort

```bash
export AMALGKIT_DATA_ROOT=/Volumes/blue/data/amalgkit
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
uv run python scripts/rna/process_species.py \
  --species apis_mellifera \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT"
```

## Evidence and resume

The producer records SQLite sample state, current metadata and quantification
hashes, and resumable download progress. The downstream runner verifies those
receipts before reusing any matrix. See the [command reference](commands.md),
[workflow stages](steps/README.md), and [project storage contract](../../../projects/hymenoptera_amalgkit/doc/01_infrastructure/02_storage_contract.md).
