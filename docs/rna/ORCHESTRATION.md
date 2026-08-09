# RNA orchestration

The current RNA execution model has one producer and one downstream checkpoint
owner for a selected data root. The producer is implemented by
`StreamingPipelineOrchestrator` and is exposed through the cohort and
single-species scripts.

## State flow

```mermaid
flowchart LR
    C[Project YAML inventory] --> P[run_all_species.py]
    P --> M[metadata]
    M --> S[select]
    S --> G[getfastq]
    G --> I[integrate]
    I --> Q[quant]
    Q --> L[SQLite progress and provenance]
    L --> R[producer lock released]
    R --> D[merge]
    D --> W[wsfilter]
    W --> F[finalize]
    F --> Y[sanity]
    Y --> E[evidence manifest and analysis]
```

## Inspect the cohort

```bash
export AMALGKIT_DATA_ROOT=/Volumes/blue/data/amalgkit
uv run python scripts/rna/validate_all_species_workflow.py
uv run python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run
```

The dry run is read-only. It resolves the same config directory, species
inventory, and data root that a real producer run will use.

## Execute the complete campaign

```bash
bash projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh
```

The project launcher acquires a lock, validates resource floors, starts the
producer, waits for it to finish, and then invokes the downstream checkpoint
runner only when the producer exits successfully. It preserves the progress
database, partial downloads, quantification receipts, and logs on interruption.

Run one species through the same producer when a bounded diagnostic is needed:

```bash
uv run python scripts/rna/process_species.py \
  --species apis_mellifera \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT"
```

Do not start this command while the full campaign lock is held.

## Downstream checkpoints

After the producer has stopped and the cohort gate passes:

```bash
bash projects/hymenoptera_amalgkit/scripts/run_all_finalization.sh \
  --data-root "$AMALGKIT_DATA_ROOT"
```

The runner executes `merge`, `wsfilter`, `finalize`, and `sanity` in order. Each
stage reuses a current valid checkpoint and otherwise writes its output through
an atomic temporary path. A failed stage leaves its inputs intact and records
the failure for a bounded retry.

## Status and evidence

```bash
uv run python projects/hymenoptera_amalgkit/scripts/report_campaign_status.py \
  --data-root "$AMALGKIT_DATA_ROOT"
uv run python projects/hymenoptera_amalgkit/scripts/generate_pipeline_report.py \
  --data-root "$AMALGKIT_DATA_ROOT"
```

Status is observational. Completion requires valid metadata, selected sample
IDs, readable abundance files, contract provenance, downstream stage receipts,
and a regenerated evidence manifest. See the [storage contract](../../projects/hymenoptera_amalgkit/doc/01_infrastructure/02_storage_contract.md)
and [analysis readiness](../../projects/hymenoptera_amalgkit/doc/manuscript/analysis_readiness.md).
