# Current RNA workflow

The current per-species plan is fixed and explicit:

```text
metadata → select → getfastq → integrate → quant
         → merge → wsfilter → finalize → sanity
```

The first five stages are run by the resumable streaming producer. The final
four are run by the project checkpoint runner after the producer lock is
released.

## Read-only plan validation

```bash
uv run python scripts/rna/validate_all_species_workflow.py
uv run python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run
```

## Production execution

```bash
export AMALGKIT_DATA_ROOT=/Volumes/external_drive/Data/amalgkit
bash projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh
```

The launcher supplies the worker, thread, quant-slot, extraction-slot, and
validation-slot budgets. Resource settings must remain bounded by the
external-volume and system-volume free-space checks.

## Checkpoint execution

```bash
bash projects/hymenoptera_amalgkit/scripts/run_all_finalization.sh \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run
```

The checkpoint runner validates the cohort database and stage-specific
provenance before reusing an output. It is safe to repeat after an interrupted
run; valid checkpoints are reused and incomplete checkpoints are regenerated
atomically.

## Completion evidence

A species is complete only when all nine stage receipts, selected-ID
reconciliation, valid abundance files, current quantification sidecars, final
matrix checks, and the regenerated project evidence manifest agree. A process
exit code or populated directory alone is insufficient.
