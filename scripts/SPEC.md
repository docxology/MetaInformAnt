# Script specification

## Scope

Executable wrappers and utilities for MetaInformAnt workflows. Domain logic
belongs in `src/metainformant`; scripts own CLI parsing, explicit paths, and
evidence-safe orchestration.

## RNA production surface

| Script | Contract |
|---|---|
| `scripts/rna/run_all_species.py` | Discover the canonical config inventory and run the shared producer |
| `scripts/rna/process_species.py` | Run one named species through the same producer |
| `scripts/rna/check_pipeline_status.py` | Inspect SQLite progress and downstream evidence |
| `scripts/rna/validate_configs.py` | Validate current Amalgkit YAML and selection rules |
| `scripts/rna/validate_all_species_workflow.py` | Read-only inventory and nine-stage plan validation |
| `projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh` | Lock-owned producer then downstream checkpoint campaign |
| `projects/hymenoptera_amalgkit/scripts/run_all_finalization.sh` | Lock-owned downstream merge, filter, finalize, and sanity checkpoints |

## Current invocation contract

```bash
export AMALGKIT_DATA_ROOT=/Volumes/external_drive/Data/amalgkit
uv run python scripts/rna/run_all_species.py \
  --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
  --data-root "$AMALGKIT_DATA_ROOT" --dry-run

bash projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh
```

All commands are resumable and data-root scoped. A producer lock, SQLite
progress database, current provenance receipts, and downstream evidence
manifest are part of the completion contract.

## Script conventions

- use `argparse` and expose `--help`;
- resolve configuration and data roots explicitly;
- write large artifacts outside Git;
- make repeated discovery, validation, and checkpoint calls idempotent;
- test real filesystem behavior and current CLI contracts.
