# RNA engine

The engine implements the current Amalgkit workflow as a bounded producer,
SQLite progress store, typed workflow planner, and hash-bound evidence layer.

## Components

| Module | Purpose |
|---|---|
| `streaming_orchestrator.py` | ENA-first metadata, acquisition, integration, and quantification |
| `progress_db.py` | Concurrent-safe SQLite sample state, exclusions, and resume queries |
| `exclusions.py` | CLI for recording and inspecting sample exclusions |
| `workflow.py` | Public re-export hub for configuration, planning, and execution |
| `workflow_core.py` | Typed configuration, path mapping, validation, and defaults |
| `workflow_planning.py` | Current nine-stage plan and parameter resolution |
| `workflow_execution.py` | Amalgkit step execution and result records |
| `workflow_steps.py` | Individual current step implementations |
| `workflow_cleanup.py` | Provenance-gated raw cleanup and disk checks |
| `provenance.py` | Hash-bound metadata, quantification, and downstream receipts |
| `species.py` | Shared config and data-root discovery |
| `pipeline.py` | Matrix and downstream table helpers |
| `discovery.py` | Read-only species and sample discovery utilities |
| `sra_extraction.py` | SRA fallback extraction helpers |

## Key interfaces

- `StreamingPipelineOrchestrator.run_all()` starts the bounded producer for a
  declared config set.
- `ProgressDB` stores the states `pending`, `downloading`, `downloaded`,
  `quantifying`, `quantified`, `quarantined`, and `failed`, plus quantification
  compatibility audit records and durable `sample_exclusions` rows
  (`permanent_drop` removes accessions from task eligibility; `re_download`
  marks a stale transfer for a fresh ENA fetch without blocking eligibility).
  Record them with `python -m metainformant.rna.engine.exclusions`.
- `plan_workflow()` resolves the fixed per-species chain:
  `metadata → select → getfastq → integrate → quant → merge → wsfilter → finalize → sanity`.
- `provenance.py` rejects missing, stale, or hash-mismatched receipts.

## Python example

```python
from pathlib import Path

from metainformant.rna.engine.progress_db import ProgressDB
from metainformant.rna.engine.streaming_orchestrator import StreamingPipelineOrchestrator

data_root = Path("/Volumes/external_drive/Data/amalgkit")
db = ProgressDB(data_root / "pipeline_progress.db")
orchestrator = StreamingPipelineOrchestrator(
    config_dir=Path("projects/hymenoptera_amalgkit/config/amalgkit"),
    log_dir=data_root / "logs",
    db_path=data_root / "pipeline_progress.db",
)
```

The project shell entrypoint supplies bounded resource budgets and owns the
producer/downstream lock boundary. See the [running guide](../../../../projects/hymenoptera_amalgkit/doc/00_setup/04_running_the_pipeline.md)
and [storage contract](../../../../projects/hymenoptera_amalgkit/doc/01_infrastructure/02_storage_contract.md).
