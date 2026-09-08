# RNA engine

The engine implements the current Amalgkit workflow as a bounded producer,
SQLite progress store, typed workflow planner, and hash-bound evidence layer.
The producer refuses to start into a broken environment: `run_all()` executes
the campaign preflight before any discovery or scheduling work.

## Components

| Module | Purpose |
|---|---|
| `streaming_orchestrator.py` | ENA-first metadata, acquisition, integration, and quantification |
| `preflight.py` | Mandatory start-of-run environment preflight (data-root write probe, amalgkit CLI resolution) |
| `progress_db.py` | Concurrent-safe SQLite sample state, exclusions, and resume queries |
| `progress_dashboard.py` | Mosaic PDF/PNG dashboard of per-species sample processing status rendered from the progress DB |
| `exclusions.py` | CLI for recording and inspecting sample exclusions |
| `workflow.py` | Public re-export hub for configuration, planning, and execution |
| `workflow_core.py` | Typed configuration, path mapping, validation, and defaults |
| `workflow_planning.py` | Current nine-stage plan and parameter resolution |
| `workflow_execution.py` | Amalgkit step execution and result records |
| `workflow_steps.py` | Individual current step implementations |
| `workflow_cleanup.py` | Provenance-gated raw cleanup and disk checks |
| `provenance.py` | Hash-bound metadata, quantification, and downstream receipts |
| `raw_cleanup.py` | Provenance-gated per-sample reclamation of raw FASTQ/SRA inputs once a current-contract quantification exists |
| `species.py` | Shared config and data-root discovery |
| `pipeline.py` | Matrix and downstream table helpers |
| `discovery.py` | Read-only species and sample discovery utilities |
| `sra_extraction.py` | SRA fallback extraction helpers |

## Key interfaces

- `StreamingPipelineOrchestrator.run_all()` starts the bounded producer for a
  declared config set. It first runs `run_campaign_preflight()`
  (data-root write probe plus bare `amalgkit` PATH resolution) and refuses to
  start on failure; the same check is available standalone via
  `python -m metainformant.rna.engine.preflight --data-root <root>`.
  A start-of-run preflight prevents the observed 2026-09-03 failure class in
  which a producer without external-volume write access failed thousands of
  tasks before its first successful write.
- `quant_sample()` bounds each quantification with
  `AMALGKIT_PIPELINE_QUANT_TIMEOUT_SECONDS` (default 7200) and, when
  `AMALGKIT_PIPELINE_QUANT_STALL_TIMEOUT_SECONDS` is set to a positive
  number of seconds (default 0, disabled), with an output-growth watchdog:
  the sample's quant output directory plus any local scratch workspace is
  size-sampled every 30 seconds, and a batch showing no growth for that
  duration is terminated early and recorded as
  `Quant stalled (no output growth for ...)` instead of holding its slot
  for the full timeout. Opt-in because kallisto can legitimately spend long
  stretches in its EM phase before writing abundance files; see the
  [performance contract](../../../../projects/hymenoptera_amalgkit/doc/01_infrastructure/04_performance_and_resume.md).
- `classify_sample_error()` (in `progress_db.py`) maps stored sample
  failure text to a durable class via the order-sensitive markers in
  `SAMPLE_ERROR_CLASSES` (first match wins): `environment_write_denied`,
  `environment_missing_tool`, `transfer_all_sources_failed`,
  `extraction_timeout`, `quantification_timeout`, `quantification_failed`,
  `quantification_exception`; text matching no marker is `unclassified`,
  and a missing/empty error is `unrecorded`. Terminal-failure audits can
  thus separate environmental damage (fully retryable) from genuine
  per-sample failures; the campaign status report surfaces the counts as
  `db_failure_classes`.
- `build_cohort_funnel()` in `metainformant.rna.analysis.cohort_accounting`
  consumes the progress DB read-only plus the per-species amalgkit config
  directory and produces a fail-closed `FunnelReport`: stages `configured`,
  `with_progress`, `quantified_runs`, `failed_runs`, `excluded_runs`,
  `pending_runs`, `active_runs` in `STAGE_NAMES` order, plus `reason_codes`
  classified with `classify_sample_error()`. Missing inputs, a DB without a
  `samples` table, or a DB above the optional `max_gb` guard raise
  `CohortFunnelError`; `FunnelReport.to_tsv()` writes a byte-deterministic
  two-column TSV and `render_funnel_lines()` renders
  `cohort_funnel_<stage>: <count>` summary lines.
- `ProgressDB` stores the states `pending`, `downloading`, `downloaded`,
  `quantifying`, `quantified`, `quarantined`, and `failed`, plus
  quantification compatibility audit records and durable `sample_exclusions`
  rows (`permanent_drop` removes accessions from task eligibility;
  `re_download` marks a stale transfer for a fresh ENA fetch without blocking
  eligibility). Record them with
  `python -m metainformant.rna.engine.exclusions`.
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
