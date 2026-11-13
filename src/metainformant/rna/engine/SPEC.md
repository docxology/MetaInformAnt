# Specification: engine

## 🎯 Scope
Current RNA/Amalgkit execution engine.

`ProgressDB(..., journal_mode=...)` accepts only `WAL` or `DELETE`, defaulting
to `AMALGKIT_PROGRESS_DB_JOURNAL_MODE` or WAL. Mode selection is validated
before directory or database creation. DELETE must actually apply rollback
journaling, FULL synchronization and zero page mmap; refusal closes the
connection.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 36 Python modules
- **Key Concepts**: Refer to Pydantic models in source. `progress_db.py`
  owns the sample state machine (`pending → downloading → downloaded →
  quantifying → quantified`, with `failed` branches and provenance-audit
  `quarantined` transitions) and the
  failure-class contract: `classify_sample_error()` maps stored
  `samples.error` text via first-match order-sensitive markers
  (`SAMPLE_ERROR_CLASSES`) to `environment_write_denied`,
  `environment_access_denied`, `environment_read_only`,
  `environment_no_space`, `environment_missing_tool`,
  `transfer_all_sources_failed`,
  `extraction_timeout`, `quantification_timeout`, `quantification_failed`,
  `quantification_exception`, `unrecorded` (no error stored), or
  `unclassified` (no marker matched); the campaign status report surfaces
  the counts as `db_failure_classes`, and
  `metainformant.rna.analysis.cohort_accounting.build_cohort_funnel()`
  consumes the DB read-only for fail-closed cohort funnel accounting with
  durable failure-class reason codes. `streaming_orchestrator.py` holds the
  `campaign_lock.py` writer flock of the data root for the whole run, runs the
  mandatory campaign preflight before discovery, resets orphaned in-flight
  rows and environment failures at startup, reconciles outputs in one
  change-proportional pass, stops admission through an environment circuit
  breaker, returns a run summary with exit code 0/75/69 (73 for a busy
  root, 70 internal) and supports the opt-in
  quant stall watchdog (`AMALGKIT_PIPELINE_QUANT_STALL_TIMEOUT_SECONDS`,
  default 0/disabled) around each quantification subprocess. In frozen
  cohort mode (`cohort.py`) discovery consumes only the hash-verified cohort
  and passes a separate task table to `amalgkit quant`; `prune_species()` is
  report-only and `migrate_to_cohort()` binds a progress store to a cohort
  without deleting rows. `quant_import.py` imports verified outputs from other
  data roots without ever replacing a destination sample.

## 🔌 API Definition

### Exports
- `fastq_compression`: strict pigz level (1–9, default 6) and thread settings;
  optional level 1 for temporary acquisition scratch. Native FASTQ validation
  gates promotion.
- `acquisition_manifest`: `load_snapshot()` reads a frozen acquisition
  snapshot (`snapshot.json`, schema `metainformant.rna.acquisition_snapshot.v1`)
  and refuses a `manifest.jsonl` whose sha256 differs; `sha256_file()` is the
  shared file digest.
- `acquisition_order`, `acquisition_size_overlay`: deterministic acquisition
  order over frozen compressed FASTQ weights (see below).
- `__init__.py`
- `discovery.py`
- `campaign_lock.py`
- `cohort_host.py`
- `exclusions.py`
- `pipeline.py`
- `preflight.py`
- `progress_dashboard.py`
- `progress_db.py`
- `provenance.py`
- `quant_snapshot.py`
- `raw_cleanup.py`
- `raw_reclamation.py`
- `species.py`
- `sra_extraction.py`
- `streaming_orchestrator.py`
- `workflow.py`
- `workflow_cleanup.py`
- `workflow_core.py`
- `workflow_execution.py`
- `workflow_planning.py`
- `workflow_steps.py`

### Portable restoration contract

`durable_quant.restore_quantification` optionally accepts `verified_config_path`
and requires a frozen reference-index binding in that mode. Local copies are
recorded under a contract-bound restoration witness; `provenance.classify_quantification`
rejects tampered witnesses, changed input checksums and paths escaping the sample.
Original quantification sidecar bytes are preserved.

## Acquisition admission contract

An optional per-orchestrator whole-sample semaphore is selected with
`AMALGKIT_PIPELINE_ACQUISITION_SLOTS` (0 through 256, zero disables it).
Its separate finite positive queue timeout defaults to 7,200 seconds via
`AMALGKIT_PIPELINE_ACQUISITION_QUEUE_TIMEOUT_SECONDS`. Invalid settings raise
before opening the progress database. Queue exhaustion raises `StagingDeadline`
without acquisition I/O; leases release on success, refusal and exceptions.
The absolute first-pass acquisition deadline begins after admission and covers
file-slot waits, metadata, transfers and validation. The quantification stage
runs after release. The host profile sets acquisition slots independently of
task workers, file connections and quantifiers. Queue and acquisition bounds
add; this is not an end-to-end task deadline or a scientific processing
exclusion rule.

## Global small-sample admission

`plan_acquisition_order(..., global_size_bands=True)` applies inclusive upper
bands of 100 MiB, 500 MiB, 1/2/5/20 GiB, followed by larger known sizes and a
final unknown-size band. Species round-robin applies within the active band.
Without either flag the order is species-fair. Duplicate identities and
malformed sizes fail; checkpoint IDs are interleaved and unique. Inputs are not
mutated. Order is deterministic under input permutation.

The producer's scheduling order (`AMALGKIT_PIPELINE_SCHEDULING_ORDER`, default
`size_bands`) maps `size_bands` to `global_size_bands=True`, `smallest_first`
to `global_size_order=True` and `species_order` to neither. Strict size order
is one globally ascending acquisition lane with task-ID tie breaks and
unknown/nonpositive weights last. Size bands and strict size order are mutually
exclusive. The weights are frozen total compressed estimates, not remaining
wire bytes. This is an admission priority; separate pools, capacity
constraints and concurrent execution do not promise globally ordered starts or
completions, and the order is never a cohort filter, partial deletion,
deadline reset, global completion barrier or promise of higher wire throughput.

`effective_fastq_bytes(task, sizes)` returns a task's frozen `fastq_bytes`;
an optional `sizes` map (passed as `plan_acquisition_order(...,
size_overrides=...)`) may fill only a zero weight, with a positive integer
below 2**63, and refuses anything else. The task records are never rewritten.
