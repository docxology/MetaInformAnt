# Specification: engine

## 🎯 Scope
Current RNA/Amalgkit execution engine.

`configured_progress_db_path(default, configured)` retains the default when
unset and otherwise requires an existing absolute nonsymlink regular file
without parent traversal. Acquisition workers use it for an explicit
`AMALGKIT_PROGRESS_DB_PATH`; continuous runtime and scope receipts record the
selected path and journal policy. `ProgressDB(..., journal_mode=...)` accepts
only `WAL` or `DELETE`, defaulting to `AMALGKIT_PROGRESS_DB_JOURNAL_MODE` or WAL.
Mode selection is validated before directory or database creation. DELETE
must actually apply rollback journaling, FULL synchronization and zero page
mmap; refusal closes the connection. Recovery migration belongs to a
separately reviewed operator and must retain all schema, rows and original
database sidecars. Changing the operational store does not reset attempts.

The opt-in local frontier persists stage observations and bounds planned plus
admitted identities separately from acquisition, extraction and quantification
slots. Raw reservations use frozen compressed-size weights and are retained
until the finished callback completes its durable checkpoint. They supplement,
not replace, free-space checks and SSD extraction/staging reservations. Queued
planned items own no native worker or acquisition deadline. STOP leaves them
unsubmitted and drains already admitted work. Restart revalidates current raw
and quantification evidence; frontier events alone never establish success.

`run_acquisition_frontier(reuse=...)` accepts an optional read-only probe returning
a verified reused success or `None`. Probes are serial and attempted once per
planned identity only when admission, raw reservation and STOP checks permit.
The manifest worker binds this to current-runtime, content-digest and frozen
configuration checks. A hit emits `reusing`, journals and checkpoints before
releasing ownership; callback failure propagates. A miss keeps its normal
acquisition path. Completed native futures are drained between individual probe
attempts, including misses, without waiting for another network completion.
Quant payload and restored-witness paths are checked for regular-file ownership
and containment before the classifier reads content. Unsafe files produce a
probe miss before the success journal. The summary's
`reused_without_acquisition` is a subset of reused work, never new output.

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 18 Python modules
- **Key Concepts**: Refer to Pydantic models in source. `progress_db.py`
  owns the sample state machine (`pending → downloading → downloaded →
  quantifying → quantified`, with `failed` branches and provenance-audit
  `quarantined` transitions) and the
  failure-class contract: `classify_sample_error()` maps stored
  `samples.error` text via first-match order-sensitive markers
  (`SAMPLE_ERROR_CLASSES`) to `environment_write_denied`,
  `environment_missing_tool`, `transfer_all_sources_failed`,
  `extraction_timeout`, `quantification_timeout`, `quantification_failed`,
  `quantification_exception`, `unrecorded` (no error stored), or
  `unclassified` (no marker matched); the campaign status report surfaces
  the counts as `db_failure_classes`, and
  `metainformant.rna.analysis.cohort_accounting.build_cohort_funnel()`
  consumes the DB read-only for fail-closed cohort funnel accounting with
  durable failure-class reason codes. `streaming_orchestrator.py` runs the
  mandatory campaign preflight before discovery and supports the opt-in
  quant stall watchdog (`AMALGKIT_PIPELINE_QUANT_STALL_TIMEOUT_SECONDS`,
  default 0/disabled) around each quantification subprocess.

## 🔌 API Definition

`campaign_status.reconcile()` partitions frozen eligible tasks into cloud stages
and recorded local SQLite stages with species/column marginals. File coverage
is a separate presence partition; transfer gaps are cloud-receipt tasks without
complete canonical local files. Duplicate/malformed observations fail closed.
`campaign_status_cli.main()` collects a bounded live SSM/S3/EC2 observation and
writes timestamped Markdown/JSON and sample/transfer TSVs without changing
producer state. It does not certify receipt contents or biological readiness.

### Exports
- `fastq_compression`: strict pigz level (1–9, default 6) and thread settings;
  optional level 1 for temporary acquisition scratch. Native FASTQ validation
  still gates promotion; worker results and AWS jobs retain the selected level.
- `acquisition_prerequisites`: native metadata batch/quantifier classification and
  executable dependency validation before acquisition or durable-store access.
  Explicit 16S titles paired with PCR selection and GENOMIC or METAGENOMIC
  source metadata remain operationally unresolved, with the normalized source
  retained in the reason; frozen identities are preserved.
  Metadata reference aliases retain the frozen index digest. Unsupported long-read
  FASTA/MMI/tool bindings remain unresolved; no backend substitution is allowed.
- `acquisition_manifest`, `acquisition_references`, `acquisition_sample`,
  `acquisition_worker`, `acquisition_worker_cli`, `acquisition_lanes`: frozen local/AWS execution;
  configurations/input hashes and exclusive root ownership checked before work;
  per-run journals retained and reused/newly quantified counts distinguished.
  `run_manifest()` invokes an optional per-sample `on_quantified` checkpoint
  after journaling and before replacement admission. Checkpoint errors propagate.
  An optional `should_stop` predicate drains submitted tasks and reports
  `stopped_early` plus `unsubmitted_task_count` without admitting the remaining queue.
- `acquisition_scheduling`: explicit rate assumptions, positive workload counts,
  aggregate stage-duration planning and remaining-time admission. New AWS batches
  fit the job target or isolate an oversized task within the maximum duration;
  unknown/unsupported work remains unresolved. Original admissions are immutable.
- `acquisition_allocation`, `acquisition_snapshot`, `acquisition_estimates`,
  `acquisition_pricing`, `acquisition_cli`: immutable disjoint lane planning,
  conflict-refusing staging and explicit gross cost/time scenarios. No estimator
  grants execution permission or substitutes for runtime cost reservations.
- `__init__.py`
- `discovery.py`
- `exclusions.py`
- `pipeline.py`
- `preflight.py`
- `progress_dashboard.py`
- `progress_db.py`
- `provenance.py`
- `raw_cleanup.py`
- `species.py`
- `sra_extraction.py`
- `streaming_orchestrator.py`
- `workflow.py`
- `workflow_cleanup.py`
- `workflow_core.py`
- `workflow_execution.py`
- `workflow_planning.py`
- `workflow_steps.py`

### Operational observation contract

`campaign_observation.capture_campaign_observation` requires explicit task
selection and journal paths; it captures immutable controls and append-only
prefixes with bounded stability retries. Producer startup totals and new
restore/provenance/abundance-bound successes have separate fields. Missing or
stale transfer measurements are unavailable. `campaign_observation_report`
writes isolated, non-overwriting bundles with artifact hashes, tables, canonical
figures and an optional Pandoc working supplement. Every bundle declares
`scientific_release_ready=false`; no cohort or biological inference is implied.
`campaign_observation_cli` owns argument parsing for both thin script adapters.

### Portable restoration contract

`durable_quant.restore_quantification` optionally accepts `verified_config_path`
and requires a frozen reference-index binding in that mode. Local copies are
recorded under a contract-bound restoration witness; `provenance.classify_quantification`
rejects tampered witnesses, changed input checksums and paths escaping the sample.
Original quantification sidecar bytes are preserved. `verify_locked_campaign`
accepts an optional configuration directory for portable full-cohort restoration.
# Acquisition admission contract

An optional per-orchestrator whole-sample semaphore is selected with
`AMALGKIT_PIPELINE_ACQUISITION_SLOTS` (0 through 256, zero disables it).
Its separate finite positive queue timeout defaults to 7,200 seconds via
`AMALGKIT_PIPELINE_ACQUISITION_QUEUE_TIMEOUT_SECONDS`. Invalid settings raise
before opening the progress database. Queue exhaustion raises `StagingDeadline`
without acquisition I/O; leases release on success, refusal and exceptions.
The existing absolute first-pass acquisition deadline begins after admission,
and still covers file-slot waits, metadata, transfers and validation. The
quantification stage runs after release. The local profile uses four acquisition
slots independently of twelve task workers, eight file connections and four
quantifiers. Queue and acquisition bounds add; this is not an end-to-end task
deadline or a scientific processing exclusion rule.

## Global small-sample admission

`plan_acquisition_order(..., global_size_bands=True)` applies inclusive upper
bands of 100 MiB, 500 MiB, 1/2/5/20 GiB, followed by larger known sizes and a
final unknown-size band. Species round-robin remains within the active band.
The default false setting preserves legacy ordering. Duplicate identities and
malformed sizes must still fail; checkpoint IDs remain interleaved and unique.
Inputs are not mutated. Order is deterministic under input permutation.

The local launcher defaults `AMALGKIT_LOCAL_ACQUISITION_ORDER` to `size_bands`;
direct continuous invocation defaults to `species_fair`. Those are the only
accepted legacy values. The explicit `smallest_first` mode sets
`global_size_order=True`, creating one globally ascending acquisition lane
with task-ID tie breaks and unknown/nonpositive weights last. Size bands and
strict size order are mutually exclusive. Checkpoint verification remains
interleaved without changing the acquisition subsequence's order. The weights
are frozen total compressed estimates, not remaining wire bytes. This is a
selection/admission priority; separate pools, capacity constraints and concurrent
execution do not promise globally ordered starts or completions.
Runtime and scope receipts record the policy. Selection is
an operational scheduling priority, never a cohort filter, partial deletion,
deadline reset, global completion barrier or promise of higher wire throughput.


## Local size-overlay contract

Schema metainformant.rna.size_overlay.v1 has exactly schema, frozen_manifest_sha256, and records. Each record has task_id, accession, expected_paired, reference_index_sha256, fastq_bytes, and files. Each file has name, url, bytes, and md5. Only zero frozen weights may be supplemented; all task IDs and records must match the hash-bound manifest, mate names/cardinality must match the declared layout, and file sizes must sum to a positive bounded integer total. The file must be absolute, regular, nonsymlink, at most 1 MiB, and unchanged during its pinned-hash read. Duplicate JSON keys/identities/files, foreign URLs, invalid checksums, missing mates, changed references, and partial path/hash configuration refuse. This optional POSIX local-frontier planning input does not mutate the manifest or alter quantification provenance, checkpoints, validation, deadlines, or scientific eligibility. Reservation release uses the originally admitted effective charge. Default behavior without an overlay is unchanged.
