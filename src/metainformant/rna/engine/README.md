# RNA engine

The engine implements the current Amalgkit workflow as a bounded producer,
SQLite progress store, typed workflow planner, and hash-bound evidence layer.
The producer refuses to start into a broken environment: `run_all()` executes
the campaign preflight before any discovery or scheduling work.

## Components

The opt-in local `acquisition_frontier.py` scheduler separates acquisition,
recovery and quantification pools. `AMALGKIT_LOCAL_FRONTIER_SIZE=48` bounds
planned plus admitted identities; `AMALGKIT_LOCAL_RAW_RESERVATION_BYTES` bounds
the sum of frozen compressed-size reservations held through the checkpoint.
These are planning weights, not a measurement or hard ceiling on physical
bytes: preexisting partials, size-estimation error and extraction/staging space
need independent filesystem capacity checks. Missing/oversized reservations
remain explicit unresolved tasks. With the option unset or zero, the existing
whole-sample scheduler remains active.

Every frontier transition is fsynced to the run's `frontier_events.jsonl` with
manifest, selection and optional runtime source hashes. The manifest/selection,
current validated files and progress database drive restart; the transition
ledger is observational and never grants permission to skip validation. STOP
prevents planned items from starting while acquired/recovering/ready items
drain through the existing result journal and checkpoint. Four acquisition
slots do not imply 48 active downloads, and four quant slots share the total
thread budget.

The local frontier can check existing current-contract quantifications while
acquisition slots are busy. One synchronous reuse check runs at a time; it
cannot download, extract or quantify. A hit retains the normal raw reservation,
result journal and checkpoint sequence, and increments reused work only. A miss
is probed once and remains ordinary bounded acquisition work. STOP and admission
checks still apply. `reused_without_acquisition` counts this path; it does not
increase newly quantified samples or measured network capacity.
All checkpoint payload paths, including abundance, run information, H5 files
and restored input witnesses, must be regular files inside their owned roots
before content reads. A symlink or nonregular payload is a clean probe miss;
the downstream checkpoint remains unchanged. Ready futures are drained between
serial probe attempts, including misses.

| Module | Purpose |
|---|---|
| `streaming_orchestrator.py` | ENA-first metadata, acquisition, integration, and quantification |
| `preflight.py` | Mandatory start-of-run environment preflight (data-root write probe, amalgkit CLI resolution) |
| `progress_db.py` | Concurrent-safe SQLite sample state, exclusions, and resume queries |
| `campaign_status.py`, `campaign_status_io.py`, `campaign_status_cli.py` | Read-only frozen-inventory cloud/local status, complete stage marginals and per-sample transfer candidates; receipt/file presence is distinct from checksum validation |
| `campaign_observation.py`, `campaign_observation_report.py`, `campaign_observation_cli.py` | Explicit manifest/phase journal capture, restore-bound new-output counts, isolated tables/figures/working supplement; operational evidence only |
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
| `fastq_compression.py` | Strict lossless pigz level and thread command settings; default level 6, optional level 1 for temporary FASTQ scratch |

## Key interfaces

- `capture_campaign_observation()` captures stable control records and append-only
  journal prefixes from explicit paths. Startup totals remain producer-reported;
  new increments require matching restore/provenance contracts and abundance hashes.
  `write_observation_bundle()` refuses live-root destinations and existing outputs.
  The thin entry point is `scripts/rna/observe_acquisition.py`; optional transfer
  observations are explicit inputs, never selected by modification time.
  [Working observation protocol](../../../../projects/hymenoptera_amalgkit/docs/manuscript/working_operational_update.md).

- `acquisition_*` implements generic frozen local/AWS manifest acquisition,
  disjoint allocations, stage-resource controls, replay journals, live regional
  price quotation and explicit cost/time scenarios. Use the public facade
  `metainformant.rna.amalgkit.acquisition` or `scripts/rna/acquisition.py`.
  [Generic acquisition guide](../../../../docs/rna/GENERIC_ACQUISITION.md).
  Library `run_manifest()` callers can checkpoint each successful or reused sample
  through `on_quantified(species, accession)` and stop new admissions through
  `should_stop()`. Callback failures propagate; stopped summaries retain the
  unsubmitted task count.
  Opt-in `AMALGKIT_PIPELINE_DEDICATED_RECOVERY=yes` uses `acquisition_lanes.py`
  to isolate SRA recovery while retaining admission, storage and deadline
  ownership through each sample's checkpoint. `AMALGKIT_NCBI_PREFETCH_FIRST=yes`
  resumes verified archives inside the sample's raw tree before SSD extraction.
  `AMALGKIT_LOCAL_ENA_SCRATCH_DIR` enables bounded fresh SSD acquisition and
  serialized verified promotion to the selected data root. The opt-in
  `AMALGKIT_PIPELINE_DOWNLOAD_FIRST_PASS_SECONDS` reports `deferred_download`
  as pending without immediate recovery; retry phases must reset it to zero.
  `AMALGKIT_PIPELINE_VALIDATE_FASTQ_RECORDS=yes` applies to ENA, recovery,
  local SRA extraction and reuse. SSD acquisition also requires record
  validation for retained inputs when the explicit flag is disabled. Promotion
  validates the prospective retained and extracted mate union before modifying
  inputs. Conflicting extracted reads
  remain in a unique `<accession>-conflict-*/output/` directory beside the
  extraction workspace; recovery refuses invalid complete canonical inputs.
  [Transfer and preservation contract](../retrieval/README.md).

- `acquisition_prerequisites.py` resolves native quantifier requirements from
  frozen batch metadata before raw acquisition. `acquisition_references.py`
  binds metadata target aliases to the existing reference checksum. Long-read
  oarfish tasks require an amended frozen reference/tool/output contract; they
  remain unresolved under the Kallisto-bound acquisition envelope. See the
  [method boundary](../../../../docs/rna/HYMENOPTERA_METHODS.md#quantifier-and-reference-prerequisites).
  Explicit 16S titles together with GENOMIC/PCR or METAGENOMIC/PCR metadata remain
  unresolved for host-transcriptome quantification, without changing frozen
  membership. `acquisition_order.py` interleaves unsealed current outputs with
  species-fair, size-ordered acquisition; duplicate task identities are refused.
  Reference preflight shares full index digests only within that invocation
  when the resolved path, device, inode, size, modification time and change time
  still match. Each distinct native index copy is hashed independently.
  An unchanged invalid-SRA validation witness rejects extraction before waiting
  for a fallback slot, so ENA recovery can proceed while another SRA is extracting.

- `acquisition_scheduling.py` models new batches from explicit transfer/extraction/
  quantification rates and records remaining-deadline refusals. Unknown sizes/rates
  do not authorize new work. See the
  [deadline planning contract](../../../../docs/rna/GENERIC_ACQUISITION.md#deadline-planning-for-new-admissions).

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
- `durable_quant.py` preserves validated per-sample outputs in local and S3 stores;
  `aws_completion.py` resumes byte-bounded missing-task partitions within a gross budget.
  See [the durable quantification guide](../../../../docs/rna/DURABLE_QUANT.md).

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

### Restored quantification inputs

Portable restoration accepts `verified_config_path` alongside the frozen reference
index hash. It preserves original provenance and creates checksum-bound local
configuration/reference witnesses for strict downstream classification. See
[durable quantification](../../../../docs/rna/DURABLE_QUANT.md).
# Whole-sample acquisition admission

`AMALGKIT_PIPELINE_ACQUISITION_SLOTS` optionally bounds simultaneous
`download_fastq` calls separately from task workers, file connections and
quantification slots. Zero retains the existing unlimited sample admission.
`AMALGKIT_PIPELINE_ACQUISITION_QUEUE_TIMEOUT_SECONDS` is a finite positive
queue allowance (default 7,200 seconds). A queue expiry defers the task without
touching its inputs. The first-pass clock starts after admission; no file
request, retry, validation, or mate resets that clock. The lease is released
before quantification and on every error path.

The prepared local profile selects four acquisition slots for eight native
file connections and two file workers per paired sample, with twelve task
workers and four quantifiers. This reduces competition between samples for
file leases; it does not increase the file cap or prove maximum bandwidth.
A task may spend up to the queue allowance plus its acquisition allowance
before quantification. STOP still drains admitted work without forced termination.

## Cross-species size bands

The local launcher selects `AMALGKIT_LOCAL_ACQUISITION_ORDER=size_bands`.
Known positive frozen compressed-size weights enter global bands ending at
100 MiB, 500 MiB, 1 GiB, 2 GiB, 5 GiB and 20 GiB; larger known sizes follow,
then unknown/nonpositive sizes. Species rotate inside each band, with smaller
samples first within a species. Checkpoint verification stays interleaved.
The general API and direct continuous runner retain `species_fair` by default.
Unknown order settings fail before runtime controls or locks are written.

Bands govern admission order, not completion barriers, eligibility or remote
bandwidth. Every selected identity appears once; held prerequisites stay held.
Frozen size weights can contain metadata fallbacks and do not subtract retained
local inputs or resumable parts. Only a new phase can use a new ordering policy;
never rewrite an active partition or discard its download progress.

