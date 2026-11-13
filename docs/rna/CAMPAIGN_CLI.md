# Campaign CLI reference: run_all_species.py and run_full_campaign.sh

Current flags for the streaming producer and its campaign launcher, verified
against `scripts/rna/run_all_species.py` and
`projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh`. The launcher is a
thin wrapper: it resolves settings, runs preflight, then runs
`run_all_species.py` (with retry passes) and finalization.

## Producer: `uv run python scripts/rna/run_all_species.py`

| Flag | Default | Meaning |
| --- | --- | --- |
| `--config-dir PATH` | repo `config/amalgkit` | Directory of `amalgkit_*.yaml` species configs |
| `--species NAME` (repeatable) | all discovered | Restrict to these species |
| `--data-root PATH` | `AMALGKIT_DATA_ROOT` or repo-local | Data root for work files and the progress DB; also sets `AMALGKIT_DATA_ROOT` |
| `--dry-run` | off | List resolved species configs without running |
| `--max-gb FLOAT` | built-in default | Max single-sample size in GB |
| `--workers N` / `--threads N` | `PIPELINE_WORKERS`/`PIPELINE_THREADS`, else the host profile | Parallel workers / total kallisto threads; leaving either out derives the whole profile from the host (below) |
| `--discovery-workers N` | 4 | Concurrent species discovery (1 = serial diagnostics) |
| `--quant-slots N` | bounded by workers/threads | Concurrent Kallisto quantifications |
| `--fastq-threads N` / `--fastq-slots N` | bounded / up to 4 | fasterq-dump fallback threads and concurrency (limits external-disk contention) |
| `--compression-threads N` | bounded from quant budget | Threads per pigz process |
| `--validation-slots N` | up to 4 | Concurrent full local FASTQ validations |
| `--max-in-flight N` | 2 × workers | Maximum submitted sample tasks |
| `--requantification-policy {preserve,version-drift,all}` | `$AMALGKIT_REQUANTIFICATION_POLICY` or `preserve` | Re-quantification trigger policy |
| `--acquisition-slots N` | host profile: 8 | Concurrent whole-sample acquisitions (0 unbounded) |
| `--ena-file-workers {1,2}` | host profile: 2 | Concurrent ENA files of one paired run |
| `--scheduling-order {size_bands,smallest_first,species_order}` | `size_bands` | Fresh-acquisition order across species (`AMALGKIT_PIPELINE_SCHEDULING_ORDER`) |
| `--fast-fail-ena` | off | No ENA retries, NCBI prefetch first, 64 KiB/s floor over 20 s (only for unset variables) |
| `--print-profile` | off | Print the resolved resource profile as JSON and exit (no lock, no progress DB) |

### Host resource profile

Without explicit `--workers` and `--threads`, `build_host_resource_profile`
derives the budgets from the host and the measured c19 run (8 cores 4P+4E,
16 GiB, USB data volume: 32.4 samples/h average, 46.5/h best 6 h): 8
acquisitions with 2 ENA files each (16 ENA file connections, 2 API
connections), quant slots = min(CPUs, 8, memory guard) with one kallisto
thread each (the memory guard budgets 256 MiB + 4 x the largest index per
run against half of RAM), 2 full FASTQ validations, 1 NCBI recovery slot,
workers = acquisitions + quant slots, `max_in_flight` = workers, and a raw
reservation of min(256 GiB, 25% of the data root's free space)
(`AMALGKIT_PIPELINE_RAW_RESERVATION_GB` overrides). Every value can be
overridden by its flag or environment variable; the resolved profile is logged
at startup.

### Scheduling and retries

Each species' tasks join the queue when its own discovery finishes, so work
starts before the slowest species is discovered; a species whose discovery
fails is skipped for that run (exit 75). The queue is ordered reusable local
inputs, resumable partials, then fresh acquisitions in global size bands with
species round-robin (ties by species and run). Every finished attempt is
classified in the progress DB (`attempts`, `last_error_class`,
`next_attempt_at`): environment failures retry at the next start, transient
failures back off 1 h, 2 h, 4 h ... (cap 24 h, `AMALGKIT_RETRY_BACKOFF_*`) and
are held after 6 attempts (`AMALGKIT_MAX_SAMPLE_ATTEMPTS`), layout mismatches
are quarantined, and kallisto exit-1 failures are held. Held rows are not
retried until `python -m metainformant.rna.engine.progress_db release-holds
--db <root>/pipeline_progress.db` (under the writer lock); `retry-status`
reports them read-only. `AMALGKIT_RETRY_WAIT_MAX_SECONDS` (default 0) lets a
producer whose only work left is in backoff wait for it instead of exiting 75.

### Requantification policy

- `preserve` (default): keep any complete, checksum-verified quantification;
  only quarantined/failed work is redone. This is what a resume sweep runs —
  quarantined samples re-enter as re-quantification candidates.
- `version-drift`: additionally rebuild quantifications whose provenance audit
  reports compatible version drift.
- `all`: rebuild every quantification regardless of audit status.

The policy can also be set through `AMALGKIT_REQUANTIFICATION_POLICY`; the
orchestrator rejects unknown values at startup.

## Launcher: `projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh`

The Hymenoptera project's one-command lane. `AMALGKIT_DATA_ROOT` defaults to
`/Volumes/external_drive/Data/amalgkit`; a run exits 69 when that directory
does not exist (`--print-config` still works). It runs the producer above in
frozen-cohort mode (`--cohort <data-root>/cohort/hymenoptera27-20261005` with the project's
committed `COHORT.json` as `--cohort-definition` and its pinned sha256), with
the host-derived profile and no shell clamps. It owns `--data-root`,
`--config-dir` and the `--cohort*` options (set them through
`AMALGKIT_CONFIG_DIR`, `AMALGKIT_COHORT_DIR`, `AMALGKIT_COHORT_DEFINITION` and
`AMALGKIT_COHORT_DEFINITION_SHA256`); any other producer flag passes through,
and `--species` also limits finalization.

| Mode or setting | Default | Meaning |
| --- | --- | --- |
| `--print-config` | — | Resolved settings, command lines and the producer's `--print-profile` JSON; no lock, no database, no writes |
| `--dry-run` | — | Preflight plus the read-only plan (cohort and configuration digests, references, progress per species); exit 0/73/69/2 |
| `--status` | — | Writer lock holder, launcher status file and per-species completion (read-only) |
| `AMALGKIT_PIPELINE_RETRY_PASSES` | 1 | Extra producer passes after exit 75 |
| `AMALGKIT_PIPELINE_RETRY_DELAY_SECONDS` | 300 | Pause before a retry pass |
| `AMALGKIT_PIPELINE_WORKERS` / `_THREADS` and the other `AMALGKIT_PIPELINE_*` profile variables | host profile | Read by the producer's profile |
| `AMALGKIT_MIN_EXTERNAL_FREE_GB` / `AMALGKIT_MIN_SYSTEM_FREE_GB` | 8 / 4 | Preflight free-space floors |
| `AMALGKIT_FAST_FAIL_ENA` | `no` | `yes` adds `--fast-fail-ena` |
| `AMALGKIT_NCBI_PREFETCH_FIRST` | `yes` | Resumable verified prefetch before fasterq fallback |
| `AMALGKIT_TEMP_ROOT` / `AMALGKIT_CACHE_ROOT` | `<data-root>/.tmp` / `.cache` | Producer `TMPDIR` / `XDG_CACHE_HOME`; uv's cache stays at its default unless `AMALGKIT_UV_CACHE_DIR` is set |
| `AMALGKIT_CAFFEINATE` | `yes` | Re-execute under `caffeinate -dims` unless already under `caffeinate` |

One launch runs preflight, then the producer; exit 75
leads to bounded retry passes, and only a drained producer (exit 0) is followed
by merge → wsfilter → finalize → sanity under
`python -m metainformant.rna.engine.campaign_lock run --kind finalization`. The
launcher itself takes no lock: it refuses to start (73) while another writer
holds `.campaign.flock` or another launcher is live in
`<data-root>/results/full_campaign.status` (verified by PID and process start
time), and keeps that status file only for humans and that check. Child
processes start with default SIGINT/SIGQUIT handling; the launcher forwards
INT/TERM/HUP after a grace period and exits 130/143/129.

The producer itself (`scripts/rna/run_all_species.py`, `process_species.py`)
holds the writer `flock` on `<data-root>/.campaign.flock` for its whole life;
the kernel releases it on any death, so there is no stale lock to remove.
Quantification imports and progress-DB migrations take the same lock, and
`python -m metainformant.rna.engine.campaign_lock status --data-root <root>`
reports the holder (exit 0 quiescent, 73 busy). Producer exit codes:

| Code | Meaning |
| --- | --- |
| 0 | Drained: every scheduled task is quantified, quarantined, held by the retry bookkeeping or permanently excluded |
| 75 | Retryable work remains (failed, pending or backed-off tasks, or a species whose discovery failed); rerun the same command (the summary's `next_attempt_due_at` says when backed-off work is due) |
| 73 | Another writer holds the data root |
| 69 | Environment failure: preflight, cohort refusal, or the environment circuit breaker (`AMALGKIT_ENV_FAILURE_TRIP`, default 5 consecutive EPERM/EACCES/EROFS/ENXIO/ENOSPC failures) |
| 70 | Internal error |

The launcher starts finalization only after exit 0, so a run that still ends
with retryable work (75) after its retry passes leaves finalization to a later
run. Each
producer run logs INFO to `<data-root>/logs/run_all_species_<utc>_<pid>.log`.

## Disk-throttle behavior

While the campaign runs, the orchestrator re-checks both volume floors. When
the external data volume drops below its gate it logs `[Disk Throttle]`
warnings, stops starting new downloads, and checks again roughly every five
minutes; it logs `Space available again` and resumes the queue automatically
when space is reclaimed. No restart is needed — reclaim regenerable caches and
wait.

## Checkpoint and resume

The orchestrator is resume-safe: SQLite progress DB (`PRAGMA quick_check` must
pass), `.part` transfer files, and cached SRA inputs are all reused on restart.
`scripts/checkpoint_full_campaign.sh` stops the launcher recorded in
`results/full_campaign.status` (or else the live writer-lock holder, such as a
bare `run_all_species.py`) and its process tree at a resume-safe checkpoint. The
resume sweep re-audits every recorded quantification under the `preserve`
policy and re-queues what it cannot certify current; that is the audit working,
not data loss. All hashed provenance artifacts are content-deterministic: the
reference alias manifest (`work/<species>/reference/reference_aliases.json`)
contains no wall-clock fields, so an unchanged reference produces a
byte-identical manifest across restarts and prior quantification survives. A
`quarantined` spike after a restart therefore indicates a real contract change
(reference rebuilt, configuration edited, schema bump) — or a pre-2026-08-30
manifest sidecar from before the determinism fix, which re-quantifies once and
never recurs.
