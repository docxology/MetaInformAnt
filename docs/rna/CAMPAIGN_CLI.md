# Campaign CLI reference: run_all_species.py and run_full_campaign.sh

Current flags for the streaming producer and its campaign launcher, verified
against `scripts/rna/run_all_species.py` and
`projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh`. The launcher is a
thin wrapper: it resolves resources and env defaults, then execs
`run_all_species.py` with the same flags.

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

`AMALGKIT_DATA_ROOT` defaults to `/Volumes/external_drive/Data/amalgkit` and
the launcher exits 2 when that directory does not exist; export it explicitly
for any other location. It accepts any producer flag above (e.g.
`--species apis_mellifera`) and adds:

| Setting | Default | Meaning |
| --- | --- | --- |
| `AMALGKIT_PIPELINE_PROFILE` | `local-throughput` | `local-throughput` (workers 16–24, threads 8, quant-slots 4, fastq-slots 2, validation-slots 4) or `compat` (workers ≤ 8, fastq-slots 1) |
| `--results-dir` / `--print-config` | — | Results directory; print resolved config and exit |
| `AMALGKIT_MIN_SYSTEM_FREE_GB` | 32 | Startup floor on system-volume free space; the run is refused below it |
| `AMALGKIT_MIN_EXTERNAL_FREE_GB` | 32 | Same floor for the external data volume |
| `AMALGKIT_RECLAIM_RAW_AFTER_QUANT` | `yes` | Reclaim validated raw FASTQ/SRA after quant |
| `AMALGKIT_NCBI_PREFETCH_FIRST` | `yes` | Resumable verified prefetch before fasterq fallback |
| Per-stage timeouts (`AMALGKIT_PIPELINE_*_TIMEOUT_SECONDS`) | download 600 s low-speed window; fasterq 7200; pigz 1800; quant 7200 | curl/fasterq/pigz/kallisto bounds |

The launcher takes the campaign lock at
`<data-root>/results/.full_campaign.lock`, refuses to run while another
producer holds it, and installs a single-campaign downstream runner
(merge → wsfilter → finalize → sanity) that starts only after the producer
exits cleanly.

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

The launcher starts finalization only after exit 0, so a run that ends with
retryable sample failures (75) leaves finalization to a later run. Each
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
`scripts/checkpoint_full_campaign.sh` stops the launcher-managed producer tree
at a resume-safe checkpoint (it reads the lock-dir PID; a bare
`run_all_species.py` process has no lock and must be stopped manually). The
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
