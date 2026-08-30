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
| `--workers N` / `--threads N` | built-in defaults | Parallel workers / total threads |
| `--discovery-workers N` | 4 | Concurrent species discovery (1 = serial diagnostics) |
| `--quant-slots N` | bounded by workers/threads | Concurrent Kallisto quantifications |
| `--fastq-threads N` / `--fastq-slots N` | bounded / up to 4 | fasterq-dump fallback threads and concurrency (limits external-disk contention) |
| `--compression-threads N` | bounded from quant budget | Threads per pigz process |
| `--validation-slots N` | up to 4 | Concurrent full local FASTQ validations |
| `--max-in-flight N` | 2 × workers | Maximum submitted sample tasks |
| `--requantification-policy {preserve,version-drift,all}` | `$AMALGKIT_REQUANTIFICATION_POLICY` or `preserve` | Re-quantification trigger policy |

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

`AMALGKIT_DATA_ROOT` must be exported before invoking the launcher; its own
default points at a missing path and it exits 2 without it. It accepts any
producer flag above (e.g. `--species apis_mellifera`) and adds:

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
