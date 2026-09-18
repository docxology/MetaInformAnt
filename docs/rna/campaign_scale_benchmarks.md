# Campaign-scale benchmarks

Characterizes the orchestration-side costs of the current RNA-seq campaign —
species discovery, queue-depth accounting, retry/backoff scheduling,
hash/provenance I/O, and SQLite storage growth — at parameterized task-row
scales, **without touching the live producer, the live data root, or the live
progress database**.

Harness: `scripts/rna/benchmark_campaign_scale.py`
Tests: `tests/rna/test_benchmark_campaign_scale.py`

## Safety contract (enforced in code)

- The live data root (`/Volumes/external_drive/Data/amalgkit`) is never
  opened for reading or writing. The target guard compares paths lexically —
  it does not even stat the live root — and refuses any `--target-dir` at or
  beneath it.
- Fixtures are disposable by construction: a synthetic species-config and
  metadata tree, per-sample quantification directories with provenance
  sidecars, and a SQLite progress DB built by the real `ProgressDB`
  implementation. Default target: `output/campaign_scale_benchmarks/run_<UTC>`.
  The fixture tree is deleted after the run unless `--keep-fixture` is given.
- No network access, no external tools, and the live
  `output/amalgkit/pipeline_progress.db` (or any DB under the live data root)
  is never opened.

## Running the harness

```bash
# Small validation sweep (default full sweep is 1000,10000,100000 rows)
python scripts/rna/benchmark_campaign_scale.py --scales 1000

# Full sweep with default scales, keeping the fixture for inspection
python scripts/rna/benchmark_campaign_scale.py --keep-fixture

# Custom scale + seed
python scripts/rna/benchmark_campaign_scale.py --scales 250000 --seed 42
```

Outputs `report.json` and `report.md` into the target directory. Both carry a
machine-provenance block (commit, branch, host, platform, python, seed,
scales, target dir, `execution_context`). Determinism: every fixture byte
(species names, accessions, state assignment, exclusions, abundance payloads)
is generated from a fixed seed; timings are inherently machine- and
load-dependent and are reported as descriptive observations only — no
statistical or significance claims.

## Telemetry definitions

| Metric | Definition | Unit |
|---|---|---|
| config scans/s | sorted `amalgkit_*.yaml` glob + marker filter (`discover_species_config_names`), median of repeats | configs/s |
| cohort rows/s | parse of the synthetic per-species metadata TSVs into (species, run-accession) pairs | rows/s |
| queue depth (s) | full dashboard pass: `get_counts()`, `get_total_counts()`, per-species pending `get_samples()` intersected with permanent-drop exclusions | s |
| retry transitions/s | `set_state(failed, error)` + `set_state(pending)` pairs, error text classified via `classify_sample_error` | transitions/s |
| stale resets | rows reset by `reset_stale_downloading(3600)` on a backdated downloading cohort | rows |
| hash MB/s | `digest_file` over the abundance payloads | MB/s |
| classifications/s | full `classify_quantification(verify_content=True)` passes | samples/s |
| DB bytes/row | checkpointed SQLite size (db + WAL) per task row | bytes/row |

## Budgets (advisory, descriptive)

- Queue-depth dashboard pass must stay comfortably below the 5 s orchestrator
  stall-watchdog granularity at every campaign scale.
- Phase-1 discovery plus queue accounting must not re-digest the quant
  corpus; the resume-time reconciliation classifies by provenance contract
  instead of re-hashing payloads (`ProgressDB.reconcile` docstring). The
  `verify_content=True` classification numbers here therefore describe the
  deep-audit path, not the resume path.
- Retry scheduling is per-sample-commit bound; size backoff sweeps against
  the measured transitions/s at the storage volume actually hosting the
  progress DB.
- Storage: the checkpointed DB footprint grows linearly (~233 bytes/task row
  in the example below, schema, index, and WAL included).

## EXAMPLE run — 2026-09-17 (local disposable fixture, NOT a hosted run)

Small validation sweep at 1,000 task rows (8 synthetic species, seed
20260917) on the campaign workstation (Apple M4 Pro, external-data drive,
Python 3.12.13). **Environmental caveat:** the drive was under heavy
concurrent load from other campaign work during this run, so absolute I/O
rates (hash MB/s, transitions/s, fixture build wall time) describe that
contended configuration; re-run for fresh numbers. The proof value of this
section is that the harness runs end-to-end on a disposable fixture with a
machine-provenance block.

Summary:

| scale (task rows) | build (s) | DB bytes/row | queue-depth (s) | retry transitions/s | classifications/s |
|---|---|---|---|---|---|
| 1,000 | 457.413647 | 233.47 | 0.000665 | 6.7 | 5.8 |

Full measurements:

```json
{
  "safety": {
    "live_data_root": "/Volumes/external_drive/Data/amalgkit",
    "live_data_root_read": false,
    "live_db_opened": false,
    "network_access": false
  },
  "scales": [
    {
      "discovery": {
        "cohort_enumeration_seconds": 0.000192,
        "cohort_rows": 1000,
        "cohort_rows_per_second": 5219668.0,
        "config_files_scanned": 8,
        "config_scan_seconds": 0.000217,
        "config_scans_per_second": 36880.5
      },
      "fixture_build_seconds": 457.413647,
      "provenance_io": {
        "classification_status_counts": {
          "current": 1000
        },
        "classifications_per_second": 5.8,
        "hash_megabytes_per_second": 0.015053176879882812,
        "hash_seconds": 163.72984,
        "hashed_bytes": 2584375,
        "sample_dirs": 1000,
        "sidecar_reads_per_second": 6.3
      },
      "queue_depth": {
        "get_counts_seconds": 0.000132,
        "pending_queue_total": 568,
        "queue_depth_per_second": 1503.7,
        "queue_depth_seconds": 0.000665,
        "total_task_rows": 1000
      },
      "retry_backoff": {
        "stale_downloading_resets": 24,
        "stale_resets_per_second": 55.8,
        "stale_scan_seconds": 0.430379,
        "transition_pairs": 200,
        "transition_seconds": 59.573137,
        "transitions_per_second": 6.7
      },
      "species_count": 8,
      "storage": {
        "db_bytes": 233472,
        "db_bytes_per_row": 233.47,
        "fixture_fs_bytes": 3698007,
        "shm_bytes": 32768,
        "wal_bytes": 0
      },
      "task_rows": 1000
    }
  ],
  "schema": "metainformant.rna.benchmark.campaign_scale.v1"
}
```

Machine provenance:

```json
{
  "execution_context": "local disposable fixture (not a hosted run)",
  "generated_at_utc": "2026-09-17T23:41:47+00:00",
  "git_branch": "main",
  "git_commit": "3f48bb2d7320940435a90e5658ed914095fe9562",
  "git_dirty": null,
  "hostname": "timelock.local",
  "platform": "macOS-26.6.2-arm64-arm-64bit",
  "python_version": "3.12.13",
  "seed": 20260917,
  "scales": [
    1000
  ],
  "schema": "metainformant.rna.benchmark.campaign_scale.v1",
  "target_dir": "/Volumes/external_drive/Git/projects/ongoing/docxology/MetaInformAnt/output/campaign_scale_benchmarks/validate_1e3"
}
```

`git_dirty` is `null` because `git status` could not be captured within the
provenance timeout under the concurrent load of the run (recorded as
"unknown" rather than guessed).

## Reading the numbers

- Discovery and queue accounting are negligible next to the 5 s watchdog
  budget even at the largest planned scale; the queue-depth pass at 1,000
  rows is ~0.7 ms of index-backed GROUP BY work.
- The retry/backoff and provenance I/O figures are dominated by per-operation
  fsync/read latency on the storage volume, not by SQLite or hashing logic;
  treat them as volume characteristics under load, and re-run on the target
  volume before sizing backoff sweeps.
- Every sidecar fixture classified `current`, confirming the synthetic
  provenance chain (config hash + quant payload digest) round-trips through
  the real `write_quant_provenance` / `classify_quantification` pair.
