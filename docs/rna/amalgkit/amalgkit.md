# Amalgkit Wrapper

The `amalgkit` module provides a Python wrapper for the Amalgkit CLI, enabling seamless integration into the MetaInformAnt pipeline. It handles command construction, execution, monitoring, and error handling.

## Key Features

- **Version Validation**: `validate_amalgkit_version()` and
  `verify_installed_amalgkit()` check the installed build against the
  0.16.108 pin (upstream commit 1d1566638abe), and the project launcher
  verifies the installed release before a real campaign.
- **Bounded per-sample concurrency**: Processes samples through a bounded
  executor. Task workers cover validation/acquisition while a separate
  semaphore caps Kallisto quantification slots, so I/O overlap does not
  increase the quantification memory budget. Each sample flows through
  `download → quant → provenance-gated cleanup` independently.
- **Streaming Output**: Long-running steps (`getfastq`, `quant`) stream stdout/stderr in real time via threaded `subprocess.Popen`.
- **ENA-First Downloads**: Downloads use the implemented `ENADownloader`
  class with SRA Toolkit fallback in the streaming orchestrator.
- **Resume Support**: `redo: no` (default) skips already-quantified samples, enabling crash recovery.
- **Reference-index cache**: `AMALGKIT_LOCAL_INDEX_CACHE_DIR` optionally copies
  only `*.idx` files to a fast local volume with a size/mtime manifest; raw
  reads and quantification outputs remain under `AMALGKIT_DATA_ROOT`.
- **Local quant-output scratch**: `AMALGKIT_LOCAL_QUANT_SCRATCH_DIR` optionally
  isolates each Kallisto output on a fast local volume and promotes a validated
  result back into the external `work/quant/<run>/` directory atomically; raw
  reads remain external.
- **Local SRA-extraction scratch**: `AMALGKIT_LOCAL_FASTERQ_SCRATCH_DIR`
  optionally places cached-SRA or bare-accession `fasterq-dump` output,
  temporary files, and compression work on a fast host volume. A cached SRA
  on another filesystem is copied there ephemerally so repeated source reads
  are local. Validated gzip files are copied through an atomic
  destination-side rename. A shared free-space reservation falls back to
  external scratch when the local reserve would be crossed; the authoritative
  source SRA/cache and canonical FASTQs remain under `AMALGKIT_DATA_ROOT`.
- **Cached-SRA integrity witnesses**: `vdb-validate` checks a local archive
  before extraction. An atomic size/mtime-bound witness prevents unchanged
  corrupt caches from consuming repeated `fasterq-dump` attempts, while a
  replaced or completed source is validated again.
- **Acquisition free-space floors**: `AMALGKIT_MIN_EXTERNAL_FREE_GB` (default
  8 GiB) protects the data-root volume and `AMALGKIT_MIN_SYSTEM_FREE_GB`
  (default 4 GiB) protects the host system volume. Set the external floor to
  zero only for isolated temporary roots in tests or controlled scratch runs.
- **Explicit timeouts**: transfer, `fasterq-dump`, compression, and quantification
  timeouts are independently controlled by the four
  `AMALGKIT_PIPELINE_*_TIMEOUT_SECONDS` variables and are recorded as failures,
  never as quantified evidence.

## Running the Pipeline

### Full Campaign

Production runs use the project launcher, which runs the producer in
frozen-cohort mode under the campaign writer lock:

```bash
export AMALGKIT_DATA_ROOT=/Volumes/external_drive/Data/amalgkit
bash projects/hymenoptera_amalgkit/scripts/run_full_campaign.sh
```

Each producer run logs to
`$AMALGKIT_DATA_ROOT/logs/run_all_species_<utc>_<pid>.log`.

### Single Species (diagnostic)

```bash
env -u VIRTUAL_ENV uv run --frozen --no-sync python scripts/rna/run_all_species.py \
    --config-dir projects/hymenoptera_amalgkit/config/amalgkit \
    --data-root "$AMALGKIT_DATA_ROOT" \
    --cohort "$AMALGKIT_DATA_ROOT/cohort/hymenoptera27-20261005" \
    --species pogonomyrmex_barbatus
```

### Monitoring Progress

```bash
# Active acquisition, quantification, and orchestration processes
pgrep -af 'curl .*fastq|fasterq-dump|kallisto quant|run_all_species.py' || true

# Current database-backed status
uv run python scripts/rna/check_pipeline_status.py \
    --data-root "$AMALGKIT_DATA_ROOT" --species pogonomyrmex_barbatus

# Idempotent downstream checkpoint chain after acquisition stops
bash projects/hymenoptera_amalgkit/scripts/run_all_finalization.sh \
    --data-root "$AMALGKIT_DATA_ROOT" --dry-run
```

## Python API

```python
from metainformant.rna.amalgkit import amalgkit

params = {
    "work_dir": "output/amalgkit/pbarbatus/work",
    "threads": 16,
    "species_list": ["Pogonomyrmex barbatus"]
}

# Run metadata step
amalgkit.metadata(params, search_string='"Pogonomyrmex barbatus"[Organism] AND RNA-Seq[Strategy]')
```

## Configuration

Parameters are managed via YAML config files. Key parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `work_dir` | — | Base working directory (required) |
| `--threads` | `PIPELINE_THREADS`, else one per quant slot of the host profile | Total quantification thread budget |
| `--workers` | `PIPELINE_WORKERS`, else acquisition + quant slots of the host profile | Active sample tasks, including validation, acquisition, and quantification |
| `--fastq-threads` | Derived | Threads for an NCBI `fasterq-dump` fallback |
| `--compression-threads` | Derived | Threads per fallback `pigz` process |
| `--validation-slots` | host profile: 2 | Concurrent full local FASTQ gzip scans; limits mounted-volume contention |
| `--max-in-flight` | host profile: workers | Submitted task window |
| `redo` | `no` | Skip already-quantified samples (enables resume) |
| `max_bp` | Amalgkit default `999999999999999`; the repository configs set `50000000` | Amalgkit `getfastq` per-sample base-pair download target; the ENA-first producer does not read it and never skips a sample for it (its size budget is `--max-gb`) |
| `clean_fastq` | `no` in the project lane | Retain reads until MetaInformAnt writes verified contract provenance; then the guarded reclamation utility deletes them |

## Download Architecture

Downloads use **direct ENA retrieval** via
`metainformant.rna.retrieval.ena_downloader.ENADownloader`:

1. Queries ENA API for FTP URLs for each SRR accession
2. Downloads pre-compressed `.fastq.gz` directly (~5 GB/sample)
3. Verifies gzip integrity (`verify_gzip_integrity()`); `download_run()` does not compare MD5 checksums
4. Falls back to SRA Toolkit (`prefetch` / `fasterq-dump`) only when the ENA path cannot complete; a local `<run>.sra` is extracted before any network request

This approach:

- Avoids SRA LITE files (metadata-only, zero reads)
- Uses ~5× less disk space than SRA cache downloads
- Runs 8 concurrent acquisitions (2 ENA files each) in the default host profile, bounded by the free-space floors

## FASTQ Lifecycle

```
existing validated FASTQ/SRA extraction → ENA HTTPS download → work/getfastq/<SRR>/<SRR>_1.fastq.gz
     ↓
kallisto quant → work/quant/<SRR>/abundance.tsv + quantification provenance sidecar
     ↓
guarded reclamation deletes the FASTQ files only after the sidecar validates
```

The launcher runs the species of the frozen cohort with their configurations
from `projects/hymenoptera_amalgkit/config/amalgkit/`; without `--cohort`, the
producer discovers runnable species from `--config-dir` and excludes template,
test, and cross-species orchestration configurations. The mounted
data root is selected with `AMALGKIT_DATA_ROOT`; an absent or partial data root
must be reported as incomplete rather than inferred to be a completed
historical run.

## See Also

- [steps/README.md](steps/README.md) — Current Amalgkit stage reference
- [monitoring.md](monitoring.md) — Pipeline monitoring
- [../ORCHESTRATION.md](../ORCHESTRATION.md) — Orchestrator scripts
