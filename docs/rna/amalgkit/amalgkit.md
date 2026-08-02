# Amalgkit Wrapper

The `amalgkit` module provides a Python wrapper for the Amalgkit CLI, enabling seamless integration into the MetaInformAnt pipeline. It handles command construction, execution, monitoring, and error handling.

## Key Features

- **Version Validation**: The project launcher verifies the exact installed
  Amalgkit release before a real campaign; do not rely on an unavailable
  helper import as a substitute for that executable check.
- **Bounded per-sample concurrency**: Processes samples through a bounded
  executor. Task workers cover validation/acquisition while a separate
  semaphore caps Kallisto quantification slots, so I/O overlap does not
  increase the quantification memory budget. Each sample flows through
  `download → quant → provenance-gated cleanup` independently.
- **Streaming Output**: Long-running steps (`getfastq`, `quant`) stream stdout/stderr in real time via threaded `subprocess.Popen`.
- **Automatic Installation**: Can automatically install/upgrade `amalgkit` via `uv` if missing or outdated.
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

### All Configured Species (Recommended)

```bash
export AMALGKIT_DATA_ROOT=/path/to/amalgkit-data
uv run python scripts/rna/run_all_species.py \
  --data-root "$AMALGKIT_DATA_ROOT" \
  --workers 4 --threads 8 \
  --fastq-threads 2 --compression-threads 1 --validation-slots 4 --max-in-flight 8 \
  > "$AMALGKIT_DATA_ROOT/results/run_all_species.log" 2>&1 &

# Monitor
tail -f "$AMALGKIT_DATA_ROOT/results/run_all_species.log"
```

### Single Species

```bash
uv run python scripts/rna/run_all_species.py \
    --config-dir config/amalgkit \
    --data-root "$AMALGKIT_DATA_ROOT" \
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
| `--threads` | 8 in the general launcher; 10 in the measured Hymenoptera profile | Total quantification thread budget |
| `--workers` | 4 in the general launcher; 16 in the measured Hymenoptera profile | Active sample tasks, including validation, acquisition, and quantification |
| `--fastq-threads` | Derived | Threads for an NCBI `fasterq-dump` fallback |
| `--compression-threads` | Derived | Threads per fallback `pigz` process |
| `--validation-slots` | 4 | Concurrent full local FASTQ gzip scans; limits mounted-volume contention |
| `--max-in-flight` | `2 × workers` | Submitted task window |
| `redo` | `no` | Skip already-quantified samples (enables resume) |
| `max_bp` | `4000000000` | Max bases per sample; larger samples are auto-skipped |
| `clean_fastq` | `no` in the project lane | Retain reads until MetaInformAnt writes current provenance; then the guarded reclamation utility deletes them |

## Download Architecture

Downloads use **direct ENA retrieval** via
`metainformant.rna.retrieval.ena_downloader.ENADownloader`:

1. Queries ENA API for FTP URLs for each SRR accession
2. Downloads pre-compressed `.fastq.gz` directly (~5 GB/sample)
3. Verifies MD5 checksum and gzip integrity
4. Falls back to SRA Toolkit (`prefetch` / `fasterq-dump`) only if the sample is absent from ENA; a local `<run>.sra` is extracted before any network request

This approach:
- Avoids SRA LITE files (metadata-only, zero reads)
- Uses ~5× less disk space than SRA cache downloads
- Supports 16+ concurrent workers without exhausting disk

## FASTQ Lifecycle

```
existing validated FASTQ/SRA extraction → ENA HTTPS download → work/getfastq/<SRR>/<SRR>_1.fastq.gz
     ↓
kallisto quant → work/quant/<SRR>/abundance.tsv
     ↓
FASTQ files deleted → abundance.tsv is canonical proof of work
```

The launcher discovers runnable species from `config/amalgkit/` and excludes
template, test, and cross-species orchestration configurations. The mounted
data root is selected with `AMALGKIT_DATA_ROOT`; an absent or partial data root
must be reported as incomplete rather than inferred to be a completed
historical run.

## See Also

- [steps/README.md](steps/README.md) — Current Amalgkit stage reference
- [monitoring.md](monitoring.md) — Pipeline monitoring
- [../../ORCHESTRATION.md](../../ORCHESTRATION.md) — Orchestrator scripts
