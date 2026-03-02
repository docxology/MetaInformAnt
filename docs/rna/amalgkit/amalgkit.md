# Amalgkit Wrapper

The `amalgkit` module provides a Python wrapper for the Amalgkit CLI, enabling seamless integration into the MetaInformAnt pipeline. It handles command construction, execution, monitoring, and error handling.

## Key Features

- **Version Validation**: Ensures the installed `amalgkit` version meets requirements (`validate_amalgkit_version`).
- **Per-Sample Concurrency**: Processes samples concurrently using `ThreadPoolExecutor`. Each sample flows through `download → quant → cleanup` independently.
- **Streaming Output**: Long-running steps (`getfastq`, `quant`) stream stdout/stderr in real time via threaded `subprocess.Popen`.
- **Automatic Installation**: Can automatically install/upgrade `amalgkit` via `uv` if missing or outdated.
- **ENA-First Downloads**: Downloads use direct ENA wget (`scripts/rna/download_ena.py`) with SRA Toolkit as a fallback only.
- **Resume Support**: `redo: no` (default) skips already-quantified samples, enabling crash recovery.

## Running the Pipeline

### All 23 Species (Recommended)

```bash
nohup python3 scripts/rna/run_all_species.py \
  > output/amalgkit/run_all_species_incremental.log 2>&1 &

# Monitor
tail -f output/amalgkit/run_all_species_incremental.log
```

### Single Species

```bash
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pbarbatus.yaml

# Or with stream output and specific steps
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pbarbatus.yaml \
    --steps getfastq quant merge
```

### Monitoring Progress

```bash
# Active download workers
ps aux | grep wget | grep -v grep | wc -l

# Report quantified counts per species
python3 scripts/rna/report_completed.py

# Single-species status
python3 scripts/rna/run_workflow.py config/amalgkit/amalgkit_pbarbatus.yaml --status
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
| `threads` | 16 | Threads per kallisto quant call |
| `num_download_workers` | 16 | Parallel ENA download workers |
| `redo` | `no` | Skip already-quantified samples (enables resume) |
| `max_bp` | `4000000000` | Max bases per sample; larger samples are auto-skipped |
| `keep_fastq` | `no` | Delete FASTQs after quantification (saves disk) |

## Download Architecture

Downloads use **ENA direct wget** via `scripts/rna/download_ena.py`:

1. Queries ENA API for FTP URLs for each SRR accession
2. Downloads pre-compressed `.fastq.gz` directly (~5 GB/sample)
3. Verifies MD5 checksum and gzip integrity
4. Falls back to SRA Toolkit (`prefetch` / `fasterq-dump`) only if the sample is absent from ENA

This approach:
- Avoids SRA LITE files (metadata-only, zero reads)
- Uses ~5× less disk space than SRA cache downloads
- Supports 16+ concurrent workers without exhausting disk

## FASTQ Lifecycle

```
ENA wget download → fastq/getfastq/<SRR>/<SRR>_1.fastq.gz
     ↓
kallisto quant → work/quant/<SRR>/abundance.tsv
     ↓
FASTQ files deleted → abundance.tsv is canonical proof of work
```

## See Also

- [steps/README.md](steps/README.md) — All 11 steps
- [monitoring.md](monitoring.md) — Pipeline monitoring
- [../../ORCHESTRATION.md](../../ORCHESTRATION.md) — Orchestrator scripts
