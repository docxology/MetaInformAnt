# Amalgkit Wrapper Documentation

The `amalgkit` module provides a Python wrapper for the Amalgkit CLI, enabling seamless integration into the MetaInformAnt pipeline. It handles command construction, execution, monitoring, and error handling.

## Key Features

- **Version Validation**: Ensures the installed `amalgkit` version meets requirements (`validate_amalgkit_version`).
- **Robust Downloading (New)**: Supports a `metainformant` backend for `getfastq` to utilize the enhanced `ena_downloader` module for reliable, widespread data retrieval with integrity checks.
- **Parallel Execution**: Splits metadata files to run `getfastq` in parallel across multiple workers.
- **Process Monitoring**: Tracks progress of long-running steps (e.g., `quant`, `getfastq`) and logs heartbeat files.
- **Automatic Installation**: Can automatically install/upgrade `amalgkit` via `uv` if missing or outdated.

## Usage

### Basic Command Execution

```python
from metainformant.rna.amalgkit import amalgkit

params = {
    "work_dir": "output/analysis",
    "threads": 8,
    "species_list": ["Homo sapiens"]
}

# Run metadata step
amalgkit.metadata(params, search_string='"Homo sapiens"[Organism] AND RNA-Seq[Strategy]')
```

### Robust Fastq Download

To use the robust `metainformant` backend for downloading (recommended):

```python
# Uses ena_downloader internally for validation and fallback
amalgkit.getfastq(params, backend="metainformant")
```

Standard `amalgkit` CLI execution:

```python
# Uses standard amalgkit getfastq command
amalgkit.getfastq(params, backend="amalgkit")
```

## Configuration

Parameters are typically managed via `AmalgkitParams` or dictionaries. Key parameters include:

- `work_dir`: Base working directory (required).
- `threads`: Number of threads/cores.
- `jobs`: Number of parallel jobs for `getfastq` (when using standard backend).
- `backend`: storage backend or execution mode (e.g., `metainformant` for downloads).
