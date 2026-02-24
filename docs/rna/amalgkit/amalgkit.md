# Amalgkit Wrapper Documentation

The `amalgkit` module provides a Python wrapper for the Amalgkit CLI, enabling seamless integration into the MetaInformAnt pipeline. It handles command construction, execution, monitoring, and error handling.

## Key Features

- **Version Validation**: Ensures the installed `amalgkit` version meets requirements (`validate_amalgkit_version`).
- **Per-Sample Concurrency**: Processes samples concurrently within chunks using `ThreadPoolExecutor`. Each sample flows through `getfastq → quant → cleanup` independently.
- **Process Monitoring**: Tracks progress of long-running steps (e.g., `quant`, `getfastq`) and logs heartbeat files.
- **Automatic Installation**: Can automatically install/upgrade `amalgkit` via `uv` if missing or outdated.
- **AWS-Preferred Downloads**: Configured to prefer AWS Open Data for reliable SRA downloads (`aws: yes`, `ncbi: no`).

## Usage

### Running the Full Pipeline

```bash
# Run all 23 species sequentially (recommended)
nohup bash scripts/rna/run_all_species.sh > output/amalgkit/run_all_species_incremental.log 2>&1 &

# Run a single species
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pbarbatus.yaml --stream --chunk-size 6
```

### Monitoring Progress

```bash
# Table with ETAs
.venv/bin/python scripts/package/generate_custom_summary.py

# Check active processes
ps -fC amalgkit | grep SRR
```

### Python API

```python
from metainformant.rna.amalgkit import amalgkit

params = {
    "work_dir": "output/analysis",
    "threads": 16,
    "species_list": ["Homo sapiens"]
}

# Run metadata step
amalgkit.metadata(params, search_string='"Homo sapiens"[Organism] AND RNA-Seq[Strategy]')
```

## Configuration

Parameters are typically managed via `AmalgkitParams` or YAML config files. Key parameters include:

- `work_dir`: Base working directory (required).
- `threads`: Number of threads/cores (default: 16, dynamically divided across concurrent samples).
- `chunk_size`: Number of samples processed concurrently per chunk (default: 6).
- `aws: yes` / `ncbi: no`: Prefer AWS Open Data for SRA downloads.
- `redo: no`: Skip already-completed samples (enables resume after crash).
- `keep_fastq: no`: Delete FASTQ files immediately after quantification.
