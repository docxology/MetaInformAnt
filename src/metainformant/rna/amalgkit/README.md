# Amalgkit Module

Python wrapper for the amalgkit RNA-seq quantification tool.

## Components

| File | Purpose |
|------|---------|
| [`amalgkit.py`](amalgkit.py) | Public facade that re-exports the CLI wrapper |
| [`_amalgkit_impl.py`](_amalgkit_impl.py) | CLI wrapper implementation: command builder, step functions, and the 0.16.108 pin checks (`validate_amalgkit_version`, `verify_installed_amalgkit`) |
| [`commands.py`](commands.py) | Canonical command names and CLI option tables for the pinned Amalgkit build |
| [`sra_environment.py`](sra_environment.py) | Campaign-local subprocess environments for SRA and FASTQ tooling |
| [`genome_prep.py`](genome_prep.py) | Reference genome download and kallisto indexing |
| [`index_prep.py`](index_prep.py) | Kallisto index preparation |
| [`metadata_filter.py`](metadata_filter.py) | Sample metadata filtering and validation |
| [`metadata_utils.py`](metadata_utils.py) | Metadata parsing utilities |
| [`tissue_normalizer.py`](tissue_normalizer.py) | Tissue label normalization via mapping files |
| [`__main__.py`](__main__.py) | CLI entry point |

## Key Functions

### amalgkit.py

- `run_amalgkit()` — Execute any amalgkit CLI step
- `build_amalgkit_command()` — CLI argument construction
- `AmalgkitParams` — Step parameters dataclass
- Step functions: `metadata()`, `select()`, `getfastq()`, `quant()`, `merge()`, `wsfilter()`, `finalize()`, `sanity()`

Monitored download, quantification, merge, filtering, and finalization steps accept the
standard `timeout` keyword. It bounds the monitored child process, writes a
timeout heartbeat with an error message, terminates the child, and returns
exit code `124`. Use `monitor_timeout` when the subprocess timeout and the
progress-monitor timeout need different values.

### genome_prep.py

- `download_genome()` - Fetch reference from NCBI
- `build_kallisto_index()` - Create transcript index

### metadata_filter.py

- `filter_rna_seq_samples()` - Filter to valid RNA-seq
- `select_qualified_samples()` - Apply quality criteria

## Usage

```python
from metainformant.rna.amalgkit import (
    run_amalgkit,
    build_amalgkit_command,
    AmalgkitParams,
)

# Build and inspect a command
cmd = build_amalgkit_command("quant", {"out_dir": "work", "threads": 4})

# Run a step
result = run_amalgkit("quant", AmalgkitParams(work_dir="work", threads=4))
```

## Related

- [engine/](../engine/) - Workflow orchestration
- [Hymenoptera configurations](../../../../projects/hymenoptera_amalgkit/config/amalgkit/) - Canonical species YAMLs and rules
