# Amalgkit Module

Python wrapper for the amalgkit RNA-seq quantification tool.

## 📦 Components

| File | Purpose |
|------|---------|
| [`amalgkit.py`](amalgkit.py) | Main CLI wrapper and command builder |
| [`genome_prep.py`](genome_prep.py) | Reference genome download and kallisto indexing |
| [`metadata_filter.py`](metadata_filter.py) | Sample metadata filtering and validation |
| [`metadata_utils.py`](metadata_utils.py) | Metadata parsing utilities |
| [`__main__.py`](__main__.py) | CLI entry point |

## 🔑 Key Functions

### amalgkit.py

- `run_amalgkit_step()` - Execute single workflow step
- `build_amalgkit_command()` - CLI argument construction
- `AmalgkitParams` - Step parameters dataclass

### genome_prep.py

- `download_genome()` - Fetch reference from NCBI
- `build_kallisto_index()` - Create transcript index

### metadata_filter.py

- `filter_rna_seq_samples()` - Filter to valid RNA-seq
- `select_qualified_samples()` - Apply quality criteria

## 🚀 Usage

```python
from metainformant.rna.amalgkit import (
    run_amalgkit_step,
    download_genome,
    build_kallisto_index
)

# Download and index genome
download_genome("GCF_000187915.1", dest_dir)
build_kallisto_index(transcripts_file, index_file)

# Run quantification
run_amalgkit_step("quant", params)
```

## 🔗 Related

- [engine/](../engine/) - Workflow orchestration
- [config/amalgkit/](../../../../config/amalgkit/) - Configurations
