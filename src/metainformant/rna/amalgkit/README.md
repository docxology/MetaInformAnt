# Amalgkit Module

Python wrapper for the amalgkit RNA-seq quantification tool.

## 📦 Components

| File | Purpose |
|------|---------|
| [`amalgkit.py`](amalgkit.py) | Main CLI wrapper and command builder |
| [`genome_prep.py`](genome_prep.py) | Reference genome download and kallisto indexing |
| [`index_prep.py`](index_prep.py) | Kallisto index preparation |
| [`metadata_filter.py`](metadata_filter.py) | Sample metadata filtering and validation |
| [`metadata_utils.py`](metadata_utils.py) | Metadata parsing utilities |
| [`tissue_normalizer.py`](tissue_normalizer.py) | Tissue label normalization via mapping files |
| [`__main__.py`](__main__.py) | CLI entry point |

## 🔑 Key Functions

### amalgkit.py

- `run_amalgkit()` — Execute any amalgkit CLI step
- `build_amalgkit_command()` — CLI argument construction
- `AmalgkitParams` — Step parameters dataclass
- Step functions: `metadata()`, `select()`, `getfastq()`, `quant()`, `merge()`, `curate()`, `sanity()`

### genome_prep.py

- `download_genome()` - Fetch reference from NCBI
- `build_kallisto_index()` - Create transcript index

### metadata_filter.py

- `filter_rna_seq_samples()` - Filter to valid RNA-seq
- `select_qualified_samples()` - Apply quality criteria

## 🚀 Usage

```python
from metainformant.rna.amalgkit import (
    run_amalgkit,
    build_amalgkit_command,
    AmalgkitParams,
)

# Build and inspect a command
cmd = build_amalgkit_command("quant", {"out_dir": "work", "threads": 4})

# Run a step
result = run_amalgkit("quant", AmalgkitParams(out_dir="work", threads=4))
```

## 🔗 Related

- [engine/](../engine/) - Workflow orchestration
- [config/amalgkit/](../../../../config/amalgkit/) - Configurations
