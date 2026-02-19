# Agent Rules: Metabolomics Module

## Purpose

Guidelines for developing and maintaining the `metabolomics` module within METAINFORMANT. This module provides tools for metabolite identification, mass spectrometry data processing, pathway mapping, and metabolite-gene integration analysis.

## Dependencies

### Required

- `metainformant.core` (I/O, config, logging, paths)

### Optional

- `pyopenms` - Mass spectrometry data processing
- `pandas` - Tabular data manipulation
- `scipy` - Statistical analysis

## Key Submodules

| Submodule | Purpose |
|-----------|---------|
| `analysis` | Metabolite identification and quantification |
| `io` | Mass spectrometry data I/O (mzML, mzXML, CSV) |
| `pathways` | Metabolic pathway mapping and enrichment |
| `visualization` | Metabolomics-specific plots and figures |

## Patterns

### Data Loading

```python
from metainformant.core import io
from metainformant.metabolomics import io as metab_io

# Load mass spec data
data = metab_io.load_mzml("data/metabolomics/sample.mzML")
```

### Configuration

```python
from metainformant.core.config import load_mapping_from_file, apply_env_overrides

config = load_mapping_from_file("config/metabolomics/default.yaml")
config = apply_env_overrides(config, prefix="METAB")
```

## Configuration

- **Prefix**: `METAB_` (e.g., `METAB_THREADS`, `METAB_WORK_DIR`)
- **Config**: `config/metabolomics/default.yaml`

## Output Paths

- `output/metabolomics/identification/` - Metabolite ID results
- `output/metabolomics/quantification/` - Quantification tables
- `output/metabolomics/pathways/` - Pathway enrichment results
- `output/metabolomics/plots/` - Visualization outputs

## Integration

- **Multi-Omics** (`multiomics`): Metabolite-gene integration
- **Networks** (`networks`): Metabolic pathway networks
- **Visualization** (`visualization`): Shared plotting utilities
- **Quality** (`quality`): QC metrics for mass spec data

## Testing

- **NO MOCKING POLICY**: All tests use real data and implementations
- **Test files**: `tests/test_metabolomics_*.py`
- Use `tmp_path` fixture for all test outputs
- Skip gracefully when optional dependencies unavailable
