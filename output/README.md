# Output Directory

This directory contains all outputs from tests and real runs by default, following METAINFORMANT repository conventions.

## Purpose

The `output/` directory is for **program-produced outputs only**:
- Analysis results and quantification data
- Downloaded genomic data and reference files
- Test artifacts and coverage reports
- Pipeline execution outputs

**NOT for**: Documentation, status reports, summaries, or scripts. Those belong in `docs/`, `scripts/`, or appropriate source directories.

## Repository Rules

According to METAINFORMANT's cursor rules:
- **All outputs from tests and real runs must go here by default**
- **Treat as ephemeral and reproducible**
- **No documentation files or status reports** - only program outputs
- **Scripts belong in `scripts/`, never in `output/`**
- **Use `metainformant.core.io` for file I/O operations**
- **Use `metainformant.core.paths` for path handling and containment checks**

## Current Status

### Active Outputs
- **`amalgkit/`** - RNA-seq analysis outputs for multiple ant species (P. barbatus, C. floridanus, M. pharaonis, S. invicta)
  - Genome downloads and metadata
  - Sample quantification (run_info.json files)
  - Workflow manifests and reports
  - Species-specific analysis results

### Documentation
- **`README.md`** (this file) - Directory purpose and usage guidelines
- **`AGENTS.md`** - AI assistance in output management architecture

All other files are program-produced outputs from real analysis runs and tests.

## Usage Guidelines

### For Developers
- **Never write outside `output/` unless explicitly requested**
- **Never write scripts or documentation to `output/`**
- **Use functional APIs that accept destination paths**
- **If no destination specified, default to appropriate `output/` subpath**
- **Clean outputs regularly to maintain repository size**

### For Users
- **Outputs are ephemeral** - regenerate as needed
- **Use `--output` or `--dest` flags to specify custom locations**
- **Check this directory for test results and analysis outputs**
- **Outputs may be removed in repository maintenance**

## File Management

### What Belongs in `output/`
- ✅ Analysis results (.json, .csv, .tsv, quantification data)
- ✅ Downloaded data files (genomes, SRA, FASTQ)
- ✅ Test artifacts and coverage reports
- ✅ Logs and execution status files
- ✅ Pipeline outputs and intermediate files

### What Does NOT Belong in `output/`
- ❌ **Documentation files** (*.md summaries/reports) → use `docs/`
- ❌ **Scripts** (*.sh, *.py, *.R) → use `scripts/`
- ❌ **Source code** → use `src/`
- ❌ **Configuration files** → use `config/`
- ❌ **Status reports** → belongs in appropriate permanent location

### Integration with Code
```python
# Example: Using core I/O utilities
from metainformant.core import io, paths

# Write outputs to appropriate location
output_path = paths.expand_path("output/analysis/results.json")
io.write_json(results, output_path)
```

## Related Documentation

- See [`docs/core/paths.md`](../docs/core/paths.md) for path handling utilities
- See [`docs/core/io.md`](../docs/core/io.md) for I/O operation guidelines
- See [`.cursorrules`](../.cursorrules) for complete output management policies

