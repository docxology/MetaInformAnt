# Ontology Output Directory

This directory contains outputs from ontology processing and analysis, following METAINFORMANT repository conventions.

## Repository Rules

According to METAINFORMANT's cursor rules:
- **All outputs from tests and real runs must go here by default**
- **Treat as ephemeral and reproducible**
- **Use `metainformant.core.io` for file I/O operations**

## Directory Contents

### Gene Ontology Outputs
- **`go_summary.json`**: Summary of GO term analysis results
- **`enrichment_results.json`**: Gene set enrichment analysis outputs
- **`semantic_similarity.json`**: Ontology-based similarity calculations

## Usage Guidelines

### For Developers
- **Never write outside `output/` unless explicitly requested**
- **Use functional APIs that accept destination paths**
- **If no destination specified, default to appropriate `output/` subpath**
- **Clean outputs regularly to maintain repository size**

### For Users
- **Outputs are ephemeral** - regenerate as needed
- **Use `--output` or `--dest` flags to specify custom locations**
- **Check this directory for ontology analysis results**
- **Outputs may be removed in repository maintenance**

## File Management

### Cleanup Recommendations
- Remove large output files after analysis
- Use `.gitignore` patterns for generated outputs
- Archive important results outside the repository
- Regenerate outputs for verification

## Related Documentation

- See `src/metainformant/ontology/README.md` for ontology module documentation
- See `src/metainformant/core/paths.md` for path handling utilities
- See `src/metainformant/core/io.md` for I/O operation guidelines
