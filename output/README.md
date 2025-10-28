# Output Directory

This directory contains all outputs from tests and real runs by default, following METAINFORMANT repository conventions.

## Repository Rules

According to METAINFORMANT's cursor rules:
- **All outputs from tests and real runs must go here by default**
- **Treat as ephemeral and reproducible**
- **Use `metainformant.core.io` for file I/O operations**
- **Use `metainformant.core.paths` for path handling and containment checks**

## Directory Structure

### Analysis Outputs

#### `amalgkit/` Subdirectory
Contains RNA-seq analysis outputs from AMALGKIT workflows:
- Complete workflow consolidation results (`CONSOLIDATION_COMPLETE.md`)
- Species-specific analysis outputs (e.g., `pbarbatus/`)
- Sample count summaries and metadata
- Expression quantification results

#### `ontology/` Subdirectory
Contains ontology-related output files:
- Gene Ontology (GO) term summaries (`go_summary.json`)
- Ontology enrichment results
- Semantic similarity calculations
- Cross-ontology mapping results

#### `processing/` Subdirectory
Contains general data processing outputs:
- Processing status and results (`processing_results.json`)
- Intermediate analysis files
- Pipeline execution logs

### Test Outputs

#### Test-specific subdirectories
- `test_*/` - Outputs from domain-specific test suites
- `coverage_html/` - HTML coverage reports
- `coverage.xml` - XML coverage data for CI/CD integration

#### Test artifacts
- `test_amalgkit/` - AMALGKIT integration test outputs
- `test_docs/` - Documentation generation test outputs
- `test_integration/` - Cross-module integration test outputs
- `test_invalid_config.yaml` - Configuration validation test outputs
- `test_performance/` - Performance benchmark results
- `test_robustness/` - Error handling and robustness test outputs
- `test_steps/` - Workflow step testing outputs
- `test_workflow/` - Complete workflow testing outputs

## Usage Guidelines

### For Developers
- **Never write outside `output/` unless explicitly requested**
- **Use functional APIs that accept destination paths**
- **If no destination specified, default to appropriate `output/` subpath**
- **Clean outputs regularly to maintain repository size**

### For Users
- **Outputs are ephemeral** - regenerate as needed
- **Use `--output` or `--dest` flags to specify custom locations**
- **Check this directory for test results and analysis outputs**
- **Outputs may be removed in repository maintenance**

## File Management

### Cleanup Recommendations
- Remove large output files after analysis
- Use `.gitignore` patterns for generated outputs
- Archive important results outside the repository
- Regenerate outputs for verification

### Integration with Code
```python
# Example: Using core I/O utilities
from metainformant.core import io, paths

# Write outputs to appropriate location
output_path = paths.expand_path("output/ontology/go_summary.json")
io.write_json(results, output_path)
```

## Related Documentation

- See `src/metainformant/core/paths.md` for path handling utilities
- See `src/metainformant/core/io.md` for I/O operation guidelines
- See main project README for output management policies
