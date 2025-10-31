# RNA Domain Documentation

This directory contains comprehensive documentation for METAINFORMANT's RNA analysis capabilities.

## Overview

The RNA domain provides tools for transcriptomic analysis, workflow orchestration, and integration with external bioinformatics tools like amalgkit.

## Documentation Files

### Core RNA Analysis
- **`index.md`**: RNA domain overview and module index
- **`workflow.md`**: Complete workflow planning and execution
- **`configs.md`**: Configuration management and species profiles
- **`steps.md`**: Individual workflow step documentation

### Amalgkit Integration
- **`amalgkit/`**: Amalgkit CLI wrapper documentation
  - **`amalgkit.md`**: Complete transcriptomic analysis pipeline
  - **`comprehensive_guide.md`**: Detailed workflow documentation
  - **`testing_coverage.md`**: Testing, validation, and production results
  - **`steps/`**: Individual step documentation (11 steps)
- **`workflow.md`**: Workflow orchestration with troubleshooting and optimizations

## Related Source Code

- See `src/metainformant/rna/` for implementation details
- See `tests/test_rna_*.py` for comprehensive test coverage
- See `src/metainformant/rna/README.md` for module-specific documentation

## Usage Examples

The RNA domain supports complete transcriptomic workflows:

```python
from metainformant.rna import workflow

# Execute complete RNA analysis workflow
cfg = AmalgkitWorkflowConfig(
    work_dir="output/amalgkit/work",
    species_list=["Homo sapiens"],
    threads=8
)
results = workflow.execute_workflow(cfg)
```

## Integration

RNA analysis integrates with:
- **DNA sequences** for genomic context
- **Protein annotations** for functional analysis
- **Statistical methods** for differential expression
- **Visualization tools** for expression plots

## Testing

Comprehensive tests ensure workflow reliability:
- Workflow execution validation
- Configuration parsing and validation
- External tool integration testing
- Error handling and recovery
- Production validation with 5 ant species (6,500+ samples)

See `amalgkit/testing_coverage.md` for production validation results.

## Contributing

When adding new RNA analysis functionality:
1. Update workflow and step documentation
2. Add comprehensive integration tests
3. Update configuration templates
4. Ensure compatibility with amalgkit ecosystem

This documentation provides complete coverage of METAINFORMANT's RNA analysis capabilities.
