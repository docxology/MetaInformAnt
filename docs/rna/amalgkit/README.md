# Amalgkit RNA Analysis Documentation

This directory contains comprehensive documentation for METAINFORMANT's integration with the amalgkit RNA-seq analysis toolkit.

## Overview

The amalgkit integration provides a complete transcriptomic analysis pipeline from raw SRA data to publication-ready results, with intelligent error recovery and progress monitoring.

## Documentation Files

### Core Amalgkit Integration
- **`amalgkit.md`**: Complete transcriptomic analysis pipeline documentation
- **`complete_success_summary.md`**: Production deployment success documentation

## Related Source Code

- See `src/metainformant/rna/` for implementation details
- See `src/metainformant/rna/amalgkit.py` for CLI wrapper implementation
- See `src/metainformant/rna/workflow.py` for workflow orchestration
- See `tests/test_rna_*.py` for comprehensive test coverage

## Usage Examples

The amalgkit integration supports complete RNA-seq workflows:

```python
from metainformant.rna import amalgkit, workflow

# Check amalgkit availability
ok, help_text = amalgkit.check_cli_available()

# Execute complete workflow
cfg = workflow.AmalgkitWorkflowConfig(
    work_dir="output/amalgkit/work",
    species_list=["Homo sapiens"],
    threads=8
)
results = workflow.execute_workflow(cfg)
```

## Integration

Amalgkit integration connects with:
- **NCBI SRA database** for raw data retrieval
- **Kallisto** for pseudoalignment quantification
- **R environment** for statistical analysis
- **Visualization tools** for publication figures

## Testing

Comprehensive tests ensure workflow reliability:
- CLI tool availability validation
- Workflow execution testing
- Error handling and recovery verification
- Integration testing with real data

## Contributing

When enhancing amalgkit integration:
1. Update workflow and pipeline documentation
2. Add comprehensive integration tests
3. Ensure compatibility with amalgkit ecosystem
4. Update success metrics and benchmarks

## Related Documentation

- See `docs/rna/README.md` for RNA domain overview
- See `docs/rna/workflow.md` for workflow orchestration details
- See `docs/rna/steps.md` for individual step documentation

This documentation provides complete coverage of METAINFORMANT's amalgkit integration capabilities.
