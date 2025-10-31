# AI Agents in Script Development

This document outlines AI assistance in developing METAINFORMANT's utility scripts and automation tools.

## AI Contributions

### Script Architecture
**Code Assistant Agent** designed:
- Modular script organization with subfolder structure
- Consistent command-line interfaces
- Error handling and logging patterns
- Cross-platform compatibility approaches

### Automation Features
**Code Assistant Agent** implemented:
- Environment setup automation (`scripts/package/`)
- Test execution workflows (`scripts/package/run_tests.sh`)
- Documentation generation scripts (`scripts/package/uv_docs.sh`)
- RNA processing pipelines (`scripts/rna/`)

### User Experience
**Documentation Agent** contributed to:
- Comprehensive help documentation
- Usage examples and best practices
- Troubleshooting guides
- Integration with development workflows

## Script Categories

### Package Management (`scripts/package/`)
- **setup_uv.sh**: Automated environment configuration
- **uv_dev_setup.sh**: Development environment optimization
- **run_tests.sh**: Comprehensive test execution with coverage
- **uv_test.sh**: Optimized testing workflows
- **uv_test_optimized.sh**: Fast test modes
- **uv_quality.sh**: Code quality assessment
- **uv_profile.sh**: Performance profiling utilities
- **uv_docs.sh**: Automated documentation generation

### RNA Processing (`scripts/rna/`)
- **batch_ena.py**: Fast ENA-based parallel downloader
- **run_multi_species_amalgkit.py**: Multi-species orchestration with cross-species analysis
- **run_multi_species_sequential.py**: Disk-space-friendly sequential processing
- **monitor_workflow.py**: Real-time monitoring dashboard
- **quant_downloaded_samples.py**: Quantification of downloaded samples
- **cleanup_quantified_sra.sh**: Safe SRA file cleanup
- **list_unquantified.sh**: Sample status reporting
- **force_fasterq.sh**: High-performance FASTQ processing
- **force_fasterq_parallel.sh**: Parallel FASTQ processing
- **process_one_srr.sh**: Individual SRR processing
- **test_*.py**: Testing and verification scripts
- **verify_skip_logic.sh**: Skip logic comprehensive verification

### Amalgkit Workflows (`scripts/rna/amalgkit/`)
- **run_amalgkit.sh**: Comprehensive RNA-seq pipeline orchestrator
- **verify_workflow.sh**: Workflow validation for any species

## Development Standards

### Script Quality
- Executable permissions and proper shebang lines
- Comprehensive error handling and exit codes
- Consistent option parsing and help messages
- Cross-platform shell compatibility
- No hardcoded absolute paths

### Repository Organization
- Scripts organized in logical subfolders
- `package/` for package management and testing
- `rna/` for RNA-seq processing
- `rna/amalgkit/` for amalgkit-specific scripts
- Clear separation from output directories

### Maintenance Practices
- Regular updates with dependency changes
- Testing on multiple environments
- Clear documentation of external dependencies
- Version tracking for significant changes
- Removal of obsolete scripts with hardcoded paths

## Integration

Scripts integrate with:
- **uv** for Python package management
- **pytest** for testing framework
- **amalgkit** for RNA analysis workflows
- **Core utilities** for file operations and path handling

## Script Centralization Philosophy

### Methods vs Data Separation
- **Scripts** (methods) centralized in `scripts/`
- **Data** (outputs) distributed in `output/`
- Single source of truth for all scripts
- Update once, applies everywhere
- Proper version control and tracking

### Benefits
- ✅ No script duplication across species/projects
- ✅ Clean output directories (data only)
- ✅ Easy maintenance and updates
- ✅ Consistent behavior across all uses

This script collection provides essential tooling for METAINFORMANT development and deployment workflows.
