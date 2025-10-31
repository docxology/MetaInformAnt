# Scripts Directory

This directory contains utility scripts for METAINFORMANT development, testing, and deployment.

## Directory Structure

```
scripts/
├── package/        # Package management and testing scripts
├── rna/           # RNA-seq processing scripts
│   └── amalgkit/  # Amalgkit-specific workflow scripts
├── README.md      # This file
└── AGENTS.md      # AI agent contribution documentation
```

## Package Scripts (`scripts/package/`)

### Environment Setup
- **`setup_uv.sh`**: One-shot uv environment setup script
  - Installs uv package manager
  - Creates virtual environment
  - Installs all dependencies from `pyproject.toml`
  - Use: `bash scripts/package/setup_uv.sh`

- **`uv_dev_setup.sh`**: Development environment setup
  - Extended setup for development workflows
  - Includes additional development tools

### Testing Scripts
- **`run_tests.sh`**: Comprehensive test execution
  - Multiple execution modes (fast, network, pattern matching)
  - Coverage reporting and performance tracking
  - Use: `scripts/package/run_tests.sh --help` for options

- **`uv_test.sh`**: uv-based testing
  - Optimized test execution with uv
  - Parallel test execution

- **`uv_test_optimized.sh`**: Optimized testing workflow
  - Fastest test execution configuration

### Profiling and Quality
- **`uv_profile.sh`**: Performance profiling script
  - Memory and CPU profiling utilities

- **`uv_quality.sh`**: Code quality checking
  - Linting, formatting, and quality metrics

### Documentation
- **`uv_docs.sh`**: Documentation generation
  - Build and update documentation

## RNA Processing Scripts (`scripts/rna/`)

### Amalgkit Workflows (`scripts/rna/amalgkit/`)
- **`run_amalgkit.sh`**: Comprehensive RNA-seq pipeline orchestrator
- **`verify_workflow.sh`**: Workflow validation for any species
- See `scripts/rna/amalgkit/README.md` for details

### Multi-Species Orchestration
- **`run_multi_species_amalgkit.py`**: Full workflow for all species with cross-species analysis
- **`run_multi_species_sequential.py`**: Disk-space-friendly sequential processing

### Download and Processing
- **`batch_ena.py`**: Fast ENA-based parallel downloader
- **`force_fasterq.sh`**: SRA FASTQ processing utilities
- **`force_fasterq_parallel.sh`**: Parallel FASTQ processing
- **`process_one_srr.sh`**: Single SRR processing

### Quantification
- **`quant_downloaded_samples.py`**: Quantify already-downloaded samples
- **`cleanup_quantified_sra.sh`**: Safe deletion of quantified SRA files
- **`list_unquantified.sh`**: Report samples needing quantification

### Monitoring and Testing
- **`monitor_workflow.py`**: Real-time workflow monitoring dashboard
- **`monitor_amalgkit_progress.sh`**: Multi-species progress monitor
- **`test_pbarbatus_workflow.py`**: P. barbatus workflow testing
- **`test_single_species.py`**: Single species workflow testing
- **`test_skip_logic.py`**: Skip logic verification
- **`verify_skip_logic.sh`**: Comprehensive skip logic verification

See `scripts/rna/README.md` for detailed RNA script documentation.

## Usage Guidelines

### Development Workflow
```bash
# Setup environment
bash scripts/package/setup_uv.sh

# Run tests
bash scripts/package/run_tests.sh

# Check code quality
bash scripts/package/uv_quality.sh

# Generate documentation
bash scripts/package/uv_docs.sh
```

### Production Usage
```bash
# Run RNA analysis workflow
bash scripts/rna/amalgkit/run_amalgkit.sh --config config/amalgkit/amalgkit_pbarbatus.yaml

# Multi-species workflow
python3 scripts/rna/run_multi_species_amalgkit.py

# Monitor workflow progress
python3 scripts/rna/monitor_workflow.py
```

## Script Standards

### Repository Organization
- **Scripts always belong in `scripts/`** - never in `output/`, `data/`, or other directories
- Scripts may read from `data/` and `config/` directories
- Scripts must write outputs to `output/` directory by default
- Use subdirectories for domain-specific scripts (e.g., `scripts/rna/`)

### Requirements
- All scripts must be executable (`chmod +x`)
- Include help documentation and usage examples
- Handle errors gracefully with informative messages
- Support common flags (`--help`, `--verbose`, `--dry-run`)
- Compatible with both bash and zsh

### Maintenance
- Update scripts when dependencies change
- Test scripts on multiple systems
- Document any external dependencies
- Version scripts with significant changes

## Contributing New Scripts

### Template
```bash
#!/bin/bash
# Script Name: Brief description
#
# Usage: script_name [options]
#
# Options:
#   -h, --help     Show this help message
#   -v, --verbose  Verbose output
#   -d, --dry-run  Show what would be done

set -euo pipefail

# Script implementation
```

### Review Process
1. Test script functionality thoroughly
2. Add comprehensive help documentation
3. Ensure cross-platform compatibility
4. Update this README with new script
5. **Place script in `scripts/` directory** (not `output/` or other locations)

## Integration

Scripts integrate with:
- **uv** for Python package management
- **pytest** for testing framework
- **amalgkit** for RNA analysis
- **Core utilities** for file operations

## Troubleshooting

### Common Issues
- **Permission denied**: Run `chmod +x script_name.sh`
- **Command not found**: Ensure dependencies are installed
- **Environment errors**: Source virtual environment first

### Getting Help
- Use `--help` flag on any script
- Check script source code for detailed comments
- See related documentation for context

These scripts provide essential tooling for METAINFORMANT development and usage workflows.
