# Scripts Directory

This directory contains utility scripts for METAINFORMANT development, testing, and deployment.

## Available Scripts

### Environment Setup
- **`setup_uv.sh`**: One-shot uv environment setup script
  - Installs uv package manager
  - Creates virtual environment
  - Installs all dependencies from `pyproject.toml`
  - Use: `bash scripts/setup_uv.sh`

- **`uv_dev_setup.sh`**: Development environment setup
  - Extended setup for development workflows
  - Includes additional development tools

### Testing Scripts
- **`run_tests.sh`**: Comprehensive test execution
  - Multiple execution modes (fast, network, pattern matching)
  - Coverage reporting and performance tracking
  - Use: `./scripts/run_tests.sh --help` for options

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

### Workflow Scripts
- **`run_amalgkit.sh`**: Amalgkit workflow execution
  - RNA-seq analysis pipeline runner

- **`process_one_srr.sh`**: Single SRR processing
  - Individual SRA accession processing

- **`force_fasterq.sh`**: FASTQ processing utilities
  - High-performance sequence file processing

- **`force_fasterq_parallel.sh`**: Parallel FASTQ processing
  - Multi-threaded sequence processing

### Documentation
- **`uv_docs.sh`**: Documentation generation
  - Build and update documentation

## Usage Guidelines

### Development Workflow
```bash
# Setup environment
bash scripts/setup_uv.sh

# Run tests
./scripts/run_tests.sh

# Check code quality
bash scripts/uv_quality.sh

# Generate documentation
bash scripts/uv_docs.sh
```

### Production Usage
```bash
# Run RNA analysis workflow
bash scripts/run_amalgkit.sh --config config/amalgkit_pbarbatus.yaml

# Process large datasets
bash scripts/force_fasterq_parallel.sh --input data/ --output results/
```

## Script Standards

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
