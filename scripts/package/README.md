# Package Management Scripts

Comprehensive environment setup, testing, and package management utilities for METAINFORMANT.

## Directory Structure

```
scripts/package/
├── setup.sh                      # Unified environment setup ⭐
├── test.sh                       # Unified test runner ⭐
├── verify.sh                     # Unified verification script ⭐
├── _common.sh                    # Shared utilities library
├── uv_docs.sh                    # Documentation generation
├── uv_quality.sh                 # Code quality checking
├── uv_profile.sh                 # Performance profiling
├── run_tests.sh                  # → test.sh (deprecated)
├── uv_test.sh                    # → test.sh (deprecated)
├── uv_test_optimized.sh          # → test.sh --mode fast (deprecated)
├── setup_uv.sh                   # → setup.sh (deprecated)
├── verify_test_deps.sh           # → verify.sh --mode deps (deprecated)
├── verify_uv_setup.sh            # → verify.sh --mode setup (deprecated)
├── uv_dev_setup.sh               # Development environment setup
├── uv_test_setup.sh              # Test environment configuration
├── fix_tmp_space.sh              # Temporary space cleanup utility
└── README.md                     # This file
```

## Core Scripts (RECOMMENDED)

### Environment Setup (`setup.sh`)

Unified environment setup script using `uv` package manager with comprehensive dependency management.

**Features:**
- Virtual environment creation with filesystem detection
- Complete dependency installation (dev, scientific, optional)
- External CLI tool setup (SRA tools, kallisto, etc.)
- Development environment configuration
- Automatic filesystem compatibility (FAT/external drive support)

**Usage:**
```bash
# Basic setup (recommended)
bash scripts/package/setup.sh

# Development environment
bash scripts/package/setup.sh --dev

# Full installation with all dependencies
bash scripts/package/setup.sh --with-all --with-deps

# Skip amalgkit installation
bash scripts/package/setup.sh --skip-amalgkit

# Set NCBI email for API access
bash scripts/package/setup.sh --ncbi-email your.email@example.com
```

**What it installs:**
- **Core**: uv package manager, virtual environment
- **Development**: pre-commit, documentation tools, testing
- **Scientific**: NumPy, SciPy, scikit-learn, pandas, etc.
- **RNA-seq**: amalgkit, kallisto, SRA tools (optional)
- **Optional**: Database clients, network libraries, scraping tools

### Test Runner (`test.sh`)

Unified test execution with multiple modes, coverage reporting, and parallel execution.

**Features:**
- Multiple test execution modes (ultra-fast to comprehensive)
- Coverage analysis and reporting
- Parallel test execution
- Network and external tool test support
- Real implementation testing (no mocks)

**Usage:**
```bash
# Fast essential tests (recommended default)
bash scripts/package/test.sh --mode fast

# Full coverage analysis
bash scripts/package/test.sh --mode coverage --coverage

# Parallel execution with coverage
bash scripts/package/test.sh --mode parallel --coverage

# Include network-dependent tests
bash scripts/package/test.sh --mode network

# Ultra-fast core functionality only
bash scripts/package/test.sh --mode ultra-fast
```

**Test Modes:**
- `ultra-fast`: Core functionality only (~5s)
- `fast`: Essential tests (~15s) **[DEFAULT]**
- `coverage`: Full coverage analysis (~1-2min)
- `coverage-fast`: Fast coverage on core modules (~20s)
- `parallel`: Parallel execution with coverage
- `network`: Include network-dependent tests (REAL API calls)
- `external`: Include external CLI tool tests
- `integration`: End-to-end integration tests
- `smoke`: Quick smoke tests for deployment

### Verification (`verify.sh`)

Unified environment and dependency verification with automatic fixes.

**Features:**
- UV setup verification
- Test dependency checking
- Environment configuration validation
- Automatic issue resolution
- Comprehensive status reporting

**Usage:**
```bash
# Verify everything
bash scripts/package/verify.sh

# Verify specific component
bash scripts/package/verify.sh --mode setup
bash scripts/package/verify.sh --mode deps

# Auto-fix issues
bash scripts/package/verify.sh --fix
```

## Development and Quality Tools

### Documentation (`uv_docs.sh`)

Automated documentation generation and serving.

**Usage:**
```bash
# Build documentation
bash scripts/package/uv_docs.sh build

# Serve documentation locally
bash scripts/package/uv_docs.sh serve

# Clean documentation build
bash scripts/package/uv_docs.sh clean
```

### Code Quality (`uv_quality.sh`)

Comprehensive code quality checking with linting, formatting, and type checking.

**Usage:**
```bash
# Run all quality checks
bash scripts/package/uv_quality.sh all

# Individual checks
bash scripts/package/uv_quality.sh format     # Black formatting
bash scripts/package/uv_quality.sh lint       # Ruff linting
bash scripts/package/uv_quality.sh typecheck  # mypy type checking
```

### Performance Profiling (`uv_profile.sh`)

CPU and memory profiling utilities for performance optimization.

**Usage:**
```bash
# CPU profiling
bash scripts/package/uv_profile.sh cpu

# Memory profiling
bash scripts/package/uv_profile.sh memory

# Benchmarking
bash scripts/package/uv_profile.sh benchmark
```

### Temporary Space Cleanup (`fix_tmp_space.sh`)

Quick utility to clean up `/tmp` space by removing pytest temp files and uv cache.

**Usage:**
```bash
bash scripts/package/fix_tmp_space.sh
```

**What it does:**
- Checks `/tmp` disk usage before and after cleanup
- Removes pytest temp directories (`/tmp/pytest-of-*`)
- Removes uv cache (`/tmp/uv-cache`)
- Shows large temp files for manual review

**Note**: Use when `/tmp` space is low and affecting workflow execution.

## Legacy Scripts (DEPRECATED)

### Migration Guide
```bash
# OLD → NEW
run_tests.sh              → test.sh --mode fast
uv_test.sh               → test.sh --mode coverage
uv_test_optimized.sh     → test.sh --mode ultra-fast
setup_uv.sh             → setup.sh
verify_test_deps.sh     → verify.sh --mode deps
verify_uv_setup.sh      → verify.sh --mode setup
```

## Shared Utilities (`_common.sh`)

Common functions used across all package management scripts:
- Filesystem detection and UV cache configuration
- Virtual environment management
- Color output and status reporting
- Error handling and logging
- Dependency checking utilities

## Integration

These scripts integrate with:
- **uv**: Modern Python package management
- **pytest**: Comprehensive testing framework
- **pre-commit**: Code quality enforcement
- **GitHub Actions**: CI/CD pipeline
- **Development workflow**: Standardized setup and testing

## Dependencies

- **Core**: bash, curl, git
- **Python**: uv package manager, Python 3.11+
- **Optional**: External bioinformatics tools (SRA Toolkit, kallisto, etc.)
- **Development**: pre-commit, mypy, ruff, black

## Key Features

✅ **Unified Interface**: Single entry points for setup, testing, verification
✅ **Filesystem Agnostic**: Automatic detection of FAT/external drives
✅ **Comprehensive Coverage**: All aspects of development environment
✅ **Real Testing**: No mocks - actual external API and tool testing
✅ **Automatic Fixes**: Self-healing verification and setup
✅ **Development Ready**: Pre-configured development environment
✅ **CI/CD Integration**: GitHub Actions compatible
✅ **Backward Compatible**: Legacy script support with migration messages

## Related Documentation

- [UV Setup Guide](../../docs/UV_SETUP.md)
- [Testing Documentation](../../docs/testing.md)
- [Development Setup](../../docs/setup.md)
- [METAINFORMANT CLI](../../docs/cli.md)
