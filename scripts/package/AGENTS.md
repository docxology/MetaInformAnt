# Agent Directives: scripts/package

## Role
Package management, build, and testing infrastructure scripts.

## Key Scripts

### Environment Setup
- `setup.sh` - Complete environment setup
- `uv_dev_setup.sh` - UV-based development setup
- `setup_uv.sh` - Initialize UV environment
- `verify_uv_setup.sh` - Verify UV installation

### Testing
- `test.sh` - Main test runner (supports modes: ultra-fast, coverage, parallel)
- `run_tests.sh` - Simple test execution
- `uv_test.sh` - Run tests via UV
- `uv_test_optimized.sh` - Optimized test execution
- `uv_test_setup.sh` - Test environment setup
- `verify_test_deps.sh` - Verify test dependencies

### Building
- `build.sh` - Build package
- `validate_build.sh` - Validate build artifacts
- `verify.sh` - Verify installation

### Code Quality
- `uv_quality.sh` - Run format, lint, typecheck
- `uv_profile.sh` - Profiling utilities
- `uv_docs.sh` - Build documentation

### Utilities
- `release.sh` - Release workflow
- `fix_tmp_space.sh` - Fix temp space issues
- `_common.sh` - Shared shell utilities
- `build_utils.sh` - Build utilities

## Usage
```bash
# Setup
bash scripts/package/setup.sh

# Testing
bash scripts/package/test.sh --mode ultra-fast

# Quality
bash scripts/package/uv_quality.sh
```
