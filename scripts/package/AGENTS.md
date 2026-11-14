# Agent Directives: scripts/package

## Role
Package management, build, and testing infrastructure scripts.

## Key Scripts

### Environment Setup
- `setup.sh` - Complete environment setup
- `uv_dev_setup.sh` - UV-based development setup
- `verify_amalgkit_runtime.py` - Read-only check of the pinned Amalgkit runtime contract

### Testing
- `test.sh` - Main test runner (supports modes: ultra-fast, coverage, parallel)

### Building
- `build.sh` - Build package
- `validate_build.sh` - Validate build artifacts
- `verify.sh` - Verify installation: UV setup (`--mode setup`) and test dependencies (`--mode deps`)

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

# Recreate an existing environment only after confirming no process depends on it
bash scripts/package/setup.sh --recreate-venv --python 3.14

# Testing
bash scripts/package/test.sh --mode ultra-fast

# Quality
bash scripts/package/uv_quality.sh
```
