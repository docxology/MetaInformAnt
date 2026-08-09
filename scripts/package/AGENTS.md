# Agent Directives: scripts/package

## Role
Package management, build, and testing infrastructure scripts.

## Key Scripts

### Environment Setup
- `setup.sh` - Complete environment setup
- `uv_dev_setup.sh` - UV-based development setup
- `verify_uv_setup.sh` - Verify UV installation
- `verify_amalgkit_runtime.py` - Read-only check of the pinned Amalgkit runtime contract

### Testing
- `test.sh` - Main test runner (supports modes: ultra-fast, coverage, parallel)
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
- `generate_cursor_skills.py` - Emit `.cursor/skills/<name>/SKILL.md` for every `AGENTS.md` in the repo; `--check` verifies complete wrapper parity, canonical link targets, front matter, and orphan absence
- `release.sh` - Release workflow
- `fix_tmp_space.sh` - Fix temp space issues
- `_common.sh` - Shared shell utilities
- `build_utils.sh` - Build utilities

## Usage
```bash
# Setup
bash scripts/package/setup.sh

# Recreate an existing environment only after confirming no process depends on it
bash scripts/package/setup.sh --recreate-venv --python 3.12

# Testing
bash scripts/package/test.sh --mode ultra-fast

# Quality
bash scripts/package/uv_quality.sh
```
