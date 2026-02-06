# GitHub Actions Workflows

This document describes the CI/CD workflows for the MetaInformAnt project.

## Overview

| Workflow | Trigger | Purpose |
|----------|---------|---------|
| [Build Package](#build-package) | Push, PR, Manual | Build and validate package |
| [Test Suite](#test-suite) | Push, PR, Manual | Run tests across platforms |
| [Release](#release) | Tag, Release, Manual | Publish to PyPI |

---

## Build Package

**File:** `build.yml`

### Triggers

- Push to `main` or `develop`
- Pull requests to `main` or `develop`
- Manual dispatch with build type option

### Jobs

| Job | Description |
|-----|-------------|
| `build` | Build wheel/sdist for Python 3.11 & 3.12 |
| `build-validation` | Install and validate package in clean env |
| `summary` | Report overall status |

### Artifacts

- `metainformant-{python-version}-{sha}` - Built distribution files

---

## Test Suite

**File:** `test.yml`

### Triggers

- Push to `main` or `develop`
- Pull requests to `main` or `develop`
- Manual dispatch with test type option

### Jobs

| Job | Description |
|-----|-------------|
| `test` | Matrix tests (ubuntu/macos × py3.11/3.12 × fast/network/external/all) |
| `test-infrastructure` | Validate test infrastructure itself |
| `quality-checks` | Black, isort, flake8, bandit, safety |
| `test-examples` | Run example scripts |
| `fat-filesystem-test` | Test on FAT-like filesystem |
| `summary` | Report overall status |

### Test Types

- `fast` - Quick core tests (~15s)
- `network` - Tests requiring internet
- `external` - Tests requiring external CLI tools (muscle, seqkit)
- `all` - Full suite with coverage

---

## Release

**File:** `release.yml`

### Triggers

- GitHub release published
- Push to `v*.*.*` tags
- Manual dispatch (test/production)

### Jobs

| Job | Description |
|-----|-------------|
| `release-validation` | Extract and validate version |
| `build-release` | Build distribution packages |
| `build-docs` | Generate Sphinx documentation |
| `publish-test` | Upload to TestPyPI (prereleases) |
| `publish-production` | Upload to PyPI (releases) |
| `deploy-docs` | Deploy to GitHub Pages |
| `summary` | Report overall status |

### Required Secrets

| Secret | Purpose |
|--------|---------|
| `TEST_PYPI_API_TOKEN` | TestPyPI API token |
| `PYPI_API_TOKEN` | PyPI API token |
| `GITHUB_TOKEN` | Auto-provided for releases |

---

## Local Testing

Run workflows locally using the package scripts:

```bash
# Build
bash scripts/package/build.sh --check

# Test (various modes)
bash scripts/package/test.sh --mode fast
bash scripts/package/test.sh --mode coverage

# Verify environment
bash scripts/package/verify.sh --mode all
```

---

## Environment Variables

| Variable | Purpose |
|----------|---------|
| `PYTHONUNBUFFERED` | Real-time output |
| `UV_CACHE_DIR` | UV package cache location |
