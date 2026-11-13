# GitHub Actions Workflows

This document describes the CI/CD workflows for the MetaInformAnt project.

## Overview

| Workflow | Trigger | Purpose |
|----------|---------|---------|
| [Build Package](#build-package) | Push, PR, Manual | Build and validate package |
| [Test Suite](#test-suite) | Push, PR, Manual | Run tests across platforms |
| [Release](#release) | Release, Manual | Publish to PyPI |

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
| `quality-checks` | Black, isort, flake8, bandit, pip-audit |
| `test-examples` | Run example scripts |
| `fat-filesystem-test` | Test on FAT-like filesystem |
| `summary` | Report overall status |

The public Test Suite checks out only the parent repository. Private project
submodules are not fetched with the repository-scoped `GITHUB_TOKEN`; their
project-specific tests report explicit skips when unavailable. This workflow
does not validate the nested Hymenoptera or BeeWAS suites. Run those suites
separately in an authenticated project checkout before publishing project changes.

Each pytest invocation inherits a five-minute per-test timeout. Test matrix and
module-lane jobs have a 45-minute limit; other jobs use shorter limits. New runs
cancel superseded runs for the same workflow and branch. A timeout fails the job
and is reported as incomplete verification.

### Test Types

- `fast` - Quick core tests (~15s)
- `network` - Tests requiring internet
- `external` - Tests requiring external CLI tools (muscle, seqkit)
- `all` - Full suite with coverage

---

## Release

**File:** `release.yml`

### Triggers

- GitHub release published (full validation, production publish for non-prereleases)
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

### Structural-variant diagnostics in full coverage

The Ubuntu full-coverage job requires the native pysam import, then runs the PCA
and translocation contract files with their own JUnit report. A failure stops the
job with the complete diagnostic trace before the broad suite consumes its budget.
These preliminary checks do not replace full coverage: after they pass, the same
complete suite, coverage floor and 45-minute job limit remain required.

The full-coverage invocation uses two xdist workers with `loadgroup` scheduling,
which honors existing shared-service groups such as the Ollama tests. OpenMP,
OpenBLAS, MKL and NumExpr each receive a one-thread budget per process. Worker
restarts are disabled so a crashed worker cannot silently resume the run. The
complete `tests/` selection, 60-second full-suite per-test timeout and 74.90%
coverage floor are unchanged; worker coverage is combined by pytest-cov.

A job cancelled at its time limit has no complete coverage result. The required
summary continues to fail for a failed or cancelled matrix job, including a live
network failure; parallel execution does not waive unavailable-service checks.
