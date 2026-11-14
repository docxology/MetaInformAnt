# GitHub Actions Workflows

This document describes the CI/CD workflows for the MetaInformAnt project.

## Overview

| Workflow | Trigger | Purpose |
|----------|---------|---------|
| [Build Package](#build-package) | Push, PR, Manual | Build and validate package |
| [Test Suite](#test-suite) | Push, PR, Manual | Run tests across platforms |
| [Release](#release) | Release, Manual | Publish to PyPI |
| Security Checks (`security.yml`) | Push and PR to `main`, Manual | Blocking Bandit, pip-audit and pattern audits |

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
| `build` | Build wheel/sdist and install the wheel on Python 3.14 and 3.15 |
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
| `test` | Fast lane on ubuntu/macos × Python 3.14/3.15; network, external and all (coverage) lanes on ubuntu × Python 3.14 |
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

### Python versions

The supported floor is Python 3.14 (`requires-python = ">=3.14"`), and the
repository's `.python-version` selects 3.14 for local uv commands. Every
single-version job runs on 3.14. The `test` job sets `UV_PYTHON` to the matrix
version so the re-sync inside `test.sh` stays on that interpreter.

The fast lane installs only the core dependencies and the `test-fast` group
(`uv sync --frozen --group test-fast`), which is the set `test.sh --mode fast`
runs against. That install works on Python 3.15. The full extras do not install
on 3.15 (ray has no cp315 wheels), so the network, external and all lanes run on
3.14 only.

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
| `install-smoke` | Install the built wheel in clean Python 3.14 and 3.15 environments |
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
