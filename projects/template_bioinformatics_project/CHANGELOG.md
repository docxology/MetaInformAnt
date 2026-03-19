# Changelog

All notable changes to this template are documented here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.0.0/). Version: [SemVer](https://semver.org/).

---

## [Unreleased]

> _(Changes staged for the next release)_

---

## [0.1.0] — 2026-03-19

### Added

- `pyproject.toml` replacing `requirements.txt` (uv-compatible, pinnable, declarative)
- `main.py` top-level CLI entry point with `--stage` and `--config` routing
- `run.sh` stage-selectable shell runner with `--clean`, `--synthetic`, and `--stage` flags
- `scripts/02_analyze_results.py` — downstream statistical analysis (summary, correlation, PCA)
- `scripts/03_visualize.py` — distribution grid + correlation heatmap with all parameters from config
- `scripts/99_create_synthetic_data.py` — reproducible synthetic test data generator
- `tests/conftest.py` and `tests/test_pipeline.py` — Zero-Mock end-to-end test suite
- `doc/` documentation suite: `index.md`, `architecture.md`, `configuration.md`, `data.md`
- `doc/stages/` per-stage technical references for all four pipeline stages
- `AGENTS.md` — AI agent responsibilities and thin-orchestration standards
- `SPEC.md` — living project specification and extension guidelines
- `CHANGELOG.md` — this file
- `config/samples.tsv` — sample manifest template
- `.github/workflows/ci.yml` — GitHub Actions CI pipeline
- Enhanced `config/default.yaml` with `metadata`, `analysis`, `visualization`, and `logging` sections
- `--force` flag on all pipeline scripts for unconditional reprocessing
- Idempotency guards in all numbered scripts

### Changed

- `scripts/01_process_data.py` — rewritten with real functional logic (missing-value filtering,
  z-score normalisation, row-count filtering, idempotency, timed logging, `--force` flag)
- `.gitignore` — extended with uv cache, dist/, egg-info/, pytest cache, coverage outputs
- `README.md` — comprehensive upgrade with quickstart, pipeline flow, directory tree, doc links

### Principles

- **Zero hardcoding**: all paths and thresholds from `config/default.yaml`
- **Immutability**: `data/raw/` is read-only
- **Idempotency**: every stage skips existing outputs; use `--force` to reprocess
- **Zero-Mock testing**: tests exercise real code against synthetic data
- **uv-first**: environment management via `uv sync`
