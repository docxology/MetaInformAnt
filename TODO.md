# Metainformant Project Master TODO

This document tracks high-level architectural goals, technical debt, and pending enhancements, particularly focusing on the industrial-scale Amalgkit integration.

## 🚀 Priority: Amalgkit & RNA-Seq Pipeline Hardening

- [x] **Standardize Index Complexity Filtering**
  - The `fix_harpegnathos_index.py` script proved essential for `Harpegnathos`.
  - **Goal**: Create a formal `metainformant.rna.index.filter` module to automatically detect and clean filtering artifacts (ncRNA, duplicates, rRNA) for *any* species.
  - **Status**: Implemented `IndexComplexityManager` in `src/metainformant/rna/amalgkit/index_prep.py`.

- [x] **Unified Download Orchestrator**
  - Currently, we have `process_apis_mellifera.py` (NCBI/fasterq-dump) and `process_apis_mellifera_ena.py` (ENA/direct).
  - **Goal**: Merge these into `StreamingPipelineOrchestrator`.
  - **Logic**: Try ENA first (fast), fallback to NCBI (comprehensive but slow).
  - **Status**: Implemented `StreamingPipeline` in `metainformant.rna.engine.orchestrator`.

- [ ] **Dynamic Resource Allocation**
  - **Challenge**: `Harpegnathos` stalled on high-complexity index; `Apis` requires high I/O.
  - **Goal**: Implement resource profiles (e.g., `--profile high_mem` vs `--profile high_io`) in `orchestrator.py`.

## 📚 Documentation & Rules

- [ ] **Entrench "Zero-Mock" Policy**
  - Ensure all new modules (longread, singlecell) have functional tests interacting with real (or minimized real) data, not mocks.
  - **verification**: Check `tests/` for any regression to mocking.

- [ ] **Architecture Documentation**
  - Update `docs/rna/README.md` to reflect the "ENA First" strategy.
  - **Consolidate Documentation**: `docs/rna/README.md` is a stale duplicate of `src/metainformant/rna/README.md`. Decide on a Single Source of Truth (likely `src/`) and symlink or redirect.
  - Document the `work/` vs `fastq/` directory structure explicitly to avoid confusion (like the symlink issue).
  - Document TUI tools (`monitor_tui.py`, `run_workflow_tui.py`).

- [x] **Script Consolidation**
  - Merge `process_apis_mellifera.py`, `_ena.py`, and `_parallel.py` into a single robust `process_species.py` with strategy flags (`--strategy ena`, `--strategy sra`).
  - **Status**: Completed in `scripts/rna/process_species.py`.

## 🛠️ Infrastructure

- [x] **Infrastructure Stability**
  - **Fix Circular Dependency**: Refactored `src/metainformant/__init__.py` to remove top-level imports.
  - **Unit Test Harness**: Setup for unblocked unit testing.

- [x] **MCP Tooling for Amalgkit**
  - Created `src/metainformant/mcp/tools/amalgkit_monitor.py` (verified via dry-run).

- [x] **Error Recovery Automation**
  - **Goal**: Auto-detect "stalled" processes (0 CPU for >1h) and kill/restart them.
  - **Mechanism**: Implemented `ProcessWatchdog` in `src/metainformant/core/utils/watchdog.py`.
  - **Status**: Integrated into `Orchestrator` with `--watchdog` flag. Verified via `tests/core/test_watchdog_integration.py`. <!-- id: 163 -->

## 🔬 Scientific Modules

- [ ] **GWAS Integration**
  - Verify `gwas` module integration points with `rna` (eQTL analysis).
  - Ensure `gwas` pipe accepts upstream `amalgkit` quantification tables directly.

---

*Last Updated: 2026-02-18*
