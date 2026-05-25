# MetaInformAnt — Master TODO

> **Version**: 0.3.0 · **Last Updated**: 2026-05-25
>
> Historical note: earlier sections preserve the March 2026 pipeline roadmap.
> Treat the status block below as the current checkout snapshot.

---

## Recent Accomplishments (v0.2.5 → v0.2.7)

-**GCP Cloud Deployment** — Full `cloud/` module with VM lifecycle management, Docker pipelines, autonomous bootstrapping, genome prep, and micromamba integration
-**Streaming Orchestrator** — `streaming_orchestrator.py` with SQLite-backed progress tracking, size-ordered scheduling, and 16-worker concurrency
-**Download Resilience** — ENA-first + NCBI `fasterq-dump` fallback, thread-safe symlink handling, concurrent SQLite locking with retry
-**Index Complexity Filtering** — `IndexComplexityManager` in `index_prep.py` for automatic ncRNA/duplicate/rRNA cleaning on any species
-**Process Watchdog** — `ProcessWatchdog` in `core/utils/watchdog.py` with `--watchdog` flag; auto-kills stalled jobs (0 CPU > 1h)
-**Pipeline Housekeeping** — Cleaned ~7k stale sentinel files (heartbeats + `.safely_removed` markers); `.gitignore` hardened
-**Script Consolidation** — All legacy scripts physically deleted; fully superseded by `run_workflow.py`
-**Dynamic Resource Profiles** — `ResourceProfile` dataclass with 4 presets (`high_mem`, `high_io`, `default`, `minimal`) and species auto-detection in `orchestrator.py`
-**MCP Amalgkit Monitor** — `mcp/tools/amalgkit_monitor.py` dry-run verified
-**GWAS P. barbatus Runner** — End-to-end real-data GWAS via BWA + BCFtools + SRA reads
-**GWAS ↔ RNA Integration** — `gwas/data/expression.py` reads amalgkit `abundance.tsv` directly; `eqtl_coloc` and `from_rna_expression` exported and functional
-**Infrastructure Fixes** — Circular dependency removal, t-SNE perplexity correction, diffusion map dense fallback, 72% import-error reduction

### Documentation (completed 2026-03-04)

-**RNA docs consolidated** — `docs/rna/README.md` → symlink to `src/metainformant/rna/README.md` (single source of truth)
-**ENA-first strategy documented** — Updated RNA README with streaming orchestrator architecture and dual-tier download strategy
-**TUI tools documented** — Added `monitor_tui.py` and `run_workflow_tui.py` sections to `docs/cli.md`
-**Directory conventions documented** — `work/` vs `fastq/` vs `output/` hierarchy with symlink explanation in `docs/cli.md`
-**Zero-Mock audit passed** — Only 1 dead mock import found and removed (`tests/rna/test_ena_downloader.py`); no actual mock usage in test suite
-**Stale configs verified** — `amalgkit_test.yaml` and `amalgkit_template.yaml` are referenced by tests/scripts and remain valid

---

## In Progress: Amalgkit Multi-Species Pipeline (GCP Cloud Run)

### Pipeline Status (GCP VM: `metainformant-pipeline`)

| Species | Samples Quantified | Status |
|---------|-------------------|--------|
| *Harpegnathos saltator* | 368 | [DONE] Done |
| *Solenopsis invicta* | 323 | [DONE] Done |
| *Ooceraea biroi* | 217 | [DONE] Done |
| *Temnothorax longispinosus* | 148 | [DONE] Done |
| *Temnothorax nylanderi* | 115 | [DONE] Done |
| *Linepithema humile* | 111 | [DONE] Done |
| *Monomorium pharaonis* | 98 | [DONE] Done |
| *Pogonomyrmex barbatus* | 78 | [DONE] Done |
| *Apis mellifera* | 57+ | Active (~5,700 total) |
| *Temnothorax americanus* | 54 | [DONE] Done |
| 10 other ant species | 182 | [DONE] Done |
| *Cardiocondyla obscurior* | 3 | Starting |
| *Atta/Camponotus* | 0 | ⏳ Metadata deployed |
| **TOTAL** | **1,751+** | **20/22 species active** |

### Pipeline Tasks

- [x] **Deploy pipeline to GCP** — Docker container with 28 workers, spot pricing
- [x] **Deploy metadata for blocked species** — Uploaded local `metadata_selected.tsv` for amellifera + 3 species
- [ ] **Complete Apis mellifera processing** — ~5,700 samples total, ~57 done, ETA ~40h
- [ ] **Process remaining 2 blocked species** — Atta, Camponotus (metadata deployed, waiting on pipeline loop)
- [ ] **Download results to local** — `bash scripts/cloud/download_results.sh output/amalgkit`
- [ ] **Post-quant analysis** — Run `merge`, `curate`, `sanity` steps locally after download
- [ ] **Cross-species expression matrix** — Generate unified multi-species expression table

---

## Upcoming: v0.3.0 Milestones

### Cloud & Scalability

- [ ] **GCP parallel multi-species processing** — Deploy multiple species concurrently on separate VMs instead of sequential local processing
- [ ] **Auto-scaling disk management** — Monitor and resize GCP disks dynamically during large pipeline runs

### Scientific Analysis

- [ ] **Cross-species comparative transcriptomics** — Downstream analysis of completed multi-species expression matrices (PCA, differential expression, phylogenetic signal)
- [ ] **Expand GWAS to additional species** — Replicate `run_pbarbatus_gwas.py` pattern for other species with public SRA/variant data

### Module Maturation

- [ ] **Long-read sequencing** (65%) — Complete PacBio/ONT assembly and error correction workflows
- [ ] **Metagenomics** (60%) — Finish taxonomic profiling and functional annotation pipelines
- [ ] **Spatial transcriptomics** (50%) — Tissue mapping and spatial statistics integration
- [ ] **Single-cell** (74%) — Complete trajectory analysis and pseudo-time workflows
- [ ] **Metabolomics** (50%) — MS data processing and pathway mapping

---

## Current Stabilization Status

- [x] **Test collection** — `uv run pytest --collect-only -q` collected 7,736 tests in the 2026-05-25 stabilization pass
- [x] **Local non-network suite** — `uv run pytest tests -q -m 'not network and not external_tool' --tb=short` passed with 7,495 passed, 71 skipped, 170 deselected
- [ ] **Runner/doc drift** — Keep `scripts/package/test.sh`, README/SPEC examples, and generated report destinations synchronized with the nested test layout

---

*Next release target: **v0.3.0** — Complete multi-species quantification + cross-species analysis*
