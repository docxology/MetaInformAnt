# Changelog

All notable changes to METAINFORMANT are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.7] - 2026-03-01

### Added

- **`scripts/rna/fetch_pending_metadata.py`** — batch metadata fetcher for SRA runs not yet in the pipeline manifest; handles rate limiting and incremental writes
- **`scripts/rna/find_missing_tissues.py`** — scans all species configs and tissue patches to report bioprojects with incomplete tissue coverage; used for targeted tissue-patch research
- **`scripts/rna/extract_bioproject_details.py`** — fetches NCBI BioProject title and SRA run counts for a list of accessions; used alongside `find_missing_tissues.py`
- **`scripts/rna/report_completed.py` / `report_completed.sh`** — reports quantified sample counts per species from disk, distinct from log-based counts; useful for progress monitoring independent of the orchestrator log
- **`scripts/rna/link_indices.sh`** — symlinks shared kallisto indices across species output dirs to avoid redundant index builds
- **`scripts/gwas/run_pbarbatus_gwas.py`** — end-to-end GWAS runner for *Pogonomyrmex barbatus*: downloads SRA reads, aligns with BWA, calls variants with BCFtools, and runs association with real phenotype data

### Changed

- **`.gitignore`** — added `*.heartbeat.json`, `**/.downloads/`, and `*.safely_removed` to prevent pipeline sentinel files from polluting git history
  - `*.heartbeat.json`: transient JSON progress files written during downloads; write-only, never read back, auto-recreated by active jobs
  - `*.safely_removed`: marker files confirming FASTQ deletion after successful kallisto quantification; redundant given that a non-empty `*_abundance.tsv` is the canonical proof of work
- **Pipeline housekeeping**: removed 3,898 existing heartbeat files and 3,055 `.safely_removed` markers from disk (∼7k stale sentinel files cleaned up)

### Pipeline Status (as of this release)

- *Apis mellifera*: **2,661 samples quantified** (of 5,715 total; ~3,900 auto-skipped for size >4B bases, pending second pass)
- *Acromyrmex echinatior*: **3 samples quantified** (pipeline just started)
- All other 21 species: queued — pipeline processes species sequentially

## [0.2.6] - 2026-02-20


### Fixed

- t-SNE perplexity auto-adjustment: threshold corrected to `(n_samples-1)/3` per sklearn
- `compute_diffusion_map`: removed `sigma=1.0` from `eigsh` that caused singular matrix crash; added diagonal regularization and dense fallback
- `compute_neighbors`: `n_pcs` now stored in neighbor params for downstream inspection
- Disk space: cleaned `.uv-cache` (476M) and ran `git gc --prune=all` (315MB recovered)

## [0.2.5] - 2026-02-20

### Changed

- Version bump to 0.2.5 reflecting accumulated test hardening
- CHANGELOG updated with v0.2.3–v0.2.4 entries
- ENA mention added to `scripts/rna/AGENTS.md`

## [0.2.4] - 2026-02-20

### Fixed

- `SingleCellData` diffusion map: `len(data)` → `data.n_obs` (TypeError fix)
- `ENADownloader` re-exported from `rna.retrieval.__init__`
- RNA test config YAML paths corrected (`blue/` → `output/`)
- `run_workflow.py` startup: deferred heavy imports for fast `--help` (<1s)
- Early-exit for no-args/`--list-configs` before heavy imports

## [0.2.3] - 2026-02-20

### Fixed

- **Protein**: 20 relative → absolute imports in `orchestration.py`; `read_taxon_ids` return type `List[str]`
- **Visualization**: 7 plot function re-exports; config self-reference; 22 analysis re-exports; phylo import fix; `WONG` palette constant
- **Quality**: 6 standalone contamination functions rewritten with per-sequence-index returns; `per_base_quality` quantiles guard for <2 data points; `generate_contamination_report` METAINFORMANT header
- **RNA**: `calculate_md5`, `clean_stagnant_file`, `verify_gzip_integrity` implemented in `ena_downloader.py`
- **Phenotype**: `extract_phenotypes_from_events`, `aggregate_temporal_phenotypes`, `map_events_to_traits` rewritten for `life_events.core.events` compatibility; `Measurement`/`MorphometricProfile` re-exported from `morphological/__init__`; backward-compat `career_changes`/`health_score` keys
- **Singlecell**: PCA `n_components` auto-clamp to `min(n_samples-1, n_features)`

## [0.2.2] - 2026-02

### Added

- Comprehensive repo-wide documentation audit and consistency fixes
- Specialized modules: longread, metagenomics, pharmacogenomics, spatial, structural_variants
- Enhanced Mermaid architecture diagrams with correct node references

### Changed

- Updated all documentation to use consistent plot count (70+)
- Synchronized version references across README, CHANGELOG, PAI

### Fixed

- Broken Mermaid node IDs in README data flow diagrams
- Stale version references (0.2.0 → 0.2.2) across all docs

## [0.2.1] - 2026-02

### Added

- Streaming RNA-seq orchestrator with size-ordered scheduling
- Tiered fallback recovery for failed pipeline steps
- Concurrency-hardened workspace initialization (mkdir race conditions)
- Clean-restart operational protocols for large-scale transcriptomic runs

### Changed

- Pipeline optimized for 16 parallel workers with multi-threaded scaling
- Improved subprocess timeout handling (30-minute cap on downstream calls)

### Fixed

- Pipeline stalls from missing subprocess timeouts
- FASTQ extraction directory structure handling
- Permission errors on external volume log writes

## [0.2.0] - 2025-12

### Added

- Complete UV package management integration (replaces pip)
- Enhanced RNA workflow with automated genome indexing
- Comprehensive GWAS visualization suite (manhattan, QQ, regional plots)
- Information theory module with syntactic and semantic measures
- Life events analysis with sequence modeling
- Single-cell analysis with Leiden clustering
- Multi-omics integration workflows
- Progress tracking with real-time TUI
- Intelligent caching for expensive computations

### Changed

- Python 3.11+ minimum requirement (previously 3.10+)
- Reorganized configuration system with environment variable overrides
- Improved error handling with context-aware messages
- Enhanced parallel processing with better thread management

### Fixed

- Import errors reduced from ~225 to ~63 (72% improvement)
- Test collection success rate improved to 87%
- Memory optimization for large datasets

### Removed

- pip-based installation (use `uv` exclusively)
- Legacy configuration formats

## [0.1.0] - 2025-06

### Added

- Initial release with core bioinformatics modules
- DNA analysis: sequences, alignment, population genetics, phylogeny
- RNA analysis: amalgkit integration, workflow orchestration
- Protein analysis: sequences, structures, AlphaFold integration
- GWAS pipeline: association testing, QC, visualization
- Network analysis: community detection, centrality measures
- Visualization: 57+ plot types
- Quality control: FASTQ analysis, contamination detection
- Machine learning: classification, regression, feature selection
- Mathematical biology: coalescent, population genetics theory

### Notes

- NO_MOCKING_POLICY established: all tests use real implementations
- UV-only package management policy established
- Output directory isolation policy established
