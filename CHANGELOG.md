# Changelog

All notable changes to METAINFORMANT are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
