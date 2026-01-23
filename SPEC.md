# METAINFORMANT Technical Specification

Comprehensive bioinformatics toolkit for multi-omic analysis. Domain-driven, modular architecture designed for performance and scientific rigor.

## Design Philosophy

- **Domain-Driven Design (DDD)**: Logic partitioned into biological domains (DNA, RNA, Protein, etc.)
- **No-Mock Testing**: All tests use real implementations; mocks are strictly prohibited
- **UV Package Management**: Dependencies managed exclusively via `uv` (never pip)
- **AI-Native Documentation**: High-fidelity function indices for AI agent navigation

## Repository Structure

```
metainformant/
├── src/metainformant/     # Source code (19 domain modules)
│   ├── core/              # Shared infrastructure (I/O, config, logging)
│   ├── dna/               # Genomic analysis, alignment, population genetics
│   ├── rna/               # Transcriptomic workflows, Amalgkit integration
│   ├── protein/           # Proteomic analysis, structure modeling
│   ├── gwas/              # Genome-wide association studies
│   ├── epigenome/         # Methylation, ChIP-seq, ATAC-seq
│   ├── networks/          # Biological networks, community detection
│   ├── multiomics/        # Multi-omic data integration
│   ├── singlecell/        # Single-cell RNA-seq analysis
│   ├── visualization/     # 57+ plot types, publication-quality output
│   ├── quality/           # QC metrics, contamination detection
│   ├── ml/                # Machine learning pipelines
│   ├── math/              # Population genetics theory, coalescent
│   ├── information/       # Information theory (entropy, MI)
│   ├── ontology/          # GO analysis, semantic similarity
│   ├── phenotype/         # Trait analysis, curation
│   ├── ecology/           # Community diversity
│   ├── simulation/        # Synthetic data generation
│   └── life_events/       # Event sequence analysis
├── scripts/               # Thin wrapper orchestrators
├── tests/                 # Pytest test suite (real implementations only)
├── docs/                  # Documentation by domain
├── config/                # YAML configuration templates
├── data/                  # Input data (read-mostly)
└── output/                # Program-generated results (ephemeral)
```

## Architectural Patterns

### 1. Sequential Failover (I/O)
Prioritize local data sources before remote acquisition (NCBI, SRA).

### 2. Thin Wrapper Orchestration
Scripts in `scripts/` are thin wrappers around core methods. Business logic resides in `src/`.

### 3. Configuration with Environment Overrides
YAML configs in `config/` can be overridden via environment variables with domain prefixes (`AK_`, `GWAS_`, `DNA_`, etc.).

### 4. Output Isolation
All program-generated results go to `output/`. Never create documentation or reports in `output/`.

## Implementation Standards

### Testing Policy
- **No Mocking**: Tests use real implementations with actual file I/O and API calls
- **Graceful Skips**: When external dependencies unavailable, skip with clear messages
- **Markers**: `@pytest.mark.network`, `@pytest.mark.external_tool`, `@pytest.mark.slow`

### Code Quality
- Python 3.11+ minimum
- Black formatting (120 char lines)
- mypy type checking (strict)
- All functions must have type hints

### Documentation
- Module README.md: User-facing documentation and examples
- docs/<domain>/: Extended guides and tutorials

## Data Flow

```
Raw Data (FASTQ/VCF) → Preprocessing & QC → Domain Analysis → Multi-Omic Integration → Visualization
         ↑                    ↑                    ↑                    ↑
    data/ directory    quality/ module      domain modules      visualization/
```

## Cross-Module Communication

Modules communicate via:
- Standard Python protocols
- Shared data structures (`pandas` DataFrames, `numpy` arrays)
- Core infrastructure (`metainformant.core`) for I/O, config, and logging
