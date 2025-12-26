# METAINFORMANT Documentation Guide

A comprehensive guide to navigating and understanding the METAINFORMANT documentation system.

## Quick Start

If you're new to METAINFORMANT, follow this path:

1. **[QUICKSTART.md](../QUICKSTART.md)** - Installation and basic usage examples
2. **[Architecture](architecture.md)** - System design and module relationships
3. **[CLI Reference](cli.md)** - Command-line interface for all modules
4. **Domain-specific documentation** - See below for your area of interest

## Documentation Organization

METAINFORMANT documentation is organized hierarchically by domain and purpose:

### Repository Root
- **[README.md](../README.md)** - Project overview, features, and quick examples
- **[QUICKSTART.md](../QUICKSTART.md)** - Fast setup and common commands
- **[AGENTS.md](../AGENTS.md)** - AI assistance in development

### Core Documentation (`docs/`)
- **[index.md](index.md)** - Complete documentation index with navigation
- **[architecture.md](architecture.md)** - System design, module dependencies, data flow
- **[cli.md](cli.md)** - Unified CLI reference for all modules
- **[setup.md](setup.md)** - Installation and environment configuration
- **[testing.md](testing.md)** - Test suite documentation and guidelines
- **[UV_SETUP.md](UV_SETUP.md)** - Package management with `uv`

### Domain Documentation

Each biological domain has its own subdirectory with consistent structure:

| Domain | Location | Primary Focus |
|--------|----------|---------------|
| **Core** | `docs/core/` | Shared utilities (I/O, config, logging, paths) |
| **DNA** | `docs/dna/` | Sequence analysis, alignment, phylogeny, population genetics |
| **RNA** | `docs/rna/` | RNA-seq workflows, amalgkit integration |
| **GWAS** | `docs/gwas/` | Genome-wide association studies |
| **Protein** | `docs/protein/` | Protein sequences, structures, databases |
| **Single-Cell** | `docs/singlecell/` | scRNA-seq preprocessing, clustering, trajectory |
| **Networks** | `docs/networks/` | Biological networks, community detection |
| **Multi-Omics** | `docs/multiomics/` | Cross-omics data integration |
| **Math** | `docs/math/` | Population genetics theory, epidemiology |
| **ML** | `docs/ml/` | Machine learning methods |
| **Visualization** | `docs/visualization/` | Plots, animations, trees |
| **Simulation** | `docs/simulation/` | Synthetic data generation |
| **Quality** | `docs/quality/` | Data quality assessment |
| **Information** | `docs/information/` | Information theory methods |
| **Life Events** | `docs/life_events/` | Life course event analysis |
| **Ontology** | `docs/ontology/` | Gene Ontology, functional annotation |
| **Phenotype** | `docs/phenotype/` | Phenotypic trait analysis |
| **Epigenome** | `docs/epigenome/` | Methylation, ChIP-seq, ATAC-seq |
| **Ecology** | `docs/ecology/` | Community analysis, diversity metrics |

### Domain Documentation Structure

Each domain directory follows a consistent pattern:

```
docs/<domain>/
├── README.md       # Domain overview and quick start
├── AGENTS.md       # AI contributions to this domain
├── index.md        # Detailed index of all submodules
└── <topic>.md      # Topic-specific documentation
```

## Finding What You Need

### By Task

| Task | Documentation |
|------|---------------|
| Install METAINFORMANT | [QUICKSTART.md](../QUICKSTART.md), [setup.md](setup.md) |
| Run CLI commands | [cli.md](cli.md) |
| Analyze DNA sequences | [docs/dna/](dna/index.md) |
| Run RNA-seq workflow | [docs/rna/workflow.md](rna/workflow.md) |
| Perform GWAS | [docs/gwas/workflow.md](gwas/workflow.md) |
| Analyze single-cell data | [docs/singlecell/](singlecell/index.md) |
| Visualize results | [docs/visualization/](visualization/index.md) |
| Run tests | [testing.md](testing.md) |
| Understand architecture | [architecture.md](architecture.md) |

### By Module

Source code documentation mirrors the `src/metainformant/` structure:

| Source Module | Source README | User Documentation |
|--------------|---------------|-------------------|
| `core/` | `src/metainformant/core/README.md` | `docs/core/` |
| `dna/` | `src/metainformant/dna/README.md` | `docs/dna/` |
| `rna/` | `src/metainformant/rna/README.md` | `docs/rna/` |
| ... | ... | ... |

## Documentation Types

### README Files
- Located in each module directory
- Provide quick overview and basic usage
- Link to detailed documentation

### AGENTS.md Files
- Document AI contributions to each module
- Track development history and design decisions
- Present in all major directories (65+ files)

### Topic Documentation
- Detailed coverage of specific functionality
- Code examples and API references
- Performance considerations

### Configuration Documentation
- YAML configuration file format
- Environment variable overrides
- Species-specific configurations

## Best Practices

### Reading Documentation
1. Start with domain `README.md` for overview
2. Check `index.md` for available topics
3. Read topic-specific files for details
4. Reference source code READMEs for API specifics

### Contributing Documentation
1. Follow existing structure and formatting
2. Include practical code examples
3. Update cross-references when adding new docs
4. Place new docs in appropriate `docs/<domain>/` directory
5. Never create documentation in `output/` directory

### Documentation Standards
- Clear, technical writing style
- Consistent markdown formatting
- Runnable code examples
- Cross-references between related topics
- Regular updates with code changes

## Cross-Reference Navigation

### Key Entry Points
- **[Documentation Index](index.md)** - Complete hierarchical navigation
- **[Architecture](architecture.md)** - Module relationships diagram
- **[CLI Reference](cli.md)** - All available commands

### Related Resources
- **[cursorrules/README.md](../cursorrules/README.md)** - Module-specific coding rules
- **[tests/README.md](../tests/README.md)** - Test suite organization
- **[scripts/README.md](../scripts/README.md)** - Workflow scripts

## Getting Help

- **Documentation Issues**: Check for broken links or outdated content
- **Code Questions**: Reference source code READMEs and docstrings
- **Workflow Help**: See domain-specific workflow documentation
- **Bug Reports**: https://github.com/q/metainformant/issues

---

*This guide provides navigation for 70+ documentation files organized across 20+ domains, ensuring comprehensive coverage of METAINFORMANT's bioinformatics capabilities.*





