# Documentation Directory

This directory contains comprehensive documentation for all METAINFORMANT modules and functionality, organized by domain and feature.

**📖 New to METAINFORMANT documentation?** Start with the **[Documentation Guide](DOCUMENTATION_GUIDE.md)** for complete navigation and best practices.

**⚠️ Package Management**: METAINFORMANT uses `uv` for all Python package management. See **[UV Setup Guide](UV_SETUP.md)** for complete setup instructions. Always use `uv venv`, `uv pip install`, `uv run`, `uv sync`, `uv add`, and `uv remove` - never use `pip` directly.

## Quick Links

- **[Documentation Guide](DOCUMENTATION_GUIDE.md)** - Complete guide to navigating all documentation
- **[Documentation Index](index.md)** - Hierarchical navigation to all docs
- **[Setup Guide](setup.md)** - Installation and environment configuration
- **[Disk Space Management](DISK_SPACE_MANAGEMENT.md)** - ⚠️ **CRITICAL**: External drive temp directory setup
- **[Testing Guide](testing.md)** - Comprehensive testing documentation
- **[Architecture](architecture.md)** - System design and diagrams
- **[CLI Reference](cli.md)** - Command-line interface
- **[AI Agents](AGENTS.md)** - AI assistance in development

## Structure

### Core Documentation (`core/`)
Documentation for shared utilities and infrastructure:
- `cache.md` - Caching mechanisms and configuration
- `config.md` - Configuration management and environment variables
- `db.md` - Database integration and connection management
- `hash.md` - Content hashing and file integrity
- `io.md` - Input/output utilities and file format support
- `logging.md` - Logging configuration and best practices
- `parallel.md` - Parallel processing and thread management
- `paths.md` - Path handling and validation utilities
- `text.md` - Text processing and normalization

### Domain Documentation
Each biological domain has its own subdirectory with detailed documentation:

#### DNA (`dna/`)
- `accessions.md` - Genome accession validation and retrieval
- `alignment.md` - Sequence alignment algorithms and methods
- `codon.md` - Codon usage analysis and genetic code
- `composition.md` - Sequence composition and GC content analysis
- `consensus.md` - Consensus sequence generation
- `distances.md` - Evolutionary distance calculations
- `fastq.md` - FASTQ format processing and quality analysis
- `index.md` - DNA domain overview and module index
- `motifs.md` - Motif discovery and pattern analysis
- `msa.md` - Multiple sequence alignment
- `mutations.md` - Mutation analysis and variant calling
- `ncbi.md` - NCBI database integration
- `phylogeny.md` - Phylogenetic tree construction and analysis
- `population.md` - Population genetics and diversity metrics
- `restriction.md` - Restriction enzyme analysis
- `sequences.md` - DNA sequence manipulation and I/O
- `transcription.md` - DNA to RNA transcription
- `translation.md` - Genetic code translation
- `variants.md` - Variant detection and analysis

#### RNA (`rna/`)
- `amalgkit/` - Amalgkit integration documentation
  - `steps/README.md` - Individual workflow step documentation
- `CONFIGURATION.md` - RNA workflow configuration
- `index.md` - RNA domain overview
- `workflow.md` - Complete workflow orchestration

#### GWAS (`gwas/`)
- `README.md` - GWAS module overview and quick start
- `workflow.md` - Step-by-step GWAS workflow guide
- `config.md` - Configuration reference and options
- `pbarbatus_config.md` - Example configuration for P. barbatus
- `verification_report.md` - Verification and validation report

#### Protein (`protein/`)
- `index.md` - Protein domain overview
- `proteomes.md` - Proteome analysis and databases

#### Math (`math/`)
- `coalescent.md` - Coalescent theory and genealogy
- `ddm.md` - Drift-diffusion models for decision making
- `epidemiology.md` - Disease modeling and epidemiology
- `index.md` - Mathematical biology overview
- `ld.md` - Linkage disequilibrium analysis
- `popgen.md` - Population genetics theory
- `price.md` - Price equation and selection analysis
- `selection.md` - Natural selection modeling

#### Machine Learning (`ml/`)
- `index.md` - Machine learning methods overview

#### Multi-Omics Integration (`multiomics/`)
- `index.md` - Multi-omics integration overview

#### Network Analysis (`networks/`)
- `community.md` - Network community detection
- `graph.md` - Graph algorithms and analysis
- `index.md` - Network analysis overview
- `pathway.md` - Pathway analysis and enrichment
- `ppi.md` - Protein-protein interaction networks
- `regulatory.md` - Gene regulatory networks

#### Specialized Domains
- `ecology/` - Ecological metadata and community analysis
- `epigenome/` - DNA methylation and chromatin analysis
- `ontology/` - Gene Ontology and functional annotation
- `phenotype/` - Phenotypic trait analysis
- `quality/` - Data quality assessment
- `simulation/` - Synthetic data generation
- `singlecell/` - Single-cell transcriptomic analysis
- `visualization/` - Plotting and animation utilities
- `information/` - Information-theoretic analysis (entropy, mutual information, semantic similarity)
- `life_events/` - Life course and event sequence analysis

## Documentation Standards

### Writing Guidelines
- Use clear, concise technical writing
- Include practical code examples
- Document both API and CLI interfaces
- Explain biological context and relevance
- Provide performance considerations

### Organization Principles
- Each module has comprehensive documentation
- Cross-references between related modules
- Consistent formatting and structure
- Regular updates with code changes

### Maintenance
- Documentation is updated with each major release
- Technical writers review for clarity and accuracy
- Community contributions are welcome and documented
- Automated checks ensure broken links are identified

## Contributing to Documentation

### Adding New Documentation
1. Follow the established format and structure
2. Include practical examples and use cases
3. Cross-reference related modules
4. Update indices and navigation

### Updating Existing Documentation
1. Keep examples current with API changes
2. Add new features and capabilities
3. Fix broken links and references
4. Improve clarity and completeness

## Navigation

- Start with domain-specific index files for overviews
- Use cross-references to explore related functionality
- Check the main project README for high-level navigation
- See `src/metainformant/*/README.md` for module-specific details

This documentation provides comprehensive coverage of all METAINFORMANT functionality, from basic usage to advanced implementation details.
