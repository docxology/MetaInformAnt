# METAINFORMANT Documentation Guide

Complete guide to navigating and using the METAINFORMANT documentation system.

## Quick Navigation

### Essential Starting Points

- **[README.md](README.md)** - Documentation structure and organization
- **[index.md](index.md)** - Complete navigation to all documentation
- **[setup.md](setup.md)** - Installation and environment configuration
- **[testing.md](testing.md)** - Comprehensive testing documentation
- **[architecture.md](architecture.md)** - System design and component relationships
- **[cli.md](cli.md)** - Command-line interface reference

### Core Documentation

**Core Utilities** - Foundational infrastructure for all modules

- **[Core Overview](core/README.md)** - Introduction to core utilities
- [cache.md](core/cache.md) - JSON caching with TTL support
- [config.md](core/config.md) - Configuration management and environment variables
- [db.md](core/db.md) - Database integration and connection management
- [hash.md](core/hash.md) - Content hashing and file integrity
- [io.md](core/io.md) - File I/O utilities and format support
- [logging.md](core/logging.md) - Logging configuration
- [parallel.md](core/parallel.md) - Parallel processing utilities
- [paths.md](core/paths.md) - Path handling and validation
- [text.md](core/text.md) - Text processing and normalization

### Domain-Specific Documentation

#### DNA Analysis
**[DNA Overview](dna/index.md)** - Comprehensive DNA analysis toolkit

- [sequences.md](dna/sequences.md) - DNA sequence manipulation
- [alignment.md](dna/alignment.md) - Pairwise sequence alignment
- [msa.md](dna/msa.md) - Multiple sequence alignment
- [phylogeny.md](dna/phylogeny.md) - Phylogenetic tree construction
- [population.md](dna/population.md) - Population genetics metrics
- [fastq.md](dna/fastq.md) - FASTQ format processing
- [codon.md](dna/codon.md) - Codon usage analysis
- [composition.md](dna/composition.md) - Sequence composition
- [distances.md](dna/distances.md) - Evolutionary distances
- [transcription.md](dna/transcription.md) - DNA to RNA transcription
- [translation.md](dna/translation.md) - Genetic code translation
- [variants.md](dna/variants.md) - Variant analysis

#### RNA Analysis
**[RNA Overview](rna/index.md)** - Transcriptomic analysis workflows

- [workflow.md](rna/workflow.md) - Complete workflow orchestration
- [configs.md](rna/configs.md) - Workflow configuration
- [steps.md](rna/steps.md) - Individual workflow steps

#### GWAS (Genome-Wide Association Studies)
**[GWAS Overview](gwas/index.md)** - Complete GWAS pipeline

- [GWAS README](gwas/README.md) - Module overview and features
- [workflow.md](gwas/workflow.md) - Step-by-step workflow guide
- [config.md](gwas/config.md) - Configuration reference
- [pbarbatus_config.md](gwas/pbarbatus_config.md) - Example configuration
- [verification_report.md](gwas/verification_report.md) - Validation report

#### Protein Analysis
**[Protein Overview](protein/index.md)** - Protein structure and sequence analysis

- [proteomes.md](protein/proteomes.md) - Proteome databases and analysis

#### Single-Cell Genomics
**[Single-Cell Overview](singlecell/index.md)** - Single-cell transcriptomics

- [preprocessing.md](singlecell/preprocessing.md) - Data preprocessing
- [dimensionality.md](singlecell/dimensionality.md) - Dimensionality reduction
- [clustering.md](singlecell/clustering.md) - Cell clustering
- [trajectory.md](singlecell/trajectory.md) - Trajectory analysis
- [visualization.md](singlecell/visualization.md) - Visualization methods
- [integration.md](singlecell/integration.md) - Data integration

#### Mathematical Biology
**[Math Overview](math/index.md)** - Theoretical and quantitative biology

- [coalescent.md](math/coalescent.md) - Coalescent theory
- [ddm.md](math/ddm.md) - Drift-diffusion models
- [dynamics.md](math/dynamics.md) - Population dynamics
- [epidemiology.md](math/epidemiology.md) - Disease modeling
- [ld.md](math/ld.md) - Linkage disequilibrium
- [popgen.md](math/popgen.md) - Population genetics theory
- [price.md](math/price.md) - Price equation
- [selection.md](math/selection.md) - Natural selection models

#### Machine Learning
**[ML Overview](ml/index.md)** - Machine learning for bioinformatics

Complete ML pipeline documentation with classification, regression, feature selection, and model validation.

#### Network Analysis
**[Networks Overview](networks/index.md)** - Biological network analysis

- [community.md](networks/community.md) - Community detection
- [graph.md](networks/graph.md) - Graph algorithms
- [pathway.md](networks/pathway.md) - Pathway analysis
- [ppi.md](networks/ppi.md) - Protein-protein interactions
- [regulatory.md](networks/regulatory.md) - Gene regulatory networks

#### Multi-Omics Integration
**[Multi-Omics Overview](multiomics/index.md)** - Cross-omic data integration

Comprehensive guide to integrating genomics, transcriptomics, proteomics, and other omic layers.

#### Quality Control
**[Quality Overview](quality/index.md)** - Data quality assessment

- [fastq.md](quality/fastq.md) - FASTQ quality analysis (comprehensive 645-line guide)

#### Visualization
**[Visualization Overview](visualization/index.md)** - Data visualization

- [trees.md](visualization/trees.md) - Phylogenetic tree plotting
- [plots.md](visualization/plots.md) - General plotting utilities
- [animations.md](visualization/animations.md) - Animation generation

#### Simulation
**[Simulation Overview](simulation/index.md)** - Synthetic data generation

- [sequences.md](simulation/sequences.md) - Sequence generators
- [rna_counts.md](simulation/rna_counts.md) - RNA count simulation
- [agents.md](simulation/agents.md) - Agent-based modeling

#### Specialized Domains

- **[Ontology](ontology/index.md)** - Gene Ontology and functional annotation
- **[Phenotype](phenotype/index.md)** - Phenotypic data curation
- **[Epigenome](epigenome/index.md)** - DNA methylation and chromatin
- **[Ecology](ecology/index.md)** - Community ecology and metadata

## Documentation by Use Case

### Getting Started
1. [setup.md](setup.md) - Environment setup with UV
2. [cli.md](cli.md) - Command-line interface basics
3. [Core README](core/README.md) - Core utilities overview
4. Domain-specific index (e.g., [DNA](dna/index.md), [RNA](rna/index.md))

### Running Workflows
1. **RNA-seq**: [RNA workflow](rna/workflow.md) → [RNA configs](rna/configs.md)
2. **GWAS**: [GWAS workflow](gwas/workflow.md) → [GWAS config](gwas/config.md)
3. **Single-cell**: [Single-cell index](singlecell/index.md) → Individual analysis docs

### API Development
1. [Core utilities](core/README.md) - Use shared infrastructure
2. Domain-specific READMEs - Module APIs and examples
3. [architecture.md](architecture.md) - System design patterns

### Testing and Quality
1. [testing.md](testing.md) - Testing policy and execution
2. [comprehensive_test_analysis.md](comprehensive_test_analysis.md) - Test coverage report
3. Module-specific test documentation

## Documentation Standards

### File Organization
- **README.md** - Module overview, quick start, features
- **index.md** - Complete navigation and reference
- **AGENTS.md** - AI assistance documentation
- **Specific docs** - Detailed guides for individual components

### Navigation Patterns
```
Top-level README → docs/index.md → Domain index → Specific documentation
                → docs/README.md (documentation structure)
```

### Cross-References
All documentation includes cross-references to:
- Related modules and functionality
- Core utilities used (see [Core Utilities](./core.md))
- Integration examples
- Testing documentation (see [Testing Guide](./testing.md))
- Orchestrator scripts (see [Scripts README](../../scripts/README.md))
- Configuration templates (see [Config README](../../config/README.md))
- CLI commands (see [CLI Reference](./cli.md))

## AI Assistance Documentation

AI contributions are documented at multiple levels:

- **Repository level**: [AGENTS.md](../AGENTS.md) - Overall AI assistance
- **Documentation level**: [docs/AGENTS.md](AGENTS.md) - Documentation development
- **Module level**: Each module has its own AGENTS.md (e.g., [core/AGENTS.md](../src/metainformant/core/AGENTS.md))
- **Domain level**: Domain-specific AI contributions (e.g., [gwas/AGENTS.md](gwas/AGENTS.md))

## Additional Resources

### Comprehensive Guides
- [comprehensive_uv_setup.md](comprehensive_uv_setup.md) - UV environment setup and optimization
- [comprehensive_test_analysis.md](comprehensive_test_analysis.md) - Complete test suite analysis

### Development Resources
- [architecture.md](architecture.md) - System architecture with diagrams
- [testing.md](testing.md) - Testing policy (NO_MOCKING_POLICY)
- [core.md](core.md) - Core utilities overview

## Finding What You Need

### By Task
- **Setup/Installation**: [setup.md](setup.md)
- **Running analyses**: Domain-specific workflows (RNA, GWAS, single-cell)
- **Writing code**: Core utilities + module READMEs
- **Testing**: [testing.md](testing.md)
- **Understanding architecture**: [architecture.md](architecture.md)

### By Domain
Start with domain index files:
- DNA → [dna/index.md](dna/index.md) - See also: [DNA Analysis Scripts](../../scripts/README.md#dna-analysis-scripts-scriptsdna)
- RNA → [rna/index.md](rna/index.md) - See also: [RNA Workflows](../../scripts/rna/README.md) and [CLI Commands](./cli.md#rna-workflows)
- GWAS → [gwas/index.md](gwas/index.md) - See also: [GWAS Scripts](../../scripts/README.md#gwas-analysis-scripts-scriptsgwas) and [Config Templates](../../config/README.md#gwas-analysis-configurations)
- Single-cell → [singlecell/index.md](singlecell/index.md) - See also: [Single-Cell Scripts](../../scripts/README.md#single-cell-genomics-scripts-scriptssinglecell) and [Config Template](../../config/README.md#single-cell-analysis-singlecell_templateyaml)
- Math → [math/index.md](math/index.md) - See also: [Math Scripts](../../scripts/README.md#mathematical-biology-scripts-scriptsmath)
- ML → [ml/index.md](ml/index.md) - See also: [ML Scripts](../../scripts/README.md#machine-learning-scripts-scriptsml)
- Networks → [networks/index.md](networks/index.md) - See also: [Network Scripts](../../scripts/README.md#network-analysis-scripts-scriptsnetworks) and [Config Template](../../config/README.md#network-analysis-networkstemplateyaml)
- Multi-omics → [multiomics/index.md](multiomics/index.md) - See also: [Multi-Omics Scripts](../../scripts/README.md#multi-omics-scripts-scriptsmultiomics) and [Config Template](../../config/README.md#multi-omics-integration-multiomicstemplateyaml)
- Quality → [quality/index.md](quality/index.md) - See also: [Quality Scripts](../../scripts/README.md#quality-control-scripts-scriptsquality)
- Visualization → [visualization/index.md](visualization/index.md) - See also: [Visualization Scripts](../../scripts/README.md#visualization-scripts-scriptvisualization)

### By Component
- **Core utilities**: Start at [core/README.md](core/README.md)
- **Specific modules**: Navigate via [index.md](index.md)
- **Examples**: Check module READMEs for quick starts

## Documentation Maintenance

### Updating Documentation
1. Follow existing structure and formatting
2. Update cross-references when adding new modules
3. Maintain consistency across similar documentation
4. Update indices when adding new files

### Documentation Principles
- **Clear**: Technical but accessible writing
- **Complete**: Cover all functionality
- **Current**: Update with code changes
- **Connected**: Cross-reference related content
- **Practical**: Include runnable examples

## Support and Contribution

### Getting Help
1. Check domain-specific documentation
2. Review examples in module READMEs
3. Consult [testing.md](testing.md) for test examples
4. See [cli.md](cli.md) for command-line usage

### Contributing Documentation
1. Follow existing patterns and structure
2. Add examples for new functionality
3. Update indices and navigation
4. Maintain AGENTS.md for AI-assisted work
5. Cross-reference related modules

---

**This guide provides comprehensive navigation through the METAINFORMANT documentation system. Start with the essential documents above, then explore domain-specific content based on your needs.**

For the most current overview, always check [docs/index.md](index.md) and [docs/README.md](README.md).

