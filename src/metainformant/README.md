# METAINFORMANT Package Documentation

This directory contains the core METAINFORMANT package, a comprehensive bioinformatics and systems biology toolkit for integrated multi-omic analysis.

## Package Overview

METAINFORMANT provides a modular, extensible framework for biological data analysis across multiple domains, from molecular sequences to complex biological systems. The package is organized into specialized modules that can be used independently or composed together for comprehensive analyses.

## Core Modules

### 🧬 Core (`core/`)
Shared utilities and infrastructure used across all domains:
- **Configuration Management**: Environment-based config loading and validation
- **I/O Operations**: Robust file I/O with support for JSON, JSONL, CSV/TSV, and compressed formats
- **Logging**: Consistent, structured logging across all modules
- **Parallel Processing**: Thread-based parallel execution with order preservation
- **Caching**: JSON-based caching with TTL support
- **Path Management**: Path expansion, resolution, and containment validation
- **Database Integration**: Optional database client helpers

### 🧬 DNA (`dna/`)
Comprehensive DNA sequence analysis and manipulation:
- **Sequence I/O**: FASTA/FASTQ file reading and writing
- **Alignment**: Pairwise and multiple sequence alignment
- **Phylogenetics**: Neighbor-joining trees and Newick format support
- **Population Genetics**: Nucleotide diversity, Tajima's D, Fst statistics
- **Genomic Data**: NCBI datasets integration and accession validation
- **Motif Analysis**: Pattern discovery and restriction enzyme analysis
- **Variant Calling**: SNP and indel detection and analysis

### 🔬 RNA (`rna/`)
Transcriptomic analysis and pipeline orchestration:
- **Amalgkit Integration**: Modular wrapper around amalgkit CLI tools
- **Workflow Management**: Complete pipeline planning and execution
- **Metadata Handling**: Transcriptomic metadata retrieval and curation
- **Quantification**: Gene expression quantification and normalization

### 🧫 Protein (`protein/`)
Proteomic analysis and sequence manipulation:
- **Proteome Retrieval**: Database integration for protein sequence data
- **Sequence Analysis**: Protein sequence manipulation and analysis
- **Structure Integration**: Support for structural data and annotations

### 🧬 Epigenome (`epigenome/`)
Epigenetic modification analysis:
- **Methylation Analysis**: DNA methylation pattern detection and analysis
- **Chromatin Analysis**: Chromatin state and modification studies
- **Track Processing**: Genomic track file processing and visualization

### 📊 Ontology (`ontology/`)
Functional annotation and semantic analysis:
- **Gene Ontology**: GO term enrichment and annotation
- **Semantic Similarity**: Ontology-based similarity measures
- **Annotation Integration**: Multi-source annotation harmonization

### 🐜 Phenotype (`phenotype/`)
Phenotypic trait analysis and curation:
- **Morphological Traits**: Shape and structural phenotype analysis
- **Behavioral Traits**: Behavioral phenotype quantification
- **Web Scraping**: Automated phenotype data collection (e.g., AntWiki)

### 🌿 Ecology (`ecology/`)
Ecological metadata and community analysis:
- **Community Composition**: Species diversity and abundance analysis
- **Environmental Metadata**: Ecological parameter integration
- **Population Dynamics**: Community structure and interaction analysis

### 📈 Math (`math/`)
Theoretical and quantitative biology:
- **Population Genetics**: Mathematical models of evolutionary processes
- **Selection Theory**: Natural selection and quantitative trait analysis
- **Price Equation**: Evolutionary change decomposition
- **Kin Selection**: Hamilton's rule and inclusive fitness
- **Decision Theory**: Drift-diffusion models for behavioral analysis

### 📊 Visualization (`visualization/`)
Unified plotting and animation framework:
- **Statistical Plots**: Line plots, heatmaps, scatter plots, pair plots
- **Phylogenetic Trees**: Tree visualization and manipulation
- **Animations**: Time-series and dynamic data animation
- **Publication Graphics**: High-quality figure generation

### 🔄 Simulation (`simulation/`)
Synthetic data generation and agent-based modeling:
- **Sequence Simulation**: DNA, RNA, and protein sequence generation
- **Expression Simulation**: Gene expression count simulation
- **Agent-Based Models**: Grid-world and multi-agent simulations
- **Evolutionary Dynamics**: Population genetic simulations

### 🔬 Single Cell (`singlecell/`)
Single-cell transcriptomic analysis:
- **Preprocessing**: Quality control and normalization
- **Dimensionality Reduction**: PCA, t-SNE, UMAP integration
- **Clustering**: Cell type identification and marker analysis
- **Trajectory Analysis**: Pseudotime and developmental trajectory inference
- **Integration**: Multi-sample batch correction

### 🧬 Quality (`quality/`)
Data quality assessment and control:
- **Sequence Quality**: FASTQ quality metrics and filtering
- **Assembly Quality**: Genome assembly validation
- **Expression Quality**: Transcriptome quality assessment

### 🕸️ Networks (`networks/`)
Biological network analysis:
- **Protein Interaction**: PPI network construction and analysis
- **Regulatory Networks**: Gene regulatory network inference
- **Pathway Analysis**: Biological pathway enrichment and visualization
- **Community Detection**: Network module identification

### 🤖 Machine Learning (`ml/`)
Statistical and machine learning methods:
- **Classification**: Supervised learning for biological prediction
- **Regression**: Continuous trait prediction models
- **Feature Selection**: Dimensionality reduction and feature importance
- **Validation**: Cross-validation and model assessment
- **Dimensionality Reduction**: PCA, t-SNE, UMAP for biological data

### 🔗 Multi-omics (`multiomics/`)
Integrated multi-omic data analysis:
- **Data Integration**: Cross-platform data harmonization
- **Joint Analysis**: Multi-omic correlation and interaction analysis
- **Systems Biology**: Network-based integration approaches

## Architecture Principles

### Modularity
Each domain is self-contained with clear interfaces, allowing independent use and testing while enabling composition across domains.

### Extensibility
New organisms, assays, and analysis methods can be added incrementally without disrupting existing functionality.

### Performance
Efficient implementations with parallel processing capabilities and memory-conscious data structures.

### Reproducibility
Consistent random seeds, deterministic algorithms, and comprehensive configuration management.

## Usage Patterns

### Individual Domain Analysis
```python
from metainformant.dna import sequences, alignment, phylogeny

# Load and analyze DNA sequences
seqs = sequences.read_fasta("data/sequences.fasta")
aln = alignment.global_pairwise(seqs.values())
tree = phylogeny.neighbor_joining(aln)
```

### Cross-Domain Integration
```python
from metainformant.dna import population
from metainformant.rna import workflow
from metainformant.protein import proteomes

# Combine genetic and expression data
genetic_diversity = population.genetic_diversity(dna_sequences)
expression_patterns = workflow.extract_expression_patterns(rna_data)
protein_annotations = proteomes.annotate_sequences(protein_sequences)
```

### Pipeline Orchestration
```python
from metainformant.core import config, logging
from metainformant.multiomics import integration

# Configure and run integrated analysis
cfg = config.load_config("analysis_config.yaml")
logger = logging.setup_logger("multiomic_analysis")
results = integration.run_integrated_pipeline(cfg)
```

## Testing

The package includes comprehensive tests for all modules:
- Unit tests for individual functions
- Integration tests for cross-module functionality
- Performance benchmarks for critical operations
- Real-data validation tests

Run tests with:
```bash
pytest tests/
```

## Contributing

### Adding New Modules
1. Create new subdirectory in `src/metainformant/`
2. Implement core functionality in module files
3. Add comprehensive tests in `tests/`
4. Update documentation in `docs/`
5. Add module to main `__init__.py`

### Code Style
- Follow PEP 8 and type hints
- Use descriptive variable and function names
- Include docstrings for all public functions
- Maintain backward compatibility when possible

## Dependencies

Core dependencies are managed through `uv` and specified in `pyproject.toml`. Optional dependencies for specific functionality are imported defensively to avoid hard failures.

## License

Licensed under the Apache License, Version 2.0. See the main project LICENSE file for details.
