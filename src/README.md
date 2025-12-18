# Source Code Directory

This directory contains the METAINFORMANT source code organized into a modular, domain-specific structure.

## Structure Overview

```
src/
└── metainformant/           # Main package directory
    ├── __init__.py         # Package initialization and exports
    ├── __main__.py         # CLI entry point (`python -m metainformant`)
    ├── core/              # Shared utilities and infrastructure
    ├── dna/               # DNA sequence analysis and genomics
    ├── rna/               # RNA transcriptomic analysis
    ├── protein/           # Protein sequence and structure analysis
    ├── epigenome/         # Epigenetic modification analysis
    ├── ontology/          # Functional annotation and ontologies
    ├── phenotype/         # Phenotypic trait analysis
    ├── ecology/           # Ecological metadata and community analysis
    ├── math/              # Mathematical and theoretical biology
    ├── gwas/              # Genome-wide association studies
    ├── information/       # Information-theoretic analysis
    ├── life_events/       # Life course event analysis
    ├── visualization/     # Plotting and animation utilities
    ├── simulation/        # Synthetic data and agent-based modeling
    ├── singlecell/        # Single-cell transcriptomic analysis
    ├── quality/           # Data quality assessment and control
    ├── networks/          # Biological network analysis
    ├── ml/                # Machine learning for biological data
    └── multiomics/        # Multi-omic data integration
```

## Module Organization

### Core Infrastructure (`core/`)
Foundational utilities used across all domains:
- Configuration management and environment integration
- Robust I/O operations for multiple file formats
- Consistent logging and error handling
- Parallel processing and performance optimization
- Path handling and validation
- Caching mechanisms for expensive operations

### Domain Modules
Each domain represents a specific area of biological analysis:

#### 🧬 DNA (`dna/`)
Comprehensive DNA sequence analysis toolkit covering:
- Sequence I/O and manipulation (FASTA/FASTQ)
- Multiple sequence alignment and phylogenetics
- Population genetics and evolutionary analysis
- Genomic data retrieval from NCBI
- Restriction enzyme analysis and motif discovery

#### 🔬 RNA (`rna/`)
Complete transcriptomic analysis pipeline featuring:
- Integration with amalgkit CLI toolkit
- Workflow orchestration and pipeline management
- Metadata retrieval and curation
- Expression quantification and analysis

#### 🧫 Protein (`protein/`)
Protein sequence and structure analysis including:
- Proteome retrieval and annotation
- Structure integration and analysis
- Functional annotation and database integration

#### 🧬 Epigenome (`epigenome/`)
Epigenetic modification analysis covering:
- DNA methylation pattern detection
- Chromatin state analysis
- Track file processing and visualization

#### 📊 Ontology (`ontology/`)
Functional annotation using biological ontologies:
- Gene Ontology (GO) enrichment analysis
- Semantic similarity calculations
- OBO format parsing and processing

#### 🐜 Phenotype (`phenotype/`)
Phenotypic trait analysis and curation:
- Morphological and behavioral phenotype data
- Web scraping for phenotype databases (AntWiki)
- Trait quantification and analysis

#### 🌿 Ecology (`ecology/`)
Ecological metadata and community analysis:
- Species diversity and abundance metrics
- Environmental parameter integration
- Community structure analysis

#### 📈 Math (`math/`)
Theoretical and quantitative biology:
- Population genetics models (Price equation, kin selection)
- Decision theory (drift-diffusion models)
- Evolutionary simulation experiments
- Epidemiological modeling

#### 📊 Visualization (`visualization/`)
Unified plotting and animation framework:
- Statistical plots and heatmaps
- Phylogenetic tree visualization
- Animation and time-series plotting

#### 🔄 Simulation (`simulation/`)
Synthetic data generation and modeling:
- DNA/RNA/protein sequence generation
- Agent-based ecosystem modeling
- Expression count simulation

#### 🔬 Single Cell (`singlecell/`)
Single-cell transcriptomic analysis pipeline:
- Quality control and preprocessing
- Dimensionality reduction (PCA, t-SNE, UMAP)
- Clustering and trajectory analysis

#### 🧬 Quality (`quality/`)
Data quality assessment and validation:
- Sequence quality metrics
- Data integrity checking
- Format validation

#### 🕸️ Networks (`networks/`)
Biological network analysis:
- Protein-protein interaction networks
- Gene regulatory network inference
- Pathway enrichment analysis

#### 🤖 Machine Learning (`ml/`)
Statistical learning for biological data:
- Classification and regression
- Feature selection and importance
- Model validation and assessment

#### 🔗 Multi-omics (`multiomics/`)
Integrated multi-omic analysis:
- Cross-platform data harmonization
- Joint statistical analysis
- Systems biology integration

#### 🧬 GWAS (`gwas/`)
Genome-wide association studies:
- Variant calling and quality control
- Association testing and statistical analysis
- Population structure analysis
- Comprehensive visualization

#### 📊 Information Theory (`information/`)
Information-theoretic analysis:
- Shannon entropy and mutual information
- Semantic similarity calculations
- Information content analysis
- Network information measures

#### 📅 Life Events (`life_events/`)
Life course and event sequence analysis:
- Event sequence modeling
- Temporal pattern prediction
- Embedding learning
- Interpretability analysis

## Development Guidelines

### Code Organization
- **Modular Design**: Each domain is self-contained with clear interfaces
- **Consistent Patterns**: Similar structure and API patterns across modules
- **Defensive Imports**: Optional dependencies handled gracefully
- **Type Safety**: Comprehensive type hints throughout

### API Design Principles
- **Functional APIs**: Accept destination paths rather than writing to fixed locations
- **Composable Operations**: Functions designed to work together
- **Error Handling**: Clear error messages with biological context
- **Performance Awareness**: Efficient algorithms for large biological datasets

### Integration Patterns
```python
# Cross-module integration example
from metainformant.dna import sequences, phylogeny
from metainformant.visualization import trees

# Load and analyze DNA sequences
seqs = sequences.read_fasta("data/sequences.fasta")
tree = phylogeny.neighbor_joining_tree(seqs)
trees.visualize_tree(tree, "output/tree.png")
```

## Documentation

Each module includes:
- **README.md** with usage examples
- **API documentation** for all public functions
- **Integration examples** showing cross-module usage
- **Performance notes** and best practices

## Testing

The source code includes tests:
- Unit tests for individual functions
- Integration tests for cross-module functionality
- Performance benchmarks for scalability
- Real data validation tests

## Contributing

### Adding New Modules
1. Create new subdirectory in `src/metainformant/`
2. Follow established patterns for `__init__.py` and module structure
3. Add tests in `tests/` directory
4. Update main package `__init__.py` exports
5. Create detailed README.md documentation

### Code Standards
- Follow PEP 8 style guidelines
- Use descriptive variable and function names
- Include docstrings for all public APIs
- Maintain backward compatibility when possible
- Add type hints for better IDE support

This source code organization provides a foundation for biological data analysis while maintaining modularity and extensibility.
