# METAINFORMANT Package Documentation

This directory contains the core METAINFORMANT package, a bioinformatics and systems biology toolkit for integrated multi-omic analysis.

## Package Overview

METAINFORMANT provides a framework for biological data analysis across multiple domains, from molecular sequences to complex biological systems. The package is organized into specialized modules that can be used independently or composed together for analyses.

## Core Modules

### 🧬 Core (`core/`)
Shared utilities and infrastructure used across all domains:

**Configuration & Data Management:**
- **Configuration Management** (`config`): Environment-based config loading and validation
- **I/O Operations** (`io`): Robust file I/O with support for JSON, JSONL, CSV/TSV, and compressed formats
- **Caching** (`cache`): JSON-based caching with TTL support
- **Database Integration** (`db`): Optional database client helpers (PostgreSQL)
- **Disk Operations** (`disk`): Disk space and file system utilities

**Path & File Management:**
- **Path Management** (`paths`): Path expansion, resolution, and containment validation
- **Hashing** (`hash`): Content and file hashing functions

**Processing & Execution:**
- **Parallel Processing** (`parallel`): Thread-based parallel execution with order preservation
- **Workflow Utilities** (`workflow`): Config-based data processing and workflow orchestration functions
- **Progress Tracking** (`progress`): Progress tracking utilities

**Development & Discovery:**
- **Discovery** (`discovery`): Symbolic mapping and context discovery utilities
- **Symbol Management** (`symbols`): Symbol indexing and cross-referencing utilities

**Utilities:**
- **Logging** (`logging`): Consistent, structured logging across all modules
- **Text Processing** (`text`): Text normalization and manipulation utilities
- **Error Handling** (`errors`): Error handling and resilience utilities
- **Validation** (`validation`): Data validation utilities

See [`core/README.md`](core/README.md) for detailed documentation.

### 🧬 DNA (`dna/`)
Comprehensive DNA sequence analysis and manipulation:
- **Sequence I/O**: FASTA/FASTQ file reading and writing
- **Alignment**: Pairwise and multiple sequence alignment
- **Phylogenetics**: Neighbor-joining trees and Newick format support
- **Population Genetics**: Nucleotide diversity, Tajima's D, Fst statistics
- **Genomic Data**: NCBI datasets integration and accession validation
- **Motif Analysis**: Pattern discovery and restriction enzyme analysis
- **Variant Calling**: SNP and indel detection and analysis

See [`dna/README.md`](dna/README.md) for detailed documentation.

**Related Modules**: RNA (transcriptomic analysis), Protein (translation), GWAS (variant analysis), Math (population genetics)

### 🔬 RNA (`rna/`)
Transcriptomic analysis and pipeline orchestration:
- **Amalgkit Integration**: Modular wrapper around amalgkit CLI tools
- **Workflow Management**: Complete pipeline planning and execution
- **Metadata Handling**: Transcriptomic metadata retrieval and curation
- **Quantification**: Gene expression quantification and normalization

See [`rna/README.md`](rna/README.md) for detailed documentation.

**Related Modules**: DNA (genomic coordinates), Single-cell (scRNA-seq), Multi-omics (transcriptomic integration), Networks (regulatory networks)

### 🧫 Protein (`protein/`)
Proteomic analysis and sequence manipulation:
- **Proteome Retrieval**: Database integration for protein sequence data
- **Sequence Analysis**: Protein sequence manipulation and analysis
- **Structure Integration**: Support for structural data and annotations
- **Structure Prediction**: AlphaFold integration and model retrieval
- **Functional Annotation**: InterPro domain annotation and UniProt integration

See [`protein/README.md`](protein/README.md) for detailed documentation.

**Related Modules**: DNA (genomic coordinates), RNA (translation), Networks (PPI networks), Ontology (functional annotation)

### 🧬 Epigenome (`epigenome/`)
Epigenetic modification analysis:
- **Methylation Analysis**: DNA methylation pattern detection and analysis
- **Chromatin Analysis**: Chromatin state and modification studies
- **Track Processing**: Genomic track file processing and visualization

See [`epigenome/README.md`](epigenome/README.md) for detailed documentation.

**Related Modules**: DNA (genomic coordinates), Networks (chromatin interactions), Multi-omics (epigenomic integration)

### 📊 Ontology (`ontology/`)
Functional annotation and semantic analysis:
- **Gene Ontology**: GO term enrichment and annotation
- **Semantic Similarity**: Ontology-based similarity measures
- **Annotation Integration**: Multi-source annotation harmonization

See [`ontology/README.md`](ontology/README.md) for detailed documentation.

**Related Modules**: Protein (functional annotation), Networks (pathway enrichment), Information Theory (semantic similarity)

### 🐜 Phenotype (`phenotype/`)
Phenotypic trait analysis and curation:
- **Morphological Traits**: Shape and structural phenotype analysis
- **Behavioral Traits**: Behavioral phenotype quantification
- **Web Scraping**: Automated phenotype data collection (e.g., AntWiki)

See [`phenotype/README.md`](phenotype/README.md) for detailed documentation.

**Related Modules**: Life Events (life course analysis), GWAS (trait association), Networks (trait networks)

### 🌿 Ecology (`ecology/`)
Ecological metadata and community analysis:
- **Community Composition**: Species diversity and abundance analysis
- **Environmental Metadata**: Ecological parameter integration
- **Population Dynamics**: Community structure and interaction analysis

See [`ecology/README.md`](ecology/README.md) for detailed documentation.

**Related Modules**: Networks (ecological interactions), DNA (population genetics), Information Theory (diversity measures)

### 📈 Math (`math/`)
Theoretical and quantitative biology:
- **Population Genetics**: Mathematical models of evolutionary processes
- **Selection Theory**: Natural selection and quantitative trait analysis
- **Price Equation**: Evolutionary change decomposition
- **Kin Selection**: Hamilton's rule and inclusive fitness
- **Decision Theory**: Drift-diffusion models for behavioral analysis

See [`math/README.md`](math/README.md) for detailed documentation.

**Related Modules**: DNA (population genetics), GWAS (statistical models), Information Theory (theoretical foundations)

### 📊 Visualization (`visualization/`)
Unified plotting and animation framework:
- **Statistical Plots**: Line plots, heatmaps, scatter plots, pair plots
- **Phylogenetic Trees**: Tree visualization and manipulation
- **Animations**: Time-series and dynamic data animation
- **Publication Graphics**: High-quality figure generation

See [`visualization/README.md`](visualization/README.md) for detailed documentation.

**Related Modules**: All modules (visualization support), GWAS (statistical plots), Networks (graph visualization)

### 🔄 Simulation (`simulation/`)
Synthetic data generation and agent-based modeling:
- **Sequence Simulation**: DNA, RNA, and protein sequence generation
- **Expression Simulation**: Gene expression count simulation
- **Agent-Based Models**: Grid-world and multi-agent simulations
- **Evolutionary Dynamics**: Population genetic simulations
- **Module-Specific Scripts**: Dedicated simulation scripts for each domain module (DNA, RNA, protein, GWAS, single-cell, etc.)

See [`simulation/README.md`](simulation/README.md) for detailed documentation. Module-specific simulation scripts are documented in `scripts/simulation/README.md`.

**Related Modules**: All modules (synthetic data generation for testing and validation)

### 🔬 Single Cell (`singlecell/`)
Single-cell transcriptomic analysis:
- **Preprocessing**: Quality control and normalization
- **Dimensionality Reduction**: PCA, t-SNE, UMAP integration
- **Clustering**: Cell type identification and marker analysis
- **Trajectory Analysis**: Pseudotime and developmental trajectory inference
- **Integration**: Multi-sample batch correction

See [`singlecell/README.md`](singlecell/README.md) for detailed documentation.

**Related Modules**: RNA (transcriptomic analysis), Multi-omics (quality control), ML (dimensionality reduction)

### 🧬 Quality (`quality/`)
Data quality assessment and control:
- **Sequence Quality**: FASTQ quality metrics and filtering
- **Assembly Quality**: Genome assembly validation
- **Expression Quality**: Transcriptome quality assessment

See [`quality/README.md`](quality/README.md) for detailed documentation.

**Related Modules**: All modules (quality assessment), DNA (sequence quality), RNA (expression quality)

### 🕸️ Networks (`networks/`)
Biological network analysis:
- **Protein Interaction**: PPI network construction and analysis
- **Regulatory Networks**: Gene regulatory network inference
- **Pathway Analysis**: Biological pathway enrichment and visualization
- **Community Detection**: Network module identification

See [`networks/README.md`](networks/README.md) for detailed documentation.

**Related Modules**: Protein (PPI networks), RNA (regulatory networks), Ontology (pathway networks), Information Theory (network entropy)

### 🤖 Machine Learning (`ml/`)
Statistical and machine learning methods:
- **Classification**: Supervised learning for biological prediction
- **Regression**: Continuous trait prediction models
- **Feature Selection**: Dimensionality reduction and feature importance
- **Validation**: Cross-validation and model assessment
- **Dimensionality Reduction**: PCA, t-SNE, UMAP for biological data

See [`ml/README.md`](ml/README.md) for detailed documentation.

**Related Modules**: All modules (feature engineering), Single-cell (dimensionality reduction), Networks (network-based features)

### 🔗 Multi-omics (`multiomics/`)
Integrated multi-omic data analysis:
- **Data Integration**: Cross-platform data harmonization
- **Joint Analysis**: Multi-omic correlation and interaction analysis
- **Systems Biology**: Network-based integration approaches

See [`multiomics/README.md`](multiomics/README.md) for detailed documentation.

**Related Modules**: DNA, RNA, Protein, Epigenome (data integration), Networks (systems biology), Visualization (integrated plots)

### 📅 Life Events (`life_events/`)
Life course and event sequence analysis:
- **Event Sequences**: Temporal event representation with domain categorization
- **Event Embeddings**: Word2Vec-style embeddings for life events
- **Sequence Prediction**: Models for predicting outcomes from event sequences
- **Workflow Integration**: End-to-end pipelines for life course analysis
- **Population Comparison**: Compare event patterns across groups

See [`life_events/README.md`](life_events/README.md) for detailed documentation.

**Related Modules**: Phenotype (life course analysis), ML (sequence prediction), Information Theory (temporal patterns)

### 📊 Information Theory (`information/`)
Information-theoretic analysis for biological data:
- **Syntactic Information**: Shannon entropy, mutual information, KL divergence
- **Semantic Information**: Information content, semantic similarity measures
- **Continuous Data**: Differential entropy and continuous information measures
- **Estimation Methods**: Bias-corrected entropy and MI estimation
- **Integration**: Cross-module integration with DNA, RNA, single-cell, and multi-omics
- **Workflows**: Batch processing and information analysis pipelines

See [`information/README.md`](information/README.md) and [`information/EXAMPLES.md`](information/EXAMPLES.md) for detailed documentation.

**Related Modules**: DNA (variant analysis), Networks (information flow), Ontology (semantic information), Multi-omics (information integration)

### 🧬 GWAS (`gwas/`)
Genome-Wide Association Studies workflow:
- **Variant Calling**: Integration with bcftools and GATK
- **Quality Control**: MAF filtering, missingness, HWE testing
- **Population Structure**: PCA and kinship matrix computation
- **Association Testing**: Linear and logistic regression models
- **Multiple Testing**: Bonferroni, FDR, and genomic control correction
- **Visualization**: Manhattan plots, Q-Q plots, regional plots, and visualizations
- **Data Acquisition**: SRA download and reference genome retrieval

See [`gwas/README.md`](gwas/README.md) for detailed documentation.

**Related Modules**: DNA (variant calling), Math (statistical models), Phenotype (trait association), Visualization (GWAS plots), Quality (QC filters)

## Architecture Principles

### Modularity
Each domain is self-contained with clear interfaces, allowing independent use and testing while enabling composition across domains.

### Extensibility
New organisms, assays, and analysis methods can be added incrementally without disrupting existing functionality.

### Performance
Efficient implementations with parallel processing capabilities and memory-conscious data structures.

### Reproducibility
Consistent random seeds, deterministic algorithms, and configuration management.

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

**Related Modules**: DNA, RNA, and Protein modules work together for multi-omic analysis. See Multi-omics module for integrated workflows.

### Multi-Omic Integration
```python
from metainformant.dna import variants
from metainformant.rna import workflow
from metainformant.epigenome import methylation
from metainformant.multiomics import integration

# Integrate genomic, transcriptomic, and epigenomic data
variants_data = variants.read_vcf("data/variants.vcf")
expression_data = workflow.load_expression("data/expression.tsv")
methylation_data = methylation.load_methylation("data/methylation.bed")

# Joint analysis
integrated_results = integration.joint_analysis(
    genomics=variants_data,
    transcriptomics=expression_data,
    epigenomics=methylation_data
)
```

**Related Modules**: DNA (variants), RNA (expression), Epigenome (methylation), Multi-omics (integration)

### Network and Systems Biology
```python
from metainformant.networks import ppi, regulatory
from metainformant.ontology import go
from metainformant.information import mutual_information

# Build regulatory network from expression data
reg_network = regulatory.infer_regulatory_network(expression_data)

# Annotate with GO terms
annotated_network = go.annotate_network(reg_network)

# Analyze information flow
info_flow = mutual_information(expression_data, network_structure)
```

**Related Modules**: Networks (PPI, regulatory), Ontology (GO annotation), Information Theory (network analysis)

### Pipeline Orchestration
```python
from metainformant.core import config, logging
from metainformant.multiomics import integration

# Configure and run integrated analysis
cfg = config.load_config("analysis_config.yaml")
logger = logging.setup_logger("multiomic_analysis")
results = integration.run_integrated_pipeline(cfg)
```

### Life Course Analysis
```python
from metainformant.life_events import EventSequence, Event, analyze_life_course
from datetime import datetime

# Create event sequences
events = [
    Event("degree", datetime(2010, 6, 1), "education"),
    Event("job_change", datetime(2015, 3, 1), "occupation"),
]
sequence = EventSequence(person_id="person_001", events=events)

# Analyze life course
results = analyze_life_course([sequence], outcomes=None)
```

## Testing

The package includes tests for all modules:
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
3. Add tests in `tests/`
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
