# METAINFORMANT

A comprehensive bioinformatics and systems biology toolkit for integrated multi-omic analysis, developed with AI assistance for enhanced code quality and documentation.

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## Overview

METAINFORMANT is a modular, extensible framework for biological data analysis spanning multiple domains: genomics, transcriptomics, proteomics, epigenomics, and systems biology. Built with modern Python practices and comprehensive documentation, it provides professional-grade tools for bioinformatics research.

### Key Features

- 🧬 **Comprehensive Multi-Omic Analysis**: DNA, RNA, protein, and epigenome analysis
- 📊 **Statistical & ML Methods**: GWAS, population genetics, machine learning pipelines
- 🔬 **Single-Cell Genomics**: Complete scRNA-seq analysis workflows
- 🌐 **Network Analysis**: Biological networks, pathways, and community detection
- 📈 **Publication-Quality Visualizations**: High-resolution plots and interactive animations
- 🔧 **Modular Architecture**: Use individual modules or compose complete workflows
- 📝 **Extensive Documentation**: 70+ README files, comprehensive guides
- ✅ **Real Implementation Testing**: No mocks - all tests use real methods

## Quick Start

### Installation

```bash
# Clone repository
git clone https://github.com/q/MetaInformAnt.git
cd MetaInformAnt

# Setup environment with uv
bash scripts/package/setup_uv.sh

# Or manual setup
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

### Quick Example

```python
from metainformant.dna import sequences, composition
from metainformant.visualization import lineplot

# Load DNA sequences
seqs = sequences.read_fasta("data/sequences.fasta")

# Analyze GC content
gc_values = [composition.gc_content(seq) for seq in seqs.values()]

# Visualize
ax = lineplot(None, gc_values)
ax.set_ylabel("GC Content")
ax.set_title("GC Content Across Sequences")
ax.figure.savefig("output/gc_content.png", dpi=300)
```

### Complete Workflow Demonstration

```bash
# Run comprehensive workflow demo
python3 scripts/core/run_demo.py

# Demonstrates:
# - Configuration management
# - Data processing
# - Visualization with informative names
# - Complete output organization
```

See `scripts/core/run_demo.py` for the workflow demonstration. Outputs are saved to `output/demo/` with workflow configuration, processed data, visualizations, and summary reports.

## Module Overview

### Core Infrastructure

- **[core/](src/metainformant/core/)** - Shared utilities (I/O, logging, configuration, parallel processing)

### Molecular Analysis

- **[dna/](src/metainformant/dna/)** - DNA sequences, alignment, phylogenetics, population genetics
- **[rna/](src/metainformant/rna/)** - RNA-seq workflows, amalgkit integration, transcriptomics
- **[protein/](src/metainformant/protein/)** - Protein sequences, structure, AlphaFold, proteomics
- **[epigenome/](src/metainformant/epigenome/)** - Methylation analysis, chromatin tracks

### Statistical & ML Methods

- **[gwas/](src/metainformant/gwas/)** - Genome-wide association studies, variant calling, visualization
- **[math/](src/metainformant/math/)** - Mathematical biology, population genetics theory, dynamics
- **[ml/](src/metainformant/ml/)** - Machine learning pipelines, classification, regression
- **[information/](src/metainformant/information/)** - Information theory methods (Shannon entropy, mutual information, semantic similarity)

### Systems Biology

- **[networks/](src/metainformant/networks/)** - Biological networks, community detection, pathways
- **[multiomics/](src/metainformant/multiomics/)** - Multi-omic data integration
- **[singlecell/](src/metainformant/singlecell/)** - Single-cell RNA-seq analysis
- **[simulation/](src/metainformant/simulation/)** - Synthetic data generation, agent-based models

### Annotation & Metadata

- **[ontology/](src/metainformant/ontology/)** - Gene Ontology, functional annotation
- **[phenotype/](src/metainformant/phenotype/)** - Phenotypic data curation
- **[ecology/](src/metainformant/ecology/)** - Ecological metadata, community analysis
- **[life_events/](src/metainformant/life_events/)** - Life course and event sequence analysis, temporal pattern prediction

### Utilities

- **[quality/](src/metainformant/quality/)** - Quality control and assessment
- **[visualization/](src/metainformant/visualization/)** - Plots, heatmaps, animations, trees

## Documentation

### Quick Links

- **[Documentation Guide](docs/DOCUMENTATION_GUIDE.md)** - Complete navigation guide
- **[Quick Start](QUICKSTART.md)** - Fast setup commands
- **[Architecture](docs/architecture.md)** - System design
- **[Testing Guide](docs/testing.md)** - Comprehensive testing documentation
- **[CLI Reference](docs/cli.md)** - Command-line interface

### Module Documentation

Each module has comprehensive documentation in `src/metainformant/<module>/README.md` and `docs/<module>/`.

## Scripts & Workflows

The [`scripts/`](scripts/) directory contains production-ready workflow orchestrators:

- **Package Management**: Setup, testing, quality control
- **RNA-seq**: Multi-species workflows, amalgkit integration
- **GWAS**: Genome-scale association studies
- **Module Orchestrators**: ✅ Complete workflow scripts for all domains (DNA, protein, networks, multiomics, single-cell, quality, simulation, visualization, epigenome, ecology, ontology, phenotype, ML, math)

See [`scripts/README.md`](scripts/README.md) for complete documentation.

### CLI Interface

All modules are accessible via the unified CLI:

```bash
# Setup and environment
uv run metainformant setup --with-amalgkit

# Domain workflows
uv run metainformant dna fetch --assembly GCF_000001405.40
uv run metainformant rna run --work-dir output/rna --threads 8
uv run metainformant gwas run --config config/gwas/gwas_template.yaml

# New module workflows
uv run metainformant networks run --input data/interactions.tsv --output output/networks
uv run metainformant multiomics run --genomics data/genomics.tsv --output output/multiomics
uv run metainformant singlecell run --input data/counts.h5ad --output output/singlecell --qc
uv run metainformant quality run --fastq data/reads.fq --output output/quality
uv run metainformant ontology run --go data/go.obo --output output/ontology
uv run metainformant ml run --features data/features.csv --output output/ml --classify

# See all available commands
uv run metainformant --help
```

See [`docs/cli.md`](docs/cli.md) for complete CLI documentation.

## Usage Examples

### DNA Analysis

```python
from metainformant.dna import alignment, population

# Pairwise alignment
align_result = alignment.global_align("ACGTACGT", "ACGTAGGT")
print(f"Score: {align_result.score}")

# Population genetics
sequences = ["ATCGATCG", "ATCGTTCG", "ATCGATCG"]
diversity = population.nucleotide_diversity(sequences)
print(f"π = {diversity:.4f}")
```

### RNA-seq Workflow

```bash
# Multi-species RNA-seq pipeline
python3 scripts/rna/run_multi_species.py

# Single species
bash scripts/rna/amalgkit/run_amalgkit.sh --config config/amalgkit/config.yaml
```

### GWAS Analysis

```python
from metainformant.gwas import run_gwas, manhattan_plot, load_gwas_config

# Load configuration and run workflow
config = load_gwas_config("config/gwas/gwas_template.yaml")
results = run_gwas(
    vcf_path="data/variants/cohort.vcf.gz",
    phenotype_path="data/phenotypes/traits.tsv",
    config={"association": {"model": "linear"}},
    output_dir="output/gwas"
)

# Visualize results
manhattan_plot(results["association_results"], output_path="output/gwas/manhattan.png")
```

### Visualization

```python
from metainformant.visualization import heatmap, animate_time_series

# Heatmap
heatmap(correlation_matrix, cmap="viridis", annot=True)

# Animation
fig, anim = animate_time_series(time_series_data)
anim.save("output/animation.gif")
```

### Network Analysis

```python
from metainformant.networks import create_network, detect_communities, centrality_measures

# Create network from interactions
network = create_network(edges, directed=False)

# Detect communities
communities = detect_communities(network)

# Calculate centrality
centrality = centrality_measures(network)
```

### Multi-Omics Integration

```python
from metainformant.multiomics import integrate_omics_data, joint_pca

# Integrate multiple omics datasets
multiomics = integrate_omics_data(
    genomics=genomics_data,
    transcriptomics=rna_data,
    proteomics=protein_data
)

# Joint dimensionality reduction
pca_result = joint_pca(multiomics)
```

### Information Theory

```python
from metainformant.information import shannon_entropy, mutual_information, information_content

# Calculate Shannon entropy
probs = [0.5, 0.3, 0.2]
entropy = shannon_entropy(probs)

# Mutual information between sequences
mi = mutual_information(sequence_x, sequence_y)

# Information content for hierarchical terms
ic = information_content(term_frequencies, "GO:0008150")
```

### Life Events Analysis

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

## Development

### Running Tests

```bash
# All tests
bash scripts/package/run_tests.sh

# Fast tests only
bash scripts/package/run_tests.sh --fast

# Specific module
pytest tests/dna/ -v
```

### Code Quality

```bash
# Check code quality
bash scripts/package/uv_quality.sh

# Run linting
ruff check src/

# Type checking
mypy src/metainformant
```

## Project Structure

```
MetaInformAnt/
├── src/metainformant/       # Main package
│   ├── core/               # Core utilities
│   ├── dna/                # DNA analysis
│   ├── rna/                # RNA analysis
│   ├── protein/            # Protein analysis
│   ├── gwas/               # GWAS analysis
│   └── ...                 # Additional modules
├── scripts/                # Workflow scripts
│   ├── package/            # Package management
│   ├── rna/                # RNA workflows
│   ├── gwas/               # GWAS workflows
│   └── ...                 # Module scripts
├── docs/                   # Documentation
├── tests/                  # Test suite
├── config/                 # Configuration files
├── output/                 # Analysis outputs
└── data/                   # Input data
```

## AI-Assisted Development

This project was developed with AI assistance (grok-code-fast-1 via Cursor) to enhance:

- Code generation and algorithm implementation
- Comprehensive documentation
- Test case generation
- Architecture design

All AI-generated content undergoes human review. See [`AGENTS.md`](AGENTS.md) for details.

## Known Limitations

### Module Completeness

Some modules have partial implementations or optional dependencies:

- **Machine Learning**: Framework exists; some methods may need completion (see [ML Documentation](docs/ml/README.md))
- **Multi-omics**: Core integration methods implemented; advanced features may require additional dependencies
- **Single-cell**: Requires `scipy`, `scanpy`, `anndata` for full functionality (see [Single-Cell Documentation](docs/singlecell/README.md))
- **Network Analysis**: Basic algorithms implemented; advanced regulatory network features may need enhancement

### GWAS Module

- **Variant Download**: Database download (dbSNP, 1000 Genomes) is a placeholder; use SRA-based workflow or provide VCF files
- **Functional Annotation**: Requires external tools (ANNOVAR, VEP, SnpEff) for variant annotation
- **Advanced Mixed Models**: Basic relatedness adjustment implemented; advanced MLM methods may require GCTA/EMMAX integration

### Test Coverage

Some modules have lower test success rates due to optional dependencies:
- **Single-cell**: Requires scientific dependencies (`scanpy`, `anndata`)
- **Multi-omics**: Framework exists, tests may skip without dependencies
- **Network Analysis**: Core tests pass; advanced features may need additional setup

See [Test Analysis](docs/test_analysis.md) for detailed coverage information.

## Best Practices

### File Naming

- ✅ Use informative names: `sample_pca_biplot_colored_by_treatment.png`
- ❌ Avoid generic names: `plot1.png`, `output.png`

### Output Organization

- All outputs in `output/` directory
- Configuration saved with results
- Visualizations in subdirectories with metadata

### No Mocking Policy

- All tests use real implementations
- No fake/mocked/stubbed methods
- Real API calls or graceful skips
- Ensures actual functionality

## Requirements

- Python 3.11+
- Optional: SRA Toolkit, kallisto (for RNA workflows)
- Optional: samtools, bcftools, bwa (for GWAS)

## Contributing

Contributions are welcome! Please:

1. Follow the existing code style
2. Add tests for new features
3. Update documentation
4. Use informative commit messages

## Citation

If you use METAINFORMANT in your research, please cite this repository:

```bibtex
@software{metainformant2025,
  author = {MetaInformAnt Development Team},
  title = {MetaInformAnt: Comprehensive Bioinformatics Toolkit},
  year = {2025},
  url = {https://github.com/q/MetaInformAnt},
  version = {0.2.0}
}
```

## License

This project is licensed under the Apache License, Version 2.0 - see [LICENSE](LICENSE) for details.

## Contact

- **Repository**: https://github.com/q/MetaInformAnt
- **Issues**: https://github.com/q/MetaInformAnt/issues
- **Documentation**: https://github.com/q/MetaInformAnt/blob/main/docs/

## Acknowledgments

- Developed with AI assistance from Cursor's Code Assistant (grok-code-fast-1)
- Built on established bioinformatics tools and libraries
- Community contributions and feedback

---

**Status**: Active Development | **Version**: 0.2.0 | **Python**: 3.11+ | **License**: Apache 2.0


