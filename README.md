# 🐜 METAINFORMANT

**Comprehensive bioinformatics toolkit for multi-omic analysis**

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Modules](https://img.shields.io/badge/modules-26-green.svg)](src/metainformant/)
[![Files](https://img.shields.io/badge/files-560-brightgreen.svg)](src/metainformant/)

---

## Overview

METAINFORMANT provides production-ready bioinformatics analysis across genomics, transcriptomics, proteomics, epigenomics, and systems biology. Built with Python 3.11+ and [`uv`](https://astral.sh/uv) for fast dependency management.

### 📊 At a Glance

| Metric | Value |
|--------|-------|
| **Modules** | 26 specialized analysis modules |
| **Python Files** | 560 implementation files |
| **Plot Types** | 70+ visualization methods |
| **Documentation** | 240+ README files |

### 🔬 Core Capabilities

| Domain | Features |
|--------|----------|
| **DNA** | Sequences, alignment, phylogenetics, population genetics, variant analysis |
| **RNA** | Amalgkit integration, SRA downloads, Kallisto quantification, industrial-scale pipelines |
| **GWAS** | Association testing, fine-mapping, visualization, complete GWAS pipelines |
| **eQTL** | Integration of GWAS variants and Amalgkit RNA-seq expression data |
| **Multi-omics** | Cross-omic integration, joint PCA, correlation analysis |
| **ML** | Classification, regression, feature selection, LLM integration |
| **Visualization** | Manhattan plots, heatmaps, networks, animations, publication-ready output |

### System Architecture

```mermaid
flowchart TB
    subgraph coreInfra["Core Infrastructure"]
        CORE["Core Utilities"]
    end

    subgraph molecular["Molecular Analysis"]
        DNA["DNA Analysis"]
        RNA["RNA Analysis"]
        PROT["Protein Analysis"]
        EPI["Epigenome Analysis"]
    end

    subgraph statsML["Statistical and ML"]
        GWAS["GWAS Analysis"]
        MATH["Mathematical Biology"]
        ML["Machine Learning"]
        INFO["Information Theory"]
    end

    subgraph systems["Systems Biology"]
        NET["Network Analysis"]
        MULTI["Multi-Omics Integration"]
        SC["Single-Cell Analysis"]
        SIM["Simulation"]
    end

    subgraph annotation["Annotation and Metadata"]
        ONT["Ontology"]
        PHEN["Phenotype Analysis"]
        ECO["Ecology"]
        LE["Life Events"]
    end

    subgraph utilities["Utilities"]
        QUAL["Quality Control"]
        VIZ["Visualization"]
    end

    subgraph specialized["Specialized Domains"]
        LR["Long-Read Sequencing"]
        METAG["Metagenomics"]
        SV["Structural Variants"]
        SPATIAL["Spatial Transcriptomics"]
        PHARMA["Pharmacogenomics"]
        METAB["Metabolomics"]
        MENU["Menu System"]
    end

    CORE --> DNA
    CORE --> RNA
    CORE --> PROT
    CORE --> EPI
    CORE --> GWAS
    CORE --> MATH
    CORE --> ML
    CORE --> INFO
    CORE --> NET
    CORE --> MULTI
    CORE --> SC
    CORE --> SIM
    CORE --> ONT
    CORE --> PHEN
    CORE --> ECO
    CORE --> LE
    CORE --> QUAL
    CORE --> VIZ
    CORE --> LR
    CORE --> METAG
    CORE --> SV
    CORE --> SPATIAL
    CORE --> PHARMA
    CORE --> METAB
    CORE --> MENU
```

### Data Flow and Integration Architecture

```mermaid
graph TD
    A["Raw Biological Data"] --> B["Data Ingestion"]
    B --> C{Data Type}

    C -->|DNA| D["DNA Module"]
    C -->|RNA| E["RNA Module"]
    C -->|Protein| F["Protein Module"]
    C -->|Epigenome| G["Epigenome Module"]
    C -->|Phenotype| H["Phenotype Module"]
    C -->|Environmental| I["Ecology Module"]

    D --> J["Quality Control"]
    E --> J
    F --> J
    G --> J
    H --> J
    I --> J

    J --> K["Core Processing"]
    K --> L{Analysis Type}

    L -->|Statistical| M["GWAS Module"]
    L -->|ML| N["ML Module"]
    L -->|Information| O["Information Module"]
    L -->|Networks| P["Networks Module"]
    L -->|Systems| Q["Multi-Omics Module"]
    L -->|Singlecell| R["Single-Cell Module"]
    L -->|Simulation| S["Simulation Module"]

    M --> T["Results Integration"]
    N --> T
    O --> T
    P --> T
    Q --> T
    R --> T
    S --> T

    T --> U["Visualization"]
    U --> V["Publication Figures"]
    V --> W["Scientific Insights"]

    subgraph "Primary Data Types"
        X["Genomic"] -.-> D
        Y["Transcriptomic"] -.-> E
        Z["Proteomic"] -.-> F
        AA["Epigenetic"] -.-> G
    end

    subgraph "Analysis Workflows"
        BB["Population Genetics"] -.-> M
        CC["Feature Selection"] -.-> N
        DD["Mutual Information"] -.-> O
        EE["Community Detection"] -.-> P
        FF["Joint PCA"] -.-> Q
        GG["Trajectory Analysis"] -.-> R
    end

    subgraph "Output Formats"
        HH["Manhattan Plots"] -.-> V
        II["Heatmaps"] -.-> V
        JJ["Network Graphs"] -.-> V
        KK["Animations"] -.-> V
    end
```

### Multi-Omic Integration Pipeline

```mermaid
graph TD
    A["Multi-Omic Datasets"] --> B["Sample Alignment"]
    B --> C["Batch Effect Correction"]

    C --> D{Integration Strategy}
    D -->|Early| E["Concatenated Matrix"]
    D -->|Late| F["Separate Models"]
    D -->|Intermediate| G["Meta-Analysis"]

    E --> H["Joint Dimensionality Reduction"]
    F --> I["Individual Analysis"]
    G --> J["Result Integration"]

    H --> K["Unified Clustering"]
    I --> L["Individual Clustering"]
    J --> M["Consensus Clustering"]

    K --> N["Functional Enrichment"]
    L --> N
    M --> N

    N --> O["Pathway Analysis"]
    O --> P["Network Construction"]

    P --> Q["Biological Interpretation"]
    Q --> R["Systems Biology Insights"]

    subgraph "Omic Layers"
        S["Genomics"] -.-> A
        T["Transcriptomics"] -.-> A
        U["Proteomics"] -.-> A
        V["Metabolomics"] -.-> A
        W["Epigenomics"] -.-> A
    end

    subgraph "Integration Methods"
        X["MOFA"] -.-> H
        Y["Joint PCA"] -.-> H
        Z["Similarity Networks"] -.-> H
    end

    subgraph "Biological Outputs"
        AA["Gene Modules"] -.-> Q
        BB["Regulatory Networks"] -.-> Q
        CC["Disease Pathways"] -.-> Q
        DD["Biomarkers"] -.-> Q
    end
```

### Quality Assurance Framework

```mermaid
graph TD
    A["Data Processing Pipeline"] --> B["Input Validation"]
    B --> C["Type Checking"]
    C --> D["Schema Validation"]

    D --> E["Processing Logic"]
    E --> F["Error Handling"]
    F --> G["Recovery Mechanisms"]

    G --> H["Output Validation"]
    H --> I["Result Verification"]
    I --> J["Quality Metrics"]

    J --> K{Acceptable Quality?}
    K -->|Yes| L["Pipeline Success"]
    K -->|No| M["Quality Issues"]

    M --> N["Diagnostic Analysis"]
    N --> O["Error Classification"]

    O --> P{Recoverable?}
    P -->|Yes| Q["Data Correction"]
    P -->|No| R["Pipeline Failure"]

    Q --> E
    L --> S["Validated Results"]
    R --> T["Error Reporting"]

    subgraph "Validation Layers"
        U["Data Integrity"] -.-> B
        V["Business Logic"] -.-> E
        W["Statistical Validity"] -.-> H
    end

    subgraph "Quality Controls"
        X["Unit Tests"] -.-> F
        Y["Integration Tests"] -.-> I
        Z["Performance Benchmarks"] -.-> J
    end

    subgraph "Error Types"
        AA["Data Errors"] -.-> O
        BB["Logic Errors"] -.-> O
        CC["System Errors"] -.-> O
        DD["External Errors"] -.-> O
    end
```

### Key Features

- **Multi-Omic Analysis**: DNA, RNA, protein, and epigenome data integration
- **Statistical & ML Methods**: GWAS, population genetics, machine learning pipelines
- **Single-Cell Genomics**: Complete scRNA-seq analysis workflows
- **Network Analysis**: Biological networks, pathways, community detection algorithms
- **Visualization Suite**: 14 specialized plotting modules with 70+ plot types and publication-quality output
- **Modular Architecture**: Individual modules or complete end-to-end workflows
- **Comprehensive Documentation**: 70+ README files with technical specifications
- **Implementation Testing**: Real methods in tests, no mocks or stubs
- **Quality Assurance**: Rigorous validation and error handling throughout
- **Performance Optimization**: Efficient algorithms for large-scale biological data

## Quick Start

### Prerequisites

- **Python 3.11+**
- **`uv`** - Fast Python package manager (**REQUIRED**)
  - Install: `curl -LsSf https://astral.sh/uv/install.sh | sh`
  - Verify: `uv --version`

### Installation

**METAINFORMANT uses `uv` for all package management. Never use `pip` directly.**

```bash
# Clone repository
git clone https://github.com/docxology/metainformant.git
cd metainformant

# Automated setup with uv (recommended - handles FAT filesystems automatically)
bash scripts/package/setup.sh

# Or manual setup with uv
curl -LsSf https://astral.sh/uv/install.sh | sh  # Install uv if needed
uv venv
source .venv/bin/activate  # or /tmp/metainformant_venv/bin/activate on FAT filesystems
uv pip install -e .
```

**Package Management**: All Python dependencies are managed via `uv`:

- Create venv: `uv venv`
- Install packages: `uv pip install -e .`
- Run commands: `uv run pytest`, `uv run metainformant --help`
- Sync dependencies: `uv sync --extra dev --extra scientific`
- Add dependencies: `uv add <package>`
- Remove dependencies: `uv remove <package>`

**Note**: Setup scripts automatically detect FAT filesystems (exFAT, FAT32) and configure UV cache and virtual environment locations accordingly. See [UV Setup Guide](docs/UV_SETUP.md) for details.

### Quick Example

```python
from metainformant.dna import sequences, composition
from metainformant.visualization import lineplot

# Load DNA sequences
seqs = sequences.read_fasta("data/sequences.fasta")

# Analyze GC content
gc_values = [sequences.gc_content(seq) for seq in seqs.values()]

# Visualize
ax = lineplot(None, gc_values)
ax.set_ylabel("GC Content")
ax.set_title("GC Content Across Sequences")
ax.figure.savefig("output/gc_content.png", dpi=300)
```

### Complete Workflow Demonstration

```bash
# Run workflow demo
python3 scripts/core/run_demo.py

# Demonstrates:
# - Configuration management and I/O operations
# - DNA sequence analysis and visualization
# - Quality control and metrics calculation
# - Real data processing with informative output names
# - Complete output organization in output/demo/ directory
```

See `scripts/core/run_demo.py` for the workflow demonstration. Outputs are saved to `output/demo/` directory with:

- Workflow configuration files
- Processed biological data (FASTA sequences, analysis results)
- Publication-quality visualizations with informative naming
- Summary reports and metadata

## Module Status Overview

### ✅ **Production-Ready Modules**

| Category | Module | Status | Key Features |
|----------|--------|--------|--------------|
| **Core** | [core/](src/metainformant/core/) | ✅ **Complete** | I/O, config, logging, parallel, cache, validation, workflow orchestration |
| **DNA** | [dna/](src/metainformant/dna/) | ✅ **Complete** | Sequences, alignment, phylogeny, population genetics, variant analysis |
| **RNA** | [rna/](src/metainformant/rna/) | ✅ **Complete & Verified** | AMALGKIT integration, workflow orchestration, expression quantification |
| **Protein** | [protein/](src/metainformant/protein/) | ✅ **Complete** | Sequences, structures, AlphaFold, UniProt, functional analysis |
| **GWAS** | [gwas/](src/metainformant/gwas/) | ✅ **Complete** | Association testing, QC, population structure, visualization |
| **Math** | [math/](src/metainformant/math/) | ✅ **Complete** | Population genetics, coalescent, selection, epidemiology |
| **Visualization** | [visualization/](src/metainformant/visualization/) | ✅ **Complete** | 70+ plot types, animations, publication-quality output |
| **Ontology** | [ontology/](src/metainformant/ontology/) | ✅ **Complete** | GO analysis, semantic similarity, functional annotation |
| **Quality** | [quality/](src/metainformant/quality/) | ✅ **Complete** | FASTQ analysis, validation, contamination detection |

### 🟡 **Functional Modules** (Partial Implementation)

| Category | Module | Status | Key Features | Coverage |
|----------|--------|--------|--------------|----------|
| **ML** | [ml/](src/metainformant/ml/) | 🟡 **Partial** | Classification, regression, feature selection | 75% |
| **Networks** | [networks/](src/metainformant/networks/) | 🟡 **Partial** | Graph algorithms, community detection | 78% |
| **Multi-Omics** | [multiomics/](src/metainformant/multiomics/) | 🟡 **Partial** | Integration, joint PCA, correlation | 72% |
| **Single-Cell** | [singlecell/](src/metainformant/singlecell/) | 🟡 **Partial** | Preprocessing, clustering, DE analysis | 74% |
| **Epigenome** | [epigenome/](src/metainformant/epigenome/) | 🟡 **Partial** | Methylation, ChIP-seq, ATAC-seq | 76% |
| **Phenotype** | [phenotype/](src/metainformant/phenotype/) | 🟡 **Partial** | AntWiki integration, trait analysis | 79% |
| **Ecology** | [ecology/](src/metainformant/ecology/) | 🟡 **Partial** | Community diversity, environmental | 77% |
| **Life Events** | [life_events/](src/metainformant/life_events/) | 🟡 **Partial** | Event sequences, embeddings | 73% |
| **Simulation** | [simulation/](src/metainformant/simulation/) | 🟡 **Partial** | Sequence simulation, ecosystems | 71% |
| **Information** | [information/](src/metainformant/information/) | 🟡 **Partial** | Entropy, mutual information | 80% |

### 🔬 **Specialized Domain Modules**

| Category | Module | Status | Key Features | Coverage |
|----------|--------|--------|--------------|----------|
| **Long-Read** | [longread/](src/metainformant/longread/) | 🟡 **Partial** | PacBio/ONT sequencing, assembly, error correction | 65% |
| **Metagenomics** | [metagenomics/](src/metainformant/metagenomics/) | 🟡 **Partial** | Taxonomic profiling, functional annotation | 60% |
| **Structural Variants** | [structural_variants/](src/metainformant/structural_variants/) | 🟡 **Partial** | SV/CNV detection, breakpoint resolution | 55% |
| **Spatial** | [spatial/](src/metainformant/spatial/) | 🟡 **Partial** | Spatial transcriptomics, tissue mapping | 50% |
| **Pharmacogenomics** | [pharmacogenomics/](src/metainformant/pharmacogenomics/) | 🟡 **Partial** | Drug-gene interactions, variant interpretation | 55% |
| **Metabolomics** | [metabolomics/](src/metainformant/metabolomics/) | 🟡 **Partial** | MS data processing, pathway mapping | 50% |
| **Menu** | [menu/](src/metainformant/menu/) | 🟡 **Partial** | Interactive CLI menu, workflow navigation | 70% |

## Module Overview

### Complete Module Reference

All modules live in [`src/metainformant/`](src/metainformant/) with documentation in each module's `README.md`.

| Module | Files | Description | Key Components | Docs |
|--------|-------|-------------|----------------|------|
| **Core Infrastructure** |||||
| [`core/`](src/metainformant/core/) | 26 | Shared utilities, I/O, logging, config, parallel processing, caching | [`io/`](src/metainformant/core/io/), [`data/`](src/metainformant/core/data/), [`execution/`](src/metainformant/core/execution/) | [README](src/metainformant/core/README.md) |
| **Molecular Analysis** |||||
| [`dna/`](src/metainformant/dna/) | 27 | DNA sequences, alignment, phylogenetics, population genetics, variants | [`sequence/`](src/metainformant/dna/sequence/), [`alignment/`](src/metainformant/dna/alignment/), [`population/`](src/metainformant/dna/population/) | [README](src/metainformant/dna/README.md) |
| [`rna/`](src/metainformant/rna/) | 29 | RNA-seq workflows, amalgkit integration, expression quantification | [`amalgkit/`](src/metainformant/rna/amalgkit/), [`engine/`](src/metainformant/rna/engine/), [`analysis/`](src/metainformant/rna/analysis/) | [README](src/metainformant/rna/README.md) |
| [`protein/`](src/metainformant/protein/) | 17 | Protein sequences, structure analysis, AlphaFold, UniProt integration | [`sequence/`](src/metainformant/protein/sequence/), [`structure/`](src/metainformant/protein/structure/), [`database/`](src/metainformant/protein/database/) | [README](src/metainformant/protein/README.md) |
| [`epigenome/`](src/metainformant/epigenome/) | 8 | Methylation analysis, ChIP-seq, ATAC-seq, chromatin accessibility | [`assays/`](src/metainformant/epigenome/assays/), [`chromatin_state/`](src/metainformant/epigenome/chromatin_state/), [`peak_calling/`](src/metainformant/epigenome/peak_calling/) | [README](src/metainformant/epigenome/README.md) |
| **Statistical & ML** |||||
| [`gwas/`](src/metainformant/gwas/) | 39 | GWAS, fine-mapping, eQTL analysis, colocalization, visualization | [`finemapping/`](src/metainformant/gwas/finemapping/), [`visualization/`](src/metainformant/gwas/visualization/), [`analysis/`](src/metainformant/gwas/analysis/) | [README](src/metainformant/gwas/README.md) |
| [`math/`](src/metainformant/math/) | 20 | Population genetics theory, coalescent, selection, epidemiology | [`population_genetics/`](src/metainformant/math/population_genetics/), [`epidemiology/`](src/metainformant/math/epidemiology/), [`evolutionary_dynamics/`](src/metainformant/math/evolutionary_dynamics/) | [README](src/metainformant/math/README.md) |
| [`ml/`](src/metainformant/ml/) | 12 | Machine learning pipelines, classification, regression, features | [`models/`](src/metainformant/ml/models/), [`features/`](src/metainformant/ml/features/), [`llm/`](src/metainformant/ml/llm/) | [README](src/metainformant/ml/README.md) |
| [`information/`](src/metainformant/information/) | 14 | Information theory, Shannon entropy, mutual information, semantic similarity | [`metrics/`](src/metainformant/information/metrics/), [`integration/`](src/metainformant/information/integration/) | [README](src/metainformant/information/README.md) |
| **Systems Biology** |||||
| [`networks/`](src/metainformant/networks/) | 9 | Biological networks, graph algorithms, community detection, pathways | [`analysis/`](src/metainformant/networks/analysis/), [`interaction/`](src/metainformant/networks/interaction/) | [README](src/metainformant/networks/README.md) |
| [`multiomics/`](src/metainformant/multiomics/) | 6 | Multi-omic integration, joint PCA, cross-omic correlation | [`analysis/`](src/metainformant/multiomics/analysis/), [`methods/`](src/metainformant/multiomics/methods/) | [README](src/metainformant/multiomics/README.md) |
| [`singlecell/`](src/metainformant/singlecell/) | 9 | scRNA-seq preprocessing, clustering, differential expression | [`data/`](src/metainformant/singlecell/data/), [`analysis/`](src/metainformant/singlecell/analysis/), [`visualization/`](src/metainformant/singlecell/visualization/) | [README](src/metainformant/singlecell/README.md) |
| [`simulation/`](src/metainformant/simulation/) | 7 | Synthetic data, agent-based models, sequence simulation, ecosystems | [`models/`](src/metainformant/simulation/models/), [`workflow/`](src/metainformant/simulation/workflow/), [`benchmark/`](src/metainformant/simulation/benchmark/) | [README](src/metainformant/simulation/README.md) |
| **Annotation & Metadata** |||||
| [`ontology/`](src/metainformant/ontology/) | 7 | Gene Ontology, functional annotation, semantic similarity | [`core/`](src/metainformant/ontology/core/), [`query/`](src/metainformant/ontology/query/), [`visualization/`](src/metainformant/ontology/visualization/) | [README](src/metainformant/ontology/README.md) |
| [`phenotype/`](src/metainformant/phenotype/) | 15 | Phenotypic data curation, AntWiki integration, trait analysis | [`analysis/`](src/metainformant/phenotype/analysis/), [`data/`](src/metainformant/phenotype/data/), [`behavior/`](src/metainformant/phenotype/behavior/) | [README](src/metainformant/phenotype/README.md) |
| [`ecology/`](src/metainformant/ecology/) | 7 | Community diversity, environmental correlations, species matrices | [`analysis/`](src/metainformant/ecology/analysis/), [`phylogenetic/`](src/metainformant/ecology/phylogenetic/), [`visualization/`](src/metainformant/ecology/visualization/) | [README](src/metainformant/ecology/README.md) |
| [`life_events/`](src/metainformant/life_events/) | 9 | Life course analysis, event sequences, temporal embeddings | [`models/`](src/metainformant/life_events/models/), [`workflow/`](src/metainformant/life_events/workflow/) | [README](src/metainformant/life_events/README.md) |
| **Utilities** |||||
| [`quality/`](src/metainformant/quality/) | 4 | FASTQ quality assessment, validation, contamination detection | [`io/`](src/metainformant/quality/io/), [`analysis/`](src/metainformant/quality/analysis/), [`reporting/`](src/metainformant/quality/reporting/) | [README](src/metainformant/quality/README.md) |
| [`visualization/`](src/metainformant/visualization/) | 22 | 70+ plot types, heatmaps, networks, animations, publication-ready | [`plots/`](src/metainformant/visualization/plots/), [`genomics/`](src/metainformant/visualization/genomics/), [`analysis/`](src/metainformant/visualization/analysis/) | [README](src/metainformant/visualization/README.md) |
| **Specialized Domains** |||||
| [`longread/`](src/metainformant/longread/) | 19 | Long-read sequencing (PacBio, ONT), assembly, error correction | [`assembly/`](src/metainformant/longread/assembly/), [`quality/`](src/metainformant/longread/quality/) | [README](src/metainformant/longread/README.md) |
| [`metagenomics/`](src/metainformant/metagenomics/) | 11 | Metagenomic analysis, taxonomic profiling, functional annotation | [`taxonomy/`](src/metainformant/metagenomics/taxonomy/), [`functional/`](src/metainformant/metagenomics/functional/) | [README](src/metainformant/metagenomics/README.md) |
| [`pharmacogenomics/`](src/metainformant/pharmacogenomics/) | 12 | Drug-gene interactions, pharmacokinetics, variant interpretation | [`interactions/`](src/metainformant/pharmacogenomics/interactions/) | [README](src/metainformant/pharmacogenomics/README.md) |
| [`spatial/`](src/metainformant/spatial/) | 11 | Spatial transcriptomics, tissue mapping, spatial statistics | [`analysis/`](src/metainformant/spatial/analysis/) | [README](src/metainformant/spatial/README.md) |
| [`structural_variants/`](src/metainformant/structural_variants/) | 9 | SV detection, CNV analysis, breakpoint resolution | [`detection/`](src/metainformant/structural_variants/detection/) | [README](src/metainformant/structural_variants/README.md) |
| [`menu/`](src/metainformant/menu/) | 4 | Interactive CLI menu system, workflow navigation | [`ui/`](src/metainformant/menu/ui/) | [README](src/metainformant/menu/README.md) |

**Total: 26 modules, 560 Python files**

## Documentation

### Quick Links

- **[Documentation Guide](docs/DOCUMENTATION_GUIDE.md)** - Complete navigation guide
- **[Quick Start](QUICKSTART.md)** - Fast setup commands
- **[Architecture](docs/architecture.md)** - System design
- **[Technical Specification](SPEC.md)** - Design standards
- **[Testing Guide](docs/testing.md)** - Comprehensive testing documentation
- **[CLI Reference](docs/cli.md)** - Command-line interface
- **[eQTL Integration](docs/eqtl/README.md)** - eQTL pipeline documentation

### Module Documentation

Each module has documentation in `src/metainformant/<module>/README.md` and `docs/<module>/`.

## Scripts & Workflows

The [`scripts/`](scripts/) directory contains production-ready workflow orchestrators:

- **Package Management**: Setup, testing, quality control
- **RNA-seq (Amalgkit)**: Multi-species workflows, amalgkit integration
- **GWAS (Variants)**: Genome-scale association studies
- **eQTL Integration**: RNA-seq + Variant cross-omics integration pipelines
- **Module Orchestrators**: ✅ Complete workflow scripts for all domains (core, DNA, RNA, protein, networks, multiomics, single-cell, quality, simulation, visualization, epigenome, ecology, ontology, phenotype, ML, math, gwas, information, life_events)

See [`scripts/README.md`](scripts/README.md) for documentation.

### CLI Interface

All modules are accessible via the unified CLI:

```bash
# Setup and environment
uv run metainformant setup --with-amalgkit

# Domain workflows
uv run metainformant dna fetch --assembly GCF_000001405.40
uv run metainformant dna align --input data/sequences.fasta --output output/dna/alignment
uv run metainformant dna variants --input data/variants.vcf --format vcf --output output/dna/variants
uv run metainformant rna run --work-dir output/rna --threads 8 --species Apis_mellifera
uv run metainformant rna run-config --config config/amalgkit/amalgkit_pbarbatus.yaml
uv run metainformant protein taxon-ids --file data/taxon_ids.txt
uv run metainformant protein rmsd-ca --pdb-a data/structure1.pdb --pdb-b data/structure2.pdb
uv run metainformant gwas run --config config/gwas/gwas_template.yaml

# Epigenome and annotation
uv run metainformant epigenome run --methylation data/methylation.tsv --output output/epigenome
uv run metainformant ontology run --go data/go.obo --output output/ontology
uv run metainformant phenotype run --input data/traits.csv --output output/phenotype
uv run metainformant ecology run --input data/species.csv --output output/ecology --diversity

# Analysis and modeling
uv run metainformant math popgen --input data/sequences.fasta --output output/math/popgen
uv run metainformant math coalescent --n-samples 10 --output output/math/coalescent
uv run metainformant information entropy --input data/seqs.fasta --output output/information
uv run metainformant simulation run --model sequences --output output/simulation

# Systems biology
uv run metainformant networks run --input data/interactions.tsv --output output/networks
uv run metainformant multiomics run --genomics data/genomics.tsv --output output/multiomics
uv run metainformant singlecell run --input data/counts.h5ad --output output/singlecell --qc
uv run metainformant quality run --fastq data/reads.fq --output output/quality --analyze-fastq
uv run metainformant ml run --features data/features.csv --output output/ml --classify
uv run metainformant visualization run --input data/matrix.csv --plot-type heatmap --output output/visualization
uv run metainformant life-events embed --input data/events.json --output output/life_events/embeddings

# See all available commands
uv run metainformant --help
```

See [`docs/cli.md`](docs/cli.md) for CLI documentation.

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

```python
from metainformant.rna import AmalgkitWorkflowConfig, plan_workflow, execute_workflow, check_cli_available

# Check if amalgkit is available
available, help_text = check_cli_available()
if not available:
    print(f"Amalgkit not available: {help_text}")

# Configure workflow
config = AmalgkitWorkflowConfig(
    work_dir="output/amalgkit/work",
    threads=8,
    species_list=["Apis_mellifera"]
)

# Plan workflow steps
steps = plan_workflow(config)
print(f"Planned {len(steps)} workflow steps")

# Execute workflow
results = execute_workflow(config)
for step, result in results.items():
    print(f"{step}: exit code {result.returncode}")
```

```bash
# End-to-end workflow for a single species (recommended)
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pbarbatus.yaml

# Check status
python3 scripts/rna/run_workflow.py --config config/amalgkit/amalgkit_pbarbatus.yaml --status

# Alternative: Bash-based orchestrator
bash scripts/rna/amalgkit/run_amalgkit.sh --config config/amalgkit/amalgkit_pbarbatus.yaml
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

### eQTL Integration Pipeline

The eQTL pipeline bridges the genomic variants from the **GWAS** pipeline with the gene expression matrices provided by the **Amalgkit (RNA)** pipeline.

```bash
# Run the pipeline leveraging real Amalgkit RNA-seq quantification data
uv run python scripts/eqtl/run_eqtl_real.py

# Or explore the logic with synthetic data
uv run python scripts/eqtl/run_eqtl_demo.py
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

### Protein Analysis

```python
from metainformant.protein import sequences, alignment, structure

# Read protein sequences
proteins = sequences.read_fasta("data/proteins.fasta")

# Pairwise alignment
align_result = alignment.global_align(proteins["seq1"], proteins["seq2"])

# Structure analysis
structure_data = structure.load_pdb("data/structure.pdb")
contacts = structure.analyze_contacts(structure_data)
```

### Epigenome Analysis

```python
from metainformant.epigenome import methylation, chipseq

# Methylation analysis
meth_data = methylation.load_bedgraph("data/methylation.bedgraph")
regions = methylation.find_dmr(meth_data, threshold=0.3)

# ChIP-seq peak calling
peaks = chipseq.call_peaks("data/chipseq.bam", "data/control.bam")
```

### Ontology Analysis

```python
from metainformant.ontology.core import go
from metainformant.ontology.query import query

# Load Gene Ontology
go_graph = go.load_obo("data/go.obo")

# Query ontology
terms = query.get_ancestors(go_graph, "GO:0008150")
similarity = query.semantic_similarity(go_graph, "GO:0008150", "GO:0008151")
```

### Phenotype Analysis

```python
from metainformant.phenotype import life_course, antwiki

# Life course analysis
traits = life_course.load_traits("data/traits.csv")
curated = life_course.curate_traits(traits)

# AntWiki integration
species_data = antwiki.fetch_species("Pogonomyrmex_barbatus")
```

### Ecology Analysis

```python
from metainformant.ecology import community, environmental

# Community analysis
species_matrix = community.load_matrix("data/species.csv")
diversity = community.calculate_diversity(species_matrix)

# Environmental data
env_data = environmental.load_data("data/environment.csv")
correlations = environmental.analyze_correlations(species_matrix, env_data)
```

### Mathematical Biology

```python
from metainformant.math import popgen, coalescent

# Population genetics
sequences = ["ATCGATCG", "ATCGTTCG", "ATCGATCG"]
fst = popgen.fst(sequences, populations=[0, 0, 1])

# Coalescent simulation
tree = coalescent.simulate_coalescent(n_samples=10, Ne=1000)
```

### Single-Cell Analysis

```python
from metainformant.singlecell import preprocessing, clustering

# Load single-cell data
adata = preprocessing.load_h5ad("data/counts.h5ad")

# Preprocessing
adata = preprocessing.filter_cells(adata, min_genes=200)
adata = preprocessing.normalize(adata)

# Clustering
clusters = clustering.leiden(adata, resolution=0.5)
```

### Quality Control

```python
from metainformant.quality import fastq, metrics

# FASTQ quality assessment
qc_report = fastq.assess_quality("data/reads.fastq")
print(f"Mean quality: {qc_report['mean_quality']}")

# General metrics
quality_score = metrics.calculate_quality(data_matrix)
```

### Machine Learning

```python
from metainformant.ml import classification, features

# Feature extraction
features = features.extract_features(data, method="pca", n_components=50)

# Classification
model = classification.train_classifier(
    X_train, y_train, method="random_forest"
)
predictions = model.predict(X_test)
```

### Simulation

```python
from metainformant.simulation import sequences, ecosystems

# Sequence simulation
sim_seqs = sequences.simulate_sequences(
    n_sequences=100, length=1000, mutation_rate=0.01
)

# Ecosystem simulation
ecosystem = ecosystems.simulate_community(
    n_species=50, interactions="random"
)
```

### Core Utilities

```python
from metainformant.core import io, paths, logging

# I/O operations
data = io.load_json("config/example.yaml")
io.dump_json(results, "output/results.json")

# Path handling
resolved = paths.expand_and_resolve("~/data/input.txt")
is_safe = paths.is_within(resolved, base_path="/safe/directory")

# Logging
logger = logging.get_logger(__name__)
logger.info("Processing data")
```

## Development

### Running Tests

```bash
# All tests
bash scripts/package/test.sh

# Fast tests only
bash scripts/package/test.sh --mode fast

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
- **Multi-omics**: Integration methods implemented; additional dependencies may be required
- **Single-cell**: Requires `scipy`, `scanpy`, `anndata` (see [Single-Cell Documentation](docs/singlecell/README.md))
- **Network Analysis**: Algorithms implemented; regulatory network features may need enhancement

### GWAS Module

- **Variant Download**: Database download (dbSNP, 1000 Genomes) is a placeholder; use SRA-based workflow or provide VCF files
- **Functional Annotation**: Requires external tools (ANNOVAR, VEP, SnpEff) for variant annotation
- **Mixed Models**: Relatedness adjustment implemented; MLM methods may require GCTA/EMMAX integration

### Test Coverage

Some modules have lower test success rates due to optional dependencies:

- **Single-cell**: Requires scientific dependencies (`scanpy`, `anndata`)
- **Multi-omics**: Framework exists, tests may skip without dependencies
- **Network Analysis**: Tests pass; features may need additional setup

See [Testing Guide](docs/testing.md) for detailed testing documentation and coverage information.

## Best Practices

### File Naming

- ✅ Use informative names: `sample_pca_biplot_colored_by_treatment.png`
- ❌ Avoid generic names: `plot1.png`, `output.png`

### Output Organization

- All outputs in `output/` directory
- Configuration saved with results
- Visualizations in subdirectories with metadata

### No Mocking Policy

- All tests use implementations
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

## Recent Improvements

### Performance Enhancements

- **Intelligent Caching**: Automatic caching for expensive computations (Tajima's constants, entropy calculations)
- **NumPy Vectorization**: Optimized mathematical operations for 10-100x performance improvements
- **Progress Tracking**: Real-time progress bars for long-running analyses
- **Memory Optimization**: Efficient algorithms for large datasets

### Enhanced Documentation

- **Comprehensive Tutorials**: End-to-end guides for DNA, RNA, GWAS, and information theory workflows
- **Method Comparison Guides**: Decision-making guides for choosing analysis algorithms
- **Extended FAQ**: Troubleshooting and usage guidance for common scenarios
- **Standardized Docstrings**: Consistent formatting with examples and DOI citations

### Testing & Reliability

- **Expanded Test Coverage**: 37+ new comprehensive tests with real implementations
- **Validation Enhancements**: Improved parameter validation and error handling
- **Cross-Platform Compatibility**: Python 3.14 support and external drive optimization
- **Integration Testing**: Verified cross-module functionality

### New Features

- **Enhanced GWAS Visualization**: Complete visualization suite for population structure, effects, and comparisons
- **Information Theory Workflows**: Batch processing with progress tracking
- **Protein Proteome Analysis**: Taxonomy ID processing and proteome utilities
- **Advanced Error Handling**: Structured error reporting with actionable guidance

## Citation

If you use METAINFORMANT in your research, please cite this repository:

```bibtex
@software{metainformant2025,
  author = {MetaInformAnt Development Team},
  title = {MetaInformAnt: Comprehensive Bioinformatics Toolkit},
  year = {2025},
  url = {https://github.com/docxology/MetaInformAnt},
  version = {0.2.6}
}
```

## License

This project is licensed under the Apache License, Version 2.0 - see [LICENSE](LICENSE) for details.

## Contact

- **Repository**: <https://github.com/docxology/MetaInformAnt>
- **Issues**: <https://github.com/docxology/MetaInformAnt/issues>
- **Documentation**: <https://github.com/docxology/MetaInformAnt/blob/main/docs/>

## Acknowledgments

- Developed with AI assistance from Cursor's Code Assistant (grok-code-fast-1)
- Built on established bioinformatics tools and libraries
- Community contributions and feedback

---

**Status**: Active Development | **Version**: 0.2.6 | **Python**: 3.11+ | **License**: Apache 2.0
