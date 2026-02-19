# Architecture

Comprehensive system architecture of METAINFORMANT, showing the modular design, data flow patterns, and integration relationships across all biological domains.

## System Overview

METAINFORMANT is a comprehensive bioinformatics toolkit designed for multi-omic data analysis. The architecture follows a modular, domain-driven design with clear separation of concerns and extensive integration capabilities.

### Core Design Principles

- **Modular Architecture**: Each biological domain is self-contained with well-defined interfaces
- **Data Flow Integration**: Seamless data exchange between modules through standardized formats
- **Quality Assurance**: Comprehensive validation and error handling throughout
- **Extensibility**: Easy addition of new analysis methods and data types
- **Performance**: Optimized algorithms for large-scale biological data processing

## Module Architecture

```mermaid
graph TB
    %% Core Infrastructure Layer
    subgraph "Core Infrastructure"
        CORE["Core Utilities<br/>config | io | logging | parallel | paths | cache | workflow | validation"]
    end

    %% Data Processing Layer
    subgraph "Data Processing"
        DNA["DNA Analysis<br/>sequences | alignment | phylogeny | population"]
        RNA["RNA Analysis<br/>amalgkit | quantification | transcriptomics | workflow"]
        PROTEIN["Protein Analysis<br/>sequences | structure | AlphaFold | analysis"]
        EPIGENOME["Epigenome Analysis<br/>methylation | ChIP-seq | ATAC-seq | chromatin"]
    end

    %% Analysis Layer
    subgraph "Analysis Methods"
        GWAS["GWAS Analysis<br/>association | QC | population | visualization"]
        MATH["Mathematical Biology<br/>population genetics | coalescent | selection"]
        ML["Machine Learning<br/>classification | regression | feature selection"]
        INFO["Information Theory<br/>entropy | MI | similarity | semantic measures"]
    end

    %% Systems Biology Layer
    subgraph "Systems Biology"
        NETWORKS["Network Analysis<br/>PPI | pathways | community detection"]
        MULTIOMICS["Multi-Omics Integration<br/>integration | harmonization | joint analysis"]
        SINGLECELL["Single-Cell Analysis<br/>preprocessing | clustering | trajectory | DE"]
        SIMULATION["Simulation<br/>sequences | ecosystems | agent-based | evolution"]
    end

    %% Knowledge Layer
    subgraph "Knowledge & Metadata"
        ONTOLOGY["Ontology<br/>GO | functional annotation | semantic similarity"]
        PHENOTYPE["Phenotype Analysis<br/>traits | life course | AntWiki | curation"]
        ECOLOGY["Ecology<br/>communities | diversity | environmental analysis"]
        LIFEEVENTS["Life Events<br/>sequences | embeddings | temporal prediction"]
    end

    subgraph "Utilities"
        QUALITY["Quality Control<br/>FASTQ | assembly | metrics | validation"]
        VISUALIZATION["Visualization<br/>plots | animations | trees | networks"]
    end

    %% Specialized Domains Layer
    subgraph "Specialized Domains"
        LONGREAD["Long-Read Sequencing<br/>PacBio | ONT | assembly | error correction"]
        METAGENOMICS["Metagenomics<br/>taxonomy | functional | community"]
        SV["Structural Variants<br/>detection | CNV | breakpoints"]
        SPATIAL["Spatial Transcriptomics<br/>tissue mapping | spatial stats"]
        PHARMA["Pharmacogenomics<br/>drug-gene | pharmacokinetics | variants"]
        METAB["Metabolomics<br/>MS processing | pathway mapping"]
        MENU["Menu System<br/>interactive CLI | workflow navigation"]
    end

    %% Data Flow
    CORE --> DNA
    CORE --> RNA
    CORE --> PROTEIN
    CORE --> EPIGENOME
    CORE --> GWAS
    CORE --> MATH
    CORE --> ML
    CORE --> INFO
    CORE --> NETWORKS
    CORE --> MULTIOMICS
    CORE --> SINGLECELL
    CORE --> SIMULATION
    CORE --> ONTOLOGY
    CORE --> PHENOTYPE
    CORE --> ECOLOGY
    CORE --> LIFEEVENTS
    CORE --> QUALITY
    CORE --> VISUALIZATION
    CORE --> LONGREAD
    CORE --> METAGENOMICS
    CORE --> SV
    CORE --> SPATIAL
    CORE --> PHARMA
    CORE --> METAB
    CORE --> MENU

    DNA --> RNA
    DNA --> PROTEIN
    DNA --> GWAS
    RNA --> SINGLECELL
    RNA --> MULTIOMICS
    PROTEIN --> NETWORKS
    PROTEIN --> ONTOLOGY
    EPIGENOME --> DNA
    EPIGENOME --> NETWORKS
    ONTOLOGY --> MULTIOMICS
    PHENOTYPE --> GWAS
    PHENOTYPE --> LIFEEVENTS
    ECOLOGY --> NETWORKS
    MATH --> GWAS
    MATH --> DNA
    INFO --> NETWORKS
    INFO --> ONTOLOGY
    NETWORKS --> MULTIOMICS
    SINGLECELL --> MULTIOMICS
    QUALITY --> ALL["All Analysis Modules"]
    VISUALIZATION --> ALL

    %% External Dependencies
    RNA -.->|"Amalgkit CLI"| RNA
    SINGLECELL -.->|"scanpy, anndata"| SINGLECELL
    PROTEIN -.->|"AlphaFold"| PROTEIN
    GWAS -.->|"bcftools, GATK"| GWAS

    %% Data Sources
    NCBI[("NCBI Genomes")] --> DNA
    SRA[("SRA Sequences")] --> RNA
    SRA --> SINGLECELL
    SRA --> GWAS
    PDB[("PDB Structures")] --> PROTEIN
    GEO[("GEO Expression")] --> RNA
    GO[("Gene Ontology")] --> ONTOLOGY

    %% Styling

    class CORE core
    class DNA,RNA,PROTEIN,EPIGENOME data
    class GWAS,MATH,ML,INFO analysis
    class NETWORKS,MULTIOMICS,SINGLECELL,SIMULATION systems
    class ONTOLOGY,PHENOTYPE,ECOLOGY,LIFEEVENTS knowledge
    class QUALITY,VISUALIZATION utility
    class LONGREAD,METAGENOMICS,SV,SPATIAL,PHARMA,METAB,MENU specialized
    class NCBI,SRA,PDB,GEO,GO external
```

### Data Flow Architecture

```mermaid
graph TD
    A["Raw Biological Data"] --> B["Ingestion Layer"]
    B --> C["Format Detection"]
    C --> D["Data Validation"]

    D --> E{Data Type}
    E -->|Genomic| F["DNA Processing Pipeline"]
    E -->|Transcriptomic| G["RNA Processing Pipeline"]
    E -->|Proteomic| H["Protein Processing Pipeline"]
    E -->|Epigenomic| I["Epigenome Processing Pipeline"]

    F --> J["Quality Control"]
    G --> J
    H --> J
    I --> J

    J --> K["Feature Extraction"]
    K --> L["Statistical Analysis"]

    L --> M{Integration Strategy}
    M -->|Single-omic| N["Individual Analysis"]
    M -->|Multi-omic| O["Joint Analysis"]

    N --> P["Domain-Specific Results"]
    O --> Q["Systems-Level Results"]

    P --> R["Visualization"]
    Q --> R

    R --> S["Publication Figures"]
    S --> T["Scientific Insights"]


    subgraph "Data Sources"
        U["NCBI"] -.-> F
        V["SRA"] -.-> G
        W["PDB"] -.-> H
        X["ENCODE"] -.-> I
    end

    subgraph "Processing Stages"
        Y["Normalization"] -.-> K
        Z["Transformation"] -.-> K
        AA["Filtering"] -.-> K
        BB["Scaling"] -.-> K
    end

    subgraph "Analysis Types"
        CC["Differential Expression"] -.-> L
        DD["Association Testing"] -.-> L
        EE["Network Inference"] -.-> L
        FF["Clustering"] -.-> L
    end
```

### Integration Framework

```mermaid
graph TD
    A["Module Interfaces"] --> B["Standardized APIs"]
    B --> C["Data Format Standards"]

    C --> D["Interoperability Layer"]
    D --> E["Cross-Module Communication"]

    E --> F{Integration Pattern}
    F -->|Pipeline| G["Sequential Processing"]
    F -->|Workflow| H["Orchestrated Execution"]
    F -->|Joint| I["Simultaneous Analysis"]

    G --> J["Result Aggregation"]
    H --> J
    I --> J

    J --> K["Unified Output"]
    K --> L["Quality Assurance"]

    L --> M["Final Results"]


    subgraph "API Standards"
        N["Function Signatures"] -.-> B
        O["Parameter Types"] -.-> B
        P["Return Formats"] -.-> B
        Q["Error Handling"] -.-> B
    end

    subgraph "Data Standards"
        R["FASTA Format"] -.-> C
        S["Count Matrices"] -.-> C
        T["Pandas DataFrames"] -.-> C
        U["NetworkX Graphs"] -.-> C
    end

    subgraph "Integration Mechanisms"
        V["Core Utilities"] -.-> D
        W["Config Files"] -.-> D
        X["Workflow Orchestration"] -.-> D
        Y["Event System"] -.-> D
    end
```

### Quality Assurance Architecture

```mermaid
graph TD
    A["Code Development"] --> B["Static Analysis"]
    B --> C["Type Checking"]
    C --> D["Linting"]

    D --> E["Unit Testing"]
    E --> F["Integration Testing"]
    F --> G["Performance Testing"]

    G --> H["Documentation"]
    H --> I["Validation"]

    I --> J{Standards Met?}
    J -->|Yes| K["Release Ready"]
    J -->|No| L["Issue Resolution"]

    L --> A

    K --> M["Deployment"]


    subgraph "Code Quality"
        N["mypy"] -.-> C
        O["ruff"] -.-> D
        P["black"] -.-> D
        Q["isort"] -.-> D
    end

    subgraph "Testing Framework"
        R["pytest"] -.-> E
        S["Real Data"] -.-> E
        T["No Mocks"] -.-> E
        U["Integration"] -.-> F
    end

    subgraph "Performance"
        V["Benchmarks"] -.-> G
        W["Memory Profiling"] -.-> G
        X["Scalability Tests"] -.-> G
    end

    subgraph "Documentation"
        Y["README Files"] -.-> H
        Z["API Docs"] -.-> H
        AA["Examples"] -.-> H
        BB["Guides"] -.-> H
    end
```

## Detailed Module Relationships

```mermaid
flowchart LR
  subgraph CLI
    A["metainformant.__main__<br/>argparse CLI"]
  end
  subgraph Core["Core Utilities"]
    C1[config]
    C2[io]
    C3[logging]
    C4[text]
    C5[parallel]
    C6[hash]
    C7[paths]
    C8[cache]
    C9[db]
  end
  subgraph DNA
    D1[sequences]
    D2[alignment]
    D3[msa]
    D4[phylogeny]
    D5[population]
    D6[entrez/ncbi]
    D7[variants]
  end
  subgraph RNA
    R1[amalgkit]
    R2[workflow]
    R3[configs]
  end
  subgraph GWAS["Genome-Wide Association"]
    G1[quality]
    G2[structure]
    G3[association]
    G4[visualization]
  end
  subgraph SingleCell["Single-Cell Genomics"]
    SC1[preprocessing]
    SC2[dimensionality]
    SC3[clustering]
    SC4[trajectory]
    SC5[visualization]
    SC6[integration]
  end
  subgraph Quality["Quality Control"]
    Q1[fastq]
  end
  subgraph Simulation
    S1[sequences]
    S2[rna]
    S3[agents]
  end
  subgraph Math
    M1[price]
    M2[selection]
    M3[ddm]
  end
  subgraph Viz
    V1[trees]
    V2[plots]
    V3[animations]
  end
  subgraph Information["Information Theory"]
    I1[entropy]
    I2[mutual-information]
    I3[profile]
  end
  subgraph LifeEvents["Life Events"]
    L1[embed]
    L2[predict]
    L3[interpret]
  end
  subgraph Protein["Protein Analysis"]
    O1["protein"]
  end
  subgraph Ontology["Ontology"]
    O2["ontology"]
  end
  subgraph Phenotype["Phenotype Analysis"]
    O3["phenotype"]
  end
  subgraph Epigenome["Epigenome Analysis"]
    O4["epigenome"]
  end
  subgraph Ecology["Ecology"]
    O5["ecology"]
  end
  subgraph ML["Machine Learning"]
    O6["ml"]
  end
  subgraph MultiOmics["Multi-Omics Integration"]
    O7["multiomics"]
  end
  subgraph Networks["Network Analysis"]
    O8["networks"]
  end
  subgraph Specialized["Specialized Domains"]
    SP1["longread"]
    SP2["metagenomics"]
    SP3["structural_variants"]
    SP4["spatial"]
    SP5["pharmacogenomics"]
    SP6["metabolomics"]
    SP7["menu"]
  end
  A --> D1 & D2 & D3 & D4 & D5
  A --> R1 & R2
  A --> G1 & G2 & G3 & G4
  A --> SC1 & SC2 & SC3 & SC4 & SC5 & SC6
  A --> Q1
  A --> S1 & S2 & S3
  A --> M1 & M2 & M3
  A --> V1 & V2 & V3
  A --> I1 & I2 & I3
  A --> L1 & L2 & L3
  A --> O1 & O2 & O3 & O4 & O5 & O6 & O7 & O8
  A --> SP1 & SP2 & SP3 & SP4 & SP5 & SP6 & SP7
  C1 -.-> A
  C2 -.-> A
  C1 -.-> D1
  C2 -.-> D1
  C1 -.-> R2
  C2 -.-> R1
  C2 -.-> SC1
  C5 -.-> SC2
  C8 -.-> SC6
  C2 -.-> Q1
  C5 -.-> R2
  C8 -.-> R2
  C1 -.-> O6
  C2 -.-> O6
  C1 -.-> O7
  C2 -.-> O7
  C1 -.-> O8
  C2 -.-> O8
```

## Project directories and conventions

- **`config/`**: Declarative configuration and options for runs. Read by `metainformant.core.config` and consumed across domains. Environment variables may override values.
- **`data/`**: Canonical datasets and local databases. Treated as read-mostly inputs and long-lived artifacts under versioned subfolders.
- **`output/`**: All run and test outputs. Ephemeral, reproducible, safe to delete. Modules must default to writing here unless a user-specified path is provided.

Guidelines:

- `core` owns config loading and path resolution; domain modules do not hardcode absolute paths.
- Prefer parameters and env variables to override locations, but default to `config/`, `data/`, and `output/` at repo root.
- Tests and CLI runs should not write outside `output/`.

## Module Relationships and Dependencies

### Core Dependencies

All modules depend on `core` utilities:

- **config**: Configuration loading and environment variable handling
- **io**: File I/O operations (JSON, CSV, TSV, compressed formats)
- **logging**: Structured logging across all modules
- **paths**: Path validation and containment checks
- **cache**: JSON-based caching for expensive operations

### Domain Module Dependencies

#### Core Module

- **Purpose**: Shared infrastructure and utilities for all modules
- **Depends on**: None (base module)
- **Used by**: All modules (infrastructure layer)
- **Key Components**:
  - **config**: Configuration loading with environment variable overrides
  - **io**: File I/O operations (JSON, CSV, TSV, gzip-aware)
  - **logging**: Structured logging with consistent formatting
  - **paths**: Path validation and containment checks
  - **cache**: JSON-based caching for expensive operations
  - **parallel**: Parallel processing utilities
  - **hash**: Hashing and checksum functions
  - **text**: Text processing utilities
  - **db**: Database utilities
- **Integrates with**: All domain modules as the foundation layer

#### DNA Module

- **Depends on**: `core` (all utilities)
- **Used by**: `gwas`, `rna` (for genomic coordinates), `information` (sequence analysis)
- **Integrates with**: `visualization` (phylogenetic trees), `math.popgen` (population genetics theory)

#### RNA Module

- **Depends on**: `core`, `dna` (for genomic coordinates)
- **Used by**: `multiomics` (transcriptomics integration)
- **Integrates with**: External `amalgkit` CLI tool

#### GWAS Module

- **Depends on**: `core`, `dna.variants`, `dna.population`, `math.popgen`, `ml.regression`
- **Used by**: `multiomics` (genomics integration)
- **Integrates with**: `visualization` (Manhattan plots, Q-Q plots), external tools (bcftools, GATK)

#### Protein Module

- **Depends on**: `core`
- **Used by**: `networks` (protein-protein interactions), `multiomics` (proteomics integration)
- **Integrates with**: External databases (UniProt, InterPro, AlphaFold, PDB)

#### Single-Cell Module

- **Depends on**: `core`, `rna` (expression data), `ml.dimensionality`
- **Used by**: `multiomics` (single-cell omics integration), `information` (single-cell information analysis)
- **Integrates with**: External tools (scanpy, anndata) when available

#### Multi-Omics Module

- **Depends on**: `core`, `dna`, `rna`, `protein`, `singlecell`, `ml.dimensionality`
- **Used by**: Workflow orchestration scripts
- **Integrates with**: All omics modules for cross-platform integration

#### Networks Module

- **Depends on**: `core`, `protein` (PPI networks), `ontology` (functional annotation)
- **Used by**: `multiomics` (network-based integration)
- **Integrates with**: `visualization` (network plots)

#### Information Theory Module

- **Depends on**: `core`
- **Used by**: All modules (can analyze any biological data)
- **Integrates with**: `dna`, `rna`, `singlecell`, `multiomics` (cross-module integration functions)

#### Life Events Module

- **Depends on**: `core`
- **Used by**: `phenotype` (life course phenotype extraction)
- **Integrates with**: `ml` (sequence models), `visualization` (event timeline plots)

#### Phenotype Module

- **Depends on**: `core`, `life_events` (optional, for life course integration)
- **Used by**: `gwas` (phenotype-genotype associations)
- **Integrates with**: `ontology` (trait functional annotation)

#### Quality Module

- **Depends on**: `core`, `dna.fastq`
- **Used by**: `rna` (FASTQ quality control), `gwas` (variant QC)
- **Integrates with**: External tools (FastQC) when available

#### Math Module

- **Depends on**: `core`
- **Used by**: `dna.population` (population genetics theory), `gwas` (statistical models)
- **Integrates with**: Domain modules for theoretical analysis

#### Visualization Module

- **Depends on**: `core` (optional matplotlib/seaborn)
- **Used by**: All modules (plotting support)
- **Integrates with**: Domain-specific visualization (phylogenetic trees, network graphs)

#### ML Module

- **Depends on**: `core`
- **Used by**: `gwas` (regression models), `singlecell` (dimensionality reduction), `life_events` (sequence models)
- **Integrates with**: Domain modules for biological data preprocessing

#### Ontology Module

- **Depends on**: `core`
- **Used by**: `networks` (functional enrichment), `phenotype` (trait annotation)
- **Integrates with**: External databases (Gene Ontology)

#### Epigenome Module

- **Depends on**: `core`, `dna` (for genomic coordinates)
- **Used by**: `multiomics` (epigenomics integration)
- **Integrates with**: `visualization` (methylation plots)

#### Ecology Module

- **Depends on**: `core`
- **Used by**: Workflow scripts
- **Integrates with**: `math` (diversity calculations)

#### Simulation Module

- **Depends on**: `core`
- **Used by**: All modules (for generating synthetic test data)
- **Key Components**:
  - **sequences**: Synthetic sequence generation (DNA, RNA, protein)
  - **rna**: RNA expression count simulation
  - **agents**: Agent-based modeling and ecosystem simulation
- **Integrates with**: All modules for testing and validation purposes

#### Long-Read Module

- **Depends on**: `core`, `dna` (for genomic coordinates), `quality` (read QC)
- **Used by**: `multiomics` (long-read integration), `structural_variants` (SV detection from long reads)
- **Integrates with**: External tools (minimap2, Flye, Canu) when available

#### Metagenomics Module

- **Depends on**: `core`, `dna` (for sequence analysis), `quality` (read QC)
- **Used by**: `ecology` (microbiome community analysis), `multiomics` (metagenomic integration)
- **Integrates with**: External tools (Kraken2, MetaPhlAn) when available

#### Structural Variants Module

- **Depends on**: `core`, `dna` (for genomic coordinates and variant analysis)
- **Used by**: `gwas` (SV-GWAS association), `multiomics` (structural variant integration)
- **Integrates with**: External tools (Delly, Manta) when available

#### Spatial Module

- **Depends on**: `core`, `singlecell` (for expression data), `visualization` (spatial plots)
- **Used by**: `multiomics` (spatial multi-omics integration)
- **Integrates with**: External tools (Squidpy, SpatialDE) when available

#### Pharmacogenomics Module

- **Depends on**: `core`, `dna.variants` (for variant interpretation), `gwas` (association data)
- **Used by**: Clinical analysis workflows
- **Integrates with**: External databases (PharmGKB, CPIC) when available

#### Metabolomics Module

- **Depends on**: `core`
- **Used by**: `multiomics` (metabolomic integration), `networks` (metabolic pathways)
- **Integrates with**: External tools (MZmine, XCMS) when available

#### Menu Module

- **Depends on**: `core` (all utilities)
- **Used by**: CLI entry point for interactive navigation
- **Integrates with**: All domain modules for workflow discovery

### Common Integration Patterns

#### Multi-Omics Integration Workflow

```
DNA (genomics) → Multi-Omics
RNA (transcriptomics) → Multi-Omics  
Protein (proteomics) → Multi-Omics
Single-Cell → Multi-Omics
→ Joint PCA/NMF/CCA
→ Visualization
```

#### Genotype-Phenotype Association

```
DNA (variants) → GWAS
Phenotype (traits) → GWAS
→ Association Testing
→ Visualization (Manhattan plots)
```

#### Functional Annotation Pipeline

```
DNA/RNA/Protein (sequences) → Ontology
Networks (modules) → Ontology
→ Enrichment Analysis
→ Visualization
```

#### Quality Control Pipeline

```
DNA (FASTQ) → Quality
RNA (FASTQ) → Quality
→ QC Metrics
→ Filtering
→ Downstream Analysis
```

### Data Flow Patterns

1. **Input → Processing → Output**: Most modules follow this pattern using `core.io` for I/O
2. **Configuration → Workflow → Results**: Workflow modules use `core.config` for configuration
3. **Cache → Compute → Cache**: Expensive operations use `core.cache` for results
4. **Logging**: All modules use `core.logging` for consistent log messages

See also: [CLI](./cli.md), [Core](./core/README.md).
