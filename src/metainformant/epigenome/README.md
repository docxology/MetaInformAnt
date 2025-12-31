# Epigenome Module

The `epigenome` module provides tools for epigenetic modification analysis, including DNA methylation, chromatin structure, and histone modifications.

## Overview

This module handles epigenetic data processing, from raw sequencing data to biological interpretation of epigenetic modifications.

### Module Architecture

```mermaid
graph TB
    subgraph "Epigenome Module"
        Methyl[methylation<br/>DNA Methylation]
        ChIP[chipseq<br/>ChIP-seq Analysis]
        ATAC[atac<br/>ATAC-seq Analysis]
        Tracks[tracks<br/>Track Processing]
        Workflow[workflow<br/>Workflow Orchestration]
    end
    
    subgraph "Input Data"
        BSSeq[Bisulfite Seq]
        ChIPData[ChIP-seq Data]
        ATACData[ATAC-seq Data]
        BedGraph[BedGraph Files]
    end
    
    subgraph "Other Modules"
        DNA_Mod[dna]
        Networks_Mod[networks]
        MultiOmics[multiomics]
    end
    
    BSSeq --> Methyl
    ChIPData --> ChIP
    ATACData --> ATAC
    BedGraph --> Tracks
    Methyl --> Workflow
    ChIP --> Workflow
    ATAC --> Workflow
    Tracks --> Workflow
    Workflow --> DNA_Mod
    Workflow --> Networks_Mod
    Workflow --> MultiOmics
```

### Epigenome Analysis Workflow

```mermaid
graph TD
    A[Raw Epigenomic Data] --> B{Data Type}
    B -->|Bisulfite-seq| C[Methylation Analysis]
    B -->|ChIP-seq| D[ChIP-seq Analysis]
    B -->|ATAC-seq| E[ATAC-seq Analysis]
    B -->|BedGraph| F[Track Processing]

    C --> G[Load CpG Data]
    G --> H[Compute Beta Values]
    H --> I[Differential Methylation]
    I --> J[mQTL Analysis]

    D --> K[Peak Calling]
    K --> L[Peak Analysis]
    L --> M[Motif Enrichment]

    E --> N[Open Chromatin Analysis]
    N --> O[Tn5 Integration Sites]
    O --> P[Accessibility Analysis]

    F --> Q[Track Loading]
    Q --> R[Genome Browser Tracks]

    J --> S[Integration]
    M --> S
    P --> S
    R --> S

    S --> T[Multi-Omics Integration]
    T --> U[Output Results]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style S fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style U fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Quality Control"
        V[Data Validation] -.-> C
        V -.-> D
        V -.-> E
        V -.-> F
    end

    subgraph "Statistical Analysis"
        W[Differential Analysis] -.-> I
        X[Enrichment Analysis] -.-> M
        Y[Peak Overlap] -.-> L
    end
```

### Methylation Analysis Pipeline

```mermaid
graph TD
    A[CpG Methylation Data] --> B[Load CpG Table]
    B --> C[Validate Data Format]
    C --> D[Compute Beta Values]
    D --> E[Beta Value Summary]

    E --> F{Differential Analysis?}
    F -->|Yes| G[Group Comparison]
    F -->|No| H[Descriptive Statistics]

    G --> I[Statistical Testing]
    I --> J[Multiple Testing Correction]
    J --> K[DMR Identification]

    H --> L[Genome-wide Patterns]
    K --> L

    L --> M{mQTL Analysis?}
    M -->|Yes| N[Genotype-Phenotype Association]
    M -->|No| O[Finalize Results]

    N --> P[mQTL Results]
    P --> O
    H --> O

    O --> Q[Output Reports]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style I fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style Q fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px
```

### ChIP-seq Analysis Pipeline

```mermaid
graph TD
    A[ChIP-seq Data] --> B[Quality Control]
    B --> C[Peak Calling]
    C --> D[Peak Filtering]

    D --> E{Peak Analysis}
    E --> F[Peak Overlap Analysis]
    E --> G[Motif Enrichment]
    E --> H[Functional Annotation]

    F --> I[Peak Set Comparison]
    G --> J[TF Binding Sites]
    H --> K[Gene Target Prediction]

    I --> L[Integration]
    J --> L
    K --> L

    L --> M[ChIP-seq Results]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style E fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style M fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px
```

## Submodules

### Methylation Analysis (`methylation.py`)
DNA methylation pattern detection and analysis.

**Key Features:**
- Methylation calling from sequencing data
- Differential methylation analysis
- Methylation quantitative trait locus (mQTL) mapping
- Methylation age prediction

**Usage:**
```python
from metainformant.epigenome import load_cpg_table, compute_beta_values, summarize_beta_by_chromosome

# Load methylation data from CPG table
cpg_df = load_cpg_table("methylation_data.tsv")

# Compute beta values (methylation levels)
cpg_df = compute_beta_values(cpg_df)

# Summarize by chromosome
summary = summarize_beta_by_chromosome(cpg_df)
```

### Chromatin Tracks (`tracks.py`)
Chromatin state and accessibility analysis.

**Key Features:**
- ATAC-seq and DNase-seq processing
- Chromatin accessibility quantification
- Peak calling and annotation
- Chromatin state segmentation

**Usage:**
```python
from metainformant.epigenome import read_bedgraph

# Load chromatin track data
bedgraph_df = read_bedgraph("chromatin_track.bedgraph")
# DataFrame columns: chromosome, start, end, value
```

## Integration with Other Modules

### With DNA Module
```python
from metainformant.dna import sequences
from metainformant.epigenome import load_cpg_table, compute_beta_values

# Epigenetic regulation of DNA sequences

# Load DNA sequences
seqs = sequences.read_fasta("sequences.fasta")

# Load and analyze methylation
cpg_df = load_cpg_table("methylation.tsv")
methylation = compute_beta_values(cpg_df)
# Analyze promoter methylation using genomic coordinates
```

### With Visualization Module
```python
from metainformant.epigenome import summarize_beta_by_chromosome
from metainformant.visualization import lineplot

# Visualize methylation patterns

# Summarize methylation by chromosome
summary = summarize_beta_by_chromosome(methylation_df)

# Visualize methylation levels
# Use chromosome/position data from summary for plotting
ax = lineplot(None, summary["mean_beta"].values)
ax.set_title("DNA Methylation Profile")
```

## Data Formats

- **Bisulfite sequencing**: Bismark, BS-Seeker output formats
- **ATAC-seq**: BED, narrowPeak, broadPeak formats
- **ChIP-seq**: BED, SAM, BAM formats
- **WGBS**: Methylation call formats

## Performance Features

- Memory-efficient processing of large epigenomic datasets
- Parallel peak calling and analysis
- Streaming processing for genome-scale data

## Testing

Comprehensive tests cover:
- Methylation calling accuracy
- Peak detection algorithms
- Integration with genomic annotations

## Dependencies

- PyBedTools for genomic interval operations
- Optional: specialized epigenomic analysis packages

## See Also

- **[AGENTS.md](AGENTS.md)**: AI agent contributions and development details for the epigenome module

## Related Modules

The Epigenome module integrates with several other METAINFORMANT modules:

- **DNA Module**: Epigenetic modifications mapped to genomic DNA sequences; CpG island identification and methylation pattern analysis
- **Networks Module**: Chromatin interaction networks, enhancer-promoter regulatory networks, and epigenetic regulation modeling
- **Multi-omics Module**: Epigenetic data integration with DNA, RNA, and protein data; regulatory analysis
- **RNA Module**: Gene expression regulation by epigenetic modifications; transcription factor binding analysis
- **Protein Module**: Histone modifications and chromatin-associated protein analysis
- **Quality Module**: Epigenetic data quality assessment and preprocessing
- **Visualization Module**: Genome browser tracks, epigenetic modification heatmaps, and chromatin state visualization
- **ML Module**: Machine learning analysis of epigenetic patterns and regulatory predictions
- **Information Module**: Information-theoretic analysis of epigenetic complexity and regulatory patterns
- **GWAS Module**: Epigenetic modifications as quantitative traits; methylation quantitative trait loci (mQTL) analysis

This module provides essential tools for epigenetic research and chromatin biology.
