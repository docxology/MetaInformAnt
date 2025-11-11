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
flowchart LR
    Start[Raw Data] --> Type{Data Type?}
    Type -->|Methylation| Methyl[Process Methylation]
    Type -->|ChIP-seq| ChIP[Process ChIP-seq]
    Type -->|ATAC-seq| ATAC[Process ATAC-seq]
    Methyl --> QC[Quality Control]
    ChIP --> QC
    ATAC --> QC
    QC --> Analyze[Analysis]
    Analyze --> Integrate[Integration]
    Integrate --> Output[Results]
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

This module provides essential tools for epigenetic research and chromatin biology.
