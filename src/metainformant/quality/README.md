# Quality Control Module

The `quality` module provides tools for assessing and ensuring data quality across various biological data types.

## Overview

This module handles quality assessment, filtering, and validation for biological datasets including sequences, expression data, and genomic information.

### Module Architecture

```mermaid
graph TB
    subgraph "Quality Module"
        FASTQ[fastq<br/>FASTQ Quality]
        Metrics[metrics<br/>Quality Metrics]
        Contamination[contamination<br/>Contamination Detection]
    end
    
    subgraph "Input Data"
        FASTQFiles[FASTQ Files]
        Expression[Expression Data]
        Sequences[Sequences]
    end
    
    subgraph "Other Modules"
        DNA_Mod[dna]
        RNA_Mod[rna]
        All[All Modules]
    end
    
    FASTQFiles --> FASTQ
    Expression --> Metrics
    Sequences --> Metrics
    FASTQ --> Contamination
    Metrics --> Contamination
    DNA_Mod --> FASTQ
    RNA_Mod --> Metrics
    All --> Metrics
```

### FASTQ Quality Assessment Pipeline

```mermaid
graph TD
    A[Raw FASTQ Files] --> B[Per-Base Quality Scores]
    B --> C[Quality Distribution Analysis]

    A --> D[Read Length Distribution]
    D --> E[Length Statistics]

    A --> F[Base Composition Analysis]
    F --> G[GC Content & Bias]

    A --> H[Adapter Content Detection]
    H --> I[Adapter Trimming Assessment]

    A --> J[Duplicate Read Analysis]
    J --> K[PCR Duplicate Estimation]

    C --> L[Quality Metrics Aggregation]
    E --> L
    G --> L
    I --> L
    K --> L

    L --> M{Overall Quality}
    M -->|High| N[Pass Quality Control]
    M -->|Medium| O[Flag for Review]
    M -->|Low| P[Fail Quality Control]

    N --> Q[Proceed to Analysis]
    O --> R[Manual Inspection]
    P --> S[Data Rejection]

    R --> T{Acceptable?}
    T -->|Yes| Q
    T -->|No| S

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style L fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style Q fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Quality Dimensions"
        U[Per-Base Quality] -.-> B
        V[Per-Sequence Quality] -.-> C
        W[Sequence Length] -.-> D
        X[Base Content] -.-> F
        Y[Sequence Duplication] -.-> J
    end

    subgraph "Assessment Criteria"
        Z[Mean Quality Score] -.-> M
        AA[Quality Score Distribution] -.-> M
        BB[Adapter Content %] -.-> M
        CC[GC Content Range] -.-> M
        DD[Duplicate Rate] -.-> M
    end
```

### Multi-Omics Data Quality Framework

```mermaid
graph TD
    A[Multi-Omics Datasets] --> B{Data Type}
    B -->|Genomic| C[Variant Quality Control]
    B -->|Transcriptomic| D[Expression Quality Control]
    B -->|Proteomic| E[Protein Quality Control]
    B -->|Epigenomic| F[Methylation Quality Control]

    C --> G[Genotype Calling Quality]
    C --> H[Missing Data Analysis]
    C --> I[Allele Frequency Checks]

    D --> J[Library Size Normalization]
    D --> K[Gene Detection Rate]
    D --> L[Expression Distribution]

    E --> M[Peptide Identification]
    E --> N[Protein Quantification]
    E --> O[Contamination Assessment]

    F --> P[Methylation Beta Values]
    F --> Q[Detection P-values]
    F --> R[Probe Performance]

    G --> S[Quality Metrics]
    H --> S
    I --> S
    J --> S
    K --> S
    L --> S
    M --> S
    N --> S
    O --> S
    P --> S
    Q --> S
    R --> S

    S --> T[Platform-Specific QC]
    S --> U[Cross-Platform QC]

    T --> V[Intra-Platform Assessment]
    U --> W[Inter-Platform Assessment]

    V --> X[Platform-Specific Filters]
    W --> Y[Harmonization Adjustments]

    X --> Z[Filtered Datasets]
    Y --> Z

    Z --> AA[Integrated Quality Report]

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style S fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style AA fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "QC Metrics Categories"
        BB[Technical Metrics] -.-> S
        CC[Biological Metrics] -.-> S
        DD[Platform-Specific] -.-> T
        EE[Cross-Platform] -.-> U
    end

    subgraph "Quality Control Actions"
        FF[Filter Low-Quality Data] -.-> X
        GG[Batch Effect Correction] -.-> Y
        HH[Normalization] -.-> Y
        II[Imputation] -.-> Y
    end
```

### Assembly Quality Assessment

```mermaid
graph TD
    A[Genome/Transcriptome Assembly] --> B[N50 Statistics]
    B --> C[Contig Length Distribution]

    A --> D[Coverage Analysis]
    D --> E[Read Mapping Statistics]

    A --> F[Completeness Assessment]
    F --> G[BUSCO Analysis]
    F --> H[Core Gene Coverage]

    A --> I[Contamination Screening]
    I --> J[Foreign Sequence Detection]

    A --> K[Structural Validation]
    K --> L[Synteny Analysis]
    K --> M[Gene Model Quality]

    C --> N[Assembly Metrics]
    E --> N
    G --> N
    J --> N
    L --> N
    M --> N

    N --> O{Assembly Quality}
    O -->|High| P[Publishable Assembly]
    O -->|Medium| Q[Improved Assembly Needed]
    O -->|Low| R[Re-assembly Required]

    P --> S[Downstream Analysis]
    Q --> T[Gap Filling]
    R --> U[Additional Sequencing]

    T --> V{Re-assembly?}
    V -->|Yes| R
    V -->|No| S

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style N fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style S fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Assembly Statistics"
        W[Total Length] -.-> B
        X[Number of Contigs] -.-> B
        Y[Largest Contig] -.-> B
        Z[GC Content] -.-> C
    end

    subgraph "Validation Methods"
        AA[Read Alignment] -.-> D
        BB[Reference Comparison] -.-> F
        CC[BLAST Searches] -.-> I
        DD[Ortholog Detection] -.-> K
    end
```

### Batch Effect Detection and Correction

```mermaid
graph TD
    A[Multi-Batch Dataset] --> B[Principal Component Analysis]
    B --> C[Batch Effect Visualization]

    A --> D[Differential Expression Analysis]
    D --> E[Batch-Associated Genes]

    A --> F[Correlation Analysis]
    F --> G[Batch Correlation Patterns]

    A --> H[Clustering Analysis]
    H --> I[Batch-Specific Clusters]

    C --> J[Batch Effect Assessment]
    E --> J
    G --> J
    I --> J

    J --> K{Batch Effects Present?}
    K -->|Strong| L[Correction Required]
    K -->|Moderate| M[Optional Correction]
    K -->|None| N[Proceed to Analysis]

    L --> O{Correction Method}
    M --> O

    O -->|ComBat| P[Empirical Bayes Framework]
    O -->|limma| Q[Linear Mixed Models]
    O -->|PEER| R[Probabilistic Estimation]
    O -->|SVA| Q

    P --> S[Corrected Data]
    Q --> S
    R --> S

    S --> T[Post-Correction Validation]
    T --> U{Correction Successful?}

    U -->|Yes| V[Analysis Ready]
    U -->|No| W[Alternative Method]
    W --> O

    N --> V

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style J fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style V fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Detection Methods"
        X[PCA Visualization] -.-> C
        Y[Heatmaps] -.-> C
        Z[ANOVA] -.-> E
        AA[Spearman Correlation] -.-> G
    end

    subgraph "Correction Approaches"
        BB[Location/Scale Adjustment] -.-> P
        CC[Regression-Based] -.-> Q
        DD[Factor Analysis] -.-> R
    end
```

### Quality Control Reporting Framework

```mermaid
graph TD
    A[Quality Assessment Results] --> B[Metric Summarization]
    B --> C[Quality Score Calculation]

    A --> D[Visualization Generation]
    D --> E[Quality Plots]
    D --> F[Summary Statistics]

    A --> G[Threshold Evaluation]
    G --> H[Pass/Fail Criteria]

    C --> I[Quality Report]
    E --> I
    F --> I
    H --> I

    I --> J{Report Format}
    J -->|HTML| K[Interactive Report]
    J -->|PDF| L[Static Report]
    J -->|JSON| M[Structured Data]
    J -->|MultiQC| N[Multi-Sample Report]

    K --> O[Web-Based QC Review]
    L --> P[Archival Documentation]
    M --> Q[Programmatic Access]
    N --> R[Batch Processing Results]

    O --> S[Quality Assurance Complete]
    P --> S
    Q --> S
    R --> S

    style A fill:#e1f5fe,stroke:#01579b,stroke-width:2px
    style I fill:#fff3e0,stroke:#e65100,stroke-width:2px
    style S fill:#e8f5e8,stroke:#1b5e20,stroke-width:2px

    subgraph "Report Components"
        T[Executive Summary] -.-> I
        U[Detailed Metrics] -.-> I
        V[Visualizations] -.-> I
        W[Recommendations] -.-> I
    end

    subgraph "Quality Metrics"
        X[Pass Rate %] -.-> C
        Y[Mean Quality Score] -.-> C
        Z[Failure Reasons] -.-> H
        AA[Improvement Suggestions] -.-> H
    end
```

## Submodules

### Sequence Quality (`fastq.py`)
FASTQ format quality analysis and filtering.

**Key Features:**
- Per-base quality score analysis
- Read trimming and filtering
- Quality report generation
- Format validation and conversion

**Usage:**
```python
from metainformant.quality import analyze_fastq_quality, per_base_quality
from metainformant.dna.fastq import iter_fastq, FastqRecord

# Quality analysis
# First, read FASTQ records
reads = list(FastqRecord(read_id, seq, qual) for read_id, seq, qual in iter_fastq("reads.fastq"))

# Analyze quality
quality_stats = analyze_fastq_quality(reads)
per_base = per_base_quality(reads)
```

### Contamination Detection (`contamination.py`)
Sequence contamination detection and removal.

**Key Features:**
- Adapter sequence detection and trimming
- Cross-species contamination screening
- Primer dimer identification
- Vector and plasmid contamination detection

**Usage:**
```python
from metainformant.quality import (
    detect_adapter_contamination,
    detect_cross_species_contamination,
    generate_contamination_report
)
from metainformant.dna.fastq import iter_fastq

# Detect adapters
reads = list(iter_fastq("reads.fastq"))
adapter_results = detect_adapter_contamination(reads)

# Screen for cross-species contamination
contamination_results = detect_cross_species_contamination(
    reads,
    reference_species="human"
)

# Generate contamination report
report = generate_contamination_report(
    adapter_results,
    cross_species_results=contamination_results
)
```

### Quality Metrics (`metrics.py`)
Comprehensive quality scoring and assessment.

**Key Features:**
- Quality score calculations
- Completeness and accuracy metrics
- Batch effect detection
- Statistical quality summaries

**Usage:**
```python
from metainformant.quality import (
    calculate_quality_metrics,
    generate_quality_report,
    calculate_gc_metrics,
    calculate_length_metrics
)

# Calculate quality metrics (from quality scores: list of lists)
quality_scores = [[33, 35, 32, 34], [34, 36, 33, 35]]
metrics = calculate_quality_metrics(quality_scores)

# Calculate GC metrics
gc_values = [0.45, 0.50, 0.48, 0.52]
gc_metrics = calculate_gc_metrics(gc_values)

# Calculate length metrics
lengths = [100, 150, 120, 140]
length_metrics = calculate_length_metrics(lengths)

# Generate quality report
report = generate_quality_report({
    "quality": metrics,
    "gc": gc_metrics,
    "length": length_metrics
})
```

### Data Validation (`validation.py`)
General data validation and quality assessment.

**Key Features:**
- Format validation
- Completeness checking
- Consistency verification
- Metadata validation

**Usage:**
```python
from metainformant.quality import generate_quality_report
from metainformant.quality import calculate_quality_metrics

# For sequence validation, use DNA module
from metainformant.dna import sequences
seqs = sequences.read_fasta("sequences.fasta")

# Generate quality reports using quality module functions
# See metrics examples above
```

## Integration with Other Modules

### With DNA Module
```python
from metainformant.dna import sequences
from metainformant.quality import analyze_fastq_quality

# Quality control for DNA sequences
seqs = sequences.read_fasta("raw_sequences.fasta")

# Use quality module functions for FASTQ quality assessment
from metainformant.dna.fastq import iter_fastq
reads = list(iter_fastq("reads.fastq"))
quality_stats = analyze_fastq_quality(reads)
```

### With RNA Module Workflows
```python
from metainformant.quality import analyze_fastq_quality, detect_contamination
from metainformant.dna.fastq import iter_fastq

# Quality control for RNA-seq FASTQ files
# Quality assessment before downstream RNA analysis
rna_reads = list(iter_fastq("rna_reads.fastq"))
quality_stats = analyze_fastq_quality(rna_reads)

# Contamination detection for RNA-seq
contamination = detect_adapter_contamination(rna_reads)

# Filter reads based on quality before RNA workflow
# Pass quality-filtered reads to RNA quantification pipeline
```

### With GWAS Workflows
```python
from metainformant.quality import analyze_fastq_quality
from metainformant.gwas.quality import apply_qc_filters, parse_vcf_full

# Quality control for variant calling input
# FASTQ quality assessment before variant calling
genomic_reads = list(iter_fastq("genomic_reads.fastq"))
quality_stats = analyze_fastq_quality(genomic_reads)

# After variant calling, apply GWAS-specific QC filters
vcf_data = parse_vcf_full("variants.vcf")
qc_variants = apply_qc_filters(vcf_data, maf_threshold=0.05, missing_threshold=0.1)

# Integrated quality control pipeline for GWAS
```

## Performance Features

- Efficient processing of large sequence files
- Memory-conscious quality assessment
- Parallel quality checking where applicable

## Testing

Quality assessment tools are tested with:
- Known good and bad datasets
- Performance benchmarking
- Edge case handling

## Dependencies

- Biopython for sequence quality analysis
- Optional: FastQC integration for detailed reports

## See Also

- **[AGENTS.md](AGENTS.md)**: AI agent contributions and development details for the quality module

## Related Modules

The Quality module integrates with all other METAINFORMANT modules:

- **DNA Module**: DNA sequence quality assessment, FASTQ preprocessing, and variant calling validation
- **RNA Module**: RNA-seq data quality control, adapter trimming, and expression data validation
- **Protein Module**: Protein structure validation and sequence quality assessment
- **Epigenome Module**: Epigenetic data quality control and preprocessing validation
- **GWAS Module**: Genotype data quality control and population structure validation
- **Multi-omics Module**: Cross-platform data quality assessment and integration validation
- **Single-cell Module**: Single-cell data quality control and preprocessing standards
- **Networks Module**: Network data validation and quality assessment
- **ML Module**: Model validation and performance quality assessment
- **Information Module**: Information quality measures and data integrity validation
- **Visualization Module**: Quality metric visualization and validation reporting
- **Simulation Module**: Synthetic data quality validation and method benchmarking
- **Math Module**: Statistical quality control and validation methods
- **Ontology Module**: Annotation quality assessment and validation
- **Phenotype Module**: Phenotype data validation and measurement quality control
- **Life Events Module**: Life course data quality assessment and validation
- **Ecology Module**: Ecological data validation and biodiversity assessment standards

This module ensures data integrity and quality across all biological data processing pipelines.
