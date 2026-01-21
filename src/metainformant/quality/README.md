# Quality Control Module

The `quality` module provides tools for assessing and ensuring data quality across various biological data types.

## Overview

This module handles quality assessment, filtering, and validation for biological datasets including sequences, expression data, and genomic information.

### Module Architecture

```mermaid
graph TB
    subgraph "Quality Module"
        FASTQfastqFastqQuality[fastq_FASTQ Quality]
        MetricsmetricsQualityMetrics[metrics_Quality Metrics]
        ContaminationcontaminationContaminationDetection[contamination_Contamination Detection]
    end
    
    subgraph "Input Data"
        FASTQFilesfastqFiles[FASTQ Files]
        ExpressionexpressionData[Expression Data]
        Sequences[Sequences]
    end
    
    subgraph "Other Modules"
        dna[dna]
        rna[rna]
        AllallModules[All Modules]
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
    ArawFastqFiles[Raw FASTQ Files] --> Bper-baseQualityScores[Per-Base Quality Scores]
    B --> CqualityDistributionAnalysis[Quality Distribution Analysis]

    A --> DreadLengthDistribution[Read Length Distribution]
    D --> ElengthStatistics[Length Statistics]

    A --> FbaseCompositionAnalysis[Base Composition Analysis]
    F --> GgcContent&Bias[GC Content & Bias]

    A --> HadapterContentDetection[Adapter Content Detection]
    H --> IadapterTrimmingAssessment[Adapter Trimming Assessment]

    A --> JduplicateReadAnalysis[Duplicate Read Analysis]
    J --> KpcrDuplicateEstimation[PCR Duplicate Estimation]

    C --> LqualityMetricsAggregation[Quality Metrics Aggregation]
    E --> L
    G --> L
    I --> L
    K --> L

    L --> M{Overall Quality}
    M -->|High| NpassQualityControl[Pass Quality Control]
    M -->|Medium| OflagForReview[Flag for Review]
    M -->|Low| PfailQualityControl[Fail Quality Control]

    N --> QproceedToAnalysis[Proceed to Analysis]
    O --> RmanualInspection[Manual Inspection]
    P --> SdataRejection[Data Rejection]

    R --> T{Acceptable?}
    T -->|Yes| Q
    T -->|No| S


    subgraph "Quality Dimensions"
        Uper-baseQuality[Per-Base Quality] -.-> B
        Vper-sequenceQuality[Per-Sequence Quality] -.-> C
        WsequenceLength[Sequence Length] -.-> D
        XbaseContent[Base Content] -.-> F
        YsequenceDuplication[Sequence Duplication] -.-> J
    end

    subgraph "Assessment Criteria"
        ZmeanQualityScore[Mean Quality Score] -.-> M
        AAqualityScoreDistribution[Quality Score Distribution] -.-> M
        BBadapterContent%[Adapter Content %] -.-> M
        CCgcContentRange[GC Content Range] -.-> M
        DDduplicateRate[Duplicate Rate] -.-> M
    end
```

### Multi-Omics Data Quality Framework

```mermaid
graph TD
    Amulti-omicsDatasets[Multi-Omics Datasets] --> B{Data Type}
    B -->|Genomic| CvariantQualityControl[Variant Quality Control]
    B -->|Transcriptomic| DexpressionQualityControl[Expression Quality Control]
    B -->|Proteomic| EproteinQualityControl[Protein Quality Control]
    B -->|Epigenomic| FmethylationQualityControl[Methylation Quality Control]

    C --> GgenotypeCallingQuality[Genotype Calling Quality]
    C --> HmissingDataAnalysis[Missing Data Analysis]
    C --> IalleleFrequencyChecks[Allele Frequency Checks]

    D --> JlibrarySizeNormalization[Library Size Normalization]
    D --> KgeneDetectionRate[Gene Detection Rate]
    D --> LexpressionDistribution[Expression Distribution]

    E --> MpeptideIdentification[Peptide Identification]
    E --> NproteinQuantification[Protein Quantification]
    E --> OcontaminationAssessment[Contamination Assessment]

    F --> PmethylationBetaValues[Methylation Beta Values]
    F --> QdetectionP-values[Detection P-values]
    F --> RprobePerformance[Probe Performance]

    G --> SqualityMetrics[Quality Metrics]
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

    S --> Tplatform-specificQc[Platform-Specific QC]
    S --> Ucross-platformQc[Cross-Platform QC]

    T --> Vintra-platformAssessment[Intra-Platform Assessment]
    U --> Winter-platformAssessment[Inter-Platform Assessment]

    V --> Xplatform-specificFilters[Platform-Specific Filters]
    W --> YharmonizationAdjustments[Harmonization Adjustments]

    X --> ZfilteredDatasets[Filtered Datasets]
    Y --> Z

    Z --> AAintegratedQualityReport[Integrated Quality Report]


    subgraph "QC Metrics Categories"
        BBtechnicalMetrics[Technical Metrics] -.-> S
        CCbiologicalMetrics[Biological Metrics] -.-> S
        DD[Platform-Specific] -.-> T
        EE[Cross-Platform] -.-> U
    end

    subgraph "Quality Control Actions"
        FFfilterLow-qualityData[Filter Low-Quality Data] -.-> X
        GGbatchEffectCorrection[Batch Effect Correction] -.-> Y
        HH[Normalization] -.-> Y
        II[Imputation] -.-> Y
    end
```

### Assembly Quality Assessment

```mermaid
graph TD
    Agenome/transcriptomeAssembly[Genome/Transcriptome Assembly] --> Bn50Statistics[N50 Statistics]
    B --> CcontigLengthDistribution[Contig Length Distribution]

    A --> DcoverageAnalysis[Coverage Analysis]
    D --> EreadMappingStatistics[Read Mapping Statistics]

    A --> FcompletenessAssessment[Completeness Assessment]
    F --> GbuscoAnalysis[BUSCO Analysis]
    F --> HcoreGeneCoverage[Core Gene Coverage]

    A --> IcontaminationScreening[Contamination Screening]
    I --> JforeignSequenceDetection[Foreign Sequence Detection]

    A --> KstructuralValidation[Structural Validation]
    K --> LsyntenyAnalysis[Synteny Analysis]
    K --> MgeneModelQuality[Gene Model Quality]

    C --> NassemblyMetrics[Assembly Metrics]
    E --> N
    G --> N
    J --> N
    L --> N
    M --> N

    N --> O{Assembly Quality}
    O -->|High| PpublishableAssembly[Publishable Assembly]
    O -->|Medium| QimprovedAssemblyNeeded[Improved Assembly Needed]
    O -->|Low| Rre-assemblyRequired[Re-assembly Required]

    P --> SdownstreamAnalysis[Downstream Analysis]
    Q --> TgapFilling[Gap Filling]
    R --> UadditionalSequencing[Additional Sequencing]

    T --> V{Re-assembly?}
    V -->|Yes| R
    V -->|No| S


    subgraph "Assembly Statistics"
        WtotalLength[Total Length] -.-> B
        XnumberOfContigs[Number of Contigs] -.-> B
        YlargestContig[Largest Contig] -.-> B
        ZgcContent[GC Content] -.-> C
    end

    subgraph "Validation Methods"
        AAreadAlignment[Read Alignment] -.-> D
        BBreferenceComparison[Reference Comparison] -.-> F
        CCblastSearches[BLAST Searches] -.-> I
        DDorthologDetection[Ortholog Detection] -.-> K
    end
```

### Batch Effect Detection and Correction

```mermaid
graph TD
    Amulti-batchDataset[Multi-Batch Dataset] --> BprincipalComponentAnalysis[Principal Component Analysis]
    B --> CbatchEffectVisualization[Batch Effect Visualization]

    A --> DdifferentialExpressionAnalysis[Differential Expression Analysis]
    D --> Ebatch-associatedGenes[Batch-Associated Genes]

    A --> FcorrelationAnalysis[Correlation Analysis]
    F --> GbatchCorrelationPatterns[Batch Correlation Patterns]

    A --> HclusteringAnalysis[Clustering Analysis]
    H --> Ibatch-specificClusters[Batch-Specific Clusters]

    C --> JbatchEffectAssessment[Batch Effect Assessment]
    E --> J
    G --> J
    I --> J

    J --> K{Batch Effects Present?}
    K -->|Strong| LcorrectionRequired[Correction Required]
    K -->|Moderate| MoptionalCorrection[Optional Correction]
    K -->|None| NproceedToAnalysis[Proceed to Analysis]

    L --> O{Correction Method}
    M --> O

    O -->|ComBat| PempiricalBayesFramework[Empirical Bayes Framework]
    O -->|limma| QlinearMixedModels[Linear Mixed Models]
    O -->|PEER| RprobabilisticEstimation[Probabilistic Estimation]
    O -->|SVA| Q

    P --> ScorrectedData[Corrected Data]
    Q --> S
    R --> S

    S --> Tpost-correctionValidation[Post-Correction Validation]
    T --> U{Correction Successful?}

    U -->|Yes| VanalysisReady[Analysis Ready]
    U -->|No| WalternativeMethod[Alternative Method]
    W --> O

    N --> V


    subgraph "Detection Methods"
        XpcaVisualization[PCA Visualization] -.-> C
        Y[Heatmaps] -.-> C
        Z[ANOVA] -.-> E
        AAspearmanCorrelation[Spearman Correlation] -.-> G
    end

    subgraph "Correction Approaches"
        BBlocation/scaleAdjustment[Location/Scale Adjustment] -.-> P
        CC[Regression-Based] -.-> Q
        DDfactorAnalysis[Factor Analysis] -.-> R
    end
```

### Quality Control Reporting Framework

```mermaid
graph TD
    AqualityAssessmentResults[Quality Assessment Results] --> BmetricSummarization[Metric Summarization]
    B --> CqualityScoreCalculation[Quality Score Calculation]

    A --> DvisualizationGeneration[Visualization Generation]
    D --> EqualityPlots[Quality Plots]
    D --> FsummaryStatistics[Summary Statistics]

    A --> GthresholdEvaluation[Threshold Evaluation]
    G --> Hpass/failCriteria[Pass/Fail Criteria]

    C --> IqualityReport[Quality Report]
    E --> I
    F --> I
    H --> I

    I --> J{Report Format}
    J -->|HTML| KinteractiveReport[Interactive Report]
    J -->|PDF| LstaticReport[Static Report]
    J -->|JSON| MstructuredData[Structured Data]
    J -->|MultiQC| Nmulti-sampleReport[Multi-Sample Report]

    K --> Oweb-basedQcReview[Web-Based QC Review]
    L --> ParchivalDocumentation[Archival Documentation]
    M --> QprogrammaticAccess[Programmatic Access]
    N --> RbatchProcessingResults[Batch Processing Results]

    O --> SqualityAssuranceComplete[Quality Assurance Complete]
    P --> S
    Q --> S
    R --> S


    subgraph "Report Components"
        TexecutiveSummary[Executive Summary] -.-> I
        UdetailedMetrics[Detailed Metrics] -.-> I
        V[Visualizations] -.-> I
        W[Recommendations] -.-> I
    end

    subgraph "Quality Metrics"
        XpassRate%[Pass Rate %] -.-> C
        YmeanQualityScore[Mean Quality Score] -.-> C
        ZfailureReasons[Failure Reasons] -.-> H
        AAimprovementSuggestions[Improvement Suggestions] -.-> H
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
