# Quality Control Module

The `quality` module provides tools for assessing and ensuring data quality across various biological data types.

## Overview

This module handles quality assessment, filtering, and validation for biological datasets including sequences, expression data, and genomic information.

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
from metainformant.quality import fastq

# Quality analysis
quality_stats = fastq.analyze_quality("reads.fastq")
filtered_reads = fastq.filter_by_quality(reads, min_score=30)
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
from metainformant.quality import contamination

# Detect adapters
adapters_found = contamination.detect_adapters("reads.fastq")
cleaned_reads = contamination.remove_adapters(reads, adapters_found)

# Screen for contamination
contamination_report = contamination.screen_contamination(
    reads,
    reference_databases=["human", "mouse", "bacteria"]
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
from metainformant.quality import metrics

# Calculate quality scores
quality_scores = metrics.calculate_quality_metrics(dataset)
overall_score = metrics.compute_overall_quality_score(quality_scores)

# Detect batch effects
batch_effects = metrics.detect_batch_effects([dataset1, dataset2])
corrected_data = metrics.correct_batch_effects(data, batch_effects)
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
from metainformant.quality import validation

# Validate datasets
is_valid = validation.validate_fasta("sequences.fasta")
report = validation.generate_quality_report(dataset)
```

## Integration with Other Modules

### With DNA Module
```python
from metainformant.dna import sequences
from metainformant.quality import fastq

# Quality control for DNA sequences
seqs = sequences.read_fasta("raw_sequences.fasta")
quality_checked = fastq.validate_sequences(seqs)
```

## Performance Features

- Efficient processing of large sequence files
- Memory-conscious quality assessment
- Parallel quality checking where applicable

## Testing

Quality assessment tools are thoroughly tested with:
- Known good and bad datasets
- Performance benchmarking
- Edge case handling

## Dependencies

- Biopython for sequence quality analysis
- Optional: FastQC integration for detailed reports

This module ensures data integrity and quality across all biological data processing pipelines.
