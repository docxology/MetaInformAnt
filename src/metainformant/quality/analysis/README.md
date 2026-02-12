# Quality Analysis

Quality scoring, contamination detection, and outlier analysis for sequencing data (FASTQ, VCF, BAM).

## Contents

| File | Purpose |
|------|---------|
| `contamination.py` | Contamination detection: RNA, vector, adapter, cross-species, mycoplasma |
| `metrics.py` | Quality scoring, outlier detection, integrity checks, batch analysis |

## Key Classes and Functions

| Symbol | Description |
|--------|-------------|
| `ContaminationDetector` | Class for multi-type contamination screening |
| `detect_rna_contamination()` | Check DNA sequences for RNA contamination |
| `detect_vector_contamination()` | Screen for vector backbone sequences |
| `detect_adapter_contamination()` | Identify adapter sequence remnants |
| `generate_contamination_report()` | Full contamination screening report |
| `calculate_quality_score()` | Composite quality score for FASTQ, VCF, or BAM |
| `detect_outliers()` | Statistical outlier detection (IQR or z-score) |
| `calculate_data_integrity_score()` | Data completeness and consistency score |
| `generate_quality_report()` | Multi-section quality assessment report |
| `batch_quality_analysis()` | Run quality analysis across multiple samples |

## Usage

```python
from metainformant.quality.analysis.metrics import calculate_quality_score, detect_outliers
from metainformant.quality.analysis.contamination import ContaminationDetector

score = calculate_quality_score(fastq_data, data_type="fastq")
outliers = detect_outliers(quality_values, method="iqr")
detector = ContaminationDetector(sequences)
```
