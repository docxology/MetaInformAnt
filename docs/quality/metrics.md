# Quality Control Metrics

Composite quality scoring, outlier detection, data integrity checks, and
specialised metric calculators for coverage, duplication, GC content, sequence
length, quality scores, and sequence complexity.

## Key Concepts

**Quality scores** are computed as a weighted sum of component scores that vary
by data type. Each component contributes to a 0--100 overall score which is
converted to a letter grade:

| Grade | Score Range | Interpretation |
|-------|------------|----------------|
| A     | 90 -- 100  | Excellent |
| B     | 80 -- 89   | Good |
| C     | 70 -- 79   | Acceptable |
| D     | 60 -- 69   | Needs improvement |
| F     | < 60       | Poor |

**FASTQ scoring components** (default weights):

| Component | Weight | Source Metric |
|-----------|--------|---------------|
| Basic quality | 40% | Mean Phred quality score |
| Per-base quality | 30% | Quality degradation across read |
| GC content | 15% | Fraction of reads in 40--60% GC |
| Adapter content | 15% | Maximum adapter percentage |

**VCF** and **BAM** data types have analogous component breakdowns tailored to
variant quality, filter pass rates, depth, mapping quality, and duplication.

## Function Reference

### calculate_quality_score

```python
def calculate_quality_score(
    data: Dict[str, Any],
    data_type: str = "fastq",
) -> Dict[str, Any]
```

Computes `overall_score`, `total_score`, `max_possible_score`, per-`components`
breakdown, and letter `grade`. Supported data types: `"fastq"`, `"vcf"`,
`"bam"`.

### detect_outliers

```python
def detect_outliers(
    data: List[float],
    method: str = "iqr",
    threshold: float = 1.5,
) -> Dict[str, Any]
```

Detects outliers using IQR (default), z-score, or modified z-score methods.
Returns `outliers`, `outlier_indices`, `method`, `threshold`,
`outlier_percentage`.

### calculate_data_integrity_score

```python
def calculate_data_integrity_score(
    data: Dict[str, Any],
    data_type: str = "fastq",
) -> Dict[str, Any]
```

Runs structural checks (has reads, quality range, length consistency) and
returns `integrity_score` (0--100), `passed_checks`, `total_checks`, and
individual check results.

### compare_quality_metrics

```python
def compare_quality_metrics(
    dataset1: Dict[str, Any],
    dataset2: Dict[str, Any],
    data_type: str = "fastq",
) -> Dict[str, Any]
```

Computes quality scores for two datasets and compares basic statistics
(mean quality, length, GC) with absolute and percentage differences.

### generate_quality_report

```python
def generate_quality_report(
    quality_data: Dict[str, Any],
    data_type: str = "fastq",
    output_path: Optional[str | Path] = None,
) -> str
```

Formatted text report combining the quality score, integrity score, component
breakdown, basic statistics, and actionable recommendations. Optionally saves
to file.

### batch_quality_analysis

```python
def batch_quality_analysis(
    file_paths: List[str | Path],
    data_type: str = "fastq",
    n_reads: Optional[int] = None,
) -> Dict[str, Any]
```

Analyses multiple files and produces per-file results plus summary statistics
(mean, median, min, max scores across files).

### calculate_coverage_metrics

```python
def calculate_coverage_metrics(
    coverage_values: List[float],
    target_coverage: float = 30.0,
) -> Dict[str, Any]
```

Mean/median/min/max/std coverage, coverage bias and uniformity (CV),
low/high coverage percentages, breadth (>= 1x), and target achievement rate.

### calculate_duplication_metrics

```python
def calculate_duplication_metrics(
    duplication_levels: Dict[int, int],
) -> Dict[str, Any]
```

Total, unique, and duplicate read counts; duplication and uniqueness rates;
mean duplication level; and per-level distribution.

### calculate_gc_metrics

```python
def calculate_gc_metrics(gc_content: List[float]) -> Dict[str, Any]
```

Mean/median/min/max/std GC, binned distribution, and GC bias (deviation from
50%).

### calculate_length_metrics

```python
def calculate_length_metrics(lengths: List[int]) -> Dict[str, Any]
```

Mean/median/min/max/std length, most common lengths, range, and percentages
of short (< 50 bp) and long (> 1000 bp) reads.

### calculate_quality_metrics

```python
def calculate_quality_metrics(
    quality_scores: List[float],
) -> Dict[str, Any]
```

Mean/median/min/max/std quality, binned distribution, low/medium/high quality
percentages, and Q20/Q30 pass rates.

### calculate_complexity_metrics

```python
def calculate_complexity_metrics(
    sequences: List[str],
) -> Dict[str, Any]
```

Dinucleotide Shannon entropy, normalised entropy, k-mer uniqueness, composite
complexity score, sequence diversity, GC complexity, per-sequence character
complexity, and low-complexity rate.

## Usage Example

```python
from metainformant.quality.analysis.metrics import (
    calculate_quality_score,
    detect_outliers,
    generate_quality_report,
    calculate_coverage_metrics,
)

# Score FASTQ data
score = calculate_quality_score(fastq_data, data_type="fastq")
print(f"Grade: {score['grade']} ({score['overall_score']:.1f}/100)")

# Detect outlier quality values
outliers = detect_outliers(quality_values, method="zscore", threshold=3.0)

# Coverage analysis
cov = calculate_coverage_metrics(coverage, target_coverage=30.0)
print(f"Breadth: {cov['coverage_breadth']:.1f}%")

# Full report
report = generate_quality_report(fastq_data, output_path="output/quality/qc.txt")
```

## Related Modules

- `metainformant.quality.analysis.contamination` -- contamination detection
- `metainformant.quality.io.fastq` -- FASTQ parsing and per-read QC
- `metainformant.quality.reporting` -- report generation and threshold checks
