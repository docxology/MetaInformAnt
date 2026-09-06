# Specification: quality

## Scope

Quality control analysis module for METAINFORMANT. Sequence quality assessment,
contamination detection, batch effect analysis, and FASTQ parsing.

## Architecture

- **Dependency Level**: Domain
- **Component Type**: Source Code

## Data Structures

- **Sub-packages**: analysis, batch, io, reporting
- **Key Concepts**: QC metrics, batch effect detection/correction, FASTQ parsing

## API Definition

### Exports — `analysis/contamination.py`

- `ContaminationDetector` — Multi-type contamination screening (microbial, cross-species, adapter, duplication)
- `detect_rna_contamination` — Detect RNA bases in DNA sequencing data
- `detect_vector_contamination` — Screen for vector backbone sequences
- `detect_adapter_contamination` — Identify adapter sequence remnants
- `detect_cross_species_contamination` — Longest-matching-substring species screening
- `detect_mycoplasma_contamination` — Mycoplasma motif/genome screening
- `detect_rrna_contamination` — rRNA sequence screening
- `generate_contamination_report` — Combined contamination report with METAINFORMANT header

### Exports — `analysis/metrics.py`

- `calculate_quality_score` — Composite weighted score for FASTQ, VCF, or BAM data
- `calculate_data_integrity_score` — FASTQ integrity checks (has reads, quality range, length consistency)
- `compare_quality_metrics` — Side-by-side comparison of two datasets' scores
- `generate_quality_report` — Multi-section quality assessment report
- `batch_quality_analysis` — FASTQ quality analysis across multiple files
- `detect_outliers` — IQR, z-score, or modified z-score outlier detection
- `calculate_coverage_metrics`, `calculate_duplication_metrics`, `calculate_gc_metrics`,
  `calculate_length_metrics`, `calculate_quality_metrics`, `calculate_complexity_metrics` —
  Statistical metric helpers

### Exports — `io/fastq.py`

- `FastqRecord` — Single FASTQ read with validation, quality decoding, GC content
- `read_fastq_records` — Streaming FASTQ reader (plain or gzip via `core.io`)
- `analyze_fastq_quality` — Complete QC metric collection for a FASTQ file
- `basic_statistics`, `per_base_quality`, `per_sequence_quality`, `sequence_length_distribution`,
  `gc_content_distribution`, `adapter_content`, `overrepresented_sequences`, `duplication_levels`,
  `n_content_per_position`, `quality_score_distribution`, `filter_reads` — QC analysis helpers

### Exports — `batch/detection.py`

- `BatchEffectReport` — Dataclass summarizing batch effect detection
- `detect_batch_effects` — PVCA-style variance decomposition + silhouette score
- `correct_batch_combat` — Empirical-Bayes (ComBat-like) batch correction

### Exports — `reporting/multiqc_integration.py`

- `default_qc_thresholds` — Default warn/fail thresholds for common QC metrics
- `check_qc_thresholds` — Evaluate metrics against thresholds
- `aggregate_sample_qc` — Cross-sample metric aggregation (mean/median/std)
- `generate_qc_report` — JSON QC report from collected metrics
- `qc_trend_analysis` — Linear trend analysis of a metric over time
