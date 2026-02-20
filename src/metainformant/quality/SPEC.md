# Specification: quality

## 🎯 Scope

Quality control analysis module for METAINFORMANT. Sequence quality assessment,
contamination detection, batch effect analysis, and FASTQ parsing.

## 🧱 Architecture

- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures

- **Sub-packages**: analysis, batch, io, reporting
- **Key Concepts**: QC metrics, batch effect detection/correction, FASTQ parsing

## 🔌 API Definition

### Exports — `analysis/contamination.py`

- `detect_barcode_hopping` — Detect contamination via barcode hopping
- `detect_cross_species` — Cross-species contamination detection
- `detect_index_swapping` — Index swap contamination detection
- `detect_pcr_duplicates` — PCR duplicate detection
- `detect_reagent_contamination` — Reagent contamination detection
- `detect_sample_mislabeling` — Sample mislabeling detection
- `generate_contamination_report` — Combined contamination report with METAINFORMANT header

### Exports — `io/fastq.py`

- `per_base_quality` — Per-base quality score distribution
- `FastqReader` — FASTQ file reader
