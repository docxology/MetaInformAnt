# Specification: assays

## 🎯 Scope
Epigenome assay types sub-package (ATAC-seq, ChIP-seq, methylation).

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 4 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `atacseq.py`
- `chipseq.py`
- `methylation.py`

### Function contracts (behavior pinned by tests)

- `chipseq.find_motifs_in_peaks(peaks, genome_fasta, motif_patterns, window_size=200)`:
  real IUPAC motif scan (no simulated hits). Scans the window centred on
  each peak summit (clipped to peak bounds) on both strands; a minus-strand
  match is an interval of the forward window whose reverse complement
  matches the motif, reported at forward genomic coordinates with
  `sequence` in the motif's own (minus-strand) orientation. Peaks whose
  chromosome is missing from the FASTA and windows without a match
  contribute no occurrences; invalid IUPAC patterns raise `ValueError`.
  Returns `total_peaks_analyzed`, per-pattern `motif_counts` (peaks with
  >= 1 match), `match_counts` (total matches), and `motif_positions`.
