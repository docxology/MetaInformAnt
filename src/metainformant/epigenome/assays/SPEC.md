# Specification: assays

## 🎯 Scope
Epigenome assay types sub-package (ATAC-seq, ChIP-seq, methylation).

## 🧱 Architecture
- **Dependency Level**: Domain
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 5 Python modules
- **Key Concepts**: `ChIPPeak`, `ATACPeak` and `MethylationSite` are plain classes whose constructors validate their arguments.

## 🔌 API Definition
### Exports
- `__init__.py`
- `_peak_format.py` (internal, not in `__all__`: the peak-file reader and writer that `atacseq.py` and `chipseq.py` share)
- `atacseq.py`
- `chipseq.py`
- `methylation.py`

### Function contracts (behavior pinned by tests)

- Peak files (`chipseq.load_chip_peaks` / `save_chip_peaks`, `atacseq.load_atac_peaks` /
  `save_atac_peaks` and the `to_bed_format` method of both peak classes): the layouts are
  defined once in `_peak_format.py`. narrowPeak has ten columns (`chrom`, `chromStart`,
  `chromEnd`, `name`, `score`, `strand`, `signalValue`, `pValue`, `qValue`, `peak`), broadPeak
  the first nine and BED the first six (three are enough when reading). `peak` is the summit as
  an offset from `chromStart`. A missing p-value, q-value or summit is written `-1`, and `-1` or
  `.` reads back as `None`. Peaks carry p- and q-values as probabilities while files hold
  `-log10` values: a column `v` is read as `10 ** -v`, a probability `p` is written as
  `-log10(p)`, and a probability of 0 is written as the `-log10` of the smallest positive float
  (about 323.306). A negative column other than `-1` skips the line. A p-value or q-value
  outside [0, 1] raises `ValueError` when a narrowPeak or broadPeak file is saved. Layout names
  are case-insensitive and an unsupported one raises `ValueError`. Saving formats every line
  before it opens the file, so a failure writes nothing and leaves an existing file alone. The
  ATAC-seq loader reads a line with at least six columns but fewer than the layout needs as BED,
  with the score as the signal value.
- `chipseq.find_overlapping_peaks(peaks1, peaks2, min_overlap=1)`: every pair on one chromosome
  that `overlaps_with` accepts, including pairs involving peaks that are nested or overlapping
  within one set; ordered by chromosome (first appearance in `peaks1`), then `peaks1` start,
  then `peaks2` start.
- `chipseq.merge_overlapping_peaks(peaks, max_distance=0)`: the merged p-value and q-value are the
  lower of the two peaks joined, or `None` when either is missing; 0.0 is a value, not a gap.
- `chipseq.calculate_peak_enrichment(peaks, genome_size, expected_peaks=None)`: with a positive
  `expected_peaks`, `enrichment_p_value` is the Poisson upper tail summed in log space (no
  overflow for any peak count) and `enrichment_log10_p_value` is its base-10 logarithm, which
  stays finite when the p-value underflows to 0.0.
- `chipseq.find_motifs_in_peaks(peaks, genome_fasta, motif_patterns, window_size=200)`:
  real IUPAC motif scan (no simulated hits). Scans the window centred on
  each peak summit (clipped to peak bounds) on both strands; a minus-strand
  match is an interval of the forward window whose reverse complement
  matches the motif, reported at forward genomic coordinates with
  `sequence` in the motif's own (minus-strand) orientation. Every occurrence
  is reported, including one that overlaps an earlier occurrence of the same
  pattern. Peaks whose chromosome is missing from the FASTA and windows
  without a match contribute no occurrences; an empty or invalid IUPAC
  pattern raises `ValueError`.
  Returns `total_peaks_analyzed`, per-pattern `motif_counts` (peaks with
  >= 1 match), `match_counts` (total matches), and `motif_positions`.
- `atacseq.find_tf_binding_sites(peaks, tf_motifs, genome_fasta=None, peak_sequences=None)`: scans
  caller-supplied sequence for IUPAC consensus motifs on both strands and reports every
  occurrence, overlapping ones included. The sequence comes from exactly one of `genome_fasta`
  (a peak on a chromosome missing from the file is not scanned) and `peak_sequences` (one
  sequence per peak, as long as the peak); neither, both, a count or length mismatch, and an
  empty or invalid motif raise `ValueError`. Returns `peaks_analyzed`, `peaks_scanned`,
  `motifs_analyzed`, per-TF `motif_counts` (peaks with a site), `site_counts`, `expected_sites`
  (from the base composition of the scanned sequence, `None` when undefined),
  `motif_enrichment` (`site_counts / expected_sites`, `None` when undefined or zero) and
  `<TF>_positions`. There is no PWM scoring.
- `methylation.load_methylation_bedgraph(path, min_coverage=1, level_scale="auto")`: column 4 is
  read as a fraction (`"fraction"`), a percentage (`"percent"`) or, with `"auto"`, as a percentage
  when any level exceeds 1 (a warning is logged). A line whose level is outside the range of the
  scale is skipped with a warning; another `level_scale` raises `ValueError`.
- `methylation.find_differentially_methylated_regions`: the per-site test is a pooled
  two-proportion z-test with Yates' continuity correction when every expected count of the
  methylated/unmethylated by condition table is at least 10, and Fisher's exact test otherwise.
  A zero standard error does not make the p-value 1, so 0% against 100% methylation is tested;
  the p-value is set to 1 only when the table has no methylated reads or no unmethylated reads
  (or a condition has no reads), where the conditions cannot differ. Regions carry
  `mean_p_value` (the mean of the per-site p-values) and `min_p_value`, and are ordered by
  chromosome name and position.
- `methylation.find_cpg_islands(sequences, window_size=200, step_size=50, min_gc=0.5,
  min_cpg_ratio=0.6, min_length=200)`: Gardiner-Garden & Frommer islands from real sequence
  (length, GC fraction, observed/expected CpG `CpG count * length / (C count * G count)`).
  Qualifying windows are merged when they overlap or touch, and a merged stretch is kept when it is
  at least `min_length` long and meets both criteria as a whole, so boundaries resolve to about
  `window_size`. `identify_cpg_islands(..., sequences=...)` calls the islands the same way and adds
  `sites` and `avg_methylation` per island; it raises `ValueError` without `sequences`.
- `methylation.calculate_methylation_entropy`: levels are binned into five equal-width bins whose
  last bin is closed on the right, so a level of 1.0 is counted.
