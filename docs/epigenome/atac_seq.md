# ATAC-seq Analysis

Assay for Transposase-Accessible Chromatin sequencing analysis. Includes peak management, peak-length-based nucleosome metrics, TSS enrichment, transcription factor motif scanning of peak sequence, and multi-condition comparisons.

## Key Concepts

**ATAC-seq** uses hyperactive Tn5 transposase to probe open chromatin. Fragment sizes reflect nucleosome organization: sub-nucleosomal fragments (<147 bp) indicate nucleosome-free regions (NFRs), while mono-nucleosomal fragments (~147--294 bp) indicate positioned nucleosomes.

**TSS enrichment** measures signal accumulation around transcription start sites, serving as a key quality control metric. High TSS enrichment indicates successful capture of regulatory elements.

**Nucleosome-free regions (NFRs)** are the primary signal in ATAC-seq and correspond to active regulatory elements including promoters, enhancers, and transcription factor binding sites.

## Data Model

### `ATACPeak`

Class representing an ATAC-seq accessible region. The constructor requires a non-empty `chromosome`, `start >= 0`, `end >= start`, and, when given, `start <= summit <= end`. The p-value and q-value are not range-checked when a peak is built; saving a narrowPeak or broadPeak file rejects a value outside 0--1.

| Field | Type | Description |
|-------|------|-------------|
| `chromosome` | `str` | Chromosome name |
| `start` | `int` | Peak start position |
| `end` | `int` | Peak end position |
| `score` | `float` | Peak score |
| `strand` | `str` | DNA strand |
| `signal_value` | `float` | Signal value |
| `p_value` | `float \| None` | P-value as a probability |
| `q_value` | `float \| None` | Q-value (FDR) as a probability |
| `summit` | `int \| None` | Absolute summit position (narrowPeak files store it as an offset from `start`) |

Read-only properties: `length`, `center` and `accessibility_score` (the same as `signal_value`). Methods: `to_bed_format(format="narrowpeak")` (one line of a peak file in the given layout, named `atac_peak_<start>_<end>`), `overlaps_with(other, min_overlap=1)` and `to_dict()`.

## Function Reference

### `load_atac_peaks(path, format="narrowpeak") -> List[ATACPeak]`

Load ATAC-seq peaks from a narrowPeak, broadPeak or BED file (read as gzip when the name ends in `.gz`). The layouts, the summit offset, the `-1` placeholders and the `-log10` conversion of the p- and q-value columns are the ones described under [Peak File Layouts](chip_seq.md#peak-file-layouts). `format` is case-insensitive; any other value raises `ValueError`, and a missing file raises `metainformant.core.utils.errors.ValidationError`.

Blank, `#`, `track` and `browser` lines are ignored. A line with at least six columns but fewer than the layout needs is read as BED, and a BED line has no signal column, so its score stands in as `signal_value`. A line with fewer than six columns (three for `format="bed"`), or with a value that does not parse (including a negative p or q column other than `-1`), is skipped with a warning. A line that parses but breaks the `ATACPeak` constraints (for example a summit outside `start`--`end`) aborts the whole load with `metainformant.core.utils.errors.IOError`.

### `save_atac_peaks(peaks, path, format="narrowpeak") -> None`

Write peaks to a text file that `load_atac_peaks` reads back, in the layout given by `format` (case-insensitive; an unsupported value raises `ValueError`). The name column is `atac_peak_<start>_<end>`. p- and q-values are written as `-log10` values, a missing p-value, q-value or summit as `-1`, and a p- or q-value outside 0--1 raises `ValueError` before anything is written. See [Peak File Layouts](chip_seq.md#peak-file-layouts).

### `calculate_atac_statistics(peaks) -> Dict`

Compute summary statistics for a peak set: total count, length distribution (mean/median/min/max/std), score distribution, and per-chromosome signal breakdown.

### `calculate_atac_specific_metrics(peaks) -> Dict`

ATAC-seq-specific QC computed from peak lengths: `nfr_peak_fraction` ([50, 150) bp), `mononucleosome_peak_fraction` ([150, 250) bp), `dinucleosome_peak_fraction` ([250, 350) bp), and periodicity scores around expected nucleosome periods (147/200/300 bp). Takes a peak list, not raw fragment sizes.

### `identify_tss_enrichment(peaks, tss_positions, window_size=2000, genome_size=None) -> Dict`

Measure peak density around annotated transcription start sites. `tss_positions` maps chromosome names to lists of TSS positions. Returns `total_tss`, `enriched_tss`, `enrichment_ratio`, `expected_ratio`, `fold_enrichment`, and `window_size`. A positive `genome_size` in bp is needed for the chance expectation `1 - exp(-N * (window_size + mean_peak_length) / genome_size)`. Without it, `expected_ratio` and `fold_enrichment` are `None`; with no peaks the fold enrichment is also `None`.

### `find_tf_binding_sites(peaks, tf_motifs, genome_fasta=None, peak_sequences=None) -> Dict`

Scan the DNA of every accessible region for transcription factor motifs on both strands, and report every occurrence, including ones that overlap each other. `tf_motifs` maps TF names to consensus motifs written in IUPAC nucleotide codes (for example `"CCGCGNGGNGGCAG"` or `"TGASTCA"`); a position matches when the base is one the code allows.

The DNA must be supplied in one of two ways, and nothing is scanned or reported without it:

- `genome_fasta`: a FASTA file covering the peaks' chromosomes. Each peak is scanned over `start` to `end` of its chromosome, and a peak whose chromosome is missing from the file is not scanned.
- `peak_sequences`: one sequence per peak, in the order of `peaks`, each exactly as long as its peak and starting at `peak.start`.

Position weight matrix (PWM) scoring is not implemented: motifs are consensus strings.

Returns a dict with:

- `peaks_analyzed` (input peaks) and `peaks_scanned` (peaks that had sequence).
- `motifs_analyzed`: number of motifs.
- `motif_counts`: per TF, the number of peaks with at least one site.
- `site_counts`: per TF, the total number of sites (a palindromic motif counts once per strand).
- `expected_sites`: per TF, the number of sites expected by chance in the scanned sequence, from its base composition with independent bases; `None` when it cannot be computed (no scanned sequence, or none long enough to hold the motif).
- `motif_enrichment`: per TF, `site_counts / expected_sites`, or `None` when `expected_sites` is missing or zero.
- `<TF>_positions`: per TF, a list of site dicts with `chromosome`, `position` (0-based forward coordinate of the site start), `end`, `strand` (`"+"` or `"-"`), `sequence` (the site read in the motif's orientation), `peak_start`, `peak_end` and `peak_score`.

Raises `ValueError` if neither or both of `genome_fasta` and `peak_sequences` are given, if `peak_sequences` does not match `peaks` in number or length, or if a motif is empty or contains a character outside the IUPAC codes.

### `calculate_chromatin_accessibility_index(peaks, genomic_regions) -> Dict`

Accessibility index for specific genomic regions given as `(chromosome, start, end)` tuples.

### `compare_atac_conditions(condition1_peaks, condition2_peaks) -> Dict`

Compare accessibility between two conditions. Returns `condition1_total`, `condition2_total`, `overlapping_peaks`, condition-only counts, and `overlap_percentage`.

## Usage Examples

```python-snippet
from metainformant.epigenome.assays.atacseq import (
    ATACPeak,
    load_atac_peaks,
    calculate_atac_statistics,
    calculate_atac_specific_metrics,
    identify_tss_enrichment,
    find_tf_binding_sites,
    compare_atac_conditions,
)

# Load ATAC-seq peaks (p- and q-values arrive as probabilities)
peaks = load_atac_peaks("atac_peaks.narrowpeak")
stats = calculate_atac_statistics(peaks)
print(f"Accessible regions: {stats['total_peaks']}")

# Nucleosome-positioning QC from peak lengths
metrics = calculate_atac_specific_metrics(peaks)
print(f"NFR peak fraction: {metrics['nfr_peak_fraction']:.2%}")
print(f"Mono-nucleosomal: {metrics['mononucleosome_peak_fraction']:.2%}")

# TSS enrichment QC (chromosome -> TSS positions)
tss_positions = {"chr1": [1000, 50000], "chr2": [25000]}
tss = identify_tss_enrichment(peaks, tss_positions, window_size=2000, genome_size=3_000_000_000)
# Replace genome_size above with the length of the assembly used for these peaks.
print(f"Fold enrichment: {tss['fold_enrichment']:.2f}")

# TF motif scan over the peak sequence (IUPAC consensus motifs, both strands);
# the DNA comes from a genome FASTA, or pass peak_sequences=[...] with one sequence per peak
tf_motifs = {"CTCF": "CCGCGNGGNGGCAG", "SOX2": "CATTGTT"}
binding = find_tf_binding_sites(peaks, tf_motifs, genome_fasta="genome.fa")
print(binding["site_counts"])       # sites per TF
print(binding["motif_enrichment"])  # observed / expected sites per TF (None when undefined)

# Condition comparison
diff = compare_atac_conditions(treatment_peaks, control_peaks)
print(f"Overlapping peaks: {diff['overlapping_peaks']}")
print(f"Overlap percentage: {diff['overlap_percentage']:.1f}%")
```

## Configuration

No dedicated configuration file or environment prefix is used by the ATAC-seq module; call the functions above directly from Python.

## Related Modules

- `metainformant.epigenome.assays.chipseq` -- histone modification peaks
- `metainformant.epigenome.peak_calling` -- de novo peak calling
- `metainformant.epigenome.chromatin_state` -- chromatin state annotation
- `metainformant.epigenome.workflow` -- integrated analysis pipelines
