# Long-Read I/O

The I/O module handles reading and conversion of long-read sequencing file formats from Oxford Nanopore Technologies (ONT) and Pacific Biosciences (PacBio) platforms. It supports FAST5/POD5 signal files, BAM alignment files, and format conversions.

The package exposes three submodules: `fast5` (signal files), `bam` (alignments) and `formats` (conversion). Import the submodule you need (`from metainformant.longread.io import bam`); the functions are not re-exported on `metainformant.longread.io` itself.

## Key Concepts

### FAST5 Format

FAST5 is an HDF5-based file format used by Oxford Nanopore for storing raw electrical signal data, basecalled sequences, and per-read metadata. Multi-read FAST5 files contain multiple reads in a single file.

### POD5 Format

POD5 is the newer ONT format replacing FAST5, offering improved read/write performance and compression. The module supports conversion from POD5 to FAST5 for compatibility.

### BAM Format

Long-read BAM files contain aligned reads with additional tags for methylation (MM/ML tags), supplementary alignments (for structural variant detection), and other long-read-specific information.

### PAF Format

Pairwise Alignment Format, a lightweight tab-separated format for read overlaps used in assembly pipelines (minimap2 output).

## Data Structures

### Fast5Read

```python
@dataclass
class Fast5Read:
    read_id: str                # Unique read identifier
    signal: Any                 # Raw ADC signal (numpy float32 array when numpy is installed)
    sequence: str | None        # Basecalled sequence (None when the file has no basecalls)
    quality_string: str | None  # Per-base quality string (None when the file has no basecalls)
    channel_id: int             # Sequencing channel number
    mux: int                    # Pore selection number: start_mux attribute (0 if absent); pore well for POD5
    start_time: int             # Start time in samples
    duration: int               # Duration in samples
    sampling_rate: float        # ADC sampling rate (Hz)
    run_id: str                 # Sequencing run identifier
    digitisation: float         # ADC digitisation value
    offset: float               # Signal offset for pA conversion
    range_value: float          # Signal range for pA conversion
    metadata: dict[str, Any]    # Additional metadata
```

### LongReadAlignment

```python
@dataclass
class LongReadAlignment:
    read_name: str              # Query read name
    query_sequence: str         # Aligned read sequence
    query_qualities: list[int]  # Per-base quality scores
    reference_name: str         # Reference contig name
    reference_start: int        # 0-based start on reference
    reference_end: int          # 0-based end on reference
    mapping_quality: int        # MAPQ score
    cigar_string: str           # CIGAR alignment string
    cigar_tuples: list[tuple[int, int]]
    is_supplementary: bool
    is_secondary: bool
    is_reverse: bool
    is_unmapped: bool
    tags: dict                  # BAM tags
    methylation_tags: dict      # Parsed MM/ML methylation data
    query_length: int
    aligned_pairs: list         # (query_pos, ref_pos) pairs
```

### AlignmentStats

```python
@dataclass
class AlignmentStats:
    total_reads: int            # Alignments passed in
    mapped_reads: int           # Alignments that are not unmapped
    unmapped_reads: int
    supplementary_count: int    # Mapped alignments flagged supplementary
    secondary_count: int        # Mapped alignments flagged secondary
    mean_mapping_quality: float # Mean MAPQ over mapped alignments
    mean_read_length: float     # Mean query length
    mean_aligned_length: float  # Mean reference-consuming length
    mean_identity: float        # Mean identity from the NM tag or the CIGAR
    total_bases: int
    total_aligned_bases: int
    coverage_depth: float       # Aligned bases / reference_length (0.0 when not given)
```

## Function Reference

### FAST5/POD5

```python-snippet
def read_fast5(filepath: str | Path) -> list[Fast5Read]
def extract_signal(fast5_data: Fast5Read) -> Any
def extract_basecalls(fast5_data: Fast5Read) -> tuple[str | None, str | None]
def get_read_metadata(fast5_data: Fast5Read) -> dict[str, Any]
```

`read_fast5` parses a `.fast5` (needs h5py) or `.pod5` (needs pod5) file, chosen by file extension, and returns a list of `Fast5Read` objects; it raises `FileNotFoundError` for a missing file and `ValueError` for any other extension. The other three functions take one parsed `Fast5Read`, not a file path. `extract_signal` converts the raw ADC trace to picoamperes, `(raw + offset) * range / digitisation`, and returns a float32 numpy array (a plain list when numpy is not installed, `None` when the read has no signal). `extract_basecalls` returns the `(sequence, quality_string)` pair, either of which may be `None`. `get_read_metadata` returns a dictionary with `read_id`, `channel_id`, `mux`, `start_time_samples`, `duration_samples`, `duration_seconds`, `sampling_rate`, `run_id`, `has_signal`, `signal_length`, `has_basecalls`, `sequence_length` and a `calibration` dict (`digitisation`, `offset`, `range`), plus whatever the read carries in its `metadata`.

### BAM

```python-snippet
def read_long_read_bam(
    filepath: str | Path, region: str | None = None,
    include_supplementary: bool = True, include_secondary: bool = False,
    min_mapping_quality: int = 0,
) -> list[LongReadAlignment]
def extract_methylation_tags(alignment: LongReadAlignment) -> dict[str, Any]
def get_supplementary_alignments(alignment: LongReadAlignment) -> list[dict[str, Any]]
def calculate_alignment_stats(
    alignments: Sequence[LongReadAlignment], reference_length: int | None = None,
) -> AlignmentStats
```

`read_long_read_bam` parses a BAM file with pysam and returns alignments; a `region` query needs a `.bai` index, and secondary alignments are skipped unless `include_secondary` is set. `extract_methylation_tags` returns one entry per base, strand and modification code as written in MM: `C+m` and `C-m` are separate keys, while `C+mh` produces `C+m` and `C+h`.

Each entry contains:

- `base` (the letter counted), `strand`, `modification`, `modified_base` (the counted letter for `+`, its complement for `-`) and `implicit`.
- `read_positions`: ascending 0-based positions in the stored `query_sequence`.
- `genomic_positions`: `(reference_position_or_None, probability_or_None)` for each listed site, in read-position order. An insertion, soft clip or unmapped base has no reference position; a missing or unusable ML tag leaves probability unknown.
- `probabilities`: usable ML values scaled by `value / 255`, in read-position order; empty without usable ML values.
- `skipped_read_positions` and `skipped_genomic_positions`: counted bases absent from the listed MM positions. Their status is unmodified when `implicit=True` (no flag or `.`) and unknown for the `?` flag.
- `n_unplaced`: listed sites that could not be placed because the read lacks sequence or a skip runs past its counted bases; a warning reports them.

MM skip counts refer to the read as sequenced. Reverse-aligned reads are reverse-complemented for counting, then coordinates are mapped back to the stored sequence. Every entry counts its written fundamental base, including `-` entries; the strand describes the modified base relative to that counted strand. Multi-code entries interleave ML values by position and code. The decoder accepts `ML` as a list, tuple, bytes, numpy array or pysam's `array('B')`, and reads draft tag names `Mm`/`Ml` when the uppercase tags are absent. Malformed MM tags and inconsistent ML lengths produce warnings rather than invented positions or probabilities. An optional `MN` tag must be a nonnegative integer matching the stored sequence length; stale or invalid lengths raise `ValueError` before either raw or cached modification calls are used.

`get_supplementary_alignments` parses the `SA` tag and returns `reference_name`, `reference_start` (0-based), `strand`, `cigar`, `mapping_quality`, `edit_distance`, `aligned_length` and `is_reverse` per split segment. `calculate_alignment_stats` returns `AlignmentStats`: `total_reads` counts all alignments, while mean fields, `total_bases` and `total_aligned_bases` cover mapped alignments.

### Format Conversion

```python-snippet
def fast5_to_fastq(fast5_path: str | Path, output_path: str | Path,
                   min_quality: float = 0.0, gzip_output: bool = False) -> int
def convert_pod5_to_fast5(pod5_path: str | Path, output_path: str | Path) -> int
def write_paf(alignments: Sequence[dict[str, Any]], output_path: str | Path) -> int
```

`fast5_to_fastq` extracts basecalls to FASTQ format with optional quality filtering (reads without basecalls are skipped, and `min_quality` is a minimum mean Phred score). `convert_pod5_to_fast5` converts POD5 to FAST5 format (needs pod5 and h5py). `write_paf` writes alignment dictionaries in PAF format; each needs `query_name`, `query_length`, `query_start`, `query_end`, `target_name`, `target_length`, `target_start`, `target_end`, and may carry `strand`, `matches`, `block_length`, `mapping_quality` and a `tags` dict. All three return the number of records written.

## Usage Examples

```python
from metainformant.longread.io import bam, fast5, formats

# Read FAST5 file
reads = fast5.read_fast5("path/to/reads.fast5")
for read in reads:
    print(f"{read.read_id}: {len(read.sequence or '')} bp, channel {read.channel_id}")

# Convert the raw signal of one read to picoamperes
signal_pa = fast5.extract_signal(reads[0])
meta = fast5.get_read_metadata(reads[0])
print(f"{meta['signal_length']} samples over {meta['duration_seconds']:.1f} s")

# Read long-read BAM file (needs pysam and an index for region queries)
alignments = bam.read_long_read_bam("path/to/aligned.bam", region="chr1:1000000-2000000")
stats = bam.calculate_alignment_stats(alignments)
print(f"Mean MAPQ: {stats.mean_mapping_quality:.1f}, Mapped: {stats.mapped_reads / stats.total_reads:.2%}")

# Extract methylation data from BAM tags
for aln in alignments:
    meth = bam.extract_methylation_tags(aln)
    for modification, calls in meth.items():
        print(f"{aln.read_name} {modification}: {len(calls['genomic_positions'])} sites")

# Convert FAST5 to FASTQ
n_reads = formats.fast5_to_fastq("reads.fast5", "output/reads.fastq.gz", gzip_output=True)
print(f"Converted {n_reads} reads")
```

## Configuration

- **Environment prefix**: `LR_` is honoured only by `workflow.pipelines.load_pipeline_config` and `utils.batch.process_batch` (`LR_BATCH_WORKERS`); the functions on this page read no environment variables
- **Optional dependencies**: h5py (FAST5), pod5 (POD5), pysam (BAM), numpy
- FAST5 reading supports both single-read and multi-read FAST5 formats
- BAM reading supports indexed BAM files with region queries; without a region the file is read in order and needs no index

## Related Modules

- `longread.quality` -- Quality metrics and filtering on loaded reads
- `longread.analysis` -- Modified base and structural variant analysis
- `longread.assembly` -- Assembly using read overlaps
- `longread.methylation` -- Signal-level methylation calling
