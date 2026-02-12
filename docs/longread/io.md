# Long-Read I/O

The I/O module handles reading and conversion of long-read sequencing file formats from Oxford Nanopore Technologies (ONT) and Pacific Biosciences (PacBio) platforms. It supports FAST5/POD5 signal files, BAM alignment files, and format conversions.

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
    signal: Any                 # Raw electrical signal (pA or ADC)
    sequence: str               # Basecalled sequence
    quality_string: str         # Per-base quality string
    channel_id: int             # Sequencing channel number
    mux: int                    # Pore selection number
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

## Function Reference

### FAST5/POD5

```python
def read_fast5(filepath: str | Path) -> list[Fast5Read]
def extract_signal(filepath: str | Path, read_id: str | None = None) -> Any
def extract_basecalls(filepath: str | Path) -> list[dict[str, Any]]
def get_read_metadata(filepath: str | Path) -> list[dict[str, Any]]
```

`read_fast5` parses a FAST5/POD5 file and returns a list of `Fast5Read` objects. `extract_signal` returns raw signal arrays. `extract_basecalls` returns basecalled sequences. `get_read_metadata` returns per-read metadata dictionaries.

### BAM

```python
def read_long_read_bam(filepath: str | Path, region: str | None = None) -> list[LongReadAlignment]
def extract_methylation_tags(alignment: LongReadAlignment) -> dict[str, Any]
def get_supplementary_alignments(filepath: str | Path, read_name: str) -> list[LongReadAlignment]
def calculate_alignment_stats(alignments: list[LongReadAlignment]) -> dict[str, Any]
```

`read_long_read_bam` parses a BAM file and returns alignments. `extract_methylation_tags` parses MM/ML tags for modified base data. `get_supplementary_alignments` retrieves split alignments for a read (useful for SV detection). `calculate_alignment_stats` computes summary statistics.

### Format Conversion

```python
def fast5_to_fastq(fast5_path: str | Path, output_path: str | Path,
                   min_quality: float = 0.0, gzip_output: bool = False) -> int
def convert_pod5_to_fast5(pod5_path: str | Path, output_path: str | Path) -> int
def write_paf(overlaps: list[dict], output_path: str | Path) -> int
```

`fast5_to_fastq` extracts basecalls to FASTQ format with optional quality filtering. `convert_pod5_to_fast5` converts POD5 to FAST5 format. `write_paf` writes overlap records in PAF format.

## Usage Examples

```python
from metainformant.longread import io

# Read FAST5 file
reads = io.read_fast5("path/to/reads.fast5")
for read in reads:
    print(f"{read.read_id}: {len(read.sequence)} bp, channel {read.channel_id}")

# Extract raw signal for a specific read
signal = io.extract_signal("path/to/reads.fast5", read_id="read_001")

# Read long-read BAM file
alignments = io.read_long_read_bam("path/to/aligned.bam", region="chr1:1000000-2000000")
stats = io.calculate_alignment_stats(alignments)
print(f"Mean MAPQ: {stats['mean_mapq']:.1f}, Mapped: {stats['mapped_fraction']:.2%}")

# Extract methylation data from BAM tags
for aln in alignments:
    meth = io.extract_methylation_tags(aln)
    if meth:
        print(f"Methylation sites: {len(meth.get('positions', []))}")

# Convert FAST5 to FASTQ
n_reads = io.fast5_to_fastq("reads.fast5", "output/reads.fastq.gz", gzip_output=True)
print(f"Converted {n_reads} reads")
```

## Configuration

- **Environment prefix**: `LR_`
- **Optional dependencies**: h5py (FAST5), pod5 (POD5), pysam (BAM), numpy
- FAST5 reading supports both single-read and multi-read FAST5 formats
- BAM reading supports indexed BAM files with region queries

## Related Modules

- `longread.quality` -- Quality metrics and filtering on loaded reads
- `longread.analysis` -- Modified base and structural variant analysis
- `longread.assembly` -- Assembly using read overlaps
- `longread.methylation` -- Signal-level methylation calling
