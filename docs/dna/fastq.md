### DNA: FASTQ

FASTQ utilities with gzip support (a `.gz` suffix is decompressed transparently). These helper functions are designed for local inspection and lightweight QC. For production-grade QC and preprocessing, see the ecosystem notes below.

#### API

- `iter_fastq(fastq_path) -> Iterator[tuple[str, str, str]]`: Stream `(read_id, sequence, quality)` from a FASTQ file (`.gz` supported); `read_id` is the first word of the header. Raises `ValueError` on an incomplete record (including a truncated final one) and does not check that sequence and quality lengths agree.
- `average_phred_by_position(fastq_path) -> dict[int, float]`: Average Phred+33 quality at each 0-based read position. Positions past the end of shorter reads are averaged over the reads that reach them. It keeps every quality value in memory.
- `gc_content(seq) -> float`: GC fraction among A/C/G/T (N-insensitive); 0.0 for an empty sequence.
- `summarize_fastq(fastq_path) -> dict`: Lightweight stats: `total_reads`, `total_bases`, `min_length`, `max_length`, `mean_length`. Streams the file.
- `FastqRecord(header, sequence, quality_string)`: Lightweight record with `quality_scores` (Phred integers), `mean_quality`, `median_quality`, `min_quality`, `max_quality`, `gc_content()` and `trim_low_quality(min_quality=20)`. Wrap rows from `iter_fastq` in it to get the first few records as objects.

The module also provides `read_fastq` (every read as `{read_id: (sequence, quality)}`), `write_fastq`, `assess_quality`, `calculate_per_base_quality`, `filter_reads` (yields the 4-line records whose mean Phred is at least `min_quality`), `trim_reads` and `convert_fastq_to_fasta`. `filter_reads` streams; `read_fastq`, `assess_quality`, `calculate_per_base_quality`, `trim_reads` and `convert_fastq_to_fasta` load every read into memory.

#### Usage

```python-snippet
from itertools import islice
from pathlib import Path

from metainformant.dna.io import fastq

# A tiny file to work with; any .fastq or .fastq.gz path works the same way
Path("reads.fastq").write_text("@r1\nACGTN\n+\nIIII#\n@r2\nGGCC\n+\nIIHH\n")

avg = fastq.average_phred_by_position("reads.fastq")  # {0: 40.0, 1: 40.0, 2: 39.5, ...}
records = [fastq.FastqRecord(*row) for row in islice(fastq.iter_fastq("reads.fastq"), 3)]
summary = fastq.summarize_fastq("reads.fastq")  # {'total_reads': 2, 'total_bases': 9, ...}
```

`iter_fastq` and `summarize_fastq` read one record at a time and never hold the whole file in memory. Both open files by suffix, so `.fastq.gz` works without extra arguments.

#### Ecosystem tools (external)

For comprehensive QC, trimming, and preprocessing, we recommend using specialized tools and then reading results alongside these helpers:

- FastQC: quick, comprehensive QC reports
- fastp: all-in-one preprocessor (QC, adapter trimming, filtering)
- R packages for QC and exploration: Rqc, qckitfastq, seqTools, ShortRead
- Python libraries: HTSeq
- Additional utilities: fqtools, NGSUtils (fastqutils), fastQ_brew

These tools provide metrics such as per-base quality, GC content distributions, adapter content, and duplication levels that complement our lightweight summaries.
