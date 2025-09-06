### DNA: FASTQ

Minimal, streaming FASTQ utilities with gzip support. These helper functions are designed for local inspection and lightweight QC. For production-grade QC and preprocessing, see the ecosystem notes below.

#### API

- `iter_fastq(path) -> Iterator[tuple[str, str, str]]`: Stream `(read_id, sequence, quality)` from a FASTQ file (`.gz` supported). Truncates mismatched seq/qual lengths to the shorter for safety; skips incomplete trailing records.
- `average_phred_by_position(path) -> list[float]`: Average Phred+33 per position across reads, truncated to the shortest encountered read length. Streams input.
- `head(path, n=5) -> list[FastqRecord]`: First `n` records as lightweight objects.
- `gc_content(seq) -> float`: GC fraction among A/C/G/T (N-insensitive).
- `summarize_fastq(path, max_reads=None) -> dict`: Lightweight stats: `num_reads`, `length_min/max/mean`, `gc_mean`, `n_content_mean`, `avg_phred_by_pos`.

#### Usage

```python
from metainformant.dna import fastq

avg = fastq.average_phred_by_position("/path/to/reads.fastq.gz")
records = fastq.head("/path/to/reads.fastq", n=3)
summary = fastq.summarize_fastq("/path/to/reads.fastq.gz", max_reads=100000)
```

All functions stream via `metainformant.core.io.open_text_auto` and never load the full file in memory.

#### Ecosystem tools (external)

For comprehensive QC, trimming, and preprocessing, we recommend using specialized tools and then reading results alongside these helpers:

- FastQC: quick, comprehensive QC reports
- fastp: all-in-one preprocessor (QC, adapter trimming, filtering)
- R packages for QC and exploration: Rqc, qckitfastq, seqTools, ShortRead
- Python libraries: HTSeq
- Additional utilities: fqtools, NGSUtils (fastqutils), fastQ_brew

These tools provide metrics such as per-base quality, GC content distributions, adapter content, and duplication levels that complement our lightweight summaries.
