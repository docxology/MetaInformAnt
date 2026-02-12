# IO

Read, write, and assess quality of FASTQ sequencing files with support for filtering, trimming, and format conversion.

## Contents

| File | Purpose |
|------|---------|
| `fastq.py` | FASTQ parsing, writing, quality assessment, read filtering, and trimming |

## Key Functions

| Function | Description |
|----------|-------------|
| `read_fastq()` | Parse a FASTQ file into a dict of (sequence, quality) tuples |
| `write_fastq()` | Write sequences and quality strings to FASTQ format |
| `assess_quality()` | Compute per-file quality metrics (mean Q, GC%, read count) |
| `filter_reads()` | Yield reads that pass a minimum average Phred quality |
| `trim_reads()` | Quality-based or fixed-length trimming of reads |
| `convert_fastq_to_fasta()` | Strip quality data and write FASTA output |
| `average_phred_by_position()` | Per-base-position average quality scores |
| `calculate_per_base_quality()` | Detailed per-position quality statistics |
| `summarize_fastq()` | Quick summary: total reads, bases, and mean length |
| `iter_fastq()` | Memory-efficient iterator yielding (header, seq, qual) tuples |
| `FastqRecord` | Dataclass representing a single FASTQ record |

## Usage

```python
from metainformant.dna.io.fastq import read_fastq, assess_quality, filter_reads

sequences = read_fastq("reads.fastq")
metrics = assess_quality("reads.fastq")
for read in filter_reads("reads.fastq", min_quality=25):
    process(read)
```
