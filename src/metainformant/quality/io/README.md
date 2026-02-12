# Quality IO

FASTQ file reading and per-read quality analysis including base quality, GC content, adapter detection, and length distributions.

## Contents

| File | Purpose |
|------|---------|
| `fastq.py` | FASTQ record parsing, per-base quality, GC content, adapter screening |

## Key Classes and Functions

| Symbol | Description |
|--------|-------------|
| `FastqRecord` | Dataclass for a single FASTQ read (header, sequence, quality) |
| `read_fastq_records()` | Iterator over FASTQ records from file path |
| `analyze_fastq_quality()` | Complete quality analysis of a FASTQ file |
| `basic_statistics()` | Read count, total bases, mean quality, mean length |
| `per_base_quality()` | Quality score distribution at each read position |
| `per_sequence_quality()` | Distribution of mean quality scores across reads |
| `sequence_length_distribution()` | Histogram of read lengths |
| `gc_content_distribution()` | Per-read GC content distribution |
| `adapter_content()` | Adapter sequence detection rates by position |
| `overrepresented_sequences()` | Identify frequently occurring subsequences |

## Usage

```python
from metainformant.quality.io.fastq import read_fastq_records, analyze_fastq_quality

quality = analyze_fastq_quality("data/sample.fastq.gz", n_reads=10000)
for record in read_fastq_records("data/sample.fastq.gz", max_records=100):
    print(record.header, len(record.sequence))
```
