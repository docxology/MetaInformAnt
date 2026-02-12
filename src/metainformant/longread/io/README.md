# IO

Long-read I/O module for reading FAST5/POD5 signal data, BAM alignment files, and converting between common long-read sequencing formats.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports fast5, bam, formats submodules |
| `fast5.py` | FAST5/POD5 file parsing for ONT raw signal and basecalls |
| `bam.py` | Long-read BAM processing with methylation tags (MM/ML) |
| `formats.py` | Format conversion between FAST5, POD5, FASTQ, and PAF |

## Key Functions

| Function | Description |
|----------|-------------|
| `fast5.read_fast5()` | Parse FAST5/POD5 files into Fast5Read records |
| `fast5.extract_signal()` | Extract raw electrical signal arrays from reads |
| `fast5.extract_basecalls()` | Extract basecalled sequences and quality strings |
| `fast5.get_read_metadata()` | Get per-read metadata (channel, duration, etc.) |
| `bam.read_long_read_bam()` | Parse BAM files into LongReadAlignment records |
| `bam.extract_methylation_tags()` | Parse MM/ML methylation tags from alignments |
| `bam.get_supplementary_alignments()` | Extract supplementary alignments for SV detection |
| `bam.calculate_alignment_stats()` | Compute alignment statistics (identity, coverage) |
| `formats.fast5_to_fastq()` | Convert FAST5 to FASTQ format |
| `formats.convert_pod5_to_fast5()` | Convert POD5 to FAST5 format |
| `formats.write_paf()` | Write alignments in PAF (minimap2) format |

## Usage

```python
from metainformant.longread.io import fast5, bam, formats

reads = fast5.read_fast5("reads.fast5")
alignments = bam.read_long_read_bam("aligned.bam")
formats.fast5_to_fastq("input.fast5", "output.fastq")
```
