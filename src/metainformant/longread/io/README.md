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
| `fast5.read_fast5()` | Parse FAST5/POD5 files into Fast5Read records; `mux` is each read's `start_mux` attribute (0 when absent) |
| `fast5.extract_signal()` | Extract raw electrical signal arrays from reads |
| `fast5.extract_basecalls()` | Extract basecalled sequences and quality strings |
| `fast5.get_read_metadata()` | Get per-read metadata (channel, mux, duration, etc.) |
| `bam.read_long_read_bam()` | Parse BAM files into LongReadAlignment records |
| `bam.extract_methylation_tags()` | Decode MM/ML base-modification tags into entries keyed `base + strand + code` (`C+m`, `C-m`, `G-m` are separate; `C+mh` gives `C+m` and `C+h`), with read and reference positions and probabilities per site; `ML` may be a list, bytes, numpy array or pysam `array('B')` |
| `bam.get_supplementary_alignments()` | Extract supplementary alignments for SV detection |
| `bam.calculate_alignment_stats()` | Compute alignment statistics (identity, coverage) |
| `formats.fast5_to_fastq()` | Convert FAST5 to FASTQ format |
| `formats.convert_pod5_to_fast5()` | Convert POD5 to FAST5 format |
| `formats.write_paf()` | Write alignments in PAF (minimap2) format |

The MM/ML decoder follows the SAM optional-tags specification, checked against htslib: the skip counts refer to the read as sequenced, so a reverse-aligned read is reverse-complemented before counting and positions are mapped back to the stored `query_sequence`; the counted letter is the one written in the entry for both strands (`G-m` counts Gs and marks the 5mC opposite them); multi-code entries carry one interleaved `ML` value per code and position; the `.` and `?` flags (and no flag, which means `.`) set `implicit`; MM without `ML` leaves the probabilities `None`. A malformed MM tag, an `ML` of the wrong length and skip counts that run past the read are reported by warnings. See [docs/longread/io.md](../../../../docs/longread/io.md).

## Usage

```python
from metainformant.longread.io import fast5, bam, formats

reads = fast5.read_fast5("reads.fast5")
alignments = bam.read_long_read_bam("aligned.bam")
formats.fast5_to_fastq("input.fast5", "output.fastq")
```
