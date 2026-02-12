# Analysis

Long-read analysis module for modified bases, structural variants, and phasing from nanopore and PacBio sequencing data.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports modified_bases, structural, phasing submodules |
| `modified_bases.py` | 5mC/6mA methylation detection from nanopore signal features |
| `structural.py` | Structural variant detection from split/supplementary alignments |
| `phasing.py` | Graph-based haplotype phasing using heterozygous variants |

## Key Functions

| Function | Description |
|----------|-------------|
| `structural.detect_sv_from_long_reads()` | Detect SVs (DEL, INS, INV, DUP, BND) from long-read alignments |
| `structural.detect_insertions()` | Detect insertion variants from CIGAR signatures |
| `structural.detect_inversions()` | Detect inversions from split-read evidence |
| `structural.phase_structural_variants()` | Assign SVs to haplotypes |
| `modified_bases.detect_methylation()` | Detect 5mC/6mA modifications from signal data |
| `modified_bases.aggregate_methylation()` | Aggregate per-read calls to per-site frequencies |
| `modified_bases.differential_methylation()` | Find differentially methylated regions |
| `phasing.phase_reads()` | Phase reads into haplotypes via heterozygous SNPs |
| `phasing.build_haplotype_blocks()` | Construct contiguous phase blocks |

## Usage

```python
from metainformant.longread.analysis import structural, modified_bases, phasing

svs = structural.detect_sv_from_long_reads(alignments, min_sv_size=50)
meth = modified_bases.detect_methylation(signal_data)
phases = phasing.phase_reads(variants, reads)
```
