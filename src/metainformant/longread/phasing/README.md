# Phasing

Read-based haplotype phasing for long-read sequencing, using heterozygous variants to assign reads to haplotypes, construct phase blocks, and perform allele-specific analysis.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports haplotyping submodule |
| `haplotyping.py` | Greedy phasing, phase block construction, haplotagging, allele-specific analysis |

## Key Functions

| Function | Description |
|----------|-------------|
| `haplotyping.phase_reads()` | Phase long reads into haplotypes using variant information; returns read assignments and phase blocks |
| `haplotyping.build_phase_blocks()` | Construct contiguous phase blocks from already-phased variants (`haplotype`, optional `allele`) |
| `haplotyping.compute_switch_errors()` | Evaluate phasing quality by computing switch error rate |
| `haplotyping.haplotag_reads()` | Tag reads with `HP` and `PS` from phase blocks |
| `haplotyping.allele_specific_analysis()` | Compare a per-read signal (expression or methylation) between the two haplotypes in a region |

Blocks from `phase_reads()` and `build_phase_blocks()` carry `variant_haplotypes`, a map from variant position to the allele (0 reference, 1 alternate) that haplotype 1 carries there. `haplotag_reads()` votes against that map, so either function's blocks can be passed to it unchanged, including blocks reloaded from JSON (JSON turns the integer position keys into text; `haplotag_reads()` reads both). A read allele that is not 0 or 1 is not a vote. See [docs/longread/phasing.md](../../../../docs/longread/phasing.md) for the voting rule and block fields.

## Usage

```python
from metainformant.longread.phasing import haplotyping

variants = [
    {"chrom": "chr1", "position": 100, "ref": "A", "alt": "G"},
    {"chrom": "chr1", "position": 400, "ref": "C", "alt": "T"},
]
reads = [
    {"read_id": "r1", "chrom": "chr1", "start": 50, "end": 500, "signal": 0.9,
     "alleles": [{"position": 100, "allele": 1}, {"position": 400, "allele": 1}]},
    {"read_id": "r2", "chrom": "chr1", "start": 50, "end": 500, "signal": 0.1,
     "alleles": [{"position": 100, "allele": 0}, {"position": 400, "allele": 0}]},
]
region = {"chrom": "chr1", "start": 50, "end": 500}

result = haplotyping.phase_reads(variants, reads, method="greedy")
tagged = haplotyping.haplotag_reads(reads, result["phase_blocks"])
ase = haplotyping.allele_specific_analysis(reads, result["haplotype_assignments"], region)  # needs numpy
```
