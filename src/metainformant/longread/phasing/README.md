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
| `haplotyping.phase_reads()` | Phase long reads into haplotypes using variant information |
| `haplotyping.build_phase_blocks()` | Construct contiguous phase blocks from phased variants |
| `haplotyping.compute_switch_errors()` | Evaluate phasing quality by computing switch error rate |
| `haplotyping.haplotag_reads()` | Tag reads with haplotype assignments for downstream analysis |
| `haplotyping.allele_specific_analysis()` | Detect allele-specific expression or methylation |

## Usage

```python
from metainformant.longread.phasing import haplotyping

result = haplotyping.phase_reads(variants, reads, method="greedy")
blocks = haplotyping.build_phase_blocks(result["phased_variants"])
tagged = haplotyping.haplotag_reads(reads, result)
ase = haplotyping.allele_specific_analysis(tagged, expression_data)
```
