# Read Phasing and Haplotype Analysis

The phasing module provides read-based haplotype phasing, phase block construction, haplotagging, switch error computation, and allele-specific analysis for long-read sequencing data. Long reads are uniquely suited for phasing because they span multiple heterozygous variants on a single molecule.

## Key Concepts

### Read-Based Phasing

Each long read that spans two or more heterozygous variants provides direct evidence for which alleles co-occur on the same haplotype. By connecting variants through overlapping reads, haplotype blocks are constructed without requiring family data or statistical imputation.

### Phase Blocks

A phase block is a contiguous set of heterozygous variants that are phased relative to each other. Block boundaries occur where no reads span between consecutive variant pairs. Longer reads and higher coverage produce longer phase blocks.

### Haplotagging

After phasing, each read is assigned to a haplotype (HP:i:1 or HP:i:2) based on the alleles it carries at heterozygous sites. This enables haplotype-resolved analyses including allele-specific expression and allele-specific methylation.

### Switch Errors

Switch errors measure phasing accuracy by comparing inferred haplotypes to truth data. A switch error occurs when the phase relationship between two consecutive heterozygous sites is inverted relative to the true haplotype.

### Allele-Specific Analysis

With phased reads, allele-specific analyses quantify differences between the two haplotype copies, such as allele-specific expression (ASE) or allele-specific methylation (ASM).

## Function Reference

### phase_reads

```python
def phase_reads(
    alignments: list[Any],
    heterozygous_sites: list[dict],
    min_mapping_quality: int = 20,
) -> dict[str, Any]
```

Phase reads using heterozygous variant alleles. Input `heterozygous_sites` are dicts with `chrom`, `pos`, `ref`, `alt`. Returns `phased_variants` (variant-to-haplotype assignments), `read_haplotypes` (read-to-haplotype assignments), and `unphased_variants`.

### build_phase_blocks

```python
def build_phase_blocks(
    phased_data: dict[str, Any],
    min_block_size: int = 2,
) -> list[dict[str, Any]]
```

Construct contiguous phase blocks from phased variant data. Returns a list of block records with `chrom`, `start`, `end`, `n_variants`, `n_reads`, and `span` (genomic distance).

### haplotag_reads

```python
def haplotag_reads(
    alignments: list[Any],
    phase_blocks: list[dict],
) -> list[dict[str, Any]]
```

Assign haplotype tags to reads based on phase block information. Returns annotated reads with `haplotype` (1 or 2), `phase_set`, and `confidence`.

### compute_switch_errors

```python
def compute_switch_errors(
    inferred_phases: dict[str, Any],
    truth_phases: dict[str, Any],
) -> dict[str, Any]
```

Compute switch error rate by comparing inferred and truth haplotypes. Returns `switch_error_rate`, `n_switches`, `n_comparisons`, `switch_positions`, and `per_block_errors`.

### allele_specific_analysis

```python
def allele_specific_analysis(
    haplotagged_reads: list[dict],
    feature_data: dict[str, Any] | None = None,
) -> dict[str, Any]
```

Perform allele-specific analysis on haplotagged reads. If `feature_data` (expression counts or methylation levels) is provided, computes per-haplotype statistics. Returns `haplotype1_stats`, `haplotype2_stats`, `imbalance_scores`, and `significant_sites`.

## Usage Examples

```python
from metainformant.longread import io, phasing

# Load long-read alignments
alignments = io.read_long_read_bam("aligned.bam")

# Define heterozygous sites (from variant calling)
het_sites = [
    {"chrom": "chr1", "pos": 10000, "ref": "A", "alt": "G"},
    {"chrom": "chr1", "pos": 10500, "ref": "C", "alt": "T"},
    {"chrom": "chr1", "pos": 11200, "ref": "G", "alt": "A"},
]

# Phase reads
phased = phasing.phase_reads(alignments, het_sites, min_mapping_quality=20)
print(f"Phased {len(phased['phased_variants'])} variants")

# Build phase blocks
blocks = phasing.build_phase_blocks(phased, min_block_size=2)
for block in blocks:
    print(f"Block: {block['chrom']}:{block['start']}-{block['end']} "
          f"({block['n_variants']} variants, {block['span']} bp span)")

# Haplotag reads
tagged = phasing.haplotag_reads(alignments, blocks)
hap1_count = sum(1 for r in tagged if r["haplotype"] == 1)
hap2_count = sum(1 for r in tagged if r["haplotype"] == 2)
print(f"Hap1: {hap1_count} reads, Hap2: {hap2_count} reads")

# Compute switch errors against truth
if truth_available:
    errors = phasing.compute_switch_errors(phased, truth_phases)
    print(f"Switch error rate: {errors['switch_error_rate']:.2%}")

# Allele-specific analysis
ase = phasing.allele_specific_analysis(tagged, feature_data=expression_counts)
for site in ase["significant_sites"][:5]:
    print(f"Position {site['pos']}: imbalance={site['imbalance']:.2f}")
```

## Configuration

- **Environment prefix**: `LR_`
- **Optional dependencies**: numpy, pysam
- Default minimum mapping quality for phasing: 20 (MAPQ)
- Minimum phase block size: 2 variants
- Phase blocks are broken when no read spans between consecutive variants

## Related Modules

- `longread.io` -- Read BAM files containing aligned reads
- `longread.methylation` -- Allele-specific methylation using phased reads
- `longread.analysis.phasing` -- Additional phasing functions in the analysis submodule
- `longread.visualization` -- `plot_phasing_blocks` for block visualization
