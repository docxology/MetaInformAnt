# Long-Read Analysis

The analysis module provides structural variant detection, modified base calling (5mC, 6mA), and haplotype phasing from long-read sequencing data. These analyses leverage the long read lengths and base modification information unique to ONT and PacBio platforms.

## Key Concepts

### Structural Variant Detection

Long reads span structural variants that are invisible to short-read sequencing. The module detects insertions, deletions, inversions, and complex rearrangements by analyzing split/supplementary alignments and CIGAR string patterns.

### Modified Base Detection

Oxford Nanopore sequencing directly detects DNA modifications (5-methylcytosine, 6-methyladenine) from the raw electrical signal. The module calls modified bases, aggregates methylation across reads, and performs differential methylation analysis between conditions.

### Haplotype Phasing

Long reads can phase heterozygous variants by linking alleles that appear on the same physical DNA molecule. The module builds haplotype blocks, tags reads by haplotype, and computes phase block statistics.

## Function Reference

### Structural Variant Detection

```python
def detect_sv_from_long_reads(
    alignments: list[Any],
    min_sv_size: int = 50,
    min_support: int = 3,
) -> list[dict[str, Any]]
```

Detect structural variants from long-read alignments. Analyzes supplementary alignments and large CIGAR operations. Returns a list of SV records with type, coordinates, size, and supporting read count.

```python
def detect_insertions(alignments: list[Any], min_size: int = 50) -> list[dict]
def detect_inversions(alignments: list[Any], min_size: int = 500) -> list[dict]
def phase_structural_variants(svs: list[dict], phase_data: dict) -> list[dict]
```

Specialized detectors for insertions and inversions. `phase_structural_variants` assigns SVs to haplotypes using phasing data.

### Modified Base Detection

```python
def detect_methylation(alignments: list[Any]) -> dict[str, Any]
def call_5mc(signal_data: Any, positions: list[int]) -> dict[str, Any]
def call_6ma(signal_data: Any, positions: list[int]) -> dict[str, Any]
def aggregate_methylation(
    methylation_calls: list[dict], region: str | None = None
) -> dict[str, Any]
def differential_methylation(
    group1_calls: list[dict], group2_calls: list[dict],
    min_coverage: int = 5,
) -> list[dict[str, Any]]
```

`detect_methylation` extracts methylation data from BAM alignment tags. `call_5mc` and `call_6ma` call specific modification types from signal data. `aggregate_methylation` combines per-read calls into per-site statistics. `differential_methylation` identifies differentially methylated positions between two groups.

### Haplotype Phasing

```python
def phase_reads(
    alignments: list[Any],
    heterozygous_sites: list[dict],
) -> dict[str, Any]
def build_haplotype_blocks(
    phased_reads: dict, min_block_size: int = 2,
) -> list[dict[str, Any]]
def tag_reads_by_haplotype(
    alignments: list[Any], phase_blocks: list[dict],
) -> list[dict[str, Any]]
def calculate_phase_block_stats(
    phase_blocks: list[dict],
) -> dict[str, Any]
```

`phase_reads` assigns reads to haplotypes based on heterozygous variant alleles. `build_haplotype_blocks` groups consecutive phased sites into contiguous blocks. `tag_reads_by_haplotype` annotates alignments with haplotype tags. `calculate_phase_block_stats` computes N50, span, and completeness metrics for phase blocks.

## Usage Examples

```python
from metainformant.longread import io, analysis

# Load alignments
alignments = io.read_long_read_bam("aligned.bam", region="chr1:1-5000000")

# Detect structural variants
svs = analysis.detect_sv_from_long_reads(alignments, min_sv_size=50, min_support=3)
for sv in svs:
    print(f"{sv['type']}: {sv['chrom']}:{sv['start']}-{sv['end']} ({sv['size']} bp)")

# Detect insertions specifically
insertions = analysis.detect_insertions(alignments, min_size=50)

# Extract methylation from BAM tags
meth_data = analysis.detect_methylation(alignments)

# Aggregate methylation across reads
per_site = analysis.aggregate_methylation(meth_data["calls"], region="chr1:1000-2000")

# Differential methylation between conditions
dmps = analysis.differential_methylation(group1_calls, group2_calls, min_coverage=5)
for dmp in dmps[:5]:
    print(f"Position {dmp['position']}: delta={dmp['methylation_difference']:.2f}")

# Phase reads using known heterozygous sites
het_sites = [{"chrom": "chr1", "pos": 1000, "ref": "A", "alt": "G"}, ...]
phased = analysis.phase_reads(alignments, het_sites)
blocks = analysis.build_haplotype_blocks(phased)
stats = analysis.calculate_phase_block_stats(blocks)
print(f"Phase block N50: {stats['n50']} variants, {len(blocks)} blocks")
```

## Configuration

- **Environment prefix**: `LR_`
- **Optional dependencies**: numpy, pysam
- Default minimum SV size is 50bp (standard for long-read SV calling)
- Minimum support of 3 reads required for an SV call by default
- Methylation aggregation requires at least 5x coverage by default

## Related Modules

- `longread.io` -- Read BAM and FAST5 files for analysis input
- `longread.methylation` -- Signal-level methylation calling module
- `longread.phasing` -- Dedicated phasing submodule with additional methods
- `longread.quality` -- Quality metrics to filter reads before analysis
- `longread.assembly` -- Assembly from long reads
