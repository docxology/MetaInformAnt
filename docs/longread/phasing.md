# Read Phasing and Haplotype Analysis

The phasing module provides read-based haplotype phasing, phase block construction, haplotagging, switch error computation, and allele-specific analysis for long-read sequencing data. Long reads are uniquely suited for phasing because they span multiple heterozygous variants on a single molecule.

The functions live in the `haplotyping` submodule (`from metainformant.longread.phasing import haplotyping`); they are not re-exported on `metainformant.longread.phasing`. They work on plain dictionaries. A second, dataclass-based implementation (`PhaseBlock` objects, `tag_reads_by_haplotype`) lives in `metainformant.longread.analysis.phasing`.

## Key Concepts

### Read-Based Phasing

Each long read that spans two or more heterozygous variants provides direct evidence for which alleles co-occur on the same haplotype. By connecting variants through overlapping reads, haplotype blocks are constructed without requiring family data or statistical imputation.

### Phase Blocks

A phase block is a set of heterozygous variants that are phased relative to each other because reads link them. `phase_reads` builds blocks as the connected components of the variant co-observation graph, so a block ends where no read links the variants beyond it; two blocks can interleave when different reads link different variant pairs. Longer reads and higher coverage produce longer phase blocks.

A block's `start` and `end` are the positions of its first and last variant (`end` is inclusive). Its `variant_haplotypes` map records, for every variant in the block, the allele (0 reference, 1 alternate) that haplotype 1 carries there. Haplotypes 1 and 2 are labels within a block, not parental origins, and the same number in two blocks need not be the same chromosome copy.

### Haplotagging

After phasing, each read is tagged with a haplotype (HP:i:1 or HP:i:2) and the phase block it belongs to (PS) from the alleles it carries at heterozygous sites. Every allele a read carries at a variant of an overlapping block is a vote, for haplotype 1 when it equals the allele haplotype 1 carries there and for haplotype 2 otherwise; the read takes the block with the most votes and that block's majority haplotype. This enables haplotype-resolved analyses including allele-specific expression and allele-specific methylation.

### Switch Errors

Switch errors measure phasing accuracy by comparing inferred haplotypes to truth data. A switch error occurs when the phase relationship between two consecutive heterozygous sites is inverted relative to the true haplotype.

### Allele-Specific Analysis

With phased reads, allele-specific analyses quantify differences between the two haplotype copies, such as allele-specific expression (ASE) or allele-specific methylation (ASM).

## Function Reference

### phase_reads

```python-snippet
def phase_reads(
    variants: list[dict[str, Any]],
    reads: list[dict[str, Any]],
    method: str = "greedy",
) -> dict[str, Any]
```

Phase reads using heterozygous variant alleles. `variants` are dicts with `chrom`, `position`, `ref`, `alt`. `reads` are dicts with `read_id`, `chrom`, `start`, `end` and `alleles`, a list of `{"position": ..., "allele": 0 or 1}` observations (0 is the reference allele). Variants are joined into a co-observation graph and bipartitioned with a greedy max-cut; reads are then assigned by majority vote, and tied reads are left unassigned. Returns `haplotype_assignments` (`read_id` to 1 or 2), `phase_blocks` (dicts with `chrom`, `start`, `end`, `n_variants`, `n_reads` and `variant_haplotypes`; connected components of at least two variants), `switch_error_rate` (the fraction of variant links that disagree with the final bipartition) and `n_phased_reads`. `variant_haplotypes` maps each variant position in a block to the allele (0 or 1) that haplotype 1 carries there, the map `haplotag_reads` votes against, so the blocks go straight into it. Allele observations outside 0/1 do not contribute votes. Any `method` other than `"greedy"` raises `ValueError`.

### build_phase_blocks

```python-snippet
def build_phase_blocks(
    phased_variants: list[dict[str, Any]],
    max_gap: int = 1000000,
) -> list[dict[str, Any]]
```

Group already-phased variants into contiguous phase blocks. `phased_variants` are dicts with `chrom`, `position` and `haplotype`, the haplotype (1 or 2) that carries the alternate allele of the heterozygous variant (VCF genotype `1|0` is 1, `0|1` is 2). An optional `allele` (0 reference or 1 alternate, default 1) names the allele that `haplotype` carries instead, and an optional `quality` feeds `phase_quality`. Variants on a chromosome that are within `max_gap` bases of the previous one share a block, and blocks of fewer than two variants are dropped. This function does not consume `phase_reads` output. Returns block dicts with `chrom`, `start`, `end`, `n_variants`, `n_reads` (`None`, since no reads are given), `phase_quality` (mean of the recorded variant qualities, or `None` when none is supplied) and `variant_haplotypes`, the position-to-allele map of haplotype 1 derived from `haplotype` and `allele`. A variant whose `haplotype` is not 1 or 2, or whose `allele` is not 0 or 1, is left out of the map. The blocks go straight into `haplotag_reads`.

### haplotag_reads

```python-snippet
def haplotag_reads(
    reads: list[dict[str, Any]],
    phase_blocks: list[dict[str, Any]],
) -> list[dict[str, Any]]
```

Assign haplotype tags to reads based on phase block information. Reads need `chrom`, `start`, `end` and `alleles`; blocks need `chrom`, `start`, `end` and a `variant_haplotypes` dict that maps each variant position to the allele (0 or 1) that haplotype 1 carries there. The blocks returned by `phase_reads` and `build_phase_blocks` carry that map, so their output can be passed in unchanged. Position keys may be integers or their JSON string forms. Invalid position keys, map alleles outside 0/1 or conflicting keys for the same position raise `ValueError`; a non-mapping allele map raises `TypeError`. Read observations outside 0/1 are ignored in phasing and tagging.

A read overlaps a block when it starts at or before the block's `end` (the position of its last variant, inclusive) and ends after the block's `start` (the read's `end` is exclusive). In every overlapping block, each allele the read carries at a block variant is a vote: the allele that haplotype 1 carries there is a vote for haplotype 1, the other biallelic observation (0 or 1) is a vote for haplotype 2. Observations outside 0/1 do not vote. The read takes the overlapping block with the most votes (the first such block, ordered by `start`, on ties, so interleaved blocks from `phase_reads` are all reachable) and the haplotype with the majority of that block's votes.

Returns copies of the reads with `HP` (0 for untagged, 1 or 2) and `PS` (index of the block in `phase_blocks`) added. A read that overlaps no block, carries no allele at a block variant, or splits its votes evenly gets `HP` 0 and `PS` 0; because 0 is also the index of the first block, test `HP` before reading `PS`. A block without `variant_haplotypes` cannot tag any read, and a warning reports how many of the blocks passed in lack it. `longread.analysis.phasing.tag_reads_by_haplotype` is the dataclass-based counterpart.

### compute_switch_errors

```python-snippet
def compute_switch_errors(
    truth_haplotypes: list[int],
    predicted_haplotypes: list[int],
) -> dict[str, Any]
```

Compute the switch error rate by comparing predicted and truth haplotype labels (1 or 2 per variant, in position order). Both lists must have the same length, otherwise `ValueError` is raised. Returns `switch_error_rate`, `n_switches` (consecutive pairs compared), `n_errors`, `longest_correct_run` and `hamming_distance` (fraction of positions that differ after allowing a global flip).

### allele_specific_analysis

```python-snippet
def allele_specific_analysis(
    reads: list[dict[str, Any]],
    haplotype_assignments: dict[str, int],
    region: dict[str, Any],
) -> dict[str, Any]
```

Compare a per-read signal (expression or methylation level) between the two haplotypes inside a region. `reads` are dicts with `read_id`, `chrom`, `start`, `end` and `signal`; `haplotype_assignments` maps `read_id` to 1 or 2 (the `phase_reads` result has one); `region` is a dict with `chrom`, `start` and `end`. Returns `haplotype_1_signal` and `haplotype_2_signal` (mean signal), `allelic_imbalance` (absolute difference divided by the sum of the two means) and `p_value` (Mann-Whitney U test, Welch t-test approximation without scipy). Without reads on both haplotypes the signals are 0.0, the imbalance is 0.0 and the p-value is 1.0. Requires numpy.

## Usage Examples

```python
from metainformant.longread.phasing import haplotyping

# Heterozygous sites (from variant calling)
het_sites = [
    {"chrom": "chr1", "position": 10000, "ref": "A", "alt": "G"},
    {"chrom": "chr1", "position": 10500, "ref": "C", "alt": "T"},
    {"chrom": "chr1", "position": 11200, "ref": "G", "alt": "A"},
]

def make_read(read_id, alleles, signal):
    """A read spanning the three sites; `alleles` gives 0 (ref) or 1 (alt) at each of them.

    "signal" is a per-read value (e.g. methylation probability) for the allele-specific test.
    """
    return {"read_id": read_id, "chrom": "chr1", "start": 9500, "end": 12000, "signal": signal,
            "alleles": [{"position": p, "allele": a} for p, a in zip((10000, 10500, 11200), alleles)]}


reads = [
    make_read("read_001", (1, 1, 0), 0.85), make_read("read_002", (1, 1, 0), 0.80),
    make_read("read_003", (0, 0, 1), 0.20), make_read("read_004", (0, 0, 1), 0.25),
]

# Phase reads
phased = haplotyping.phase_reads(het_sites, reads)
print(f"Phased {phased['n_phased_reads']} reads, switch error rate {phased['switch_error_rate']:.2%}")
for block in phased["phase_blocks"]:
    print(f"Block: {block['chrom']}:{block['start']}-{block['end']} "
          f"({block['n_variants']} variants, {block['n_reads']} reads)")
    print(f"  allele carried by haplotype 1: {block['variant_haplotypes']}")

# Haplotag the reads with the blocks phase_reads returned: they already carry variant_haplotypes
tagged = haplotyping.haplotag_reads(reads, phased["phase_blocks"])
hap1_count = sum(1 for r in tagged if r["HP"] == 1)
hap2_count = sum(1 for r in tagged if r["HP"] == 2)
print(f"Hap1: {hap1_count} reads, Hap2: {hap2_count} reads")

# Group variants that are already phased (e.g. read from a phased VCF) into blocks. "haplotype" is
# the haplotype that carries the alternate allele (VCF 1|0 is 1, 0|1 is 2)
phased_variants = [
    {"chrom": "chr1", "position": 10000, "haplotype": 1, "quality": 0.9},
    {"chrom": "chr1", "position": 10500, "haplotype": 1, "quality": 0.8},
    {"chrom": "chr1", "position": 11200, "haplotype": 2, "quality": 0.85},
]
blocks = haplotyping.build_phase_blocks(phased_variants, max_gap=1_000_000)
print(blocks[0]["variant_haplotypes"])  # {10000: 1, 10500: 1, 11200: 0}: the allele haplotype 1 carries
tagged_from_vcf = haplotyping.haplotag_reads(reads, blocks)  # these blocks tag reads as they are

# Compute switch errors against a truth set (one haplotype label per variant, in position order)
errors = haplotyping.compute_switch_errors([1, 1, 2, 2, 1], [1, 1, 2, 1, 1])
print(f"Switch error rate: {errors['switch_error_rate']:.2%}")

# Allele-specific analysis of the per-read signal inside a region
ase = haplotyping.allele_specific_analysis(
    reads, phased["haplotype_assignments"], {"chrom": "chr1", "start": 9500, "end": 12000},
)
print(f"Haplotype 1: {ase['haplotype_1_signal']:.2f}, haplotype 2: {ase['haplotype_2_signal']:.2f}, "
      f"imbalance {ase['allelic_imbalance']:.2f} (p={ase['p_value']:.3g})")
```

## Configuration

- **Environment prefix**: `LR_` is honoured only by `workflow.pipelines.load_pipeline_config` and `utils.batch.process_batch` (`LR_BATCH_WORKERS`); the functions on this page read no environment variables
- **Optional dependencies**: numpy (required by `allele_specific_analysis`), scipy (Mann-Whitney U test; a Welch t-test approximation is used without it)
- The only phasing method is `"greedy"`
- Minimum phase block size: 2 variants
- `phase_reads` forms blocks from the reads that link variants (connected components, so blocks can interleave); `build_phase_blocks` breaks them where the gap between consecutive variants exceeds `max_gap` (1,000,000 bp by default)
- Blocks from both functions carry `variant_haplotypes`, which `haplotag_reads` needs; a block's `end` is the position of its last variant (inclusive)

## Related Modules

- `longread.io` -- Read BAM files containing aligned reads
- `longread.methylation` -- Per-read methylation calls; a per-read methylation value can serve as the `signal` for `allele_specific_analysis`
- `longread.analysis.phasing` -- The dataclass-based phasing functions (`phase_reads` that also calls alleles from `aligned_pairs`, `tag_reads_by_haplotype`) in the analysis submodule
- `longread.visualization` -- `plot_phasing_blocks` for block visualization
