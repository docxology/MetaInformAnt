# Read Phasing and Haplotype Analysis

The phasing module provides read-based haplotype phasing, phase block construction, haplotagging, switch error computation, and allele-specific analysis for long-read sequencing data. Long reads are uniquely suited for phasing because they span multiple heterozygous variants on a single molecule.

The functions live in the `haplotyping` submodule (`from metainformant.longread.phasing import haplotyping`); they are not re-exported on `metainformant.longread.phasing`. They work on plain dictionaries. A second, dataclass-based implementation (`PhaseBlock` objects, `tag_reads_by_haplotype`) lives in `metainformant.longread.analysis.phasing`.

## Key Concepts

### Read-Based Phasing

Each long read that spans two or more heterozygous variants provides direct evidence for which alleles co-occur on the same haplotype. By connecting variants through overlapping reads, haplotype blocks are constructed without requiring family data or statistical imputation.

### Phase Blocks

A phase block is a set of heterozygous variants that are phased relative to each other because reads link them. `phase_reads` builds blocks as the connected components of the variant co-observation graph, so a block ends where no read links the variants beyond it; two blocks can interleave when different reads link different variant pairs. Longer reads and higher coverage produce longer phase blocks.

### Haplotagging

After phasing, each read is assigned to a haplotype (HP:i:1 or HP:i:2) based on the alleles it carries at heterozygous sites. This enables haplotype-resolved analyses including allele-specific expression and allele-specific methylation.

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

Phase reads using heterozygous variant alleles. `variants` are dicts with `chrom`, `position`, `ref`, `alt`. `reads` are dicts with `read_id`, `chrom`, `start`, `end` and `alleles`, a list of `{"position": ..., "allele": 0 or 1}` observations (0 is the reference allele). Variants are joined into a co-observation graph and bipartitioned with a greedy max-cut; reads are then assigned by majority vote, and tied reads are left unassigned. Returns `haplotype_assignments` (`read_id` to 1 or 2), `phase_blocks` (dicts with `chrom`, `start`, `end`, `n_variants`, `n_reads`; connected components of at least two variants), `switch_error_rate` (the fraction of variant links that disagree with the final bipartition) and `n_phased_reads`. Any `method` other than `"greedy"` raises `ValueError`.

### build_phase_blocks

```python-snippet
def build_phase_blocks(
    phased_variants: list[dict[str, Any]],
    max_gap: int = 1000000,
) -> list[dict[str, Any]]
```

Group already-phased variants into contiguous phase blocks. `phased_variants` are dicts with `chrom` and `position` (they normally also carry `haplotype`, which is not used for grouping) and an optional `quality`; variants on a chromosome that are within `max_gap` bases of the previous one share a block, and blocks of fewer than two variants are dropped. This function does not consume `phase_reads` output. Returns block dicts with `chrom`, `start`, `end`, `n_variants`, `n_reads` (always 0 here, since no reads are given) and `phase_quality` (mean variant quality, 1.0 when absent).

### haplotag_reads

```python-snippet
def haplotag_reads(
    reads: list[dict[str, Any]],
    phase_blocks: list[dict[str, Any]],
) -> list[dict[str, Any]]
```

Assign haplotype tags to reads based on phase block information. Reads need `chrom`, `start`, `end` and `alleles`; blocks need `chrom`, `start`, `end` and a `variant_haplotypes` dict that maps each variant position to the allele (0 or 1) that haplotype 1 carries there. A read takes the first block it overlaps and votes over the alleles it shares with that map. Returns copies of the reads with `HP` (0 for untagged, 1 or 2) and `PS` (index of the block in `phase_blocks`) added. Blocks returned by `phase_reads` and `build_phase_blocks` do not contain `variant_haplotypes`, so every read tagged against them keeps `HP` 0 until you add the map; `longread.analysis.phasing.tag_reads_by_haplotype` is the end-to-end alternative.

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

# Reads with the allele (0 = ref, 1 = alt) they carry at each site they span;
# "signal" is a per-read value (e.g. methylation probability) for the allele-specific test
reads = [
    {"read_id": "read_001", "chrom": "chr1", "start": 9500, "end": 12000, "signal": 0.85,
     "alleles": [{"position": 10000, "allele": 1}, {"position": 10500, "allele": 1},
                 {"position": 11200, "allele": 0}]},
    # ...
]

# Phase reads
phased = haplotyping.phase_reads(het_sites, reads)
print(f"Phased {phased['n_phased_reads']} reads, switch error rate {phased['switch_error_rate']:.2%}")
for block in phased["phase_blocks"]:
    print(f"Block: {block['chrom']}:{block['start']}-{block['end']} "
          f"({block['n_variants']} variants, {block['n_reads']} reads)")

# Group variants that are already phased (e.g. read from a phased VCF) into blocks
phased_variants = [
    {"chrom": "chr1", "position": 10000, "haplotype": 1, "quality": 0.9},
    {"chrom": "chr1", "position": 10500, "haplotype": 1, "quality": 0.8},
]
blocks = haplotyping.build_phase_blocks(phased_variants, max_gap=1_000_000)

# Haplotag reads: the blocks must say which allele haplotype 1 carries at each position
tag_blocks = [dict(block, variant_haplotypes={10000: 1, 10500: 1, 11200: 0})
              for block in phased["phase_blocks"]]
tagged = haplotyping.haplotag_reads(reads, tag_blocks)
hap1_count = sum(1 for r in tagged if r["HP"] == 1)
hap2_count = sum(1 for r in tagged if r["HP"] == 2)
print(f"Hap1: {hap1_count} reads, Hap2: {hap2_count} reads")

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

## Related Modules

- `longread.io` -- Read BAM files containing aligned reads
- `longread.methylation` -- Per-read methylation calls; a per-read methylation value can serve as the `signal` for `allele_specific_analysis`
- `longread.analysis.phasing` -- Additional phasing functions in the analysis submodule
- `longread.visualization` -- `plot_phasing_blocks` for block visualization
