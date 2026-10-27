# Long-Read Analysis

The analysis module provides structural variant detection, modified base calling (5mC, 6mA), and haplotype phasing from long-read sequencing data. These analyses leverage the long read lengths and base modification information unique to ONT and PacBio platforms.

The package exposes three submodules, `structural`, `modified_bases` and `phasing`. Import the submodule you need (`from metainformant.longread.analysis import structural`); the functions are not re-exported on `metainformant.longread.analysis` itself.

## Key Concepts

### Structural Variant Detection

Long reads span structural variants that are invisible to short-read sequencing. The module detects deletions (`DEL`), insertions (`INS`), inversions (`INV`), tandem duplications (`DUP`) and breakends (`BND`) by analyzing large CIGAR operations and split/supplementary alignments (the `SA` tag).

### Modified Base Detection

Oxford Nanopore sequencing directly detects DNA modifications (5-methylcytosine, 6-methyladenine) from the raw electrical signal. The module calls modified bases, aggregates methylation across reads, and performs differential methylation analysis between conditions.

### Haplotype Phasing

Long reads can phase heterozygous variants by linking alleles that appear on the same physical DNA molecule. The module builds haplotype blocks, tags reads by haplotype, and computes phase block statistics.

## Function Reference

### Structural Variant Detection

All detectors live in `metainformant.longread.analysis.structural`.

```python-snippet
def detect_sv_from_long_reads(
    alignments: Sequence[Any],
    min_size: int = 50,
    min_support: int = 2,
    min_mapping_quality: int = 20,
) -> list[StructuralVariant]
```

Detect structural variants from long-read alignments. `alignments` are `LongReadAlignment` objects (see `longread.io`) or dicts with the same field names. The detector combines CIGAR-based calls (large `I`/`D` operations) with split-read calls from `SA` tags, clusters calls of the same type within 500 bp, and keeps clusters supported by at least `min_support` distinct reads. Alignments below `min_mapping_quality` are skipped, including dicts or objects that carry no mapping quality (it defaults to 0). Returns `StructuralVariant` objects sorted by chromosome and position, each with `sv_type`, `chromosome`, `start`, `end`, `size`, `supporting_reads`, `read_names`, `quality`, `genotype`, `haplotype` and `info`; the dataclass also has `inserted_sequence`, `mate_chromosome` and `mate_position`, but the clustering step does not fill them, so they stay empty on the calls this function returns.

```python-snippet
def detect_insertions(alignments, min_size: int = 50, min_support: int = 2) -> list[StructuralVariant]
def detect_inversions(alignments, min_size: int = 50, min_support: int = 2) -> list[StructuralVariant]
def phase_structural_variants(
    variants: Sequence[StructuralVariant], haplotype_tags: dict[str, int],
) -> list[StructuralVariant]
```

`detect_insertions` and `detect_inversions` run the general detector and keep only `INS` or `INV` calls. `phase_structural_variants` assigns each SV to haplotype 1 or 2 by majority vote over its supporting reads (an SV with no phased supporting read, or a tie, gets haplotype 0 and genotype `./.`); `haplotype_tags` maps read name to haplotype (1 or 2), and the returned copies carry `haplotype`, `genotype` and `info["hp1_support"]`, `info["hp2_support"]`, `info["phase_confidence"]`.

### Modified Base Detection

These functions live in `metainformant.longread.analysis.modified_bases`.

```python-snippet
def detect_methylation(
    signal_data: Sequence[float] | Any, model: str = "cpg", threshold: float = 0.5,
) -> list[MethylationCall]
def call_5mc(signal_features: Sequence[dict[str, Any]]) -> list[MethylationCall]
def call_6ma(signal_features: Sequence[dict[str, Any]]) -> list[MethylationCall]
def aggregate_methylation(
    calls: Sequence[MethylationCall], regions: Sequence[dict[str, Any]], min_coverage: int = 3,
) -> list[RegionMethylation]
def differential_methylation(
    sample1: Sequence[MethylationCall], sample2: Sequence[MethylationCall],
    regions: Sequence[dict[str, Any]] | None = None,
    min_coverage: int = 5, min_difference: float = 0.2, alpha: float = 0.05,
) -> list[DifferentialMethylationResult]
```

`detect_methylation` segments a raw nanopore current trace (pA values, requires numpy) into events and scores each one; `model` is `"cpg"`, `"6ma"` or `"dam"`. The `position` of each returned call is the event's start sample in the trace, not a genomic coordinate. `call_5mc` and `call_6ma` score pre-extracted per-site features and keep genomic coordinates: each dict may carry `chromosome`, `position`, `strand`, `context`, `context_shift`, `dwell_time` and `current_std`. Both return one `MethylationCall` per input site (`chromosome`, `position`, `strand`, `modification_type`, `probability`, `coverage`, `modified_count`, `context`, `read_id`).

`aggregate_methylation` groups calls into the supplied regions (dicts with `chromosome`, `start`, `end` and an optional `name`) and returns one `RegionMethylation` per region (`mean_methylation`, `median_methylation`, `num_cpgs`, `mean_coverage`, `methylation_levels`); a site needs at least `min_coverage` calls to count. `differential_methylation` compares two call sets site by site (Fisher's exact test, or a chi-squared approximation when scipy is not installed; Benjamini-Hochberg adjustment) and returns `DifferentialMethylationResult` objects with `chromosome`, `position`, `sample1_methylation`, `sample2_methylation`, `difference` (sample 2 minus sample 1), `p_value`, `adjusted_p_value`, `is_significant` and per-sample coverage. To read methylation that the basecaller already wrote into the BAM `MM`/`ML` tags, use `longread.io.bam.extract_methylation_tags` instead.

### Haplotype Phasing

These functions live in `metainformant.longread.analysis.phasing`.

```python-snippet
def phase_reads(
    variants: Sequence[Variant | dict[str, Any]], reads: Sequence[dict[str, Any]],
    min_variant_quality: int = 20, min_base_quality: int = 10,
) -> tuple[list[PhaseBlock], dict[str, int]]
def build_haplotype_blocks(
    phased_variants: Sequence[dict[str, Any]], max_gap: int = 300000,
) -> list[PhaseBlock]
def tag_reads_by_haplotype(
    reads: Sequence[dict[str, Any]], phase_blocks: Sequence[PhaseBlock],
) -> list[dict[str, Any]]
def calculate_phase_block_stats(blocks: Sequence[PhaseBlock]) -> PhaseBlockStats
```

`phase_reads` takes heterozygous variants (`Variant` objects or dicts with `chromosome`, `position`, `ref_allele`, `alt_allele`) and read dicts that carry a `read_name` and a `variants` list of observed alleles (`position`, `allele` as 0/1, `quality`, and an optional `chromosome` that defaults to the read's `reference_name`). `min_variant_quality` is accepted but not used; `min_base_quality` filters the observed alleles. It returns the phase blocks (`PhaseBlock` objects built from connected components of co-observed variants) together with a `read_name -> haplotype (1 or 2)` dict; tied reads are left out. `build_haplotype_blocks` does not consume `phase_reads` output: it groups already-phased variant dicts (`chromosome`, `position`, `haplotype`, `allele`, optional `phase_set`) into blocks, splitting at gaps larger than `max_gap`. `tag_reads_by_haplotype` adds `HP` (0, 1 or 2) and `PS` (the block id) to copies of the read dicts, using each read's `reference_name`, `reference_start`, `reference_end` and `variants`. `calculate_phase_block_stats` returns a `PhaseBlockStats` with `num_blocks`, `num_phased_variants`, `n50`, `max_block_size`, `mean_block_size`, `median_block_size`, `switch_error_rate` and `total_phased_span`; sizes and N50 are in base pairs.

## Usage Examples

```python
from metainformant.longread.analysis import modified_bases, phasing, structural
from metainformant.longread.io import bam

# Load alignments (needs pysam and an indexed BAM)
alignments = bam.read_long_read_bam("aligned.bam", region="chr1:1-5000000")

# Detect structural variants
svs = structural.detect_sv_from_long_reads(alignments, min_size=50, min_support=3)
for sv in svs:
    print(f"{sv.sv_type}: {sv.chromosome}:{sv.start}-{sv.end} ({sv.size} bp, {sv.supporting_reads} reads)")

# Detect insertions specifically
insertions = structural.detect_insertions(alignments, min_size=50)

# Call 5mC from per-site signal features (one dict per read and CpG site)
tumor_features = [
    {"chromosome": "chr1", "position": 1500, "strand": "+",
     "context_shift": 3.0, "dwell_time": 1.5, "current_std": 1.2},
    # ...
]
tumor_calls = modified_bases.call_5mc(tumor_features)
normal_calls = modified_bases.call_5mc(normal_features)  # same feature layout, from the normal sample

# Aggregate calls over regions of interest
regions = [{"chromosome": "chr1", "start": 1000, "end": 2000, "name": "promoter"}]
for region in modified_bases.aggregate_methylation(tumor_calls, regions, min_coverage=3):
    print(f"{region.name}: {region.mean_methylation:.2f} over {region.num_cpgs} CpGs")

# Differential methylation between conditions
dmps = modified_bases.differential_methylation(tumor_calls, normal_calls, min_coverage=5)
for dmp in dmps[:5]:
    print(f"{dmp.chromosome}:{dmp.position} delta={dmp.difference:.2f} adj_p={dmp.adjusted_p_value:.3g}")

# Phase reads using known heterozygous sites
het_sites = [
    {"chromosome": "chr1", "position": 1000, "ref_allele": "A", "alt_allele": "G"},
    {"chromosome": "chr1", "position": 1500, "ref_allele": "C", "alt_allele": "T"},
]
reads = [
    {"read_name": "read1", "reference_name": "chr1", "reference_start": 900, "reference_end": 2000,
     "variants": [{"position": 1000, "allele": 1, "quality": 30},
                  {"position": 1500, "allele": 1, "quality": 30}]},
    # ...
]
blocks, read_haplotypes = phasing.phase_reads(het_sites, reads)
stats = phasing.calculate_phase_block_stats(blocks)
print(f"Phase block N50: {stats.n50} bp, {stats.num_blocks} blocks")
tagged_reads = phasing.tag_reads_by_haplotype(reads, blocks)

# Assign structural variants to haplotypes using the read-level assignments
phased_svs = structural.phase_structural_variants(svs, read_haplotypes)
```

## Configuration

- **Environment prefix**: `LR_` is honoured only by `workflow.pipelines.load_pipeline_config` (`LR_THREADS`, `LR_WORK_DIR`, `LR_LOG_DIR`) and by `utils.batch.process_batch` (`LR_BATCH_WORKERS`). The functions on this page read no environment variables; every threshold is an argument
- **Optional dependencies**: numpy (signal-level calling), scipy (statistical tests), pysam (reading BAM files through `longread.io`)
- Default minimum SV size is 50bp (standard for long-read SV calling)
- Minimum support of 2 distinct reads is required for an SV call by default (the SV pipeline configuration from `workflow.pipelines.get_sv_pipeline_config` raises this to 3)
- `aggregate_methylation` requires 3 calls per site by default; `differential_methylation` requires 5 in each sample

## Related Modules

- `longread.io` -- Read BAM and FAST5 files for analysis input
- `longread.methylation` -- Signal-level methylation calling module
- `longread.phasing` -- Dedicated phasing submodule with additional functions (dict-based phasing, switch errors, allele-specific analysis)
- `longread.quality` -- Quality metrics to filter reads before analysis
- `longread.assembly` -- Assembly from long reads
