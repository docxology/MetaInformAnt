# Structural Variants: Detection

CNV detection, structural variant calling, and breakpoint refinement.

## Concept Overview

The detection submodule identifies genomic structural variants through three
complementary approaches:

- **CNV Detection** -- Copy number variants from read depth using circular binary segmentation
- **SV Calling** -- Structural variants from split reads and discordant read pairs
- **Breakpoint Refinement** -- Base-pair resolution breakpoint determination, with strand, microhomology and inserted sequence per breakpoint

These methods are designed to work together: depth-based CNV detection finds
large copy number changes, SV calling identifies balanced rearrangements
(inversions, translocations), and breakpoint refinement provides precise
coordinates for all variant types.

## CNV Detection

### Core Types

**CNVSegment** -- A contiguous genomic region with consistent copy number.

| Field | Type | Description |
|-------|------|-------------|
| `chrom` | `str` | Chromosome |
| `start` | `int` | Start position (0-based, inclusive) |
| `end` | `int` | End position (0-based, exclusive) |
| `mean_log2ratio` | `float` | Mean log2 ratio for segment |
| `n_bins` | `int` | Number of coverage bins |
| `state` | `str` | `"NEUTRAL"`, `"DEL"`, `"DUP"`, `"AMP"`, `"HOMODEL"` |
| `cn` | `int` | Estimated copy number |
| `confidence` | `float` | Call confidence score |

**CNVResult** -- Result for one chromosome: `segments`, `log2_ratios`, `bin_size`, `method` and `parameters`.

### detect_cnv_from_depth

Primary entry point for CNV detection from read depth data. It returns a dictionary of `CNVResult` objects keyed by chromosome.

```python
import numpy as np

from metainformant.structural_variants.detection.cnv import detect_cnv_from_depth

results = detect_cnv_from_depth(
    depth_data={"chr1": np.array([100, 102, 98, 50, 48, 52, 99, 101])},
    window_size=1000,
    method="segmentation",       # CBS algorithm
    significance=0.01,
    ploidy=2,
    min_segment_bins=3,
    merge_distance=1000
)

for chrom, cnv_result in results.items():
    for segment in cnv_result.segments:
        print(f"{segment.chrom}:{segment.start}-{segment.end} "
              f"CN={segment.cn} ({segment.state})")
```

**Input.** Each value is the read depth of one bin, finite and non-negative (zero is a bin with no coverage). Each chromosome is divided by the median of its positive bins before segmentation. A negative or non-finite value raises `ValueError` naming the chromosome and the number of offending values, because negative values usually mean that log2 ratios were passed by mistake; segment log2 ratios with `segment_coverage` and `call_cnv_states` instead. Chromosomes with fewer than `min_segment_bins` bins are skipped before that check.

**Cost.** The CBS scan is vectorised with NumPy but still grows with the square of the number of bins per chromosome, so pick a `window_size` that leaves at most a few tens of thousands of bins per chromosome (see [PERFORMANCE.md](PERFORMANCE.md)).

### Supporting Functions

| Function | Signature | Purpose |
|----------|-----------|---------|
| `segment_coverage` | `(coverage_array, significance=0.01, max_iterations=100, min_width=2)` | CBS segmentation of log2 ratios; `min_width` must be at least 1 |
| `call_cnv_states` | `(segments, ploidy=2, del_threshold=-0.3, dup_threshold=0.3, amp_threshold=1.0, homodel_threshold=-1.5)` | Assign CN states |
| `merge_adjacent_segments` | `(segments, max_gap=1000)` | Merge nearby segments with same state |
| `calculate_log2_ratio` | `(tumor_depth, normal_depth, gc_content=None, pseudocount=1.0)` | Compute log2 ratio with optional GC correction |

## SV Calling

### Core Types

**SVType** enum: `DEL`, `DUP`, `INV`, `TRA`, `INS`, `BND`, `UNKNOWN`

**StructuralVariant** -- A detected SV with type, coordinates, and genotype.

| Field | Description |
|-------|-------------|
| `chrom`, `start`, `end` | Position (0-based). Calls other than translocations occupy `chrom:start-end` |
| `sv_type` | `SVType` member |
| `size` | Breakpoint distance in bp; 0 means "not measured" (always 0 for TRA and INS) |
| `quality`, `genotype`, `filter_status` | Quality score, genotype string, `"PASS"` or a reason |
| `evidence` | `SVEvidence`: `split_reads`, `discordant_pairs`, `breakpoint1`/`breakpoint2`, `chrom1`/`chrom2`, strands, evidence read names (sorted) |
| `chrom2`, `pos2` | Second breakpoint of a translocation; `""` and `None` for every other call |
| `info` | Additional fields |

A translocation (`TRA`) joins two loci on different chromosomes, so it has no interval: `chrom` and `start` give the first breakpoint, `chrom2` and `pos2` the second, and `end` repeats `start`. An insertion (`INS`) has two breakpoints less than 50 bp apart and its `size` is 0, because split-read evidence does not measure the inserted length.

**SVEvidence** -- Supporting evidence from reads (split reads, discordant pairs).

### call_structural_variants

Detect SVs from alignment evidence.

```python
from metainformant.structural_variants.detection.sv_calling import (
    call_structural_variants,
    InsertSizeStats,
)

insert_stats = InsertSizeStats(mean=300.0, std=50.0, median=290.0)

variants = call_structural_variants(
    alignments=alignment_records,
    min_support=3,
    insert_size_stats=insert_stats,
    min_mapq=20,
    min_sv_size=50
)

for sv in variants:
    if sv.sv_type.name == "TRA":
        where = f"{sv.chrom}:{sv.start} -> {sv.chrom2}:{sv.pos2}"
    else:
        where = f"{sv.chrom}:{sv.start}-{sv.end}"
    support = sv.evidence.split_reads + sv.evidence.discordant_pairs
    print(f"{sv.sv_type.name}: {where} support={support} genotype={sv.genotype}")
```

`min_sv_size` filters DEL, DUP, INV and UNKNOWN calls on their breakpoint distance. Translocations have no length and insertions do not have a measured one, so TRA and INS calls are exempt from it. Genotypes are assigned by `genotype_sv` with its default `min_mapq=20`, whatever `min_mapq` is passed here.

### Supporting Functions

| Function | Signature | Purpose |
|----------|-----------|---------|
| `detect_split_reads` | `(reads, min_clip=20)` | Find soft-clipped reads indicating breakpoints |
| `detect_discordant_pairs` | `(pairs, insert_size_stats, n_std=4.0)` | Find anomalous insert sizes/orientations |
| `classify_sv_type` | `(evidence)` | Determine SV type from evidence pattern |
| `genotype_sv` | `(variant, reads, min_mapq=20)` | Genotype: 0/0, 0/1, 1/1 (`./.` when no read counts); `reads` is a read list or a `ReadIndex` |
| `ReadIndex` | `(reads)` | Per-chromosome index of reads by start position; build it once and pass it to `genotype_sv` for many variants |

`genotype_sv` counts, as reference evidence, the reads that span `start - 100 .. start + 100` on the variant's chromosome and are not among the call's evidence reads. A `ReadIndex` gives the same genotype as the plain read list but replaces a pass over every read per variant with a binary search.

## Breakpoint Refinement

### Core Types

**Breakpoint** -- Single breakpoint with position, strand, confidence, microhomology, inserted sequence, and supporting reads.

**BreakpointPair** -- Two linked breakpoints defining an SV.

### refine_breakpoints

Improve breakpoint resolution using local read evidence. `variants` are dictionaries: `dataclasses.asdict` of a `StructuralVariant` works as is, translocations included.

```python
from dataclasses import asdict

from metainformant.structural_variants.detection.breakpoints import refine_breakpoints

refined = refine_breakpoints(
    variants=[asdict(sv) for sv in sv_calls],
    reads=local_reads,
    window=200,
    min_clip=5
)

for bp_pair in refined:
    print(f"{bp_pair.sv_type}: {bp_pair.bp1.chrom}:{bp_pair.bp1.position} - "
          f"{bp_pair.bp2.chrom}:{bp_pair.bp2.position} "
          f"confidence={bp_pair.confidence:.2f}")
    for bp in (bp_pair.bp1, bp_pair.bp2):
        print(f"  {bp.position}: strand={bp.strand} support={bp.support} "
              f"microhomology={bp.microhomology!r} inserted={bp.inserted_sequence!r}")
```

The first breakpoint of a variant dictionary is `breakpoint1`, else `start`, on `chrom1` or `chrom`; the second is `breakpoint2`, else `pos2`, else `end`, on `chrom2` (absent, `None` or `""` means the first chromosome). Read records carry `chrom`, `pos`, `cigar` (a string or a list of `(op_code, length)` tuples) and `name`, and optionally `is_reverse`, `seq` and `sa_tag`.

A read is located by the reference coordinate of its soft clip, not by where its alignment starts, so a long read that spans a breakpoint from far away counts when its clip lies within `window`. The reads that clip at the consensus position are the breakpoint's supporting reads, and both breakpoints are analysed from their own:

- `strand` is the majority `is_reverse` (`'-'` when more reads are reverse-complemented than not; `'+'` when no read reports it, or on a tie).
- `microhomology` is the longest `detect_microhomology` result over the window of each read's `seq` centred on that read's clip junction. It is a mirror-match heuristic on the read sequence, not true microhomology, which would need the reference sequence.
- `inserted_sequence` holds the read bases between the primary alignment and the read's single supplementary alignment. It is filled only when a supporting read has `seq`, `is_reverse`, `cigar` and an `sa_tag` with exactly one entry located at the other breakpoint, and both CIGARs add up to the length of `seq`; otherwise it is empty.

A breakpoint with no soft clip in its window keeps its input position with `support` 0.

### Supporting Functions

| Function | Signature | Purpose |
|----------|-----------|---------|
| `detect_microhomology` | `(breakpoint, flanking_seq, max_homology=50)` | Mirror-match run of identical bases either side of the midpoint of a sequence |
| `cluster_breakpoints` | `(breakpoints, max_distance=100)` | Group nearby breakpoints |
| `calculate_breakpoint_confidence` | `(evidence)` | Score based on read support and consistency |

## End-to-End Example

```python
from dataclasses import asdict

import numpy as np

from metainformant.structural_variants.detection.cnv import (
    detect_cnv_from_depth,
    calculate_log2_ratio,
)
from metainformant.structural_variants.detection.sv_calling import (
    call_structural_variants,
)
from metainformant.structural_variants.detection.breakpoints import (
    refine_breakpoints,
)

# Step 1: CNV detection. detect_cnv_from_depth takes read depths (finite, >= 0), not log2
# ratios, so pass 2**ratio for a tumor/normal comparison
log2r = calculate_log2_ratio(tumor_depth, normal_depth, gc_content=gc)
cnv_results = detect_cnv_from_depth({"chr1": np.exp2(log2r)}, significance=0.01)

# Step 2: SV calling from alignments
sv_calls = call_structural_variants(alignments, min_support=3)

# Step 3: Breakpoint refinement (variants are dicts)
refined = refine_breakpoints([asdict(sv) for sv in sv_calls], reads, window=200)

n_segments = sum(len(result.segments) for result in cnv_results.values())
print(f"{n_segments} CNV segments, {len(sv_calls)} SV calls, {len(refined)} refined breakpoint pairs")
```

## Configuration

Detection parameters are function arguments. Two have an `SV_` environment override that replaces the matching argument when the function runs:

| Parameter | Default | Env Variable |
|-----------|---------|-------------|
| `min_support` (`call_structural_variants`) | `3` | `SV_MIN_SUPPORT` |
| `significance` (`detect_cnv_from_depth`) | `0.01` | `SV_CNV_SIGNIFICANCE` |

`window_size`, `min_mapq`, `min_sv_size` and the rest have no environment variable. See [CONFIGURATION.md](CONFIGURATION.md) for every argument.

## See Also

- [Annotation](annotation.md) -- Annotate detected variants with gene/functional impact
- [Filtering](filtering.md) -- Quality filtering and multi-caller merging
- [Population](population.md) -- Population-scale genotyping and association
