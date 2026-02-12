# Structural Variants: Detection

CNV detection, structural variant calling, and breakpoint refinement.

## Concept Overview

The detection submodule identifies genomic structural variants through three
complementary approaches:

- **CNV Detection** -- Copy number variants from read depth using circular binary segmentation
- **SV Calling** -- Structural variants from split reads and discordant read pairs
- **Breakpoint Refinement** -- Base-pair resolution breakpoint determination with microhomology detection

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
| `start` | `int` | Start position |
| `end` | `int` | End position |
| `mean_log2ratio` | `float` | Mean log2 ratio for segment |
| `n_bins` | `int` | Number of coverage bins |
| `state` | `str` | `"neutral"`, `"deletion"`, `"duplication"`, `"amplification"`, `"homozygous_deletion"` |
| `cn` | `int` | Estimated copy number |
| `confidence` | `float` | Call confidence score |

**CNVResult** -- Complete result set with segments and parameters.

### detect_cnv_from_depth

Primary entry point for CNV detection from read depth data.

```python
from metainformant.structural_variants.detection.cnv import detect_cnv_from_depth

result = detect_cnv_from_depth(
    depth_data={"chr1": np.array([100, 102, 98, 50, 48, 52, 99, 101])},
    window_size=1000,
    method="segmentation",       # CBS algorithm
    significance=0.01,
    ploidy=2,
    min_segment_bins=3,
    merge_distance=1000
)

for segment in result.segments:
    print(f"{segment.chrom}:{segment.start}-{segment.end} "
          f"CN={segment.cn} ({segment.state})")
```

### Supporting Functions

| Function | Signature | Purpose |
|----------|-----------|---------|
| `segment_coverage` | `(coverage_array, significance=0.01, max_iterations=100, min_width=2)` | CBS segmentation |
| `call_cnv_states` | `(segments, ploidy=2, del_threshold=-0.3, dup_threshold=0.3, amp_threshold=1.0, homodel_threshold=-1.5)` | Assign CN states |
| `merge_adjacent_segments` | `(segments, max_gap=1000)` | Merge nearby segments with same state |
| `calculate_log2_ratio` | `(tumor_depth, normal_depth, gc_content=None, pseudocount=1.0)` | Compute log2 ratio with optional GC correction |

## SV Calling

### Core Types

**SVType** enum: `DEL`, `DUP`, `INV`, `TRA`, `INS`, `BND`, `UNKNOWN`

**StructuralVariant** -- A detected SV with type, coordinates, and genotype.

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
    print(f"{sv.sv_type.name}: {sv.chrom}:{sv.start}-{sv.end} "
          f"support={sv.support}")
```

### Supporting Functions

| Function | Signature | Purpose |
|----------|-----------|---------|
| `detect_split_reads` | `(reads, min_clip=20)` | Find soft-clipped reads indicating breakpoints |
| `detect_discordant_pairs` | `(pairs, insert_size_stats, n_std=4.0)` | Find anomalous insert sizes/orientations |
| `classify_sv_type` | `(evidence)` | Determine SV type from evidence pattern |
| `genotype_sv` | `(variant, reads, min_mapq=20)` | Genotype: 0/0, 0/1, or 1/1 |

## Breakpoint Refinement

### Core Types

**Breakpoint** -- Single breakpoint with position, strand, confidence, microhomology, and supporting reads.

**BreakpointPair** -- Two linked breakpoints defining an SV.

### refine_breakpoints

Improve breakpoint resolution using local read evidence.

```python
from metainformant.structural_variants.detection.breakpoints import refine_breakpoints

refined = refine_breakpoints(
    variants=sv_calls,
    reads=local_reads,
    window=200,
    min_clip=5
)

for bp_pair in refined:
    print(f"{bp_pair.sv_type}: {bp_pair.bp1.chrom}:{bp_pair.bp1.position} - "
          f"{bp_pair.bp2.chrom}:{bp_pair.bp2.position} "
          f"confidence={bp_pair.confidence:.2f}")
```

### Supporting Functions

| Function | Signature | Purpose |
|----------|-----------|---------|
| `detect_microhomology` | `(breakpoint, flanking_seq, max_homology=50)` | Identify sequence homology at breakpoints |
| `cluster_breakpoints` | `(breakpoints, max_distance=100)` | Group nearby breakpoints |
| `calculate_breakpoint_confidence` | `(evidence)` | Score based on read support and consistency |

## End-to-End Example

```python
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

# Step 1: CNV detection from depth
log2r = calculate_log2_ratio(tumor_depth, normal_depth, gc_content=gc)
cnv_result = detect_cnv_from_depth({"chr1": log2r}, significance=0.01)

# Step 2: SV calling from alignments
sv_calls = call_structural_variants(alignments, min_support=3)

# Step 3: Breakpoint refinement
refined = refine_breakpoints(sv_calls, reads, window=200)

# Combine results
all_variants = cnv_result.segments + refined
print(f"Total variants: {len(all_variants)}")
```

## Configuration

Detection parameters use the `SV_` environment prefix:

| Parameter | Default | Env Variable |
|-----------|---------|-------------|
| `window_size` | `1000` | `SV_WINDOW_SIZE` |
| `significance` | `0.01` | `SV_SIGNIFICANCE` |
| `min_support` | `3` | `SV_MIN_SUPPORT` |
| `min_mapq` | `20` | `SV_MIN_MAPQ` |
| `min_sv_size` | `50` | `SV_MIN_SV_SIZE` |

## See Also

- [Annotation](annotation.md) -- Annotate detected variants with gene/functional impact
- [Filtering](filtering.md) -- Quality filtering and multi-caller merging
- [Population](population.md) -- Population-scale genotyping and association
