# Getting Started with Structural Variant Analysis

Comprehensive guide for detecting, annotating, filtering, and visualizing structural variants (CNVs, deletions, duplications, inversions, translocations, insertions) from sequencing data.

## Quick Start (15 minutes)

### 1. Installation

```bash
cd /path/to/MetaInformAnt
uv venv
source .venv/bin/activate
uv pip install -e .
uv pip install numpy scipy pysam matplotlib seaborn
```

**System Requirements:**

| Component | Version | Purpose |
|-----------|---------|---------|
| Python | 3.11+ | Core runtime |
| NumPy | ≥1.24 | Numerical operations |
| SciPy | ≥1.10 | Statistics |
| pysam | ≥0.22 | BAM/CRAM I/O |

### 2. CNV Detection from BAM

**Input:** Sorted, indexed BAM/CRAM file.

```python-snippet
import pysam
import numpy as np
from metainformant.structural_variants.detection.cnv import detect_cnv_from_depth

# Compute per-chromosome depth in 1 kb windows
window_size = 1000
depth_data = {}
bam = pysam.AlignmentFile("sample.bam", "rb")

for chrom in bam.references:
    n_bins = bam.lengths[bam.get_tid(chrom)] // window_size + 1
    depths = np.zeros(n_bins, dtype=np.float64)
    for read in bam.fetch(chrom):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        bin_idx = read.reference_start // window_size
        if bin_idx < n_bins:
            depths[bin_idx] += 1.0
    depth_data[chrom] = depths
bam.close()

# Detect CNVs using Circular Binary Segmentation (CBS)
result = detect_cnv_from_depth(
    depth_data=depth_data,
    window_size=1000,
    method="segmentation",
    significance=0.01,
    ploidy=2,
    min_segment_bins=3,
    merge_distance=1000,
)

# Inspect results
from metainformant.core import io
total = sum(len(r.segments) for r in result.values())
print(f"Detected {total} CNV segments")

for chrom, cnv_result in result.items():
    for seg in cnv_result.segments:
        print(f"  {seg.chrom}:{seg.start:,}-{seg.end:,}  "
              f"state={seg.state} CN={seg.cn} conf={seg.confidence:.2f}")

# Save
output = {c: [{"chrom": s.chrom, "start": s.start, "end": s.end,
               "state": s.state, "cn": s.cn, "confidence": s.confidence}
              for s in r.segments] for c, r in result.items()}
io.dump_json(output, "cnv_calls.json")
```

**Returns:** `CNVResult` with `segments: list[CNVSegment]` per chromosome.

### 3. SV Calling from Alignments

Detect deletions, duplications, inversions, translocations from split-read and discordant pair evidence.

```python
from metainformant.structural_variants.detection.sv_calling import call_structural_variants

# Build alignment record list from BAM
alignments = []
bam = pysam.AlignmentFile("sample.bam", "rb")
for read in bam.fetch():
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        continue
    rec = {
        "name": read.query_name,
        "chrom": bam.get_reference_name(read.reference_id),
        "pos": read.reference_start,
        "mapq": read.mapping_quality,
        "cigar": read.cigarstring,
        "read_length": read.query_length,
        "is_reverse": read.is_reverse,
    }
    if read.is_paired and not read.mate_is_unmapped:
        rec["mate_chrom"] = bam.get_reference_name(read.mate_reference_id)
        rec["mate_pos"] = read.next_reference_start
        rec["mate_is_reverse"] = read.mate_is_reverse
        rec["insert_size"] = abs(read.template_length) if read.template_length else 0
    if read.has_tag("SA"):
        rec["sa_tag"] = read.get_tag("SA")
    if read.query_sequence:
        rec["seq"] = read.query_sequence
    alignments.append(rec)
bam.close()

# Call SVs
variants = call_structural_variants(
    alignments=alignments,
    min_support=3,           # Minimum supporting reads
    min_mapq=20,             # Minimum mapping quality
    min_sv_size=50,          # Minimum size in bp
)

print(f"Detected {len(variants)} SVs")
from collections import Counter
for sv_type, count in Counter(v.sv_type.value for v in variants).most_common():
    print(f"  {sv_type}: {count}")

for sv in variants:
    print(f"  [{sv.sv_type.value}] {sv.chrom}:{sv.start:,}-{sv.end:,}  "
          f"size={sv.size:,}bp qual={sv.quality:.1f} "
          f"support={sv.evidence.split_reads}+{sv.evidence.discordant_pairs}")
```

**Returns:** `list[StructuralVariant]` with fields:
- `chrom`, `start`, `end`, `sv_type` (enum: DEL/DUP/INV/TRA/INS/BND/UNKNOWN)
- `size`, `quality`, `genotype`, `evidence: SVEvidence`

### 4. Breakpoint Refinement

Refine approximate breakpoints to single-base resolution using local split-read evidence:

```python
from metainformant.structural_variants.detection.breakpoints import refine_breakpoints

# Fetch local reads around breakpoints (window ~200 bp)
local_reads = []
bam = pysam.AlignmentFile("sample.bam", "rb")
for sv in variants:
    for chrom in {sv.chrom, getattr(sv, "chrom2", sv.chrom)}:
        for read in bam.fetch(chrom, max(0, sv.start - 500), sv.end + 500):
            if read.is_unmapped or read.is_secondary:
                continue
            local_reads.append({
                "chrom": chrom,
                "pos": read.reference_start,
                "cigar": read.cigarstring,
                "seq": read.query_sequence or "",
                "name": read.query_name,
                "is_reverse": read.is_reverse,
            })
bam.close()

# Refine
refined = refine_breakpoints(
    variants=[v.__dict__ for v in variants],
    reads=local_reads,
    window=200,
    min_clip=5,
)

for bp_pair in refined:
    print(f"  {bp_pair.sv_type}: {bp_pair.bp1.chrom}:{bp_pair.bp1.position} → "
          f"{bp_pair.bp2.chrom}:{bp_pair.bp2.position}  "
          f"conf={bp_pair.confidence:.2f}")
    if bp_pair.bp1.microhomology:
        print(f"    Microhomology: {bp_pair.bp1.microhomology}")
```

**Returns:** `list[BreakpointPair]` with `bp1, bp2: Breakpoint` objects including:
- `position`, `strand`, `confidence`, `support`, `microhomology`, `read_names`

### 5. Gene & Regulatory Overlap

```python
from metainformant.structural_variants.annotation.overlap import (
    annotate_gene_overlap,
    GenomicInterval,
)
from metainformant.core import io

# Load gene annotations
gene_db = io.load_json("genes.json")  # List of {"chrom", "start", "end", "name", ...}
# Or construct:
# gene_db = [GenomicInterval("chr1", 11873, 14409, "DDX11L1", "gene", "+", {})]

# Load SV calls
variants = io.load_json("sv_calls.json")

# Annotate gene overlaps
annotated = annotate_gene_overlap(variants, gene_db)

for v in annotated:
    if v["n_genes_affected"] > 0:
        print(f"SVs at {v['chrom']}:{v['start']} overlaps genes: {v['overlapping_genes']}")
        for overlap in v["gene_overlaps"]:
            print(f"  {overlap.feature.name}: {overlap.overlap_fraction_variant:.1%} of variant")

# Regulatory elements (enhancers, promoters)
reg_db = io.load_json("enhancers.json")
annotated = annotate_regulatory_overlap(annotated, reg_db)
```

**Returns:** Variant dicts with added fields:
- `overlapping_genes: list[str]`
- `gene_overlaps: list[OverlapResult]`
- `n_genes_affected: int`

`OverlapResult` includes: `overlap_bp`, `overlap_fraction_variant`, `overlap_fraction_feature`, `relationship`.

### 6. Functional Impact Prediction

```python
from metainformant.structural_variants.annotation.functional_impact import (
    predict_functional_impact,
    assess_dosage_sensitivity,
    predict_tad_disruption,
    score_pathogenicity,
)

# Load dosage databases (ClinGen, gnomAD)
hi_db = io.load_json("haploinsufficiency_scores.json")   # {gene: HI_score ∈ [0,1]}
pli_db = io.load_json("gnomad_pli.json")                # {gene: pLI ∈ [0,1]}
tad_bounds = io.load_json("tad_boundaries.json")

# Predict impact for all annotated variants
impacts = []
for variant in annotated:
    impact = predict_functional_impact(
        variant=variant,
        gene_annotations=gene_db,
        haploinsufficiency_db=hi_db,
        tad_boundaries=tad_bounds,
    )
    impacts.append(impact)

for imp in impacts:
    print(f"\n{imp.variant_id}: {imp.impact_level} ({imp.impact_type})")
    print(f"  Pathogenicity score: {imp.pathogenicity_score:.3f}")
    print(f"  Affected genes: {', '.join(imp.affected_genes)}")
    if imp.dosage_sensitive_genes:
        print(f"  Dosage-sensitive: {', '.join(imp.dosage_sensitive_genes)}")
    if imp.tad_disrupted:
        print(f"  ⚠ TAD disrupted ({imp.details.get('tad_prediction',{}).get('n_boundaries',0)} boundaries)")
```

**Impact levels:** `HIGH`, `MODERATE`, `LOW`, `MODIFIER`

**Impact types:** `gene_disruption`, `dosage_loss`, `dosage_gain`, `gene_fusion_candidate`, `tad_disruption`, `gene_duplication`, `gene_inversion`, `large_intergenic`, `benign`.

**Returns:** `FunctionalImpact` with:
- `impact_level`, `impact_type`, `affected_genes`, `dosage_sensitive_genes`
- `tad_disrupted: bool`
- `pathogenicity_score: float` ∈ [0, 1]

### 7. Quality Filtering

```python
from metainformant.structural_variants.filtering.quality_filter import (
    filter_by_quality,
    filter_by_size,
    filter_by_frequency,
    apply_blacklist,
)

variants = io.load_json("sv_calls.json")

# Quality: Phred score, read support, MAPQ
filtered, q_stats = filter_by_quality(
    variants,
    min_qual=20.0,
    min_support=3,
    min_mapq=30.0,
)
print(f"Quality: {q_stats.input_count}→{q_stats.output_count} ({q_stats.pass_rate:.1%})")

# Size
filtered, s_stats = filter_by_size(
    filtered,
    min_size=50,
    max_size=50_000_000,
)
print(f"Size: {s_stats.input_count}→{s_stats.output_count}")

# Population frequency (remove common variants)
pop_db = io.load_json("gnomad_sv.json")
filtered, f_stats = filter_by_frequency(
    filtered,
    population_db=pop_db,
    max_af=0.01,
    match_window=200,
    min_reciprocal_overlap=0.5,
)
print(f"Frequency: removed {f_stats.filtered_count} common variants")

# Blacklist (segmental dups, artifactual regions)
blacklist = io.load_bed("encode_blacklist.bed")
filtered, b_stats = apply_blacklist(filtered, blacklist_regions=blacklist)
print(f"Blacklist: removed {b_stats.filtered_count}")

print(f"\nFinal: {len(filtered)} high-confidence variants")
io.dump_json(filtered, "filtered_svs.json")
```

Each filter returns `(variants, FilterStats)` where `FilterStats` has:
`input_count`, `output_count`, `filtered_count`, `filter_name`, `parameters`, `pass_rate`.

### 8. Multi-Caller Merging

Merge SV callsets from Manta, Delly, Lumpy using reciprocal overlap:

```python
from metainformant.structural_variants.filtering.merge import merge_callsets

callsets = {
    "Manta": io.load_json("manta.json"),
    "Delly": io.load_json("delly.json"),
    "Lumpy": io.load_json("lumpy.json"),
}

merged, stats = merge_callsets(
    callsets=callsets,
    min_overlap=0.5,      # 50% reciprocal overlap
    type_match=True,       # SV types must match
)

print(f"Merged: {stats.n_input_variants} → {stats.n_output_variants}")
print(f"  Multi-caller (≥2): {stats.n_multi_caller}")
print(f"  Single-caller: {stats.n_single_caller}")

for caller, count in stats.caller_counts.items():
    print(f"  {caller} contributed: {count}")

for mv in merged:
    if mv.n_callers >= 2:
        print(f"  CONSENSUS [{mv.n_callers}x]: {mv.chrom}:{mv.start}-{mv.end} "
              f"type={mv.sv_type} qual={mv.consensus_quality:.1f}")
```

**Returns:** `list[MergedVariant]` + `MergeStats`.

`MergedVariant` includes: `chrom, start, end, sv_type, n_callers, caller_names, consensus_quality, support_variants`.

Alternative: `survivor_merge()` for SURVIVOR-style distance-based merging:
```python
merged = survivor_merge(
    vcf_files=["manta.vcf", "delly.vcf", "lumpy.vcf"],
    max_distance=1000,
    min_callers=2,
)
```

### 9. Population-Scale Analysis

```python
from metainformant.structural_variants.population.sv_population import (
    genotype_sv_population,
    sv_allele_frequency,
    sv_population_structure,
    sv_association_test,
)

# Genotype a set of reference SVs across a cohort
samples = ["sample1", "sample2", "sample3", ...]  # Sample IDs
sv_calls = io.load_json("reference_svs.json")

genotype_result = genotype_sv_population(
    sv_calls=sv_calls,
    samples=samples,
    method="depth",  # or "split"
)
gt_matrix = genotype_result["genotype_matrix"]  # n_svs × n_samples

# Population allele frequencies (optionally stratified)
pop_labels = ["EUR", "EUR", "AFR", "EAS", ...]  # Per-sample ancestry
freq = sv_allele_frequency(gt_matrix, sample_labels=pop_labels)

print(f"Mean overall AF: {np.mean(freq['frequencies']):.4f}")
print(f"Polymorphic SVs (MAF>0): {freq['n_polymorphic']}")
for pop, afs in freq.get("population_frequencies", {}).items():
    print(f"  {pop} mean AF: {np.mean(afs):.4f}")

# PCA for population structure
pca = sv_population_structure(gt_matrix, n_components=10)
print(f"PC1 variance: {pca['explained_variance'][0]:.2%}")
```

### 10. Visualization

```python-snippet
from metainformant.structural_variants.visualization.plots import (
    plot_circos,
    plot_coverage_track,
    plot_size_distribution,
)

chrom_sizes = {"chr1": 248956422, "chr2": 242193529, ...}

# Circos plot — whole-genome chord diagram
fig_circos = plot_circos(
    variants=filtered,
    chromosomes=chrom_sizes,
    output_path="sv_circos.png",
    title="Structural Variants — Circos View",
)
# Saves PNG; returns matplotlib Figure

# Coverage track — CNV profile per chromosome
from cnv_results import result  # From step 1
for chrom in ["chr1", "chr5"]:
    fig = plot_coverage_track(
        coverage=result[chrom].log2_ratios,
        variants=[v for v in filtered if v["chrom"] == chrom],
        region=(chrom, 0, chrom_sizes[chrom]),
        output_path=f"coverage_{chrom}.png",
        bin_size=1000,
    )

# Size distribution histogram
fig_hist = plot_size_distribution(
    variants=filtered,
    output_path="sv_size_distribution.png",
)
```

### 11. Complete Pipeline (Single Command)

```python-snippet
from metainformant.structural_variants.pipeline import run_full_pipeline

results = run_full_pipeline(
    bam_path="sample.bam",
    output_dir="output/",
    reference="hg38",
    min_support=3,
    min_quality=20,
    max_population_af=0.01,
    run_annotation=True,
    run_filtering=True,
    run_visualization=True,
)
print(f"Pipeline complete: {len(results['filtered_svs'])} high-confidence variants")
```

### 12. Key Return Types

All detection functions return dataclass objects:

```python
from dataclasses import asdict

# CNVResult
result: CNVResult
result.segments        # list[CNVSegment]
result.log2_ratios     # list[float] | None
result.bin_size        # int
result.method          # str
result.parameters      # dict

# CNVSegment
seg: CNVSegment
seg.chrom, start, end  # genomic coordinates
seg.mean_log2ratio     # float
seg.state              # 'DEL' | 'DUP' | 'NEUTRAL' | 'AMP' | 'HOMODEL'
seg.cn                 # estimated integer copy number
seg.confidence         # [0, 1]

# StructuralVariant
sv: StructuralVariant
sv.chrom, start, end
sv.sv_type             # SVType enum
sv.size, quality, genotype
sv.evidence            # SVEvidence object

# BreakpointPair
bp_pair: BreakpointPair
bp_pair.bp1, bp2       # Breakpoint objects
bp_pair.bp1.position, strand, confidence, microhomology
```

---

## Configuration

All structural variant parameters respect the `SV_` environment variable prefix:

```bash
export SV_MIN_SUPPORT=3          # Minimum supporting reads
export SV_MIN_QUAL=20.0          # Minimum quality score
export SV_MIN_SV_SIZE=50         # Minimum SV size (bp)
export SV_CNV_SIGNIFICANCE=0.01  # CBS p-value threshold
```

Or via Python config:

```python
from metainformant import config
config.set("structural_variants.min_support", 5)
config.set("structural_variants.min_quality", 30.0)
```

See **[CONFIGURATION.md](CONFIGURATION.md)** for complete settings.

---

## Expected Outputs

All steps produce:

- **JSON**: Machine-readable results (Unix-friendly)
- **VCF**: Standard variant format (compatible with IGV, bcftools)
- **PNG/PDF**: Publication-quality figures (Circos, coverage, histograms)
- **CSV**: Spreadsheet tables for downstream analysis

Console logs (INFO-level) show real-time progress. Enable DEBUG for per-variant details:
```bash
export CORE_LOG_LEVEL=DEBUG
```

---

## Common Workflows

### Workflow A: Somatic CNV Calling (Tumor vs. Normal)

```python
# Load matched normal depth
normal_depth = compute_depth("normal.bam", window_size=1000)
tumor_depth = compute_depth("tumor.bam", window_size=1000)

# Compute log2 ratio
log2_ratios = calculate_log2_ratio(tumor_depth, normal_depth, gc_content=gc_bias)

# Call CNVs
cnv_result = detect_cnv_from_depth(
    depth_data=log2_ratios,
    window_size=1000,
    ploidy=2,
    significance=0.001,  # More stringent for somatic
)
```

### Workflow B: Germline SV Discovery

```python
# Single sample, high sensitivity
variants = call_structural_variants(
    alignments=alignments,
    min_support=2,          # Lower for germline heterozygotes
    min_mapq=30,
    min_sv_size=50,
)
# Annotate, filter, visualize
```

### Workflow C: Cohort Population Analysis

```python
# Reference SV panel (e.g., from 1000G)
reference_svs = io.load_json("sv_reference_panel.json")

# Genotype all samples
all_gt = {}
for sample_id, bam_path in sample_bams.items():
    alignments = extract_alignments(bam_path)
    gt = genotype_sv_population(reference_svs, [sample_id], method="depth")
    all_gt[sample_id] = gt["genotype_matrix"][0]  # One row per SV

# Build genotype matrix: n_svs × n_samples
gt_matrix = np.column_stack([all_gt[s] for s in sample_order])

# PCA, association, etc.
```

---

## What's Next?

- **Deep Dive**: Read `detection.md`, `annotation.md`, `filtering.md`, `population.md`
- **Architecture**: See `ARCHITECTURE.md` for data flow and design patterns
- **Real Examples**: Check `EXAMPLES.md` for clinical, research, and large-cohort scenarios
- **Configuration**: Fine-tune performance in `CONFIGURATION.md`
- **Performance**: Scale to >10K samples using `PERFORMANCE.md` guidelines
- **Debug**: Diagnose issues in `TROUBLESHOOTING.md`
