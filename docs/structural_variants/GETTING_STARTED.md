# Getting Started with Structural Variant Analysis

Comprehensive guide for detecting, annotating, filtering, and visualizing structural variants (CNVs, deletions, duplications, inversions, translocations, insertions) from sequencing data.

## Quick Start (15 minutes)

### 1. Installation

```bash
cd /path/to/MetaInformAnt
uv venv
source .venv/bin/activate
uv pip install -e .
uv pip install pysam   # only needed for the BAM/CRAM-reading steps below
```

NumPy, SciPy, pandas, matplotlib and seaborn are core dependencies of the package, so `uv pip install -e .` already provides them.

**System Requirements:**

| Component | Version | Purpose |
|-----------|---------|---------|
| Python | 3.11+ | Core runtime |
| NumPy | ≥1.24 | Numerical operations |
| SciPy | ≥1.11 | Installed with the package; the structural variant functions do not call it |
| pysam | ≥0.21 (`external-tools` extra) | BAM/CRAM I/O |

### 2. CNV Detection from BAM

**Input:** Sorted, indexed BAM/CRAM file.

```python-snippet
import pysam
import numpy as np
from metainformant.structural_variants.detection.cnv import detect_cnv_from_depth

# Per chromosome, count the reads starting in each 1 kb window (a simple depth proxy)
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
total = sum(s.state != "NEUTRAL" for r in result.values() for s in r.segments)
print(f"Detected {total} gain/loss segments")

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

**Returns:** `dict[str, CNVResult]` keyed by chromosome; each `CNVResult` holds `segments: list[CNVSegment]` (plus `log2_ratios`, `bin_size`, `method`, `parameters`). `segments` includes the `NEUTRAL` stretches as well as gains and losses, so filter on `seg.state != "NEUTRAL"` to keep only the calls. Chromosomes with fewer than `min_segment_bins` bins are skipped. The CBS step grows at least quadratically with the number of bins per chromosome, so 1 kb windows are only practical for small chromosomes or regions; use a coarser `window_size` for whole human chromosomes (see [PERFORMANCE.md](PERFORMANCE.md)).

### 3. SV Calling from Alignments

Detect deletions, duplications, inversions, translocations from split-read and discordant pair evidence.

```python
import pysam
from metainformant.structural_variants.detection.sv_calling import call_structural_variants

# Build alignment record list from BAM (these records are held in memory, roughly 0.6-0.7 KB each;
# for a whole genome, process one chromosome or region at a time)
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

# The annotation, filtering and plotting functions below take dicts, so convert the calls
# and save them; the next steps load this file.
from dataclasses import asdict
from metainformant.core import io

sv_dicts = []
for sv in variants:
    record = asdict(sv)
    record["sv_type"] = sv.sv_type.value   # plain string instead of the SVType enum
    sv_dicts.append(record)
io.dump_json(sv_dicts, "sv_calls.json")
```

**Returns:** `list[StructuralVariant]` with fields:
- `chrom`, `start`, `end`, `sv_type` (enum: DEL/DUP/INV/TRA/INS/BND/UNKNOWN)
- `size`, `quality`, `genotype`, `evidence: SVEvidence`, `filter_status`, `chrom2` (`""` for intra-chromosomal calls), `info`

`evidence` holds `split_reads`, `discordant_pairs`, `breakpoint1`/`breakpoint2` and `chrom1`/`chrom2`; for inter-chromosomal (`TRA`) calls take the two breakpoints from there rather than from `start`/`end`. `INS` calls have their two breakpoints less than 50 bp apart, so the default `min_sv_size=50` removes them (lower it to keep insertions; see Example 15 in [EXAMPLES.md](EXAMPLES.md)).

### 4. Breakpoint Refinement

Refine approximate breakpoints to single-base resolution using local split-read evidence:

```python
from metainformant.structural_variants.detection.breakpoints import refine_breakpoints

# Fetch local reads around each breakpoint. sv.evidence holds both breakpoints with their own
# chromosome (sv.chrom2 is "" for intra-chromosomal calls).
local_reads = []
bam = pysam.AlignmentFile("sample.bam", "rb")
for sv in variants:
    ev = sv.evidence
    for chrom, pos in ((ev.chrom1, ev.breakpoint1), (ev.chrom2 or ev.chrom1, ev.breakpoint2)):
        for read in bam.fetch(chrom, max(0, pos - 500), pos + 500):
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

# Refine: refine_breakpoints() takes variant dicts (chrom1/breakpoint1, chrom2/breakpoint2, sv_type)
refined = refine_breakpoints(
    variants=[
        {"chrom1": v.evidence.chrom1, "breakpoint1": v.evidence.breakpoint1,
         "chrom2": v.evidence.chrom2 or v.evidence.chrom1, "breakpoint2": v.evidence.breakpoint2,
         "sv_type": v.sv_type.value}
        for v in variants
    ],
    reads=local_reads,
    window=200,   # clips within 200 bp of a breakpoint, from reads starting within 700 bp
    min_clip=5,
)

for bp_pair in refined:
    print(f"  {bp_pair.sv_type}: {bp_pair.bp1.chrom}:{bp_pair.bp1.position} → "
          f"{bp_pair.bp2.chrom}:{bp_pair.bp2.position}  "
          f"conf={bp_pair.confidence:.2f}")
    if bp_pair.bp1.microhomology:
        print(f"    Microhomology: {bp_pair.bp1.microhomology}")
```

**Returns:** `list[BreakpointPair]` (`bp1`, `bp2`, `sv_type`, `confidence`) with `Breakpoint` objects including:
- `position`, `strand`, `confidence`, `support`, `microhomology`, `read_names`

A breakpoint with `support == 0` found no soft-clipped read and keeps its input position. `microhomology` is only set on `bp1`, from the first read that carries a `seq` on that chromosome, assuming the breakpoint sits in the middle of that sequence, so treat it as meaningful only for breakpoint-centred flanking sequences. `Breakpoint.inserted_sequence` exists but is never filled in.

### 5. Gene & Regulatory Overlap

```python
from metainformant.structural_variants.annotation.overlap import (
    annotate_gene_overlap,
    annotate_regulatory_overlap,
    GenomicInterval,
)
from metainformant.core import io

# Load gene annotations
gene_db = io.load_json("genes.json")  # List of {"chrom", "start", "end", "name", ...}
# Or construct (annotate_gene_overlap accepts these, but predict_functional_impact in step 6 needs the dicts):
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

`OverlapResult` includes: `feature` (the `GenomicInterval`), `overlap_bp`, `overlap_fraction_variant`, `overlap_fraction_feature`, `relationship` (`"overlap"`, `"contains"` when the SV spans the feature, `"contained_in"` when it lies inside it).

`annotate_regulatory_overlap` adds `overlapping_regulatory: list[str]`, `regulatory_overlaps: list[OverlapResult]`, `regulatory_types: set[str]` and `n_regulatory_affected: int`. Both annotators keep the original keys of each variant dict; the `OverlapResult` objects and the `regulatory_types` set are not JSON-serializable, so drop `gene_overlaps` / `regulatory_overlaps` and convert `regulatory_types` to a list (or drop it) before `io.dump_json`.

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

# predict_functional_impact() only takes the haploinsufficiency scores. For pLI (and
# triplosensitivity) look the gene up directly: HI >= 0.5, pLI >= 0.9 or (with loeuf_db) LOEUF < 0.35
# counts as haploinsufficient.
ds = assess_dosage_sensitivity("GENE1", haploinsufficiency_db=hi_db, pli_db=pli_db)
print(ds.is_haploinsufficient, ds.haploinsufficiency_score, ds.pli_score)
```

**Impact levels:** `HIGH`, `MODERATE`, `LOW`, `MODIFIER`

**Impact types:** `dosage_loss`, `dosage_gain`, `gene_disruption`, `gene_fusion_candidate`, `gene_duplication`, `gene_inversion`, `gene_insertion`, `tad_disruption`, `large_intergenic`, `intergenic`, `benign`. A gene counts as dosage-sensitive for a deletion when its HI score is at least 0.5; because no triplosensitivity database can be passed, `dosage_gain` is not currently produced and duplications over genes come back as `gene_duplication`.

**Returns:** `FunctionalImpact` with:
- `variant_id` (the variant's `id` key, or `chrom:start-end`), `impact_level`, `impact_type`, `affected_genes`, `dosage_sensitive_genes`
- `tad_disrupted: bool`
- `pathogenicity_score: float` ∈ [0, 1]
- `details: dict` (`sv_type`, `size`, `n_overlapping_genes`, and `tad_prediction` when TAD boundaries were given)

### 7. Quality Filtering

```python
from metainformant.structural_variants.filtering.quality_filter import (
    filter_by_quality,
    filter_by_size,
    filter_by_frequency,
    apply_blacklist,
)

variants = io.load_json("sv_calls.json")

# Quality: Phred score and read support. min_mapq only filters dicts that carry a "mapq" key
# (the calls saved in step 3 do not, so it removes nothing from them)
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
pop_db = io.load_json("gnomad_sv.json")   # [{"chrom", "start", "end", "sv_type", "af"}, ...]
filtered, f_stats = filter_by_frequency(
    filtered,
    population_db=pop_db,
    max_af=0.01,
    match_window=200,
    min_reciprocal_overlap=0.5,
)
print(f"Frequency: removed {f_stats.filtered_count} common variants")

# Blacklist (segmental dups, artifactual regions): (chrom, start, end) tuples or dicts
blacklist = [(chrom, int(start), int(end)) for chrom, start, end, *_ in io.read_tsv("encode_blacklist.bed")]
filtered, b_stats = apply_blacklist(filtered, blacklist_regions=blacklist)
print(f"Blacklist: removed {b_stats.filtered_count}")

print(f"\nFinal: {len(filtered)} high-confidence variants")
io.dump_json(filtered, "filtered_svs.json")
```

Each filter returns `(variants, FilterStats)` where `FilterStats` has:
`input_count`, `output_count`, `filtered_count`, `filter_name`, `parameters`, `pass_rate` (a property). It holds counts only, not the removed variants. All four filters take variant dicts (not `StructuralVariant` objects). `filter_by_quality` takes the supporting-read count from a `support`/`n_support` key or, failing that, from the `evidence` entry (`split_reads` + `discordant_pairs`), and it honours the `SV_MIN_QUAL` and `SV_MIN_SUPPORT` environment variables, which override `min_qual` and `min_support`. `filter_by_frequency` accepts the population database as a list of dicts (as above) or a dict keyed by `"chrom:start-end"`; see Example 6 in [EXAMPLES.md](EXAMPLES.md).

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

**Returns:** `(list[MergedVariant], MergeStats)`.

`MergedVariant` includes: `chrom, start, end, sv_type, size, n_callers, caller_names, consensus_quality, genotype, support_variants`. `MergeStats` has `n_input_callsets`, `n_input_variants`, `n_output_variants`, `n_multi_caller`, `n_single_caller`, `caller_counts`. Matching is transitive: if A matches B and B matches C, all three merge even when A and C do not meet `min_overlap` (`survivor_merge` does the same with `max_distance`).

Alternative: `survivor_merge()` for SURVIVOR-style distance-based merging. It returns the same `(merged, stats)` tuple, treats each VCF file as one caller (`caller_0`, `caller_1`, ...) and keeps only variants seen by at least `min_callers` callers. Unlike `merge_callsets`, it leaves `MergedVariant.genotype` at `"./."` and derives `consensus_quality` from the number of callers only:
```python
from metainformant.structural_variants.filtering.merge import survivor_merge

merged, stats = survivor_merge(
    vcf_files=["manta.vcf", "delly.vcf", "lumpy.vcf"],
    max_distance=1000,
    min_callers=2,
)
```

### 9. Population-Scale Analysis

```python
import numpy as np
from metainformant.structural_variants.population.sv_population import (
    genotype_sv_population,
    sv_allele_frequency,
    sv_population_structure,
    sv_association_test,
)

# Genotype a set of reference SVs across a cohort
samples = ["sample1", "sample2", "sample3", ...]  # Sample IDs
# Each SV dict carries the per-sample evidence that is converted to genotypes:
# {"chrom", "start", "end", "sv_type",
#  "samples": {sample_id: {"depth_ratio": 0.5, "split_reads": 6, "discordant_pairs": 4, "total_reads": 30}}}
sv_calls = io.load_json("reference_svs.json")

genotype_result = genotype_sv_population(
    sv_calls=sv_calls,
    samples=samples,
    method="depth",  # "depth": depth_ratio; "split": split_reads / discordant_pairs / total_reads
)
gt_matrix = genotype_result["genotype_matrix"]  # n_svs × n_samples (list of lists; 0/1/2 alt alleles)

# Population allele frequencies (optionally stratified)
pop_labels = ["EUR", "EUR", "AFR", "EAS", ...]  # Per-sample ancestry
freq = sv_allele_frequency(gt_matrix, sample_labels=pop_labels)

print(f"Mean overall AF: {np.mean(freq['frequencies']):.4f}")
print(f"Polymorphic SVs (MAF>0): {freq['n_polymorphic']}")
for pop, afs in freq.get("population_frequencies", {}).items():
    print(f"  {pop} mean AF: {np.mean(afs):.4f}")

# PCA for population structure (pca["pcs"] is n_samples × n_components)
pca = sv_population_structure(gt_matrix, n_components=10)
print(f"PC1 variance: {pca['variance_explained'][0]:.2%}")

# Association of one SV with a phenotype, adjusted for the top PCs
# (phenotypes: one value per sample, in the same order as `samples`)
rec = sv_association_test(gt_matrix[0], phenotypes, covariates=pca["pcs"][:, :5])
print(rec["beta"], rec["p_value"], rec["odds_ratio"])   # odds_ratio is None for a non-binary phenotype
```

**Returns:**
- `genotype_sv_population`: `genotype_matrix`, `quality_scores` (both n_svs × n_samples) and `allele_frequencies` (one per SV). The function does not read BAM files; SVs or samples without evidence default to reference genotypes.
- `sv_allele_frequency`: `frequencies`, `maf`, `n_polymorphic`, `sfs`, and `population_frequencies` when `sample_labels` is given.
- `sv_population_structure`: `pcs`, `eigenvalues`, `variance_explained`, `loadings`.
- `sv_association_test`: `beta`, `se`, `p_value`, `odds_ratio` (binary 0/1 phenotypes only), `n_samples`.

### 10. Visualization

```python-snippet
from metainformant.structural_variants.visualization.plots import (
    plot_circos,
    plot_coverage_track,
    plot_sv_size_distribution,
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

# Coverage track — read depth per chromosome, with the SVs drawn underneath
# (depth_data: the {chrom: per-bin depth} dict built in step 2)
for chrom in ["chr1", "chr5"]:
    fig = plot_coverage_track(
        coverage=depth_data[chrom],
        variants=[v for v in filtered if v["chrom"] == chrom],
        region=(chrom, 0, chrom_sizes[chrom]),
        output_path=f"coverage_{chrom}.png",
        bin_size=1000,
    )

# Size distribution histogram
fig_hist = plot_sv_size_distribution(
    variants=filtered,
    output_path="sv_size_distribution.png",
)
```

### 11. Key Return Types

The detection functions return dataclass objects (`detect_cnv_from_depth` returns a dict of them keyed by chromosome):

```python-snippet
from dataclasses import asdict
from metainformant.structural_variants.detection.breakpoints import BreakpointPair
from metainformant.structural_variants.detection.cnv import CNVResult, CNVSegment
from metainformant.structural_variants.detection.sv_calling import StructuralVariant

# CNVResult
result: CNVResult
result.segments        # list[CNVSegment]
result.log2_ratios     # list[float] | None
result.bin_size        # int
result.method          # str
result.parameters      # dict

# CNVSegment
seg: CNVSegment
seg.chrom, seg.start, seg.end  # genomic coordinates (0-based start, exclusive end)
seg.mean_log2ratio     # float
seg.state              # 'DEL' | 'DUP' | 'NEUTRAL' | 'AMP' | 'HOMODEL'
seg.cn                 # estimated integer copy number
seg.confidence         # [0, 1]

# StructuralVariant
sv: StructuralVariant
sv.chrom, sv.start, sv.end
sv.sv_type             # SVType enum
sv.size, sv.quality, sv.genotype
sv.evidence            # SVEvidence object
asdict(sv)             # plain nested dict; sv_type stays an SVType member (use sv.sv_type.value)

# BreakpointPair
bp_pair: BreakpointPair
bp_pair.bp1, bp_pair.bp2       # Breakpoint objects
bp_pair.bp1.position, bp_pair.bp1.strand, bp_pair.bp1.confidence, bp_pair.bp1.microhomology
```

---

## Configuration

Behaviour is controlled by function arguments. The only environment variables the module reads are listed below; each is read when the function is called and **overrides** the matching argument, and a non-numeric value raises `ValueError`:

| Variable | Read by | Overrides |
|----------|---------|-----------|
| `SV_MIN_SUPPORT` (int) | `call_structural_variants`, `filter_by_quality` | `min_support` |
| `SV_MIN_QUAL` (float) | `filter_by_quality` | `min_qual` |
| `SV_CNV_SIGNIFICANCE` (float) | `detect_cnv_from_depth` | `significance` |
| `SV_CNV_DEL_THRESHOLD`, `SV_CNV_DUP_THRESHOLD` (float) | `call_cnv_states` only (`detect_cnv_from_depth` ignores them) | `del_threshold`, `dup_threshold` |

```bash
export SV_MIN_SUPPORT=3          # Minimum supporting reads
export SV_MIN_QUAL=20.0          # Minimum quality score (filter_by_quality)
export SV_CNV_SIGNIFICANCE=0.01  # CBS p-value threshold
```

Everything else (`min_mapq`, `min_sv_size`, window sizes, filter limits, ...) is a plain function argument:

```python
variants = call_structural_variants(alignments, min_support=5, min_mapq=30, min_sv_size=100)
filtered, _ = filter_by_quality(sv_dicts, min_qual=30.0, min_support=5)
```

There is no YAML or `config.get`/`config.set` configuration for this module. See **[CONFIGURATION.md](CONFIGURATION.md)** for the full reference.

---

## Expected Outputs

The detection, annotation, filtering and population functions return Python objects and write no files; the plot functions save the figure to `output_path`. You can produce:

- **JSON**: Machine-readable results via `core.io.dump_json` (convert dataclasses to dicts first)
- **VCF**: Written with the DNA module's `write_vcf` (there is no SV-specific writer; see Example 14 in [EXAMPLES.md](EXAMPLES.md))
- **PNG/PDF**: Figures (Circos, coverage, histograms); the format follows the output file extension
- **CSV**: Spreadsheet tables for downstream analysis, for example `pd.DataFrame(sv_dicts).to_csv(...)` (no structural variant function writes CSV)

The modules log stage summaries at INFO (for example `Called N structural variants from M evidence items`) and details at DEBUG, but those lines are only shown once logging is configured (warnings are printed regardless). Call `configure_logging_from_env()` at the start of your script; it reads `CORE_LOG_LEVEL` (default `INFO`):
```python
from metainformant.core.utils import logging

logging.configure_logging_from_env()
```
```bash
export CORE_LOG_LEVEL=DEBUG   # per-step details; only takes effect through the call above
```

---

## Common Workflows

### Workflow A: Somatic CNV Calling (Tumor vs. Normal)

```python
import numpy as np
from metainformant.structural_variants.detection.cnv import calculate_log2_ratio, detect_cnv_from_depth

# Load matched normal depth: {chrom: per-bin depth array}, e.g. built as in step 2
normal_depth = compute_depth("normal.bam", window_size=1000)   # your helper
tumor_depth = compute_depth("tumor.bam", window_size=1000)

# Compute log2 ratio per chromosome (gc_bias: {chrom: GC fraction per bin}, optional)
log2_ratios = {
    chrom: calculate_log2_ratio(tumor_depth[chrom], normal_depth[chrom], gc_content=gc_bias[chrom])
    for chrom in tumor_depth
}

# Call CNVs. detect_cnv_from_depth() takes depth-like positive values and computes
# log2(value / median) itself, so pass the linear ratio (2**log2), not the log2 values.
cnv_result = detect_cnv_from_depth(
    depth_data={chrom: np.exp2(lr) for chrom, lr in log2_ratios.items()},
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
import numpy as np

# Reference SV panel (e.g., from 1000G). Each SV record carries the per-sample evidence under
# "samples": {sample_id: {"depth_ratio": ..., "split_reads": ..., "discordant_pairs": ...}}.
# Computing those values from the BAMs is up to you; genotype_sv_population() opens no files.
reference_svs = io.load_json("sv_reference_panel.json")

# Genotype all samples in one call
gt = genotype_sv_population(reference_svs, sample_order, method="depth")

# Genotype matrix: n_svs × n_samples
gt_matrix = np.asarray(gt["genotype_matrix"])

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
