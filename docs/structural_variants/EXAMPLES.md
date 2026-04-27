# Structural Variants: Real-World Examples

Practical, end-to-end analysis scenarios demonstrating the full capabilities of the structural variant module, from basic CNV calling to multi-omics integration and large-cohort population genetics.

---

## Table of Contents

1. [Example 1: Minimal CNV Calling from BAM](#1-minimal-cnv-calling-from-bam)
2. [Example 2: Germline SV Discovery in WGS](#2-germline-sv-discovery-in-wgs-whole-genome-sequencing)
3. [Example 3: Somatic CNV Detection (Tumor vs Normal)](#3-somatic-cnv-detection-tumor-vs-normal)
4. [Example 4: Multi-Caller Consensus Building](#4-multi-caller-consensus-building)
5. [Example 5: Clinical Interpretation — Gene Overlap + Functional Impact](#5-clinical-interpretation--gene-overlap--functional-impact)
6. [Example 6: Population Frequency Comparison with gnomAD](#6-population-frequency-comparison-with-gnomad)
7. [Example 7: Visualizing Complex Rearrangements — Circos + Coverage](#7-visualizing-complex-rearrangements--circos--coverage)
8. [Example 8: Cohort SV Genotyping and Association Testing](#8-cohort-sv-genotyping-and-association-testing)
9. [Example 9: Long-Read Breakpoint Refinement](#9-long-read-breakpoint-refinement)
10. [Example 10: Batch Processing Pipeline for 1000+ Samples](#10-batch-processing-pipeline-for-1000-samples)
11. [Example 11: Integration with DNA Module for Variant Effect Prediction](#11-integration-with-dna-module-for-variant-effect-prediction)
12. [Example 12: Multi-Omics Correlation: SVs and Gene Expression](#12-multi-omics-correlation-svs-and-gene-expression)

---

## Example 1: Minimal CNV Calling from BAM

**Goal:** Detect copy-number variants from a tumor BAM using depth only.

**Input:** `tumor.bam` (sorted, indexed).

**Workflow:**
1. Count read depth in 1 kb windows across each chromosome
2. Compute log2 ratio (self-normalized against median)
3. CBS segmentation
4. Call DEL/DUP/NEUTRAL

**Code:**

```python
import pysam
import numpy as np
from metainformant.structural_variants.detection.cnv import detect_cnv_from_depth
from metainformant.core import io, logging

logging.setup_logging(level="INFO")

bam = pysam.AlignmentFile("tumor.bam", "rb")
window_size = 1000
depth_data = {}

for chrom in bam.references:
    n_bins = max(1, bam.lengths[bam.get_tid(chrom)] // window_size)
    depths = np.zeros(n_bins, dtype=np.float64)
    for read in bam.fetch(chrom):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        bin_idx = read.reference_start // window_size
        if bin_idx < n_bins:
            depths[bin_idx] += 1.0
    depth_data[chrom] = depths
bam.close()

results = detect_cnv_from_depth(
    depth_data=depth_data,
    window_size=window_size,
    method="segmentation",
    significance=0.01,
    ploidy=2,
)

total = sum(len(r.segments) for r in results.values())
print(f"Detected {total} CNV segments")

for chrom, cnv_result in results.items():
    for seg in cnv_result.segments:
        if seg.state != "NEUTRAL":
            print(f"  {seg.chrom}:{seg.start:,}-{seg.end:,}  "
                  f"{seg.state}  CN={seg.cn}  log2={seg.mean_log2ratio:.3f}")

io.dump_json(
    {c: [{"chrom": s.chrom, "start": s.start, "end": s.end,
          "state": s.state, "cn": s.cn, "conf": s.confidence}
         for s in r.segments] for c, r in results.items()},
    "cnv_calls.json"
)
```

**Output:**
```
Detected 47 CNV segments
  chr1:1,000-25,000  DEL   CN=1  log2=-0.612
  chr5:100,000-150,000  DUP  CN=3  log2=0.587
  chr10:45,000-78,000  DEL  CN=1  log2=-0.721
  ...
```

---

## Example 2: Germline SV Discovery in WGS (Whole-Genome Sequencing)

**Goal:** Detect germline SVs (DEL, DUP, INV, TRA) from a single WGS BAM (30×).

**Input:** `sample.bam` (paired-end, ~800M reads).

**Code:**

```python
import pysam
from metainformant.structural_variants.detection.sv_calling import (
    call_structural_variants,
)
from metainformant.structural_variants.detection.breakpoints import refine_breakpoints
from metainformant.core import logging, io

logging.setup_logging(level="INFO")

# ── Extract alignment records (streaming) ───────────────────────────────────
bam = pysam.AlignmentFile("sample.bam", "rb")
alignments = []
read_count = 0

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
    read_count += 1
    if read_count % 100_000 == 0:
        logging.get_logger(__name__).info(f"  Processed {read_count:,} reads...")

bam.close()
logging.info(f"Total aligned reads: {read_count:,}")

# ── Call SVs ─────────────────────────────────────────────────────────────────
variants = call_structural_variants(
    alignments=alignments,
    min_support=3,
    min_mapq=20,
    min_sv_size=50,
)

from collections import Counter
type_counts = Counter(v.sv_type.value for v in variants)
logging.info(f"SV type breakdown: {dict(type_counts)}")

# Print high-confidence calls (quality > 50)
high_conf = [v for v in variants if v.quality >= 50]
logging.info(f"High-confidence (qual≥50): {len(high_conf)}/{len(variants)}")

# ── Save to VCF ─────────────────────────────────────────────────────────────
from metainformant.structural_variants.io import write_vcf
write_vcf(variants, "sv_calls.vcf")
logging.info("Saved VCF: sv_calls.vcf")

# ── Refine breakpoints (optional but recommended) ───────────────────────────
# Re-fetch reads near variant breakpoints
bam = pysam.AlignmentFile("sample.bam", "rb")
local_reads = []
for sv in variants:
    chroms = {sv.chrom, getattr(sv, "chrom2", sv.chrom)}
    for chrom in chroms:
        window = 300
        for read in bam.fetch(chrom, max(0, sv.start - window), sv.end + window):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
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

refined = refine_breakpoints(
    variants=[v.__dict__ for v in variants],
    reads=local_reads,
    window=200,
    min_clip=5,
)

# Save refined coordinates
io.dump_json([{
    "sv_type": bp.sv_type,
    "chrom1": bp.bp1.chrom, "pos1": bp.bp1.position,
    "chrom2": bp.bp2.chrom, "pos2": bp.bp2.position,
    "confidence": bp.confidence,
    "microhomology": bp.bp1.microhomology or bp.bp2.microhomology,
} for bp in refined], "refined_breakpoints.json")
```

---

## Example 3: Somatic CNV Detection (Tumor vs Normal)

**Goal:** Identify somatic CNVs (present in tumor, absent in normal) from matched normal.

**Input:** `tumor.bam`, `normal.bam`.

**Approach:** Compute log2(tumor/normal) depth ratios with GC correction.

```python
import pysam
import numpy as np
from metainformant.structural_variants.detection.cnv import (
    calculate_log2_ratio,
    detect_cnv_from_depth,
)

def compute_depth_GC(bam_path, window_size=1000):
    """Compute depth and GC content per window."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    depth_array = {}
    gc_array = {}

    for chrom in bam.references:
        chrom_len = bam.lengths[bam.get_tid(chrom)]
        n_bins = chrom_len // window_size
        depths = np.zeros(n_bins, dtype=np.float64)
        gc_vals = np.zeros(n_bins, dtype=np.float64)

        for read in bam.fetch(chrom):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            bin_idx = read.reference_start // window_size
            if bin_idx < n_bins:
                depths[bin_idx] += 1.0
                # Compute GC from read sequence if available
                if read.query_sequence:
                    start_in_window = read.reference_start % window_size
                    window_seq_offset = (bin_idx * window_size) + start_in_window
                    # Simplified: you would extract the portion overlapping window
                    seq_window = read.query_sequence[:100]
                    gc = (seq_window.count('G') + seq_window.count('C')) / len(seq_window) if seq_window else 0.5
                    gc_vals[bin_idx] = gc

        depth_array[chrom] = depths
        gc_array[chrom] = gc_vals

    bam.close()
    return depth_array, gc_array

# ── Compute depth for tumor and normal ──────────────────────────────────────
tumor_depth, tumor_gc = compute_depth_GC("tumor.bam")
normal_depth, normal_gc = compute_depth_GC("normal.bam")

# Keep only chromosomes present in both
common_chroms = set(tumor_depth) & set(normal_depth)

# Compute log2 ratio with GC correction
log2_ratios = {}
for chrom in common_chroms:
    td = tumor_depth[chrom]
    nd = normal_depth[chrom]
    gc = tumor_gc[chrom]  # Use tumor GC; could average
    # Truncate to same length
    n = min(len(td), len(nd), len(gc))
    td = td[:n]
    nd = nd[:n]
    gc = gc[:n]
    lr = calculate_log2_ratio(td, nd, gc_content=gc, pseudocount=1.0)
    log2_ratios[chrom] = lr

# ── Call CNVs ────────────────────────────────────────────────────────────────
cnv_results = detect_cnv_from_depth(
    depth_data=log2_ratios,
    window_size=1000,
    method="segmentation",
    significance=0.001,  # More stringent for somatic
    ploidy=2,
    min_segment_bins=5,
)

# Filter: keep only large enough events, moderate+ confidence
somatic_cnvs = []
for chrom, cnv_result in cnv_results.items():
    for seg in cnv_result.segments:
        size = (seg.end - seg.start)
        # Consider deletions and duplications with high confidence
        if seg.state in ("DEL", "DUP") and seg.confidence >= 0.7 and size >= 10000:
            somatic_cnvs.append(seg)
            print(f"  {seg.chrom}:{seg.start:,}-{seg.end:,}  {seg.state}  "
                  f"CN={seg.cn}  size={size/1000:.1f}kb  conf={seg.confidence:.2f}")

print(f"\nSomatic CNV calls: {len(somatic_cnvs)}")
```

**Interpretation tips:**
- Somatic DEL usually means tumor loss (CN=0 or 1). Somatic DUP means gain (CN=3+).
- Germline CNVs (present in both tumor and normal) cancel out in log2 ratio; they should be NEUTRAL in somatic difference.
- Use a matched normal whenever possible to control for germline CNVs and GC bias.

---

## Example 4: Multi-Caller Consensus Building

**Goal:** Combine SV calls from Manta, Delly, Lumpy into a high-confidence consensus.

**Input:** Three VCF files (or JSON after parsing).

```python
import pysam
from metainformant.structural_variants.filtering.merge import merge_callsets
from metainformant.core import io

# Read each VCF with pysam and extract simple dicts
def read_vcf_simple(vcf_path):
    """Convert VCF records to simplified dict format."""
    vcf = pysam.VariantFile(vcf_path, "r")
    variants = []
    for rec in vcf:
        variants.append({
            "chrom": rec.chrom,
            "start": rec.pos,              # 1-based in VCF
            "end": rec.info.get("END", rec.pos + len(rec.alts[0])),
            "sv_type": rec.info.get("SVTYPE", "UNKNOWN"),
            "size": rec.info.get("SVLEN", [0])[0] if "SVLEN" in rec.info else 0,
            "quality": rec.qual or 0.0,
            "genotype": rec.samples[0]["GT"] if rec.samples else "./.",
        })
    vcf.close()
    return variants

manta_calls = read_vcf_simple("manta.vcf")
delly_calls = read_vcf_simple("delly.vcf")
lumpy_calls = read_vcf_simple("lumpy.vcf")

callsets = {"Manta": manta_calls, "Delly": delly_calls, "Lumpy": lumpy_calls}

# ── Merge ────────────────────────────────────────────────────────────────────
merged, stats = merge_callsets(
    callsets=callsets,
    min_overlap=0.5,   # 50% reciprocal overlap
    type_match=True,   # Must be same SV type
)

print(f"Input: {stats.n_input_variants} total calls")
print(f"  Manta: {stats.caller_counts['Manta']}")
print(f"  Delly: {stats.caller_counts['Delly']}")
print(f"  Lumpy: {stats.caller_counts['Lumpy']}")
print(f"\nOutput: {stats.n_output_variants} consensus variants")
print(f"  Multi-caller (≥2): {stats.n_multi_caller}")
print(f"  Single-caller: {stats.n_single_caller}")

# Consensus requires 2+ callers for high confidence
consensus = [mv for mv in merged if mv.n_callers >= 2]
print(f"\nHigh-confidence consensus (≥2 callers): {len(consensus)}")

for mv in consensus[:10]:
    print(f"  [{mv.sv_type}] {mv.chrom}:{mv.start:,}-{mv.end:,}  "
          f"callers={mv.n_callers}  qual={mv.consensus_quality:.1f}")

# Save consensus to VCF
from metainformant.structural_variants.io import write_vcf_from_merged
write_vcf_from_merged(consensus, "consensus.vcf")
```

**Survivor-style merging** (distance-based) for callers without reciprocal overlap:

```python
from metainformant.structural_variants.filtering.merge import survivor_merge

merged = survivor_merge(
    vcf_files=["manta.vcf", "delly.vcf", "lumpy.vcf"],
    max_distance=1000,       # Max bp between breakpoints
    min_callers=2,           # Must be in at least 2 callers
    type_match=True,
    strand_match=False,
)
print(f"Survivor merge: {len(merged)} consensus calls")
```

---

## Example 5: Clinical Interpretation — Gene Overlap + Functional Impact

**Goal:** Prioritize SVs for clinical review based on gene overlap, dosage sensitivity, and TAD disruption.

**Input:** Filtered SV calls, gene annotations, ClinGen dosage database.

```python
from metainformant.structural_variants.annotation.overlap import annotate_gene_overlap
from metainformant.structural_variants.annotation.functional_impact import (
    predict_functional_impact,
)
from metainformant.core import io

# Load data
svs = io.load_json("filtered_svs.json")
gene_db = io.load_bed("genes.bed")  # Must have 'name' column
hi_db = io.load_json("clinGen_haploinsufficiency.json")   # {gene: score (0-1)}
ts_db = io.load_json("clinGen_triplosensitivity.json")
pli_db = io.load_json("gnomad_pli.json")   # {gene: pLI}
tad_bounds = io.load_json("tad_boundaries.json")

# ── Step 1: Gene overlap ────────────────────────────────────────────────────
annotated = annotate_gene_overlap(svs, gene_db)
print(f"Annotated {len(annotated)} SVs; "
      f"{sum(1 for v in annotated if v['n_genes_affected']>0)} overlap genes")

# ── Step 2: Functional impact ───────────────────────────────────────────────
priority_svs = []
for variant in annotated:
    impact = predict_functional_impact(
        variant=variant,
        gene_annotations=gene_db,
        haploinsufficiency_db=hi_db,
        triplosensitivity_db=ts_db,
        pli_db=pli_db,
        tad_boundaries=tad_bounds,
    )

    # Clinical prioritization criteria
    is_high_impact = (
        impact.impact_level == "HIGH" and
        impact.pathogenicity_score >= 0.7 and
        (impact.tad_disrupted or len(impact.dosage_sensitive_genes) > 0)
    )

    if is_high_impact:
        priority_svs.append({
            "variant": variant,
            "impact": impact,
        })
        print(f"\nHIGH-PRIORITY SV:")
        print(f"  {impact.variant_id}: {impact.impact_type}")
        print(f"  Pathogenicity: {impact.pathogenicity_score:.3f}")
        print(f"  Genes: {', '.join(impact.affected_genes)}")
        if impact.dosage_sensitive_genes:
            print(f"  Dosage-sensitive: {', '.join(impact.dosage_sensitive_genes)}")
        if impact.tad_disrupted:
            n_broken = impact.details.get("tad_prediction", {}).get("n_boundaries", 0)
            print(f"  ⚠ TAD boundary disruption: {n_broken} boundaries")

print(f"\n{len(priority_svs)} variants meet clinical review criteria")

# Save prioritized list
io.dump_json(priority_svs, "clinical_priority.json")
```

**Sample output:**
```
HIGH-PRIORITY SV:
  chr5:150,234,100-150,456,789: HOMODEL — dosage_loss
  Pathogenicity: 0.89
  Genes: CTNNA2
  Dosage-sensitive: CTNNA2 (HI=0.92, pLI=0.97)
  ⚠ TAD boundary disruption: 1 boundaries
```

---

## Example 6: Population Frequency Comparison with gnomAD

**Goal:** Filter SVs that are rare (<1% AF) in gnomAD-SV population database.

**Input:** SV calls, gnomAD-SV frequency database (`gnomad_sv.json`).

```python
from metainformant.structural_variants.filtering.quality_filter import filter_by_frequency
from metainformant.core import io

svs = io.load_json("sv_calls.json")
gnomad = io.load_json("gnomad_sv.json")  # Structure: {chrom: [{pos, end, sv_type, af}]}

# ── Filter ───────────────────────────────────────────────────────────────────
filtered, stats = filter_by_frequency(
    variants=svs,
    population_db=gnomad,
    max_af=0.01,                  # Keep AF < 1%
    match_window=200,             # Position match within 200 bp
    min_reciprocal_overlap=0.5,   # 50% reciprocal overlap required
)

print(f"Filtered: removed {stats.filtered_count} common variants")
print(f"Remaining: {stats.output_count} rare variants")

# Show removed
if stats.filtered_count > 0:
    print("\nCommon variants removed (AF > 1%):")
    for v in stats.filtered_variants[:5]:
        match = v.get("pop_freq_match", {})
        print(f"  {v['chrom']}:{v['start']} matched gnomAD AF={match.get('af', '?'):.4f}")

# Save rare variants
io.dump_json(filtered, "rare_svs.json")
```

**gnomAD-SV format** (auto-converted internally):
```json
{
  "chr1": [
    {"pos": 123456, "end": 123789, "sv_type": "DEL", "af": 0.003},
    {"pos": 987654, "end": 988123, "sv_type": "DUP", "af": 0.012}
  ]
}
```

**Per-population frequency filtering:**

```python
# Filter using population-specific thresholds
filtered, _ = filter_by_frequency(
    variants=svs,
    population_db=gnomad,
    max_af=0.001,  # More stringent
    pop_filter="EUR",  # Only check European AF
)
```

---

## Example 7: Visualizing Complex Rearrangements — Circos + Coverage

**Goal:** Produce publication-quality figures for a sample with multiple rearrangements.

**Input:** SV calls + CNV segments.

```python
from metainformant.structural_variants.visualization.plots import (
    plot_circos,
    plot_coverage_track,
    plot_size_distribution,
)
from metainformant.structural_variants.detection.cnv import CNVResult

# Chromosome sizes for hg38
chrom_sizes = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
    "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
    "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345, "chr17": 83257441, "chr18": 80373285,
    "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
    "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
}

# ── Circos plot ──────────────────────────────────────────────────────────────
fig = plot_circos(
    variants=filtered_svs,
    chromosomes=chrom_sizes,
    output_path="circos.png",
    title="Sample123 — Structural Variants (n=152)",
    figsize=(14, 14),
    show_labels=True,
)
print("Circos saved: circos.png")

# ── Coverage tracks for rearranged chromosomes ──────────────────────────────
from metainformant.structural_variants.detection.cnv import CNVResult

cnv_result: dict[str, CNVResult] = ...  # load from earlier

for chrom in ["chr5", "chr8", "chr17"]:  # Chromosomes with interesting SVs
    if chrom not in cnv_result:
        continue
    cov = cnv_result[chrom].log2_ratios
    bin_size = cnv_result[chrom].bin_size
    svs_on_chrom = [v for v in filtered_svs if v["chrom"] == chrom]

    fig = plot_coverage_track(
        coverage=cov,
        variants=svs_on_chrom,
        region=(chrom, 0, chrom_sizes[chrom]),
        output_path=f"coverage_{chrom}.png",
        title=f"Copy-Number Profile — {chrom}",
        bin_size=bin_size,
    )
    print(f"Coverage track: coverage_{chrom}.png")

# ── Size distribution ───────────────────────────────────────────────────────
fig = plot_size_distribution(
    variants=filtered_svs,
    output_path="size_distribution.png",
    log_scale=True,
)
```

**Result:** `circos.png` — circular genome with colored chords connecting breakpoints; `coverage_chr5.png` — CNV profile with SV rectangles overlay.

---

## Example 8: Cohort SV Genotyping and Association Testing

**Goal:** Genotype 100 known SVs across 500 samples and test for phenotype association.

**Input:**
- `reference_svs.json`: 100 variant definitions
- `samples.tsv`: sample IDs and BAM paths
- `phenotypes.tsv`: sample IDs + phenotype (binary case/control)

```python
import numpy as np
from metainformant.structural_variants.population.sv_population import (
    genotype_sv_population,
    sv_allele_frequency,
    sv_association_test,
    sv_population_structure,
)

# ── Load cohort ──────────────────────────────────────────────────────────────
import pandas as pd
samples_df = pd.read_csv("samples.tsv", sep="\t")   # columns: sample_id, bam_path, population
pheno_df = pd.read_csv("phenotypes.tsv", sep="\t")   # columns: sample_id, phenotype (0/1)

sample_ids = samples_df["sample_id"].tolist()
bam_paths = dict(zip(samples_df["sample_id"], samples_df["bam_path"]))
populations = samples_df["population"].tolist()

# Reference SV panel (100 SVs)
ref_svs = io.load_json("reference_svs.json")

# ── Genotype across all samples ──────────────────────────────────────────────
genotype_result = genotype_sv_population(
    sv_calls=ref_svs,
    samples=sample_ids,
    method="depth",  # Faster; "split" more precise but slower
)

gt_matrix = np.array(genotype_result["genotype_matrix"])  # shape: (100, 500)
print(f"Genotype matrix: {gt_matrix.shape[0]} SVs × {gt_matrix.shape[1]} samples")

# ── Allele frequencies ───────────────────────────────────────────────────────
freq = sv_allele_frequency(
    genotype_matrix=gt_matrix,
    sample_labels=populations,
)

print("Overall allele frequency stats:")
print(f"  Mean AF: {np.mean(freq['frequencies']):.4f}")
print(f"  Polysorphic SVs: {freq['n_polymorphic']} / {len(ref_svs)}")

for pop in set(populations):
    pop_idx = [i for i, p in enumerate(populations) if p == pop]
    pop_af = freq["population_frequencies"][pop]
    print(f"  {pop} (n={len(pop_idx)}): mean AF={np.mean(pop_af):.4f}")

# ── Population structure (PCA) ───────────────────────────────────────────────
# Adjust for ancestry in association tests
pca = sv_population_structure(gt_matrix, n_components=10)
print(f"\nPCA — PC1 variance explained: {pca['explained_variance'][0]:.2%}")
print(f"         PC2 variance explained: {pca['explained_variance'][1]:.2%}")

# ── Association testing ──────────────────────────────────────────────────────
phenotypes = pheno_df["phenotype"].values  # 0 or 1

results = []
for i, sv in enumerate(ref_svs):
    rec = sv_association_test(
        genotypes=gt_matrix[i],
        phenotypes=phenotypes,
        covariates=pca["components"][:, :5],  # Top 5 PCs as covariates
    )
    results.append({
        "sv_id": i,
        "chrom": sv["chrom"],
        "start": sv["start"],
        "p_value": rec["p_value"],
        "beta": rec["beta"],
        "odds_ratio": rec["odds_ratio"],
    })

# Multiple testing correction
from metainformant.metabolomics.pathways.enrichment import benjamini_hochberg
pvals = [r["p_value"] for r in results]
qvals = benjamini_hochberg(pvals)

for i, r in enumerate(results):
    r["q_value"] = qvals[i]

sig = [r for r in results if r["q_value"] < 0.05]
print(f"\nAssociation hits (FDR < 0.05): {len(sig)}")
for hit in sig:
    print(f"  SV {hit['sv_id']} @ {hit['chrom']}:{hit['start']:,}")
    print(f"    OR={hit['odds_ratio']:.3f}, p={hit['p_value']:.2e}, q={hit['q_value']:.4g}")

io.dump_json(results, "association_results.json")
```

---

## Example 9: Long-Read Breakpoint Refinement (PacBio / ONT)

**Goal:** Use long-read sequences to refine SV breakpoints to single-nucleotide precision and detect insertion sequences.

Long reads often span the entire SV; parse their alignments to extract the precise breakpoint and inserted sequence.

```python
from metainformant.structural_variants.detection.breakpoints import refine_breakpoints
from metainformant.longread import io as lr_io  # placeholder module

# Load long-read alignments (BAM from minimap2 / pbmm2)
longreads = lr_io.read_bam("longreads.bam")

# Extract reads near pre-called SV breakpoints (from short-read)
sv_calls = io.load_json("shortread_svs.json")

# Refine using long reads
refined = refine_breakpoints(
    variants=sv_calls,
    reads=longreads,
    window=1000,        # Long reads can span larger windows
    min_clip=10,
)

for bp_pair in refined:
    if bp_pair.bp1.inserted_sequence:
        print(f"Insertion at {bp_pair.bp1.chrom}:{bp_pair.bp1.position}: "
              f"sequence='{bp_pair.bp1.inserted_sequence}'")
    if bp_pair.bp1.microhomology:
        print(f"Microhomology: {bp_pair.bp1.microhomology}")

# Output precise breakpoints for VCF conversion
io.dump_json([{
    "chrom": bp.bp1.chrom,
    "pos": bp.bp1.position,
    "chrom2": bp.bp2.chrom,
    "pos2": bp.bp2.position,
    "sv_type": bp.sv_type,
    "inserted_seq": bp.bp1.inserted_sequence,
    "confidence": bp.confidence,
} for bp in refined], "precise_breakpoints.json")
```

---

## Example 10: Batch Processing Pipeline for 1000+ Samples

**Goal:** Scalable processing of a large cohort with cluster-friendly chunking.

**Strategy:** Chromosome-wise parallelization + per-sample parallelization.

```python
from metainformant import parallel, config
from metainformant.structural_variants import detection, annotation, filtering
from pathlib import Path
import pandas as pd

config.set("structural_variants.parallel", True)
config.set("structural_variants.max_workers", 8)

# List of all sample BAMs
bam_dir = Path("/data/cohort/bams/")
bam_files = sorted(bam_dir.glob("*.bam"))
sample_ids = [b.stem for b in bam_files]

def process_sample(bam_path):
    """Full SV pipeline for one sample."""
    sample_id = bam_path.stem

    # CNV detection
    from metainformant.structural_variants.detection.cnv import detect_cnv_from_depth
    depth_data = compute_depth(bam_path)  # Your helper
    cnv = detect_cnv_from_depth(depth_data, window_size=1000)

    # SV calling
    from metainformant.structural_variants.detection.sv_calling import call_structural_variants
    alns = extract_alignments(bam_path)  # Your helper
    svs = call_structural_variants(alns, min_support=3)

    # Merge CNV segments + SVs into single callset (optional)
    all_variants = svs + cnv_to_svs(cnv)  # Convert CNVSegments to StructuralVariant-like

    # Filter
    from metainformant.structural_variants.filtering.quality_filter import filter_by_quality
    filtered, _ = filter_by_quality(all_variants, min_qual=20, min_support=3)

    # Save per-sample results
    out = Path(f"results/per_sample/{sample_id}.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    io.dump_json([v.__dict__ for v in filtered], out)

    return sample_id, len(filtered)

# ── Parallel map (per-sample) ─────────────────────────────────────────────────
results = parallel.map(
    process_sample,
    bam_files,
    workers=12,
    chunk_size=50,  # Batches of 50 samples per worker
)

print("Per-sample SV counts:")
for sample_id, count in results:
    print(f"  {sample_id}: {count} SVs")

# ── Cross-sample merging (optional) ─────────────────────────────────────────
# If you want a cohort-level consensus SV set:
from metainformant.structural_variants.filtering.merge import merge_callsets

all_samples = {}
for bam_path in bam_files:
    sid = bam_path.stem
    all_samples[sid] = io.load_json(f"results/per_sample/{sid}.json")

merged_cohort, stats = merge_callsets(
    callsets=all_samples,
    min_overlap=0.5,
)
print(f"Cohort merged SVs: {stats.n_output_variants}")

# ── Logging summary ──────────────────────────────────────────────────────────
import logging
summary = [
    f"Processed {len(bam_files)} samples",
    f"Total SVs detected: {sum(count for _, count in results)}",
    f"Average per sample: {np.mean([c for _, c in results]):.1f}",
]
logging.info("\n".join(summary))
```

---

## Example 11: Integration with DNA Module for Variant Effect Prediction

**Goal:** Use `metainformant.dna` gene and transcript models to compute precise variant effects (exon deletion, fusion).

```python
from metainformant.dna import variants as dna_variants
from metainformant.structural_variants.annotation.overlap import annotate_gene_overlap
from metainformant.structural_variants.annotation.functional_impact import predict_functional_impact

# Load GENCODE GTF via dna module
gene_db = dna_variants.load_gtf("gencode.v44.annotation.gtf")
# Returns dict of transcript objects with exon structures

# Or load pre-processed gene annotations
gene_intervals = dna_variants.get_gene_intervals(gene_db)  # List[GenomicInterval]

# Annotate SV with gene overlap using the DNA module's gene model (includes exons)
svs = io.load_json("sv_calls.json")

annotated = annotate_gene_overlap(svs, gene_intervals)

for v in annotated:
    if v["n_genes_affected"] > 0:
        for overlap in v["gene_overlaps"]:
            gene = overlap.feature.name
            relationship = overlap.relationship  # "contains", "within", "overlap"
            frac = overlap.overlap_fraction_feature  # Fraction of gene overlapped

            # Does this delete an entire exon?
            if relationship == "contains" and frac > 0.5:
                print(f"  SV likely disrupts gene {gene} (covers {frac:.1%})")

# Functional impact now can use the enriched gene annotation
```

---

## Example 12: Multi-Omics Correlation: SVs and Gene Expression

**Goal:** Test whether a CNV correlates with expression of overlapped genes (eQTL-like).

```python
import numpy as np
from scipy.stats import pearsonr
from metainformant.structural_variants.detection.cnv import CNVResult
from metainformant.structural_variants.annotation.overlap import annotate_gene_overlap

# Load matched data
cnv_results: dict[str, CNVResult] = io.load_pickle("cnv_results.pkl")  # Per-sample
gene_expression = io.load_csv("gene_expression.tsv")  # genes × samples
gene_db = io.load_bed("genes.bed")

# Get per-sample CN state for each gene
def gene_cn_state(cnv_results, gene):
    """Return CN state for a gene across samples."""
    chrom = gene.chrom
    gene_mid = (gene.start + gene.end) // 2
    states = []
    for sample_cnv in cnv_results.values():
        if chrom not in sample_cnv:
            states.append(2)  # Neutral if not found
            continue
        for seg in sample_cnv[chrom].segments:
            if seg.start <= gene_mid < seg.end:
                states.append(seg.cn)
                break
        else:
            states.append(2)
    return np.array(states)

# Genes overlapped by any CNV
all_gene_names = set()
for cnv in cnv_results.values():
    # Build simplified gene overlaps...
    pass  # In practice overlap all genes

# Correlation for one gene
sample_ids = list(cnv_results.keys())
gene = "BRCA1"
cn = gene_cn_state(cnv_results, gene_db[gene])
expr = gene_expression.loc[gene, sample_ids].values

r, p = pearsonr(cn, expr)
print(f"{gene}: r={r:.3f}, p={p:.2e}")

# Loop over all genes affected by SVs
correlations = []
for gene_name in all_gene_names:
    cn = gene_cn_state(cnv_results, gene_db[gene_name])
    expr = gene_expression.loc[gene_name].values
    r, p = pearsonr(cn, expr)
    correlations.append({"gene": gene_name, "r": r, "p": p})

# Multiple testing correction
from metainformant.metabolomics.pathways.enrichment import benjamini_hochberg
pvals = [c["p"] for c in correlations]
qvals = benjamini_hochberg(pvals)
for i, c in enumerate(correlations):
    c["q"] = qvals[i]

sig_cis = [c for c in correlations if c["q"] < 0.05]
print(f"CNV-eQTLs: {len(sig_cis)} significant genes")
for hit in sig_cis:
    print(f"  {hit['gene']}: r={hit['r']:.3f}, q={hit['q']:.4g}")
```

---

## Example 13: Integration with RNA-Seq Module

Correlate SV breakpoints with splicing abnormalities:

```python
from metainformant.rna import differential_splicing
from metainformant.structural_variants.detection.breakpoints import Breakpoint

# Find inversions or genes with breakpoints within introns
svs = io.load_json("filtered_svs.json")
intragenic = [v for v in svs if v["n_genes_affected"] > 0 and v.get("is_intronic")]

# For each intragenic SV, extract breakpoint position and gene
for sv in intragenic:
    gene = sv["overlapping_genes"][0]
    bp_pos = sv["breakpoint1"]  # or refine for precision

    # Query RNA-seq splicing around breakpoint
    splicing_results = differential_splicing(
        gene=gene,
        condition="tumor_vs_normal",
        region=(bp_pos - 5000, bp_pos + 5000),
    )
    # Look for novel splice junctions involving the breakpoint
    for junction in splicing_results["abnormal_junctions"]:
        if abs(junction.start - bp_pos) < 50 or abs(junction.end - bp_pos) < 50:
            print(f"SV breakpoint in {gene} co-locates with abnormal splice junction")
```

---

## Example 14: Export to Standard Formats

**VCF Export:**

```python
from metainformant.structural_variants.io import write_vcf

# From list of StructuralVariant objects
variants = io.load_json("sv_calls.json")  # or object list
write_vcf(variants, "output.vcf")
```

**BED Export:**

```python
def svs_to_bed(svs, path):
    with open(path, "w") as f:
        for sv in svs:
            name = f"{sv['sv_type']}_{sv['chrom']}:{sv['start']}"
            f.write(f"{sv['chrom']}\t{sv['start']}\t{sv['end']}\t{name}\t0\t+\n")

svs_to_bed(filtered, "sv_calls.bed")
```

**TSV/CSV Export:**

```python
import pandas as pd
rows = []
for sv in variants:
    rows.append({
        "chrom": sv.chrom,
        "start": sv.start,
        "end": sv.end,
        "type": sv.sv_type,
        "size": sv.size,
        "qual": sv.quality,
        "gt": sv.genotype,
    })
df = pd.DataFrame(rows)
df.to_csv("sv_table.csv", index=False)
```

---

## Example 15: Custom SV Type: Mobile Element Insertions (MEIs)

**Goal:** Detect Alu/L1/SVA insertions using split-read + 3' flanking evidence.

```python
# Extend sv_calling.classify_sv_type to identify MEIs
def classify_mei(evidence, split_reads, flanking_sequence=None):
    """Heuristic: many split reads with poly-A tail or target site duplication."""
    if evidence.split_reads >= 5:
        # Check for poly-A tail in clipped sequence
        if flanking_sequence and flanking_sequence[-20:].count('A') > 15:
            return "MEI_ALU"  # or LINE1, SVA based on length/signature
    return None

# Filter to MEIs
mei_calls = [sv for sv in variants if sv.sv_type == "INS" and sv.size > 200]
```

---

## What's Next?

- **[ARCHITECTURE.md](ARCHITECTURE.md)** — Understand how these examples fit into the system design
- **[PERFORMANCE.md](PERFORMANCE.md)** — Optimize for your dataset size
- **[TROUBLESHOOTING.md](TROUBLESHOOTING.md)** — Diagnose issues when things go wrong
- **[detection.md](detection.md)**, **[annotation.md](annotation.md)**, **[filtering.md](filtering.md)**, **[population.md](population.md)** — Deep-dive into each submodule
