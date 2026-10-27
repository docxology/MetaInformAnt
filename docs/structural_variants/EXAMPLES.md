# Structural Variants: Real-World Examples

Practical, end-to-end analysis scenarios demonstrating the main capabilities of the structural variant module, from basic CNV calling to multi-omics integration and large-cohort population genetics.

---

## Table of Contents

1. [Example 1: Minimal CNV Calling from BAM](#example-1-minimal-cnv-calling-from-bam)
2. [Example 2: Germline SV Discovery in WGS](#example-2-germline-sv-discovery-in-wgs-whole-genome-sequencing)
3. [Example 3: Somatic CNV Detection (Tumor vs Normal)](#example-3-somatic-cnv-detection-tumor-vs-normal)
4. [Example 4: Multi-Caller Consensus Building](#example-4-multi-caller-consensus-building)
5. [Example 5: Clinical Interpretation — Gene Overlap + Functional Impact](#example-5-clinical-interpretation--gene-overlap--functional-impact)
6. [Example 6: Population Frequency Comparison with gnomAD](#example-6-population-frequency-comparison-with-gnomad)
7. [Example 7: Visualizing Complex Rearrangements — Circos + Coverage](#example-7-visualizing-complex-rearrangements--circos--coverage)
8. [Example 8: Cohort SV Genotyping and Association Testing](#example-8-cohort-sv-genotyping-and-association-testing)
9. [Example 9: Long-Read Breakpoint Refinement](#example-9-long-read-breakpoint-refinement-pacbio--ont)
10. [Example 10: Batch Processing Pipeline for 1000+ Samples](#example-10-batch-processing-pipeline-for-1000-samples)
11. [Example 11: Integration with DNA Module for Variant Effect Prediction](#example-11-integration-with-dna-module-for-variant-effect-prediction)
12. [Example 12: Multi-Omics Correlation: SVs and Gene Expression](#example-12-multi-omics-correlation-svs-and-gene-expression)
13. [Example 13: Integration with RNA-Seq Module](#example-13-integration-with-rna-seq-module)
14. [Example 14: Export to Standard Formats](#example-14-export-to-standard-formats)
15. [Example 15: Custom SV Type: Mobile Element Insertions (MEIs)](#example-15-custom-sv-type-mobile-element-insertions-meis)

Conventions used throughout: the annotation, filtering and plotting functions take SV **dicts** (`chrom`, `start`, `end`, `sv_type`, ...), not `StructuralVariant` objects, so convert calls with `dataclasses.asdict(sv)` before passing them on (see Example 10). Coordinates are 0-based. There is no YAML or `config.get`/`config.set` configuration: behaviour is controlled by function arguments, plus the few `SV_*` environment variables listed in [CONFIGURATION.md](CONFIGURATION.md) (an environment variable overrides the matching argument).

---

## Example 1: Minimal CNV Calling from BAM

**Goal:** Detect copy-number variants from a tumor BAM using depth only.

**Input:** `tumor.bam` (sorted, indexed).

**Workflow:**
1. Count read depth in 1 kb windows across each chromosome
2. Compute log2 ratio (self-normalized against median)
3. CBS segmentation
4. Call DEL/DUP/NEUTRAL (plus HOMODEL and AMP at extreme ratios)

**Code:**

```python-snippet
import pysam
import numpy as np
from metainformant.structural_variants.detection.cnv import detect_cnv_from_depth
from metainformant.core import io
from metainformant.core.utils import logging

logging.configure_logging_from_env()  # show the modules' INFO logs (CORE_LOG_LEVEL=DEBUG for more)

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

**Output** (from a synthetic depth profile with one 30 kb loss and one 30 kb gain on chr1; the count includes NEUTRAL segments, which are not printed):
```
Detected 6 CNV segments
  chr1:100,000-130,000  DEL  CN=1  log2=-1.003
  chr1:200,000-230,000  DUP  CN=3  log2=0.592
```

`detect_cnv_from_depth` returns a `dict[str, CNVResult]` keyed by chromosome. Its input must be depths, finite and non-negative (zero is a bin with no coverage); a negative or non-finite value raises `ValueError` naming the chromosome. The CBS step is vectorised but still grows with the square of the number of bins per chromosome (see [PERFORMANCE.md](PERFORMANCE.md)): a chromosome of 64,000 bins takes a few seconds and one of 250,000 bins about a minute, so for whole human chromosomes use a coarser `window_size` (tens to hundreds of kb) rather than 1 kb.

---

## Example 2: Germline SV Discovery in WGS (Whole-Genome Sequencing)

**Goal:** Detect germline SVs (DEL, DUP, INV, TRA) from a single WGS BAM (30×).

**Input:** `sample.bam` (paired-end).

Alignment records are held in memory as Python dicts (roughly 0.6-0.7 KB each), so for a whole 30× genome (~800M reads) run the extraction and calling per chromosome or region (`bam.fetch(chrom)`) rather than over the full BAM.

**Code:**

```python-snippet
import pysam
from metainformant.structural_variants.detection.sv_calling import (
    call_structural_variants,
)
from metainformant.structural_variants.detection.breakpoints import refine_breakpoints
from metainformant.core import io
from metainformant.core.utils import logging

logging.configure_logging_from_env()  # show the modules' INFO logs (CORE_LOG_LEVEL=DEBUG for more)
logger = logging.get_logger("sv_example")

# ── Extract alignment records ───────────────────────────────────────────────
# call_structural_variants() reads only these keys: name, chrom, pos, mapq, cigar, read_length,
# is_reverse, mate_chrom, mate_pos, mate_is_reverse, insert_size, sa_tag (no read sequence needed).
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
    alignments.append(rec)
    read_count += 1
    if read_count % 100_000 == 0:
        logger.info(f"  Processed {read_count:,} reads...")

bam.close()
logger.info(f"Total aligned reads: {read_count:,}")

# ── Call SVs ─────────────────────────────────────────────────────────────────
variants = call_structural_variants(
    alignments=alignments,
    min_support=3,
    min_mapq=20,
    min_sv_size=50,
)

from collections import Counter
type_counts = Counter(v.sv_type.value for v in variants)
logger.info(f"SV type breakdown: {dict(type_counts)}")

# Count high-confidence calls (quality >= 50)
high_conf = [v for v in variants if v.quality >= 50]
logger.info(f"High-confidence (qual≥50): {len(high_conf)}/{len(variants)}")

# ── Save to VCF ─────────────────────────────────────────────────────────────
# svs_to_vcf_data() builds the VCF dict for the DNA module's writer (defined in Example 14)
from metainformant.dna.variation.variants import write_vcf
write_vcf(svs_to_vcf_data(variants), "sv_calls.vcf")
logger.info("Saved VCF: sv_calls.vcf")

# ── Refine breakpoints (optional but recommended) ───────────────────────────
# A call's first breakpoint is (chrom, start); its second is (chrom2, pos2) for a translocation
# and (chrom, end) for every other call. Re-fetch the reads around each one.
bam = pysam.AlignmentFile("sample.bam", "rb")
window = 300
local_reads = []
for sv in variants:
    second = (sv.chrom2, sv.pos2) if sv.pos2 is not None else (sv.chrom, sv.end)
    for chrom, pos in ((sv.chrom, sv.start), second):
        for read in bam.fetch(chrom, max(0, pos - window), pos + window):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            local_reads.append({
                "chrom": chrom,
                "pos": read.reference_start,
                "cigar": read.cigarstring,
                "seq": read.query_sequence or "",
                "name": read.query_name,
                "is_reverse": read.is_reverse,
                "sa_tag": read.get_tag("SA") if read.has_tag("SA") else "",
            })
bam.close()

# refine_breakpoints() takes variant dicts: the first breakpoint is (chrom, start), the second is
# breakpoint2, else pos2, else end, on chrom2 ("" or None means chrom)
refined = refine_breakpoints(
    variants=[
        {"chrom": v.chrom, "start": v.start, "end": v.end,
         "chrom2": v.chrom2, "pos2": v.pos2, "sv_type": v.sv_type.value}
        for v in variants
    ],
    reads=local_reads,
    window=200,
    min_clip=5,
)

# Save refined coordinates
io.dump_json([{
    "sv_type": bp.sv_type,
    "chrom1": bp.bp1.chrom, "pos1": bp.bp1.position,
    "chrom2": bp.bp2.chrom, "pos2": bp.bp2.position,
    "strand1": bp.bp1.strand, "strand2": bp.bp2.strand,
    "confidence": bp.confidence,
    "microhomology": [bp.bp1.microhomology, bp.bp2.microhomology],
    "inserted_sequence": bp.bp1.inserted_sequence or bp.bp2.inserted_sequence,
} for bp in refined], "refined_breakpoints.json")
```

`refine_breakpoints` snaps each breakpoint to the most frequent soft-clip position within `window` bp, locating each read by the reference coordinate of its soft clip. The reads that clip at the chosen position are the breakpoint's supporting reads, and every field below is computed for both breakpoints from their own supporting reads:

- `strand` is the majority `is_reverse` (`'+'` when no read reports it, or on a tie).
- `microhomology` is the longest `detect_microhomology` result over the window of each read's `seq` centred on that read's clip junction. It is a mirror-match heuristic on the read sequence, not true microhomology, which would need the reference sequence.
- `inserted_sequence` is the read bases between the primary alignment and a single supplementary alignment. It is filled only when a supporting read has `seq`, `is_reverse`, `cigar` and an `sa_tag` with one entry located at the other breakpoint; otherwise it is empty.

---

## Example 3: Somatic CNV Detection (Tumor vs Normal)

**Goal:** Identify somatic CNVs (present in tumor, absent in normal) from matched normal.

**Input:** `tumor.bam`, `normal.bam`.

**Approach:** Compute log2(tumor/normal) depth ratios with GC correction (`calculate_log2_ratio`), convert back to a linear ratio, and segment it with `detect_cnv_from_depth`.

```python-snippet
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
                # Placeholder GC estimate: the GC of the latest read in the bin stands in for the
                # window GC. For real analyses compute GC per window from the reference FASTA.
                if read.query_sequence:
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
# detect_cnv_from_depth() expects depth-like values (finite and >= 0) and takes log2(value / median)
# itself, so hand it the linear tumor/normal ratio (2**log2), not the log2 ratio
# (a negative value raises ValueError naming the chromosome).
ratios = {chrom: np.exp2(lr) for chrom, lr in log2_ratios.items()}

cnv_results = detect_cnv_from_depth(
    depth_data=ratios,
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
- Somatic DEL means tumor loss to CN=1 (diploid). Homozygous loss (CN=0) is called HOMODEL and very high gain AMP; the `state in ("DEL", "DUP")` check above leaves both out. Somatic DUP means gain (CN=3+).
- Germline CNVs (present in both tumor and normal) cancel out in log2 ratio; they should be NEUTRAL in somatic difference.
- Use a matched normal whenever possible to control for germline CNVs and GC bias.

---

## Example 4: Multi-Caller Consensus Building

**Goal:** Combine SV calls from Manta, Delly, Lumpy into a high-confidence consensus.

**Input:** Three VCF files (or JSON after parsing).

```python-snippet
import pysam
from metainformant.structural_variants.filtering.merge import merge_callsets
from metainformant.core import io

# Read each VCF with pysam and extract simple dicts
def read_vcf_simple(vcf_path):
    """Convert VCF records to simplified dict format."""
    vcf = pysam.VariantFile(vcf_path, "r")
    variants = []
    for rec in vcf:
        gt = rec.samples[0]["GT"] if len(rec.samples) else None   # allele-index tuple, e.g. (0, 1)
        variants.append({
            "chrom": rec.chrom,
            "start": rec.start,            # pysam: 0-based, like the rest of the module
            "end": rec.stop,               # INFO/END when present
            "sv_type": rec.info.get("SVTYPE", "UNKNOWN"),
            "size": abs(rec.stop - rec.start),
            "quality": rec.qual or 0.0,
            "genotype": "/".join("." if a is None else str(a) for a in gt) if gt else "./.",
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

# Save consensus to VCF (svs_to_vcf_data() is defined in Example 14)
from metainformant.dna.variation.variants import write_vcf
write_vcf(svs_to_vcf_data(consensus), "consensus.vcf")
```

**Survivor-style merging** (distance-based) for callers without reciprocal overlap:

```python
from metainformant.structural_variants.filtering.merge import survivor_merge

merged, stats = survivor_merge(
    vcf_files=["manta.vcf", "delly.vcf", "lumpy.vcf"],   # one file per caller
    max_distance=1000,       # Max bp between breakpoints
    min_callers=2,           # Must be in at least 2 callers
    type_match=True,
    strand_match=False,
)
print(f"Survivor merge: {len(merged)} consensus calls from {stats.n_input_variants} input calls")
```

Like `merge_callsets`, `survivor_merge` returns `(list[MergedVariant], MergeStats)`. File paths are parsed with pysam when it is installed (otherwise with a basic text parser) and the callers are named `caller_0`, `caller_1`, ... in file order; it also accepts a flat list of variant dicts, each tagged with a `_caller` key. The merged `genotype` is the majority vote among the calls' non-missing genotypes, read from each dict's `genotype` key or from the first sample of a VCF (a phased `0|1` becomes `0/1`); a tie goes to the genotype with the larger summed call quality, then to the one that sorts first, and `"./."` means that no call has a genotype. The `consensus_quality` is the mean call quality plus `10 * n_callers * log10(n_callers + 1)`: two callers with calls of quality 80 (genotype `0/1`) and 40 (genotype `1/1`) give `0/1` and 69.54. `merge_callsets` uses `100 * log10(n_callers + 1)` as its bonus and breaks a genotype tie in favour of the genotype seen first.

---

## Example 5: Clinical Interpretation — Gene Overlap + Functional Impact

**Goal:** Prioritize SVs for clinical review based on gene overlap, dosage sensitivity, and TAD disruption.

**Input:** Filtered SV calls, gene annotations, ClinGen dosage database.

```python
from dataclasses import asdict

from metainformant.structural_variants.annotation.overlap import annotate_gene_overlap
from metainformant.structural_variants.annotation.functional_impact import (
    assess_dosage_sensitivity,
    predict_functional_impact,
)
from metainformant.core import io

# Load data
svs = io.load_json("filtered_svs.json")          # SV dicts: chrom, start, end, sv_type, ...
gene_db = [                                      # BED columns: chrom, start, end, name
    {"chrom": chrom, "start": int(start), "end": int(end), "name": name}
    for chrom, start, end, name, *_ in io.read_tsv("genes.bed")
]
hi_db = io.load_json("clinGen_haploinsufficiency.json")   # {gene: score (0-1)}
ts_db = io.load_json("clinGen_triplosensitivity.json")
pli_db = io.load_json("gnomad_pli.json")   # {gene: pLI}
tad_bounds = io.load_json("tad_boundaries.json")   # [{"chrom", "start", "end", "genes"?}, ...]

# ── Step 1: Gene overlap ────────────────────────────────────────────────────
annotated = annotate_gene_overlap(svs, gene_db)
print(f"Annotated {len(annotated)} SVs; "
      f"{sum(1 for v in annotated if v['n_genes_affected']>0)} overlap genes")

# ── Step 2: Functional impact ───────────────────────────────────────────────
priority_svs = []
for variant in annotated:
    # Every database you pass feeds the dosage flag: a deletion over a haploinsufficient gene
    # (HI >= 0.5, pLI >= 0.9 or LOEUF < 0.35) is HIGH / dosage_loss, and a duplication over a
    # triplosensitive gene (TS >= 0.5) is HIGH / dosage_gain.
    impact = predict_functional_impact(
        variant=variant,
        gene_annotations=gene_db,
        haploinsufficiency_db=hi_db,
        tad_boundaries=tad_bounds,
        triplosensitivity_db=ts_db,
        pli_db=pli_db,
    )

    # Clinical prioritization criteria
    is_high_impact = (
        impact.impact_level == "HIGH" and
        impact.pathogenicity_score >= 0.7 and
        (impact.tad_disrupted or len(impact.dosage_sensitive_genes) > 0)
    )

    if is_high_impact:
        # dump_json() needs plain JSON types: drop the OverlapResult objects, flatten the dataclass
        priority_svs.append({
            "variant": {key: val for key, val in variant.items() if key != "gene_overlaps"},
            "impact": asdict(impact),
        })
        print(f"\nHIGH-PRIORITY SV:")
        print(f"  {impact.variant_id}: {impact.impact_type}")
        print(f"  Pathogenicity: {impact.pathogenicity_score:.3f}")
        print(f"  Genes: {', '.join(impact.affected_genes)}")
        for gene in impact.dosage_sensitive_genes:
            # Look the flagged gene up again to print its scores (HI and pLI are printed)
            ds = assess_dosage_sensitivity(
                gene,
                haploinsufficiency_db=hi_db,
                triplosensitivity_db=ts_db,
                pli_db=pli_db,
            )
            print(f"  Dosage-sensitive: {gene} "
                  f"(HI={ds.haploinsufficiency_score:.2f}, pLI={ds.pli_score:.2f})")
        if impact.tad_disrupted:
            n_broken = impact.details.get("tad_prediction", {}).get("n_boundaries", 0)
            print(f"  ⚠ TAD boundary disruption: {n_broken} boundaries")

print(f"\n{len(priority_svs)} variants meet clinical review criteria")

# Save prioritized list
io.dump_json(priority_svs, "clinical_priority.json")
```

**Sample output** (synthetic inputs: a 222 kb deletion over three genes, one of them with HI=0.92, with one TAD boundary inside the deletion, plus a small duplication over a benign gene):
```
Annotated 2 SVs; 2 overlap genes

HIGH-PRIORITY SV:
  chr5:150234100-150456789: dosage_loss
  Pathogenicity: 0.825
  Genes: GENE1, GENE2, GENE3
  Dosage-sensitive: GENE1 (HI=0.92, pLI=0.97)
  ⚠ TAD boundary disruption: 1 boundaries

1 variants meet clinical review criteria
```

The `0.7` pathogenicity cut-off is demanding. The score is a weighted sum whose size term is a sigmoid centred at 100 kb, so a 2 kb deletion over a single HI=0.92 gene scores about 0.53 and a 222 kb one about 0.67; this example only clears the bar because it also spans three genes and a TAD boundary.

---

## Example 6: Population Frequency Comparison with gnomAD

**Goal:** Filter SVs that are rare (AF ≤ 1%) in a gnomAD-SV-style population database.

**Input:** SV calls, population frequency database (`gnomad_sv.json`).

```python
from metainformant.structural_variants.filtering.quality_filter import filter_by_frequency
from metainformant.core import io

svs = io.load_json("sv_calls.json")
gnomad = io.load_json("gnomad_sv.json")  # [{"chrom", "start", "end", "sv_type", "af"}, ...]

# ── Filter ───────────────────────────────────────────────────────────────────
filtered, stats = filter_by_frequency(
    variants=svs,
    population_db=gnomad,
    max_af=0.01,                  # Keep AF <= 1%
    match_window=200,             # Only database entries within 200 bp of the call are candidates
    min_reciprocal_overlap=0.5,   # 50% reciprocal overlap required
)

print(f"Filtered: removed {stats.filtered_count} common variants")
print(f"Remaining: {stats.output_count} rare variants")

# FilterStats does not list the removed calls; diff the input and output lists to see what was removed
def key(v):
    return (v["chrom"], v["start"], v["end"], v["sv_type"])

kept = {key(v) for v in filtered}
removed = [v for v in svs if key(v) not in kept]
if removed:
    print("\nCommon variants removed (AF > 1%):")
    for v in removed[:5]:
        print(f"  {v['chrom']}:{v['start']}-{v['end']} ({v['sv_type']})")

# Kept calls that matched a database entry carry its frequency as "population_af"
for v in filtered:
    if "population_af" in v:
        print(f"  kept {v['chrom']}:{v['start']} (matched AF={v['population_af']:.4f})")

# Save rare variants
io.dump_json(filtered, "rare_svs.json")
```

**Population database format.** `population_db` is either a list of dicts, as above (`chrom`, `start`, `end`, optional `sv_type`, and `af` or `allele_frequency`):
```json
[
  {"chrom": "chr1", "start": 123456, "end": 123789, "sv_type": "DEL", "af": 0.003},
  {"chrom": "chr1", "start": 987654, "end": 988123, "sv_type": "DUP", "af": 0.012}
]
```
or a dict mapping `"chrom:start-end"` (optionally `"chrom:start-end-SVTYPE"`) to the allele frequency, e.g. `{"chr1:123456-123789-DEL": 0.003}`. A call matches an entry when the SV types agree (if the entry gives one) and both reciprocal-overlap fractions reach `min_reciprocal_overlap`; the highest matching AF is used. Calls without a match are kept, and a call that already has an `af` (or `allele_frequency`) key is compared directly without a lookup. Every frequency must lie in [0, 1]. A database in a layout that is not supported, such as a dict keyed by chromosome, raises `ValueError` naming the supported layouts when none of its entries can be read; when only some entries are unreadable they are skipped and counted in a logged warning, and the rest of the database is used.

**Per-population frequency filtering:**

There is no population argument. The `af` value of each database entry is what gets compared with `max_af`, so put the population-specific frequency in that field:

```python
eur_db = [{**entry, "af": entry["af_eur"]} for entry in gnomad]  # entries that carry per-population AFs
filtered, _ = filter_by_frequency(
    variants=svs,
    population_db=eur_db,
    max_af=0.001,  # More stringent
)
```

---

## Example 7: Visualizing Complex Rearrangements — Circos + Coverage

**Goal:** Produce publication-quality figures for a sample with multiple rearrangements.

**Input:** SV calls (a list of dicts, e.g. the output of the filtering steps) + per-bin read depth.

```python-snippet
from metainformant.structural_variants.visualization.plots import (
    plot_circos,
    plot_coverage_track,
    plot_sv_size_distribution,
)

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
# depth_data: {chrom: array of per-bin read depth}, e.g. the dict built in Example 1.
# The track is labelled "Read Depth", so pass depths, not CNVResult.log2_ratios.
bin_size = 1000   # the window_size used to build depth_data

for chrom in ["chr5", "chr8", "chr17"]:  # Chromosomes with interesting SVs
    if chrom not in depth_data:
        continue
    svs_on_chrom = [v for v in filtered_svs if v["chrom"] == chrom]

    fig = plot_coverage_track(
        coverage=depth_data[chrom],
        variants=svs_on_chrom,
        region=(chrom, 0, chrom_sizes[chrom]),
        output_path=f"coverage_{chrom}.png",
        title=f"Read Depth — {chrom}",
        bin_size=bin_size,
    )
    print(f"Coverage track: coverage_{chrom}.png")

# ── Size distribution ───────────────────────────────────────────────────────
fig = plot_sv_size_distribution(
    variants=filtered_svs,
    output_path="size_distribution.png",
    log_scale=True,
)
```

**Result:** `circos.png` — circular genome with colored chords connecting breakpoints; `coverage_chr5.png` — read-depth track with SV rectangles overlay. A chord runs from `chrom:start` to the second breakpoint: `chrom2:pos2` for a translocation (the call dicts from `call_structural_variants` carry both), and `chrom:end` for every other call. A hand-built dict with `chrom2` set to another chromosome and no `pos2` is read the way VCF writes it: `end` is the position on `chrom2`. For a genome-wide view of CNV segments use `plot_cnv_profile(segments, chromosomes, output_path=...)`, where each segment is a dict with `chrom`, `start`, `end`, `mean_log2ratio` and `state`.

---

## Example 8: Cohort SV Genotyping and Association Testing

**Goal:** Genotype 100 known SVs across 500 samples and test for phenotype association.

**Input:**
- `reference_svs.json`: 100 variant definitions (`chrom`, `start`, `end`, `sv_type`) with the per-sample evidence to genotype from, under a `samples` key: `{sample_id: {"depth_ratio": ..., "split_reads": ..., "discordant_pairs": ..., "total_reads": ...}}`. `genotype_sv_population` does not read BAM files; the evidence values come from your own depth/split-read counting. Missing evidence defaults to a reference call (`depth_ratio` 1.0, no supporting reads).
- `samples.tsv`: sample IDs and population labels
- `phenotypes.tsv`: sample IDs + phenotype (binary case/control)

```python
import numpy as np
import pandas as pd
from metainformant.core import io
from metainformant.structural_variants.population.sv_population import (
    genotype_sv_population,
    sv_allele_frequency,
    sv_association_test,
    sv_population_structure,
)

# ── Load cohort ──────────────────────────────────────────────────────────────
samples_df = pd.read_csv("samples.tsv", sep="\t")   # columns: sample_id, population
pheno_df = pd.read_csv("phenotypes.tsv", sep="\t")   # columns: sample_id, phenotype (0/1)

sample_ids = samples_df["sample_id"].tolist()
populations = samples_df["population"].tolist()

# Reference SV panel (100 SVs) with per-sample evidence (see Input above)
ref_svs = io.load_json("reference_svs.json")

# ── Genotype across all samples ──────────────────────────────────────────────
genotype_result = genotype_sv_population(
    sv_calls=ref_svs,
    samples=sample_ids,
    method="depth",  # "depth" uses depth_ratio; "split" uses split_reads / discordant_pairs / total_reads
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
print(f"  Polymorphic SVs: {freq['n_polymorphic']} / {len(ref_svs)}")

for pop in set(populations):
    pop_idx = [i for i, p in enumerate(populations) if p == pop]
    pop_af = freq["population_frequencies"][pop]
    print(f"  {pop} (n={len(pop_idx)}): mean AF={np.mean(pop_af):.4f}")

# ── Population structure (PCA) ───────────────────────────────────────────────
# Adjust for ancestry in association tests
pca = sv_population_structure(gt_matrix, n_components=10)
# variance_explained is each component's share of the total variance (the sum of all eigenvalues of
# the covariance matrix), so the 10 returned values sum to less than 1 here (100 SVs give 100
# components); they sum to 1 only when every component is returned
# (n_components >= min(number of SVs, number of samples))
print(f"\nPCA — PC1 variance explained: {pca['variance_explained'][0]:.2%}")
print(f"         PC2 variance explained: {pca['variance_explained'][1]:.2%}")

# ── Association testing ──────────────────────────────────────────────────────
# Align phenotypes to the sample order used for genotyping
phenotypes = pheno_df.set_index("sample_id").loc[sample_ids, "phenotype"].to_numpy()  # 0 or 1

results = []
for i, sv in enumerate(ref_svs):
    rec = sv_association_test(
        genotypes=gt_matrix[i],
        phenotypes=phenotypes,
        covariates=pca["pcs"][:, :5],  # Top 5 PCs (n_samples × 5) as covariates
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

**Goal:** Use long-read alignments to tighten the breakpoints of SVs that were called from short reads.

`refine_breakpoints` snaps each breakpoint to the most frequent soft-clip position among the reads it is given, so long reads that clip at the breakpoint give a base-pair-level estimate. It locates a read by the reference coordinate of its soft clip, so a long read that starts far upstream of the breakpoint still counts. It works on plain dict records, so the `LongReadAlignment` objects returned by the long-read BAM reader have to be converted first.

```python
from metainformant.core import io
from metainformant.longread.io.bam import read_long_read_bam   # requires pysam
from metainformant.structural_variants.detection.breakpoints import refine_breakpoints

# Load long-read alignments (BAM from minimap2 / pbmm2) as LongReadAlignment objects
alignments = read_long_read_bam("longreads.bam", include_supplementary=True, min_mapping_quality=20)

# Convert to the record format refine_breakpoints() reads: name, chrom, pos and cigar, plus
# is_reverse (breakpoint strand), seq (microhomology and inserted sequence) and sa_tag (inserted sequence).
longreads = [
    {"name": a.read_name, "chrom": a.reference_name, "pos": a.reference_start, "cigar": a.cigar_string,
     "is_reverse": a.is_reverse, "seq": a.query_sequence, "sa_tag": a.tags.get("SA", "")}
    for a in alignments
    if not a.is_unmapped
]

# SV calls from short reads: dicts with chrom, start, end, sv_type (plus chrom2 and pos2 for translocations)
sv_calls = io.load_json("shortread_svs.json")

# Refine using long reads
refined = refine_breakpoints(
    variants=sv_calls,
    reads=longreads,
    window=1000,        # soft clips within 1000 bp of a breakpoint count, wherever their read starts
    min_clip=10,
)

for bp_pair in refined:
    print(f"{bp_pair.sv_type}: {bp_pair.bp1.chrom}:{bp_pair.bp1.position} "
          f"(support {bp_pair.bp1.support}) -> "
          f"{bp_pair.bp2.chrom}:{bp_pair.bp2.position} (support {bp_pair.bp2.support})  "
          f"confidence={bp_pair.confidence:.2f}")

# Output refined breakpoints for VCF conversion
io.dump_json([{
    "chrom": bp.bp1.chrom,
    "pos": bp.bp1.position,
    "chrom2": bp.bp2.chrom,
    "pos2": bp.bp2.position,
    "sv_type": bp.sv_type,
    "support": [bp.bp1.support, bp.bp2.support],
    "confidence": bp.confidence,
} for bp in refined], "precise_breakpoints.json")
```

Limits to keep in mind with long reads:

- Only soft clips (`S` in the CIGAR) are evidence. Hard-clipped supplementary alignments contribute no breakpoint, and no `inserted_sequence` can be derived from a record whose hard clips make its CIGAR longer than its stored `seq`; a soft-clipped record with its full `seq` and a one-entry `sa_tag` at the other breakpoint does give one. For long-read insertion calls use `metainformant.longread.analysis.structural.detect_sv_from_long_reads`.
- A read is considered when one of its soft clips lies within `window` bp of the approximate breakpoint, wherever its alignment starts, so a long read that starts tens of kb upstream and clips at a right-hand breakpoint counts. A very wide `window` lets a neighbouring clip cluster win the "most frequent position" vote. Check `bp.support` (0 means the breakpoint was left at its input position).
- `microhomology` is a mirror-match heuristic on the read bases around each supporting read's clip junction, not true microhomology, which would need the reference sequence.

---

## Example 10: Batch Processing Pipeline for 1000+ Samples

**Goal:** Scalable processing of a large cohort, one process-pool task per sample.

**Strategy:** Per-sample parallelization with `parallel.process_map` (a process pool).

```python-snippet
from dataclasses import asdict
from pathlib import Path

import numpy as np

from metainformant.core import io
from metainformant.core.execution import parallel
from metainformant.core.utils import logging
from metainformant.structural_variants.detection.cnv import detect_cnv_from_depth
from metainformant.structural_variants.detection.sv_calling import call_structural_variants
from metainformant.structural_variants.filtering.merge import merge_callsets
from metainformant.structural_variants.filtering.quality_filter import filter_by_quality

# List of all sample BAMs
bam_dir = Path("/data/cohort/bams/")
bam_files = sorted(bam_dir.glob("*.bam"))


def process_sample(bam_path):
    """Full SV pipeline for one sample (a top-level function, so the process pool can pickle it)."""
    sample_id = bam_path.stem

    # CNV detection
    depth_data = compute_depth(bam_path)  # Your helper: {chrom: per-bin depth array}
    cnv = detect_cnv_from_depth(depth_data, window_size=1000)
    cnv_calls = [asdict(seg) for res in cnv.values() for seg in res.segments
                 if seg.state != "NEUTRAL"]

    # SV calling (the filters work on dicts, so convert the StructuralVariant objects)
    alns = extract_alignments(bam_path)  # Your helper (see Example 2)
    svs = []
    for sv in call_structural_variants(alns, min_support=3):
        record = asdict(sv)
        record["sv_type"] = sv.sv_type.value   # plain string instead of the SVType enum
        svs.append(record)

    # Filter (support is read from the nested evidence counts)
    filtered, _ = filter_by_quality(svs, min_qual=20, min_support=3)

    # Save per-sample results (dump_json creates the directory)
    io.dump_json({"svs": filtered, "cnv": cnv_calls}, Path(f"results/per_sample/{sample_id}.json"))

    return sample_id, len(filtered)

# ── Parallel map (per-sample) ─────────────────────────────────────────────────
# On macOS/Windows, put everything below under `if __name__ == "__main__":`
# (workers are spawned and re-import the script).
results = parallel.process_map(process_sample, bam_files, max_workers=12)

print("Per-sample SV counts:")
for sample_id, count in results:
    print(f"  {sample_id}: {count} SVs")

# ── Cross-sample merging (optional) ─────────────────────────────────────────
# Treat each sample as a "caller" to get a cohort-level consensus SV set
# (MergedVariant.n_callers is then the number of samples carrying the SV):
all_samples = {}
for bam_path in bam_files:
    sid = bam_path.stem
    all_samples[sid] = io.load_json(f"results/per_sample/{sid}.json")["svs"]

merged_cohort, stats = merge_callsets(
    callsets=all_samples,
    min_overlap=0.5,
)
print(f"Cohort merged SVs: {stats.n_output_variants}")

# ── Logging summary ──────────────────────────────────────────────────────────
logging.configure_logging_from_env()
logger = logging.get_logger("sv_cohort")
summary = [
    f"Processed {len(bam_files)} samples",
    f"Total SVs detected: {sum(count for _, count in results)}",
    f"Average per sample: {np.mean([c for _, c in results]):.1f}",
]
logger.info("\n".join(summary))
```

`process_map` takes no chunking argument (only `max_workers`, `timeout`, `ordered`) and submits one task per item. Each worker builds its own in-memory alignment list, so peak memory is roughly the number of workers times the memory of one sample.

---

## Example 11: Integration with DNA Module for Variant Effect Prediction

**Goal:** Annotate SVs against gene models from a GENCODE GTF, classify how each SV hits each gene (whole-gene deletion, partial overlap, SV inside a gene), and bring SV calls in through the DNA module's VCF parser.

The DNA module provides the VCF reader and writer (`parse_vcf`, `write_vcf`; see Example 14) but no GTF or gene-model loader, so the GTF is read with a short helper here. `annotate_gene_overlap` accepts gene dicts (`chrom`, `start`, `end`, `name`) or `GenomicInterval` objects; load `feature="exon"` instead of `"gene"` to get exon-level overlaps from the same call (each exon is then reported as its own entry).

```python-snippet
from metainformant.dna.variation.variants import parse_vcf
from metainformant.structural_variants.annotation.overlap import annotate_gene_overlap


def load_gtf_features(path, feature="gene"):
    """Read one feature type from a GTF into 0-based, half-open dicts."""
    records = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != feature:
                continue
            attrs = dict(item.strip().split(" ", 1) for item in fields[8].split(";") if item.strip())
            records.append({
                "chrom": fields[0],
                "start": int(fields[3]) - 1,      # GTF is 1-based, inclusive
                "end": int(fields[4]),
                "name": attrs.get("gene_name", "").strip('"'),
                "feature_type": feature,
                "strand": fields[6],
            })
    return records


gene_db = load_gtf_features("gencode.v44.annotation.gtf")

# SV calls from a VCF (e.g. the file written in Example 2) through the DNA module's parser
vcf = parse_vcf("sv_calls.vcf")
svs = [
    {"chrom": rec["chrom"], "start": rec["pos"] - 1, "end": rec["info"].get("END", rec["pos"]),
     "sv_type": rec["info"].get("SVTYPE", "UNKNOWN"), "quality": rec["qual"] or 0.0}
    for rec in vcf["variants"]
]

annotated = annotate_gene_overlap(svs, gene_db)

for v in annotated:
    for overlap in v["gene_overlaps"]:
        gene = overlap.feature.name
        relationship = overlap.relationship  # "contains", "contained_in" or "overlap"
        frac = overlap.overlap_fraction_feature  # Fraction of the gene covered by the SV

        if frac > 0.5:
            print(f"  {v['sv_type']} likely disrupts gene {gene} ({relationship}, covers {frac:.1%})")

# "annotated" can be passed on to predict_functional_impact (Example 5)
```

`"contains"` means the SV spans the whole feature; `"contained_in"` means the SV lies inside the feature.

---

## Example 12: Multi-Omics Correlation: SVs and Gene Expression

**Goal:** Test whether a CNV correlates with expression of overlapped genes (eQTL-like).

```python
import pickle
from dataclasses import asdict

import numpy as np
from scipy.stats import pearsonr

from metainformant.core import io
from metainformant.structural_variants.annotation.overlap import GenomicInterval, annotate_gene_overlap

# Load matched data
with open("cnv_results.pkl", "rb") as fh:
    cnv_results = pickle.load(fh)  # {sample_id: {chrom: CNVResult}}, i.e. detect_cnv_from_depth() per sample
gene_expression = io.read_csv("gene_expression.tsv", sep="\t", index_col=0)  # genes × samples (DataFrame)
gene_db = {  # gene name -> interval; BED columns: chrom, start, end, name
    name: GenomicInterval(chrom, int(start), int(end), name=name, feature_type="gene")
    for chrom, start, end, name, *_ in io.read_tsv("genes.bed")
}

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

# Genes overlapped by any non-neutral CNV segment in any sample
all_gene_names = set()
for sample_cnv in cnv_results.values():
    altered = [asdict(seg) for res in sample_cnv.values() for seg in res.segments
               if seg.state != "NEUTRAL"]
    for annotated in annotate_gene_overlap(altered, list(gene_db.values())):
        all_gene_names.update(annotated["overlapping_genes"])

# Correlation for one gene
sample_ids = list(cnv_results.keys())
gene = "BRCA1"
cn = gene_cn_state(cnv_results, gene_db[gene])
expr = gene_expression.loc[gene, sample_ids].values

r, p = pearsonr(cn, expr)
print(f"{gene}: r={r:.3f}, p={p:.2e}")

# Loop over all genes affected by SVs
correlations = []
for gene_name in sorted(all_gene_names):
    if gene_name not in gene_expression.index:
        continue
    cn = gene_cn_state(cnv_results, gene_db[gene_name])
    if np.ptp(cn) == 0:   # same copy number in every sample: nothing to correlate
        continue
    expr = gene_expression.loc[gene_name, sample_ids].values
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

Correlate SV breakpoints with splicing abnormalities. The RNA-seq splicing module has `find_novel_junctions` (detected junctions that are not near any annotated junction) and `differential_splicing` (a PSI comparison between two groups of samples). Neither takes a gene, condition or region argument, so the SV-to-junction matching is done here:

```python
from metainformant.core import io
from metainformant.rna.splicing.splice_analysis import find_novel_junctions

# SVs annotated with gene overlaps (annotate_gene_overlap output, see Example 5, saved as JSON
# without the OverlapResult objects: overlapping_genes and n_genes_affected are what is used here)
svs = io.load_json("filtered_svs.json")
intragenic = [v for v in svs if v.get("n_genes_affected", 0) > 0]

# RNA-seq junctions ({"chrom", "start", "end", "strand", "read_count"}, e.g. from
# rna.splicing.splice_sites.detect_splice_junctions) and the annotated junction set
junctions = io.load_json("rnaseq_junctions.json")
known_junctions = io.load_json("annotated_junctions.json")   # {"chrom", "start", "end"}

# Junctions whose start and end are not both within 10 bp of an annotated junction
novel = find_novel_junctions(junctions, known_junctions, max_distance=10)

# For each intragenic SV, look for novel junctions whose start or end lies within 50 bp of one of its breakpoints
for sv in intragenic:
    genes = ", ".join(sv["overlapping_genes"])
    for bp_pos in (sv["start"], sv["end"]):  # or refined breakpoints (Example 2)
        for junction in novel:
            if junction["chrom"] != sv["chrom"]:
                continue
            if abs(junction["start"] - bp_pos) < 50 or abs(junction["end"] - bp_pos) < 50:
                print(f"SV breakpoint {sv['chrom']}:{bp_pos} in {genes} co-locates with a novel junction "
                      f"{junction['start']}-{junction['end']} ({junction['read_count']} reads)")
```

To test whether a junction is used differently in samples that carry the SV, compute a PSI per sample (`compute_psi`) and compare the two groups; `delta_psi` is the carrier mean minus the non-carrier mean:

```python
from metainformant.rna.splicing.splice_analysis import differential_splicing

noncarrier_psi = [0.81, 0.78, 0.84, 0.79]   # PSI of the junction in samples without the SV
carrier_psi = [0.62, 0.55, 0.66]            # PSI in samples carrying the SV
result = differential_splicing(noncarrier_psi, carrier_psi)  # method="empirical_bayes" by default
print(f"delta PSI={result['delta_psi']:+.2f}, p={result['p_value']:.2g}, significant={result['significant']}")
```

---

## Example 14: Export to Standard Formats

**VCF Export:**

The structural variant module has no VCF writer of its own; use the generic `write_vcf` from the DNA module, which takes a parsed-VCF dictionary (`metadata`, `samples`, `variants`). Symbolic alleles (`<DEL>`, `<DUP>`, ...) with `SVTYPE`, `END` and `SVLEN` in INFO describe the SVs:

```python-snippet
from metainformant.dna.variation.variants import write_vcf

def svs_to_vcf_data(svs):
    """Build the VCF dict expected by write_vcf from StructuralVariant or MergedVariant objects.

    Inter-chromosomal calls (TRA/BND) are skipped: they need BND records with mate
    coordinates (chrom2 and pos2 of a TRA call), which this minimal helper does not write.
    """
    records = []
    for sv in svs:
        sv_type = getattr(sv.sv_type, "value", sv.sv_type)
        if sv_type in ("TRA", "BND") or getattr(sv, "chrom2", ""):
            continue
        records.append({
            "chrom": sv.chrom,
            "pos": sv.start + 1,                      # VCF is 1-based; SV starts are 0-based
            "id": f"{sv_type}_{sv.chrom}_{sv.start}",
            "ref": "N",
            "alt": [f"<{sv_type}>"],
            "qual": getattr(sv, "quality", getattr(sv, "consensus_quality", 0.0)),
            "filter": ["PASS"],
            # size 0 means "not measured" (insertions called from split reads): leave SVLEN out
            "info": {"SVTYPE": sv_type, "END": sv.end,
                     **({"SVLEN": -sv.size if sv_type == "DEL" else sv.size} if sv.size else {})},
        })
    return {"metadata": {"fileformat": "VCFv4.2"}, "samples": [], "variants": records}

# From a list of StructuralVariant objects (e.g. the output of call_structural_variants)
write_vcf(svs_to_vcf_data(variants), "output.vcf")
```

`write_vcf` writes each `metadata` entry as a `##key=value` header line and nothing else, so the file carries no `##INFO` or `##contig` definitions unless you add them. The output can be read back with `dna.variation.variants.parse_vcf` (Example 11).

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
        "chrom2": sv.chrom2,   # second breakpoint of a translocation ("" otherwise)
        "pos2": sv.pos2,       # None otherwise
        "type": sv.sv_type,
        "size": sv.size,       # 0 = not measured (translocations and insertions)
        "qual": sv.quality,
        "gt": sv.genotype,
    })
df = pd.DataFrame(rows)
df.to_csv("sv_table.csv", index=False)
```

---

## Example 15: Custom SV Type: Mobile Element Insertions (MEIs)

**Goal:** Flag candidate Alu/L1/SVA insertions among the insertion calls, using split-read counts and a poly-A tail check.

`classify_sv_type` has no MEI class, so the MEI logic is a separate heuristic applied to the `INS` calls. An `INS` call is a split-read cluster whose two breakpoints nearly coincide (< 50 bp apart). Its `size` is 0, meaning "not measured": split-read evidence does not measure the inserted length, and the distance between the breakpoints is not that length. `INS` calls are exempt from `min_sv_size`, so the default `min_sv_size=50` keeps them and no special setting is needed.

```python
from metainformant.structural_variants.detection.sv_calling import (
    SVEvidence,
    SVType,
    call_structural_variants,
)

def classify_mei(evidence: SVEvidence, flanking_sequence=None):
    """Heuristic: many split reads plus a poly-A tail in the clipped sequence."""
    if evidence.split_reads >= 5:
        # Check for poly-A tail in clipped sequence
        if flanking_sequence and flanking_sequence[-20:].count('A') > 15:
            return "MEI_ALU"  # or LINE1, SVA based on length/signature
    return None

# Keep insertion calls (alignments: read records as in Example 2)
variants = call_structural_variants(alignments, min_support=5)
insertions = [sv for sv in variants if sv.sv_type == SVType.INS]

# Apply the heuristic; clipped_sequences is {breakpoint position: soft-clipped read sequence},
# which you extract from the reads yourself
mei_calls = [sv for sv in insertions
             if classify_mei(sv.evidence, clipped_sequences.get(sv.start))]
```

---

## What's Next?

- **[ARCHITECTURE.md](ARCHITECTURE.md)** — Understand how these examples fit into the system design
- **[PERFORMANCE.md](PERFORMANCE.md)** — Optimize for your dataset size
- **[TROUBLESHOOTING.md](TROUBLESHOOTING.md)** — Diagnose issues when things go wrong
- **[detection.md](detection.md)**, **[annotation.md](annotation.md)**, **[filtering.md](filtering.md)**, **[population.md](population.md)** — Deep-dive into each submodule
