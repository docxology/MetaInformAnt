# Structural Variants Module Architecture

Comprehensive overview of the structural variant analysis module, including component responsibilities, data flow patterns, integration points, and design rationale.

## Table of Contents

1. [Module Overview](#1-module-overview)
2. [Directory Structure](#2-directory-structure)
3. [Component Responsibilities](#3-component-responsibilities)
4. [Data Flow](#4-data-flow)
5. [Key Data Structures](#5-key-data-structures)
6. [Integration Points](#6-integration-points)
7. [Design Patterns & Rationale](#7-design-patterns--rationale)
8. [Algorithm Details](#8-algorithm-details)
9. [Configuration Architecture](#9-configuration-architecture)
10. [Extensibility](#10-extensibility)
11. [Performance Characteristics](#11-performance-characteristics)
12. [Design Decisions & Trade-offs](#12-design-decisions--trade-offs)
13. [Testing Philosophy](#13-testing-philosophy)
14. [Related Documentation](#14-related-documentation)

---

## 1. Module Overview

The `structural_variants` module provides analysis of genomic structural variants (SVs) from alignment records and read-depth arrays. It works on in-memory Python structures (alignment dicts, depth arrays, variant dicts); reading BAM or BED files is left to the caller (only `survivor_merge()` reads VCF paths itself), and the `plot_*` functions write figure files. It covers five functional domains:

| Domain | Capabilities | Primary Entry Points |
|--------|-------------|---------------------|
| **Detection** | CNV (CBS), SV (split/discordant reads), breakpoint refinement | `cnv.detect_cnv_from_depth()`, `sv_calling.call_structural_variants()`, `breakpoints.refine_breakpoints()` |
| **Annotation** | Gene/regulatory overlap, functional impact, dosage sensitivity, TAD disruption | `overlap.annotate_gene_overlap()`, `functional_impact.predict_functional_impact()` |
| **Filtering** | Quality/size/frequency filters, blacklist removal, multi-caller merging | `quality_filter.filter_by_quality()`, `merge.merge_callsets()` |
| **Population** | Cohort genotyping, allele freq, association, PCA, LD | `sv_population.genotype_sv_population()`, `sv_allele_frequency()`, `sv_association_test()` |
| **Visualization** | Circos plots, coverage tracks, size histograms | `plots.plot_circos()`, `plots.plot_coverage_track()` |

**SV types (`SVType`):** `DEL` (deletion), `DUP` (duplication), `INV` (inversion), `TRA` (translocation), `INS` (insertion), `BND` (breakend), `UNKNOWN` (evidence too weak to classify). **CNV states (`CNVSegment.state`):** `DEL`, `DUP`, `NEUTRAL`, `AMP` (amplification), `HOMODEL` (homozygous deletion).

---

## 2. Directory Structure

```
src/metainformant/structural_variants/
├── __init__.py                      # Public API: exports all submodules
├── AGENTS.md                        # Agent directives
│
├── detection/                       # SV detection algorithms
│   ├── __init__.py                  # Exports: cnv, sv_calling, breakpoints
│   ├── cnv.py                       # ~750 lines: CBS, CNV state assignment, log2 ratio
│   │   ├── CNVSegment (dataclass)
│   │   ├── CNVResult (dataclass)
│   │   ├── detect_cnv_from_depth()              # Main entry
│   │   ├── segment_coverage()                   # CBS core algorithm
│   │   ├── _cbs_recursive(), _cbs_pvalue()      # Recursive segmentation, Gumbel p-value
│   │   ├── call_cnv_states()                    # State classification
│   │   ├── calculate_log2_ratio()               # Tumor/normal ratio, optional GC correction
│   │   ├── merge_adjacent_segments()
│   │   └── _assign_cnv_state(), _calculate_segment_confidence(), _gc_correct(), _gc_correct_py(), _center_on_median()
│   │
│   ├── sv_calling.py                # ~800 lines: split/discordant evidence
│   │   ├── SVType (enum: DEL, DUP, INV, TRA, INS, BND, UNKNOWN)
│   │   ├── StructuralVariant (dataclass)
│   │   ├── SVEvidence (dataclass)
│   │   ├── InsertSizeStats (dataclass)
│   │   ├── call_structural_variants()           # Main entry
│   │   ├── detect_split_reads()                 # Soft-clip evidence
│   │   ├── detect_discordant_pairs()            # Insert size/orientation
│   │   ├── classify_sv_type()                   # Evidence → SVType
│   │   ├── genotype_sv()                        # 0/0, 0/1, 1/1 (./. without reads)
│   │   ├── _estimate_insert_size()
│   │   ├── _parse_cigar_clips(), _parse_cigar_string(), _parse_sa_tag()
│   │   └── _cluster_evidence(), _merge_evidence(), _calculate_sv_quality()
│   │
│   └── breakpoints.py               # ~510 lines: refinement, microhomology
│       ├── Breakpoint (dataclass)
│       ├── BreakpointPair (dataclass)
│       ├── refine_breakpoints()                # Consensus position
│       ├── detect_microhomology()               # Junction sequence homology
│       ├── cluster_breakpoints()
│       ├── calculate_breakpoint_confidence()
│       └── _collect_clip_positions(), _parse_clips_from_cigar(), _consensus_position(), _position_std()
│
├── annotation/                      # Variant effect prediction
│   ├── __init__.py
│   ├── overlap.py                   # ~520 lines: interval overlap
│   │   ├── GenomicInterval (dataclass)
│   │   ├── OverlapResult (dataclass)
│   │   ├── IntervalIndex (class: sorted-start binary search)
│   │   ├── annotate_gene_overlap()             # Main gene overlap
│   │   ├── annotate_regulatory_overlap()       # Enhancers/promoters
│   │   ├── find_nearest_gene()
│   │   ├── calculate_overlap_fraction()
│   │   └── _normalize_intervals(), _compute_overlap_bp(), _classify_relationship()
│   │
│   └── functional_impact.py         # ~570 lines: impact prediction
│       ├── FunctionalImpact (dataclass)
│       ├── DosageSensitivity (dataclass)
│       ├── TADPrediction (dataclass)
│       ├── predict_functional_impact()         # Main entry (all-in-one)
│       ├── assess_dosage_sensitivity()         # HI/TS/pLI/LOEUF from DBs
│       ├── predict_tad_disruption()            # Boundary crossing
│       ├── score_pathogenicity()               # Composite [0,1] score
│       ├── _find_overlapping_genes()
│       ├── _count_coding_genes()
│       └── _classify_impact()                  # Impact level/type
│
├── filtering/                       # Quality control & merging
│   ├── __init__.py                  # Exports: quality_filter, merge
│   ├── quality_filter.py            # ~480 lines: filters
│   │   ├── FilterStats (dataclass)
│   │   ├── filter_by_quality()                # Qual/support/MAPQ/GQ
│   │   ├── filter_by_size()                   # Size range
│   │   ├── filter_by_frequency()              # Population AF
│   │   ├── apply_blacklist()                  # Known bad regions
│   │   └── _build_population_lookup(), _lookup_population_frequency()
│   │
│   └── merge.py                     # ~715 lines: multi-caller consensus
│       ├── MergedVariant (dataclass)
│       ├── MergeStats (dataclass)
│       ├── merge_callsets()                   # Reciprocal overlap
│       ├── survivor_merge()                   # Distance-based
│       ├── calculate_reciprocal_overlap()
│       ├── deduplicate_variants()
│       ├── _union(), _find() (union-find)
│       ├── _median_position(), _majority_sv_type(), _consensus_genotype()
│       └── _load_vcf_variants(), _parse_vcf_basic()   # VCF input for survivor_merge
│
├── population/                      # Cohort-scale analysis
│   ├── __init__.py
│   └── sv_population.py             # ~840 lines: genotyping, stats, PCA
│       ├── genotype_sv_population()           # Genotype matrix
│       ├── sv_allele_frequency()              # AF, MAF, SFS, per-pop
│       ├── sv_association_test()              # Linear/logistic regression
│       ├── sv_population_structure()          # PCA
│       ├── sv_ld_analysis()                   # LD with SNPs
│       ├── merge_sv_callsets()                # Cross-sample merging
│       ├── _genotype_by_depth(), _genotype_by_split()
│       ├── _linear_regression(), _logistic_regression()
│       └── _should_merge(), _build_consensus(), _normal_cdf_np()
│
└── visualization/                   # Plotting
    ├── __init__.py
    └── plots.py                     # ~920 lines: Circos, coverage, histograms
        ├── plot_circos()                      # Genome-wide chord diagram
        ├── plot_coverage_track()              # Read-depth track with SV overlay
        ├── plot_sv_size_distribution()
        ├── plot_sv_type_summary()
        ├── plot_breakpoint_detail()
        ├── plot_cnv_profile()                 # Genome-wide CNV segments
        ├── _draw_chord(), _draw_chromosome_arcs()
        ├── _pos_to_angle(), _chrom_sort_key(), _check_matplotlib()
        └── Color schemes: SV_COLORS, CNV_COLORS
```

**Total (~6,200 lines, ~200 KB):** 15 Python files implementing 42 public functions, 15 dataclasses, 1 enum (`SVType`) and the `IntervalIndex` class.

---

## 3. Component Responsibilities

### 3.1 Detection Subpackage

**Purpose:** Convert alignment records (dicts built from a BAM) or depth arrays into preliminary SV calls.

#### `cnv.py`

| Function | Purpose | Input | Output |
|----------|---------|-------|--------|
| `detect_cnv_from_depth()` | Entry point for CNV detection | `dict[chrom, depth_array]` (per-bin depths, zeros allowed; each chromosome is divided by the median of its positive bins) | `dict[chrom, CNVResult]` |
| `segment_coverage()` | CBS algorithm | 1D array (log2 ratios) | `list[(start, end, mean)]` |
| `call_cnv_states()` | Assign HOMODEL/DEL/NEUTRAL/DUP/AMP | segments or CNVSegments | `list[CNVSegment]` with `state`, `cn` |
| `calculate_log2_ratio()` | Normalize tumor vs normal depth | tumor_depth, normal_depth, GC | log2 ratio array |
| `merge_adjacent_segments()` | Fuse same-state neighbors | list of CNVSegment | merged list |

**Key Data:**

- **`CNVResult`**: Container per chromosome — `segments`, `log2_ratios`, `bin_size`, `method`, `parameters`.
- **`CNVSegment`**: Single CNV region — `chrom`, `start`, `end`, `mean_log2ratio`, `n_bins`, `state`, `cn`, `confidence` (0–1).

#### `sv_calling.py`

| Function | Purpose | Key Inputs | Key Outputs |
|----------|---------|------------|-------------|
| `call_structural_variants()` | Main SV detection | `alignments: list[dict]`, `min_support`, `insert_size_stats` | `list[StructuralVariant]` |
| `detect_split_reads()` | Soft-clip breakpoint evidence | `reads`, `min_clip=20` | `list[SVEvidence]` |
| `detect_discordant_pairs()` | Insert size/orientation anomalies | `pairs`, `insert_size_stats`, `n_std=4.0` | `list[SVEvidence]` |
| `classify_sv_type()` | Evidence → DEL/DUP/INV/TRA/INS, or UNKNOWN | `SVEvidence` | `SVType` enum |
| `genotype_sv()` | 0/0, 0/1, 1/1 | `StructuralVariant`, `reads` | genotype string |

**Evidence gathering:**
- **Split-read**: CIGAR soft clips (`S` operation) + SA tag (supplementary alignments) → breakpoint pairs. A clipped read without an SA tag yields a single-position breakpoint.
- **Discordant pairs**: Abnormal insert size (> 4σ from mean) OR inter-chromosomal OR same-strand orientation OR reversed (RF) orientation.

**Clustering:** `_cluster_evidence()` sorts the evidence and puts each item into the first existing cluster (newest first) whose first item has the same chromosome pair and both breakpoints within 500 bp; otherwise it starts a new cluster. `_merge_evidence()` computes consensus breakpoints (median), the most common strands, and totals support.

**Genotyping:** Compares supporting reads (split reads + discordant pairs) to reads that span a ±100 bp window around the variant start (`mapq` at least `min_mapq`, not among the evidence reads; a read spans `pos` to `pos + read_length`, with `read_length` defaulting to 150). Allele-fraction thresholds: <0.15 → `0/0`, 0.15–0.85 → `0/1`, ≥0.85 → `1/1`; `./.` when no read counts at all. Every call rescans all reads.

#### `breakpoints.py`

| Function | Purpose |
|----------|---------|
| `refine_breakpoints()` | Mode of the soft-clip positions within `window` of each approximate breakpoint, from the reads you pass in (a read is considered only if its alignment start is within `window + 500` bp of the breakpoint) |
| `detect_microhomology()` | Identify microhomology (NHEJ/MMEJ signature) |
| `cluster_breakpoints()` | Single-linkage nearby breakpoint grouping |
| `calculate_breakpoint_confidence()` | Support count + positional consistency + clip fraction + mean MAPQ |

---

### 3.2 Annotation Subpackage

#### `overlap.py`

| Function | Purpose |
|----------|---------|
| `annotate_gene_overlap()` | Overlap SVs with gene annotations (dicts with `chrom`, `start`, `end`, `name`, or `GenomicInterval` objects) |
| `annotate_regulatory_overlap()` | Overlap with enhancers/promoters/CTCF |
| `find_nearest_gene()` | Nearest gene within `max_distance` |
| `IntervalIndex` | Sorted-start interval lookup (`query_overlap`, `query_nearest`) |

**`IntervalIndex` design:** Build once from a `GenomicInterval` list. Groups by chromosome, sorts by start, builds start-position arrays for binary search. `query_overlap(chrom, start, end)` uses `bisect_left(starts, end)` to find the intervals with `start < end`, then scans all of them and keeps those with `interval.end > start`. The scan is linear in the number of intervals that begin before the query end, so cost grows with the query's position on the chromosome. `annotate_gene_overlap()` builds a new index on every call.

#### `functional_impact.py`

| Function | Purpose |
|----------|---------|
| `predict_functional_impact()` | Unified impact assessment orchestrator |
| `assess_dosage_sensitivity()` | Look up HI/TS/pLI/LOEUF scores |
| `predict_tad_disruption()` | Check SV breakpoints vs TAD boundary midpoints |
| `score_pathogenicity()` | Weighted evidence model → [0,1] |

**Impact classification hierarchy:**
1. **HIGH**: Dosage-sensitive gene hit (HI gene for DEL/HOMODEL, TS gene for DUP) OR gene-disrupting DEL/HOMODEL OR TRA overlapping a gene (fusion candidate).
2. **MODERATE**: DUP, INV or INS involving genes, or TAD disruption.
3. **LOW**: Large (>100 kb) intergenic, or intergenic >10 kb.
4. **MODIFIER**: Small benign-appearing variants.

`predict_functional_impact()` hands only the haploinsufficiency database to `assess_dosage_sensitivity()`, so it never evaluates triplosensitivity (the TS-driven `dosage_gain` type is not produced through it).

**Pathogenicity scoring weights:**
- Coding genes (0.15 each, capped at 1) plus dosage-sensitive genes (0.3 each): 30%
- Size (sigmoid, midpoint 100 kb): 15%
- TAD disruption: 10%
- SV type (HOMODEL 0.9 > TRA/BND 0.8 > DEL 0.7 > AMP 0.6 > DUP 0.5 > INV 0.4 > INS/UNKNOWN 0.3): 15%
- Impact level (HIGH > MODERATE > LOW > MODIFIER): 15%
- Population frequency (rare = higher score, used only when `population_frequency` is supplied to `score_pathogenicity()`; otherwise it scores the maximum): 15%

---

### 3.3 Filtering Subpackage

#### `quality_filter.py`

| Function | Filters | Key Parameters |
|----------|---------|---------------|
| `filter_by_quality()` | Quality score, support count, MAPQ, GQ | `min_qual`, `min_support`, `min_mapq`, `min_gq` |
| `filter_by_size()` | Size min/max (TRA/BND always kept) | `min_size`, `max_size` |
| `filter_by_frequency()` | Population AF (gnomAD, 1000G) | `max_af`, `match_window`, `min_reciprocal_overlap` |
| `apply_blacklist()` | Remove variants overlapping blacklisted regions | `blacklist_regions`, `min_overlap_fraction` |

All return `(filtered_variants, FilterStats)`. The databases (population frequencies, blacklist regions) are passed in as Python objects; the module reads no reference files.

**Frequency matching:** A database entry matches when it lies within `match_window` (±200 bp by default) of the variant and the reciprocal overlap reaches `min_reciprocal_overlap` in both directions. The window is a pre-filter: with the default 0.5, the overlap test decides, and the window only matters when `min_reciprocal_overlap` is 0.

#### `merge.py`

| Function | Method | Description |
|----------|--------|-------------|
| `merge_callsets()` | Reciprocal overlap | Union-find over variant pairs from different callers whose `reciprocal_overlap ≥ min_overlap` (same type when `type_match`). Clusters → consensus position (median), consensus type (majority vote). |
| `survivor_merge()` | Distance-based | SURVIVOR approach: merge if both breakpoints are within `max_distance`; keeps clusters seen by at least `min_callers` callers. Accepts VCF paths or lists of variant dicts. |

**Reciprocal overlap:** `min(overlap_len/variant1_len, overlap_len/variant2_len)` over the variants' `[start, end)` intervals; variants on different chromosomes score 0.

---

### 3.4 Population Subpackage

| Function | Purpose | Returns |
|----------|---------|---------|
| `genotype_sv_population()` | Genotype known SVs across samples from per-sample evidence stored in each call (`samples[sample]` with `depth_ratio`, `split_reads`, `discordant_pairs`) | `genotype_matrix`, `quality_scores`, `allele_frequencies` |
| `sv_allele_frequency()` | AF, MAF, SFS, optionally per-population | Dict with `frequencies`, `maf`, `n_polymorphic`, `sfs`, `population_frequencies` |
| `sv_association_test()` | Linear (continuous) or logistic (binary) regression | `beta`, `se`, `p_value`, `odds_ratio` (binary only), `n_samples` |
| `sv_population_structure()` | PCA from genotype matrix | `pcs`, `eigenvalues`, `variance_explained`, `loadings` |
| `sv_ld_analysis()` | LD (r²) between SV and SNP genotypes | `ld_matrix`, `tag_snps`, `max_r2_per_sv` |
| `merge_sv_callsets()` | Combine calls across samples into unified reference | list of merged dicts (`n_samples`, `samples`, `support_count`, `merged_calls`) |

**Genotyping methods:**
- **`depth`**: Uses the observed/expected depth ratio (a linear ratio, not log2). For DEL: ratio <0.25 → 1/1, <0.75 → 0/1, otherwise 0/0. For DUP: ratio >1.75 → 1/1, >1.25 → 0/1, otherwise 0/0. Other types call 0/1 when the ratio is <0.5.
- **`split`**: Uses the fraction `(split_reads + discordant_pairs) / total_reads` (`total_reads` defaults to 30): ≥0.85 → 1/1, ≥0.15 → 0/1, otherwise 0/0.

---

### 3.5 Visualization Subpackage

| Function | Plot Type | Output |
|----------|-----------|--------|
| `plot_circos()` | Circular genome layout with SV chords | Saved to `output_path` (format from the extension), returns `Figure` |
| `plot_coverage_track()` | Read depth per bin (not log2 ratios) + SV overlay | Same |
| `plot_sv_size_distribution()` | Histogram by SV type | Same |
| `plot_sv_type_summary()` | Bar chart and pie chart of SV type counts | Same |
| `plot_breakpoint_detail()` | Three panels over `[start - flanking, end + flanking]`: read pileup (first 100 reads), coverage, split/discordant evidence | Same |
| `plot_cnv_profile()` | Genome-wide CNV segments | Same |

**Color scheme:** `SV_COLORS`: DEL(#E74C3C red), DUP(#3498DB blue), INV(#2ECC71 green), TRA and BND(#9B59B6 purple), INS(#F39C12 orange), plus UNKNOWN, NEUTRAL, AMP and HOMODEL entries; `CNV_COLORS` covers the five CNV states.

---

## 4. Data Flow

### 4.1 End-to-End: Single-Sample SV Discovery

```mermaid
flowchart TD
    A["Alignment records and depth arrays<br/>(dicts and arrays you build from BAM/CRAM)"] --> B{Detection Path}

    B --> C["CNV Detection<br/>detect_cnv_from_depth()"]
    B --> D["SV Calling<br/>call_structural_variants()"]
    B --> E["Breakpoint Refinement<br/>refine_breakpoints()"]

    C --> F["Collect calls<br/>(as dicts)"]
    D --> F
    E --> F

    F --> G["Annotation<br/>gene/regulatory overlap<br/>functional impact"]
    G --> H["Filtering<br/>quality/size/frequency/blacklist"]
    H --> I["Output<br/>JSON/VCF + Figures"]
```

**Detailed sequence:**

1. **CNV Detection** (optional but recommended):
   - Input: read depth per window (1000 bp default) per chromosome, as positive values
   - Each chromosome is normalized to its own median and converted to log2 ratios; for tumor vs normal, form the ratio with `calculate_log2_ratio()` first and pass `2**ratio` as depth
   - CBS segmentation → segments with constant mean log2
   - State calling (DEL/DUP/NEUTRAL/AMP/HOMODEL) via thresholds
   - Output: `CNVResult` per chromosome with `CNVSegment` list

2. **SV Calling**:
   - Parse CIGAR strings (soft-clip detection) and SA tags
   - Detect discordant pairs (abnormal insert size/orientation/chromosome)
   - Cluster evidence by genomic proximity (≤500 bp)
   - Merge cluster evidence (median breakpoints, total support)
   - Classify SV type from strand patterns and distance
   - Genotype via allele fraction from supporting vs reference reads
   - Output: `list[StructuralVariant]`

3. **Breakpoint Refinement** (optional, post-processing):
   - For each SV, use the reads you pass in (typically those within a 200–500 bp window of each breakpoint)
   - Collect soft-clip positions
   - Mode (most frequent position) as consensus
   - Detect microhomology from the first read sequence on the chromosome
   - Output: `list[BreakpointPair]` with `Breakpoint` objects

4. **Annotation**:
   - Build an `IntervalIndex` from gene/regulatory entries (parse BED/GTF yourself into dicts or `GenomicInterval`)
   - `query_overlap()` per variant → `OverlapResult` list
   - Dosage sensitivity lookup (haploinsufficiency scores through `predict_functional_impact()`; HI/TS/pLI/LOEUF through `assess_dosage_sensitivity()`)
   - TAD boundary proximity check (within 10 kb window)
   - Composite pathogenicity score (6 evidence lines)
   - Output: `FunctionalImpact` per variant

5. **Filtering**:
   - Quality: `quality ≥ 20`, `support ≥ 3` by default; mean MAPQ and genotype quality are optional (`min_mapq`, `min_gq`, off by default)
   - Size: `size ≥ 50 bp` by default, with an optional `max_size` (TRA/BND exempt from size)
   - Frequency: Remove variants with AF above `max_af` (1% by default) in the population database you supply
   - Blacklist: Remove overlaps with the blacklist regions you supply (for example the ENCODE blacklist or segmental duplications)
   - Output: Filtered variant list

6. **Multi-Caller Merge** (if multiple callers):
   - Tag each variant with caller name
   - Reciprocal overlap graph → connected components = consensus clusters (two variants from the same caller are never compared directly, but a call from another caller can bridge them into one cluster)
   - Consensus position = median start/end; consensus type = majority vote
   - Consensus quality increases with `n_callers`
   - Output: `list[MergedVariant]`

7. **Population Analysis** (cohort):
   - Genotype each SV per sample from the per-sample evidence in the calls (depth or split-read method)
   - Build an `n_SVs × n_samples` genotype matrix (0, 1 or 2 alternate alleles)
   - Compute AF, MAF, SFS, per-population stratification
   - Association testing (linear or logistic regression with covariates)
   - PCA for population structure
   - Output: Dict of statistics

8. **Visualization**:
   - Circos: chords connecting breakpoints across the chromosomes you pass in, colored by type
   - Coverage: CNV log2 ratio track with SV rectangles overlay
   - Histograms: size distribution per SV type
   - Output: PNG/PDF files

---

## 5. Key Data Structures

### Core Dataclasses

All major types are `@dataclass` for clarity and type safety:

```python
@dataclass
class CNVSegment:
    chrom: str
    start: int                      # 0-based inclusive
    end: int                        # 0-based exclusive
    mean_log2ratio: float
    n_bins: int
    state: str = "NEUTRAL"          # DEL/DUP/NEUTRAL/AMP/HOMODEL
    cn: int = 2                     # integer copy number
    confidence: float = 0.0         # [0, 1]

@dataclass
class StructuralVariant:
    chrom: str
    start: int
    end: int
    sv_type: SVType                 # enum: DEL, DUP, INV, TRA, INS, BND, UNKNOWN
    size: int = 0
    quality: float = 0.0
    evidence: SVEvidence = field(default_factory=SVEvidence)
    genotype: str = "./."
    filter_status: str = "PASS"
    chrom2: str = ""                # For translocations
    info: dict[str, Any] = field(default_factory=dict)

@dataclass
class Breakpoint:
    chrom: str
    position: int
    strand: str = "+"
    confidence: float = 0.0
    support: int = 0
    microhomology: str = ""
    inserted_sequence: str = ""
    read_names: list[str] = field(default_factory=list)

@dataclass
class FunctionalImpact:
    variant_id: str = ""
    impact_level: str = "MODIFIER"  # HIGH/MODERATE/LOW/MODIFIER
    impact_type: str = ""           # gene_disruption, dosage_loss, tad_disruption, ...
    affected_genes: list[str] = field(default_factory=list)
    dosage_sensitive_genes: list[str] = field(default_factory=list)
    tad_disrupted: bool = False
    pathogenicity_score: float = 0.0
    details: dict[str, Any] = field(default_factory=dict)
```

---

## 6. Integration Points

### 6.1 Within METAINFORANT

| Module | Integration | Usage |
|--------|-------------|-------|
| `core.utils.logging` | Logging | Every module calls `get_logger(__name__)`; `configure_logging_from_env()` reads `CORE_LOG_LEVEL` |
| `core.io` | File I/O (caller side) | The analysis functions do no file I/O (`survivor_merge()` reads VCF paths and the `plot_*` functions write figures); the examples use `io.load_json()`, `io.dump_json()` and `io.read_tsv()` for inputs and results |
| `core.utils.config` | Configuration | Only for loading a settings file of your own (`load_mapping_from_file`); the five `SV_` overrides are read with `os.environ` inside the detection and filtering functions |
| `dna` | VCF input/output | `dna.variation.variants.parse_vcf()` and `write_vcf()` (see [EXAMPLES.md](EXAMPLES.md)) |
| `longread` | Long-read alignments | `longread.io.bam.read_long_read_bam()` returns `LongReadAlignment` objects; convert them to the dict records used here |

The module imports nothing from other METAINFORMANT packages except `core.utils.logging`. Correlating SVs with RNA expression or other omic layers is not implemented here.

### 6.2 External Tools

- **pysam**: Optional, imported only by `filtering.merge` when `survivor_merge()` is given VCF paths (a basic text parser is used when pysam is missing). Reading BAM files and indexing them is left to the caller.
- **SURVIVOR**: Optional external merger; SURVIVOR-style merging provided in pure Python.
- **BEDTools**: Partially reimplemented in `annotation.overlap` via `IntervalIndex`.

---

## 7. Design Patterns & Rationale

### 7.1 Detection: Split-Read + Discordant Pair Hybrid

**Why both?**
- Split-read: Precise breakpoint (single-base), but only for reads crossing SV boundary (sensitivity limited by mappability). Best for SVs ≥50 bp.
- Discordant pairs: Sensitive to larger SVs (insert size deviation), but breakpoint less precise (≈insert size / 2). Best for deletions/duplications >1 kb.
- **Combined**: Complementarity increases recall; evidence clustering yields both type and coordinates.

### 7.2 Circular Binary Segmentation (CBS)

**Why CBS over simple thresholding?**
- Robust to noise: Model-based vs hard log2 ratio cutoffs.
- Statistically grounded: p-values from the asymptotic Gumbel approximation of Olshen et al. (no permutation test is run).
- Adaptive segment boundaries: No arbitrary windowing.
- **References:** Olshen et al. (2004) Biostatistics.

**Algorithm (recursive):**
1. Compute the best single split-point statistic and the best circular arc statistic (arcs of at most half the segment) over the segment
2. Take the more extreme of the two; if its Gumbel-approximated `p ≤ significance`, split
3. Recurse on the resulting pieces until no significant change-point remains (or `max_iterations` depth)
4. Merge adjacent segments with the same CN state within `merge_distance`

### 7.3 Interval Index Using Sorted-Start Binary Search

**Alternative considered:** Interval tree (O(log n + k)). **Chosen:** Sorted-start + linear scan.

**Rationale:**
- Simple implementation, lower memory than a balanced tree.
- Binary search gives the candidates with `start < query_end`; the scan then checks `end > query_start` for every one of them.
- That scan is linear in the number of intervals before the query end, so it suits gene-scale databases (tens of thousands of features) and gets slower for queries late on a chromosome with a very large database (see [PERFORMANCE.md](PERFORMANCE.md)).

### 7.4 Composite Pathogenicity Score

**Why weighted sum vs ML model?**
- Interpretable: Each weight corresponds to evidence line.
- Fixed constants in `score_pathogenicity()`, not fitted to any dataset.
- Extensible: easy to add new evidence (population frequency, constraint, etc.).
- Intended for triage (HIGH/LOW/MODIFIER), not clinical classification.

### 7.5 Union-Find for Multi-Caller Merging

**Why not all-pairs comparison?**  Comparing every pair across 3+ callers is O(n²).

**Graph clustering:** Variants are sorted by chromosome and start; each variant is compared only with the variants that follow it until the next start lies more than 10 kb past its end (`merge_callsets()`) or farther than `max_distance` plus its length (`survivor_merge()`). Matching pairs are joined with union-find (~O(n α(n))), and consensus values derive from the component members.

---

## 8. Algorithm Details

### 8.1 CBS P-value Approximation

The CBS statistic's null distribution under permutation converges to the distribution of the maximum of a standardized Brownian bridge. Olshen et al. (2004) approximates this via a Gumbel distribution:

```
a_n = sqrt(2 * log(n))
b_n = 2 * log(n) + 0.5*log(log(n)) - 0.5*log(π)
z = (T * a_n - b_n)
p = 1 - exp(-exp(-z))   # Gumbel survival function
```

Used by `_cbs_recursive()` (called from `segment_coverage()`) through `_cbs_pvalue()`. Segments with fewer than 4 points always get `p = 1`.

### 8.2 SV Type Classification Logic

From `classify_sv_type()` evidence patterns, tested in this order (for the last two rows the pair is first normalized to genomic order, upstream read first):

| Condition | SV Type |
|-----------|---------|
| `chrom2` set and different from `chrom1` | `TRA` |
| split-read support and breakpoint distance < 50 | `INS` |
| `strand1 == strand2` (both `+` or both `-`) | `INV` |
| upstream `-`, downstream `+` (reads point away from each other) | `DUP` |
| upstream `+`, downstream `-` (reads face each other) with split-read support or distance ≥ 50 | `DEL` |
| Anything else, including weak or missing strand information | `UNKNOWN` |

The classifier never defaults to `DEL` from weak evidence.

### 8.3 Dosage Sensitivity Classification

- **Haploinsufficient**: HI ≥ 0.5 OR pLI ≥ 0.9 OR LOEUF < 0.35.
- **Triplosensitive**: TS ≥ 0.5.

Genes with these properties are flagged when a variant changes copy number (DEL/HOMODEL for HI, DUP for TS). Through `predict_functional_impact()` only the HI database is consulted (see 3.2).

### 8.4 TAD Disruption Detection

A TAD boundary (midpoint) is disrupted if:
- **DEL**: Breakpoint(s) within `boundary_window` (10 kb default) of midpoint OR deletion spans midpoint.
- **INV/TRA/BND**: Either breakpoint within `boundary_window`.
- **DUP**: Tandem duplication spans midpoint (may re-loop boundary).
- **Other types**: The variant extended by `boundary_window` covers the midpoint.

Disruption score = `min(1.0, n_boundaries × 0.4 × type_weight)` where type_weight(DEL)=1.0, TRA=BND=0.95, INV=0.9, DUP=0.7, INS=0.3, other=0.5.

---

## 9. Configuration Architecture

Module parameters are keyword arguments with defaults in each function signature. Five `SV_` environment variables override the matching argument when the function is called; nothing is read from a config file and there is no `structural_variants.*` key namespace.

**Environment overrides:**

| Parameter | Default | Env Var | Used In |
|-----------|---------|---------|---------|
| `min_support` | 3 | `SV_MIN_SUPPORT` | SV calling, `filter_by_quality()` |
| `min_qual` | 20.0 | `SV_MIN_QUAL` | `filter_by_quality()` |
| `significance` | 0.01 | `SV_CNV_SIGNIFICANCE` | `detect_cnv_from_depth()` (CBS segmentation) |
| `del_threshold` | -0.3 | `SV_CNV_DEL_THRESHOLD` | `call_cnv_states()` |
| `dup_threshold` | 0.3 | `SV_CNV_DUP_THRESHOLD` | `call_cnv_states()` |

Every other parameter (for example `min_sv_size`, `max_af`, `min_overlap`, `boundary_window`) is a function argument only.

**Access:**
```python
import os

os.environ["SV_MIN_SUPPORT"] = "5"   # read on each call to call_structural_variants() and filter_by_quality()
```

See [CONFIGURATION.md](CONFIGURATION.md) for the full argument reference.

---

## 10. Extensibility

### Adding a New SV Type

1. Add enum member in `sv_calling.py`: `MYTYPE = "MYTYPE"`
2. Extend `classify_sv_type()` logic
3. Add to `SV_COLORS` in `plots.py`
4. Update `_classify_impact()` in `functional_impact.py` (and the weight tables in `score_pathogenicity()` and `predict_tad_disruption()`)

### Adding a New Filter

```python
def filter_by_annotation(variants, annotation_db):
    # Implementation
    return filtered, FilterStats(...)

# Define it in filtering/quality_filter.py (or a new module that filtering/__init__.py imports).
# Callers import it from the module path:
# from metainformant.structural_variants.filtering.quality_filter import filter_by_annotation
```

### Supporting a New File Format

The module has no I/O layer: callers convert BAM, VCF or BED into the dict records the functions take (see [GETTING_STARTED.md](GETTING_STARTED.md)). The one exception is `filtering.merge.survivor_merge()`, which reads VCF paths itself.

---

## 11. Performance Characteristics

### Time Complexity

| Operation | Complexity | Notes |
|-----------|------------|-------|
| CBS segmentation | **Super-linear (at least O(n²) per recursion level)** | Per chromosome; n = number of bins; pure-Python scans over split points and arcs |
| Split-read detection | **O(m)** | m = total reads; parse CIGAR once |
| Discordant pair detection | **O(m)** | Single pass over reads that carry mate information |
| Evidence clustering | **O(e log e) sort + O(e·c) cluster scan** | e = evidence items, c = clusters |
| Genotyping (`genotype_sv`) | **O(v · m)** | Rescans all reads for each called variant |
| Gene overlap (IntervalIndex) | **O(log g) + linear scan** | g = genes; scan covers every gene starting before the query end |
| Multi-caller merge | **O(v²)** worst case for v variants per chromosome | Sorted sweep with an early break on distance |
| Population PCA | **Dense eigendecomposition** | `numpy.linalg.eigh` on the smaller of the sample × sample and SV × SV covariance matrices |

### Memory Layout

- **Alignments**: an alignment dict with 11 keys and no sequence took about 0.6-0.7 KB (measured with `tracemalloc`), so 100M reads need about 60-70 GB before any sequence strings; process per chromosome or region.
- **CNV**: depth and log2 ratios are float64 arrays (8 bytes per bin); the result stores `log2_ratios` as a Python list. Whole-genome 1000 bp bins would be about 3.1 million bins, about 25 MB per float64 array.
- **Genotype matrix**: `genotype_sv_population()` returns nested Python lists (genotypes and quality scores), which cost far more per entry than a `uint8` array; convert with `numpy.asarray(..., dtype=np.uint8)` before downstream steps.

### Scalability

- **Single sample**: Work in slices (per chromosome or region) because the detection and annotation functions take in-memory lists and have no streaming mode; see [PERFORMANCE.md](PERFORMANCE.md) for measured scaling on synthetic data.
- **Cohort**: Parallelize across samples (embarrassingly parallel). For very large genotype matrices, the dense eigendecomposition in `sv_population_structure()` becomes the heavy step; a randomized SVD (for example `sklearn.decomposition.TruncatedSVD`) is an alternative you would run yourself.

---

## 12. Design Decisions & Trade-offs

### 12.1 Dictionary-Based Read Representation

`call_structural_variants()` accepts `list[dict[str, Any]]` rather than pysam `AlignedSegment` objects.

**Pros:** Testable with synthetic data, no pysam dependency at function level, serializable.
**Cons:** Slight overhead creating dictionaries; consumer must transform BAM→dict.

**Pattern:** The functions take dict records only; transform BAM→dict in your own code at the system boundary (see [GETTING_STARTED.md](GETTING_STARTED.md)).

### 12.2 Environment Variable Override Pattern

Five thresholds (`min_support`, `min_qual`, CNV `significance`, CNV `del_threshold` and `dup_threshold`) can be overridden with `SV_` prefixed environment variables; all other thresholds are plain arguments.

**Rationale:** Deploy-time configuration (Docker/K8s) for the thresholds most often tuned per run, without code changes; consistent with other METAINFORMANT modules.

### 12.3 No External Binary Dependencies

Although tools like SURVIVOR exist as binaries, module implements core algorithms in pure Python.

**Rationale:** Portability, reproducibility, no compilation, easier debug. Can be slower but acceptable for most cohort sizes (<10K samples). External callers can still be wrapped separately.

### 12.4 Dataclass Return Types

Detection, annotation, filtering and merging return simple dataclasses (`StructuralVariant`, `CNVResult`, `FunctionalImpact`, `MergedVariant`, `FilterStats`, ...); the population functions return plain dicts. No complex inheritance.

**Rationale:** Interoperability, transparency, avoids hidden state, easy to document. Dataclasses convert to dicts with `dataclasses.asdict()`.

---

## 13. Testing Philosophy

- **Real implementations**: Tests in `tests/structural_variants/` run the algorithms themselves on synthetic inputs under the repository's real-implementation policy (`docs/REAL_IMPLEMENTATION_POLICY.md`); `monkeypatch` is used only to set the `SV_` environment variables.
- **Synthetic datasets**: Fixtures build depth arrays and alignment dicts with injected deletions, duplications, split reads and discordant pairs of known type and size. No BAM files are used.
- **Fallback parity**: `test_numpy_fallback_parity.py` checks that the numpy and pure-Python code paths agree.
- **Targeted regressions**: classification direction (`test_sv_classification_direction.py`), consensus medians (`test_sv_consensus_median.py`), and depth tests for detection, merging, filtering and functional impact.

---

## 14. Related Documentation

- **[GETTING_STARTED.md](GETTING_STARTED.md)** — Usage examples and quick tutorial
- **[CONFIGURATION.md](CONFIGURATION.md)** — All settings and environment variables
- **[EXAMPLES.md](EXAMPLES.md)** — Real-world analysis scenarios
- **[PERFORMANCE.md](PERFORMANCE.md)** — Optimization and scaling
- **[TROUBLESHOOTING.md](TROUBLESHOOTING.md)** — Error diagnosis
- **[detection.md](detection.md)** — Detection module deep-dive
- **[annotation.md](annotation.md)** — Annotation deep-dive
- **[filtering.md](filtering.md)** — Filtering and merging
- **[population.md](population.md)** — Population genetics workflows
