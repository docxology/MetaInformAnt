# Structural Variants Module Architecture

Comprehensive overview of the structural variant analysis module, including component responsibilities, data flow patterns, integration points, and design rationale.

## Table of Contents

1. [Module Overview](#module-overview)
2. [Directory Structure](#directory-structure)
3. [Component Responsibilities](#component-responsibilities)
4. [Data Flow](#data-flow)
5. [Key Data Structures](#key-data-structures)
6. [Integration Points](#integration-points)
7. [Design Patterns & Rationale](#design-patterns--rationale)
8. [Algorithm Details](#algorithm-details)
9. [Configuration Architecture](#configuration-architecture)
10. [Extensibility](#extensibility)

---

## 1. Module Overview

The `structural_variants` module provides end-to-end analysis of genomic structural variants (SVs) from high-throughput sequencing data. It covers five functional domains:

| Domain | Capabilities | Primary Entry Points |
|--------|-------------|---------------------|
| **Detection** | CNV (CBS), SV (split/discordant reads), breakpoint refinement | `cnv.detect_cnv_from_depth()`, `sv_calling.call_structural_variants()`, `breakpoints.refine_breakpoints()` |
| **Annotation** | Gene/regulatory overlap, functional impact, dosage sensitivity, TAD disruption | `overlap.annotate_gene_overlap()`, `functional_impact.predict_functional_impact()` |
| **Filtering** | Quality/size/frequency filters, blacklist removal, multi-caller merging | `quality_filter.filter_by_quality()`, `merge.merge_callsets()` |
| **Population** | Cohort genotyping, allele freq, association, PCA, LD | `sv_population.genotype_sv_population()`, `sv_allele_frequency()`, `sv_association_test()` |
| **Visualization** | Circos plots, coverage tracks, size histograms | `plots.plot_circos()`, `plots.plot_coverage_track()` |

**SV Types Supported:** `DEL` (deletion), `DUP` (duplication), `INV` (inversion), `TRA` (translocation), `INS` (insertion), `BND` (breakend), `HOMODEL` (homozygous deletion), `AMP` (amplification).

---

## 2. Directory Structure

```
src/metainformant/structural_variants/
├── __init__.py                      # Public API: exports all submodules
├── AGENTS.md                        # Agent directives
│
├── detection/                       # SV detection algorithms
│   ├── __init__.py                  # Exports: cnv, sv_calling, breakpoints
│   ├── cnv.py                       # ~700 lines: CBS, CNV state assignment, log2 ratio
│   │   ├── CNVSegment (dataclass)
│   │   ├── CNVResult (dataclass)
│   │   ├── detect_cnv_from_depth()              # Main entry
│   │   ├── segment_coverage()                   # CBS core algorithm
│   │   ├── _cbs_recursive()                    # Recursive segmentation
│   │   ├── call_cnv_states()                    # State classification
│   │   ├── calculate_log2_ratio()               # With GC correction
│   │   └── merge_adjacent_segments()
│   │
│   ├── sv_calling.py                # ~785 lines: split/discordant evidence
│   │   ├── SVType (enum: DEL, DUP, INV, TRA, INS, BND)
│   │   ├── StructuralVariant (dataclass)
│   │   ├── SVEvidence (dataclass)
│   │   ├── InsertSizeStats (dataclass)
│   │   ├── call_structural_variants()           # Main entry
│   │   ├── detect_split_reads()                 # Soft-clip evidence
│   │   ├── detect_discordant_pairs()            # Insert size/orientation
│   │   ├── classify_sv_type()                   # Evidence → SVType
│   │   ├── genotype_sv()                        # 0/0, 0/1, 1/1
│   │   ├── _estimate_insert_size()
│   │   ├── _parse_cigar_clips(), _parse_sa_tag()
│   │   └── _cluster_evidence(), _merge_evidence()
│   │
│   └── breakpoints.py               # ~517 lines: refinement, microhomology
│       ├── Breakpoint (dataclass)
│       ├── BreakpointPair (dataclass)
│       ├── refine_breakpoints()                # Consensus position
│       ├── detect_microhomology()               # Junction sequence homology
│       ├── cluster_breakpoints()
│       ├── calculate_breakpoint_confidence()
│       └── _collect_clip_positions(), _consensus_position()
│
├── annotation/                      # Variant effect prediction
│   ├── __init__.py
│   ├── overlap.py                   # ~519 lines: interval overlap
│   │   ├── GenomicInterval (dataclass)
│   │   ├── OverlapResult (dataclass)
│   │   ├── IntervalIndex (class: sorted-start binary search)
│   │   ├── annotate_gene_overlap()             # Main gene overlap
│   │   ├── annotate_regulatory_overlap()       # Enhancers/promoters
│   │   ├── calculate_overlap_fraction()
│   │   ├── _compute_overlap_bp()
│   │   └── _classify_relationship()
│   │
│   └── functional_impact.py         # ~572 lines: impact prediction
│       ├── FunctionalImpact (dataclass)
│       ├── DosageSensitivity (dataclass)
│       ├── TADPrediction (dataclass)
│       ├── predict_functional_impact()         # Main entry (all-in-one)
│       ├── assess_dosage_sensitivity()         # HI/TS from DB
│       ├── predict_tad_disruption()            # Boundary crossing
│       ├── score_pathogenicity()               # Composite [0,1] score
│       ├── _find_overlapping_genes()
│       ├── _count_coding_genes()
│       └── _classify_impact()                  # Impact level/type
│
├── filtering/                       # Quality control & merging
│   ├── __init__.py
│   ├── quality_filter.py            # ~481 lines: filters
│   │   ├── FilterStats (dataclass)
│   │   ├── filter_by_quality()                # Qual/support/MAPQ
│   │   ├── filter_by_size()                   # Size range
│   │   ├── filter_by_frequency()              # Population AF
│   │   ├── apply_blacklist()                  # Known bad regions
│   │   └── _load_blacklist_bed()
│   │
│   └── merge.py                     # ~719 lines: multi-caller consensus
│       ├── MergedVariant (dataclass)
│       ├── MergeStats (dataclass)
│       ├── merge_callsets()                   # Reciprocal overlap
│       ├── survivor_merge()                   # Distance-based
│       ├── calculate_reciprocal_overlap()
│       ├── deduplicate_variants()
│       └── _union(), _find() (union-find)
│
├── population/                      # Cohort-scale analysis
│   ├── __init__.py
│   └── sv_population.py             # ~841 lines: genotyping, stats, PCA
│       ├── genotype_sv_population()           # Genotype matrix
│       ├── sv_allele_frequency()              # AF, MAF, SFS, per-pop
│       ├── sv_association_test()              # Linear/logistic regression
│       ├── sv_population_structure()          # PCA
│       ├── sv_ld_analysis()                   # LD with SNPs
│       ├── merge_sv_callsets()                # Cross-sample merging
│       ├── _genotype_by_depth()
│       ├── _genotype_by_split()
│       ├── _linear_regression(), _logistic_regression()
│       └── _pca_eig()
│
└── visualization/                   # Plotting
    ├── __init__.py
    └── plots.py                     # ~925 lines: Circos, coverage, histograms
        ├── plot_circos()                      # Genome-wide chord diagram
        ├── plot_coverage_track()              # CNV track with SV overlay
        ├── plot_size_distribution()
        ├── plot_type_summary()
        ├── _draw_chord(), _draw_chromosome_arcs()
        ├── _pos_to_angle()
        └── Color schemes: SV_COLORS, CNV_COLORS
```

**Total (~4,300 LOC, 430 KB):** 15 Python files implementing 50+ public functions, 20+ dataclasses.

---

## 3. Component Responsibilities

### 3.1 Detection Subpackage

**Purpose:** Convert raw sequencing alignments (BAM) or depth arrays into preliminary SV calls.

#### `cnv.py`

| Function | Purpose | Input | Output |
|----------|---------|-------|--------|
| `detect_cnv_from_depth()` | Entry point for CNV detection | `dict[chrom, depth_array]` | `dict[chrom, CNVResult]` |
| `segment_coverage()` | CBS algorithm | 1D array (log2 ratios) | `list[(start, end, mean)]` |
| `call_cnv_states()` | Assign DEL/DUP/NEUTRAL | segments or CNVSegments | `list[CNVSegment]` with `state`, `cn` |
| `calculate_log2_ratio()` | Normalize tumor vs normal depth | tumor_depth, normal_depth, GC | log2 ratio array |
| `merge_adjacent_segments()` | Fuse same-state neighbors | list of CNVSegment | merged list |

**Key Data:**

- **`CNVResult`**: Container per chromosome — `segments`, `log2_ratios`, `bin_size`, `method`, `parameters`.
- **`CNVSegment`**: Single CNV region — `chrom`, `start`, `end`, `mean_log2ratio`, `state`, `cn`, `confidence` (0–1).

#### `sv_calling.py`

| Function | Purpose | Key Inputs | Key Outputs |
|----------|---------|------------|-------------|
| `call_structural_variants()` | Main SV detection | `alignments: list[dict]`, `min_support`, `insert_stats` | `list[StructuralVariant]` |
| `detect_split_reads()` | Soft-clip breakpoint evidence | `reads`, `min_clip=20` | `list[SVEvidence]` |
| `detect_discordant_pairs()` | Insert size/orientation anomalies | `reads`, `insert_stats`, `n_std=4.0` | `list[SVEvidence]` |
| `classify_sv_type()` | Evidence → DEL/DUP/INV/TRA/INS | `SVEvidence` | `SVType` enum |
| `genotype_sv()` | 0/0, 0/1, 1/1 | `StructuralVariant`, `all_reads` | genotype string |

**Evidence gathering:**
- **Split-read**: CIGAR soft clips (`S` operation) + SA tag (supplementary alignments) → breakpoint pairs.
- **Discordant pairs**: Abnormal insert size (> 4σ from mean) OR inter-chromosomal OR same-strand orientation.

**Clustering:** `_cluster_evidence()` groups evidence within 500 bp using single-linkage. `_merge_evidence()` computes consensus breakpoints (median) and totals support.

**Genotyping:** Compares supporting reads to reference-spanning reads within ±100 bp window. AF thresholds: <0.15 → `0/0`, 0.15–0.85 → `0/1`, ≥0.85 → `1/1`.

#### `breakpoints.py`

| Function | Purpose |
|----------|---------|
| `refine_breakpoints()` | Consensus breakpoint positions from local split-reads |
| `detect_microhomology()` | Identify microhomology (NHEJ/MMEJ signature) |
| `cluster_breakpoints()` | Single-linkage nearby breakpoint grouping |
| `calculate_breakpoint_confidence()` | Support count + positional consistency + clip fraction |

---

### 3.2 Annotation Subpackage

#### `overlap.py`

| Function | Purpose |
|----------|---------|
| `annotate_gene_overlap()` | Overlap SVs with gene annotations (GTF/BED) |
| `annotate_regulatory_overlap()` | Overlap with enhancers/promoters/CTCF |
| `IntervalIndex` | Binary-search interval query (O(log n) per query) |

**`IntervalIndex` design:** Build once from `GenomicInterval` list. Groups by chromosome, sorts by start, builds start-position arrays for binary search. `query_overlap(chrom, start, end)` uses `bisect_left(starts, end)` to get candidate intervals with `start < end`, then filters by `interval.end > start`. Efficient for millions of intervals.

#### `functional_impact.py`

| Function | Purpose |
|----------|---------|
| `predict_functional_impact()` | Unified impact assessment orchestrator |
| `assess_dosage_sensitivity()` | Look up HI/TS/pLI/LOEUF scores |
| `predict_tad_disruption()` | Check SV breakpoints vs TAD boundary midpoints |
| `score_pathogenicity()` | Weighted evidence model → [0,1] |

**Impact classification hierarchy:**
1. **HIGH**: Dosage-sensitive gene hit (HI for DEL, TS for DUP) OR gene-disrupting DEL/HOMODEL OR TRA with gene fusion candidate.
2. **MODERATE**: DUP involving genes, INV disrupting genes, INS in gene, TAD disruption.
3. **LOW**: Large (>100 kb) intergenic, or intergenic >10 kb.
4. **MODIFIER**: Small benign-appearing variants.

**Pathogenicity scoring weights:**
- Gene count + dosage sensitivity: 30%
- Size (sigmoid, midpoint 100 kb): 15%
- TAD disruption: 10%
- SV type (DEL/TRA > INV > DUP/AMP > INS): 15%
- Impact level (HIGH > MODERATE): 15%
- Population frequency (rare = higher score): 15%

---

### 3.3 Filtering Subpackage

#### `quality_filter.py`

| Function | Filters | Key Parameters |
|----------|---------|---------------|
| `filter_by_quality()` | Quality score, support count, MAPQ, GQ | `min_qual`, `min_support`, `min_mapq` |
| `filter_by_size()` | Size min/max | `min_size`, `max_size` |
| `filter_by_frequency()` | Population AF (gnomAD, 1000G) | `max_af`, `match_window`, `min_reciprocal_overlap` |
| `apply_blacklist()` | Lift-over blacklisted regions | `blacklist_regions` |

All return `(filtered_variants, FilterStats)`.

**Frequency matching:** Uses window-based position matching (±200 bp default) plus reciprocal overlap to account for coordinate differences between callers/databases.

#### `merge.py`

| Function | Method | Description |
|----------|--------|-------------|
| `merge_callsets()` | Reciprocal overlap | Graph-based: edges between variants with `reciprocal_overlap ≥ min_overlap`. Union-find to cluster → consensus position (median), consensus type (majority vote). |
| `survivor_merge()` | Distance-based | SURVIVOR approach: merge if both breakpoints within `max_distance` |

**Reciprocal overlap:** `min(overlap_len/variant1_len, overlap_len/variant2_len)`. For DEL/DUP uses variant interval; for TRA uses breakpoint coordinates.

---

### 3.4 Population Subpackage

| Function | Purpose | Returns |
|----------|---------|---------|
| `genotype_sv_population()` | Genotype known SVs across samples | `genotype_matrix`, `quality_scores`, `allele_frequencies` |
| `sv_allele_frequency()` | AF, MAF, SFS, optionally per-population | Dict with `frequencies`, `maf`, `sfs`, `population_frequencies` |
| `sv_association_test()` | Linear (continuous) or logistic (binary) regression | `beta`, `se`, `p_value`, `odds_ratio` |
| `sv_population_structure()` | PCA from genotype matrix | `components`, `explained_variance`, `loadings` |
| `sv_ld_analysis()` | LD (r², D') between SV and SNP genotypes | `r2_matrix`, `d_prime_matrix`, `tag_snps` |
| `merge_sv_callsets()` | Combine calls across samples into unified reference | merged callset |

**Genotyping methods:**
- **`depth`**: Compares read depth ratio to expected (log2 ratio). Thresholds: <0.3 → 0/0, 0.3–1.7 → 0/1, ≥1.7 → 1/1 for DEL; inverse for DUP.
- **`split`**: Uses split-read + discordant pair supporting counts vs. reference-spanning reads.

---

### 3.5 Visualization Subpackage

| Function | Plot Type | Output |
|----------|-----------|--------|
| `plot_circos()` | Circular genome layout with SV chords | PNG/PDF, returns `Figure` |
| `plot_coverage_track()` | CNV profile + SV overlay | PNG/PDF |
| `plot_size_distribution()` | Histogram by SV type | PNG/PDF |
| `plot_type_summary()` | Bar chart of SV type counts | PNG/PDF |

**Color scheme:** `SV_COLORS`: DEL(#E74C3C red), DUP(#3498DB blue), INV(#2ECC71 green), TRA(#9B59B6 purple), INS(#F39C12 orange).

---

## 4. Data Flow

### 4.1 End-to-End: Single-Sample SV Discovery

```mermaid
flowchart TD
    A[BAM/CRAM<br/>alignments] --> B{Detection Path}
    
    B --> C[CNV Detection<br/>detect_cnv_from_depth()]
    B --> D[SV Calling<br/>call_structural_variants()]
    B --> E[Breakpoint Refinement<br/>refine_breakpoints()]
    
    C --> F[Aggregate<br/>all SV calls]
    D --> F
    E --> F
    
    F --> G[Annotation<br/>gene/regulatory overlap<br/>functional impact]
    G --> H[Filtering<br/>quality/size/frequency/blacklist]
    H --> I[Output<br/>JSON/VCF + Figures]
```

**Detailed sequence:**

1. **CNV Detection** (optional but recommended):
   - Input: read depth per window (1000 bp default) per chromosome
   - Compute log2 ratio (tumor vs normal, or self-normalized)
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
   - For each SV, fetch local reads (200–500 bp window)
   - Collect soft-clip positions
   - Mode (most frequent position) as consensus
   - Detect microhomology from read sequence around junction
   - Output: `list[BreakpointPair]` with `Breakpoint` objects

4. **Annotation**:
   - Build `IntervalIndex` from gene/regulatory BED/GTF
   - `query_overlap()` per variant → `OverlapResult` list
   - Dosage sensitivity lookup (ClinGen HI/TS, gnomAD pLI/LOEUF)
   - TAD boundary proximity check (within 10 kb window)
   - Composite pathogenicity score (6 evidence lines)
   - Output: `FunctionalImpact` per variant

5. **Filtering**:
   - Quality: `quality ≥ 20`, `support ≥ 3`, `MAPQ ≥ 30`
   - Size: `50 bp ≤ size ≤ 50 Mb` (TRA/BND exempt from size)
   - Frequency: Remove variants with AF > 1% in gnomAD/1000G
   - Blacklist: Remove overlaps with ENCODE DAC/segmental duplications
   - Output: Filtered variant list

6. **Multi-Caller Merge** (if multiple callers):
   - Tag each variant with caller name
   - Reciprocal overlap graph → connected components = consensus clusters
   - Consensus position = median start/end; consensus type = majority vote
   - Consensus quality increases with `n_callers`
   - Output: `list[MergedVariant]`

7. **Population Analysis** (cohort):
   - Genotype each SV per sample (depth or split-read method)
   - Build `n_SVs × n_samples` genotype matrix (0/0, 0/1, 1/1)
   - Compute AF, MAF, SFS, per-population stratification
   - Association testing (linear or logistic regression with covariates)
   - PCA for population structure
   - Output: Dict of statistics

8. **Visualization**:
   - Circos: chords connecting breakpoints across 24 chromosomes, colored by type
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
    sv_type: SVType                 # enum: DEL, DUP, INV, TRA, ...
    size: int = 0
    quality: float = 0.0
    evidence: SVEvidence
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
    variant_id: str
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

### 6.1 Within METAINFORMANT

| Module | Integration | Usage |
|--------|-------------|-------|
| `core.io` | File I/O | `io.dump_json()`, `io.load_json()`, `io.load_bed()` — always use, never `open()` |
| `core.config` | Configuration | Environment variable overrides via `SV_` prefix |
| `core.logging` | Logging | `get_logger(__name__)`; structured, level-controlled |
| `dna` | Gene annotations, VCF | `dna.io.load_gtf()`, `dna.variants.parse_vcf()` |
| `visualization` | Shared plotting utils | Color palettes, figure management |
| `multiomics` | Integration across omic layers | Future: correlate SVs with RNA expression |

### 6.2 External Tools

- **Samtools** (samtools): Implicit via `pysam`; BAM indexing assumed.
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
- Statistically grounded: Permutation-derived p-values via Gumbel approximation.
- Adaptive segment boundaries: No arbitrary windowing.
- **References:** Olshen et al. (2004) Biostatistics.

**Algorithm (recursive):**
1. Compute circular binary segmentation statistic across all arc segments
2. Find most significant change-point; if `p ≤ significance`, split
3. Recurse on left and right segments until no more change-points
4. Merge adjacent segments with same CN state within `merge_distance`

### 7.3 Interval Index Using Sorted-Start Binary Search

**Alternative considered:** Interval tree (O(log n + k)). **Chosen:** Sorted-start + linear scan of candidates.

**Rationale:**
- Gene databases typically millions of features; chromosome-wise queries are fast with binary search to get candidates with `start < query_end`.
- Then linear scan over small candidate set (typically <100 for 10 Mb query) checking `end > query_start`.
- Simpler implementation, lower memory than balanced tree.

### 7.4 Composite Pathogenicity Score

**Why weighted sum vs ML model?**
- Interpretable: Each weight corresponds to evidence line.
- Calibrated againstClinGen variant interpretations.
- Extensible: easy to add new evidence (population frequency, constraint, etc.).
- Sufficient accuracy for triage (HIGH/LOW/MODIFIER).

### 7.5 Union-Find for Multi-Caller Merging

**Why not all-pairs comparison?**  O(n²) too slow for thousands of calls across 3+ callers.

**Graph clustering:** Build sparse graph where edges = reciprocal overlap ≥ threshold. Union-find merges connected components efficiently (~O(n α(n))). Consensus derives from component member statistics.

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

Used in `segment_coverage()` via `_cbs_pvalue()`.

### 8.2 SV Type Classification Logic

From `classify_sv_type()` evidence patterns:

| Condition | SV Type |
|-----------|---------|
| `chrom1 != chrom2` | `TRA` |
| distance < 50 and split_reads > 0 | `INS` |
| `strand1 == strand2` | `INV` |
| `bp1 ≤ bp2` and `strand1="+"` and `strand2="-"' | `DEL` (if distance large) |
| `bp1 ≤ bp2` and `strand1="-"` and `strand2="+"` | `DUP` |
| Opposite orientation patterns depending on bp order | `DEL`/`DUP` |
| distance > 1000 with no other info | `DEL` (default) |

### 8.3 Dosage Sensitivity Classification

- **Haploinsufficient**: HI ≥ 0.5 OR pLI ≥ 0.9 OR LOEUF < 0.35 (ClinGen/gnomAD thresholds).
- **Triplosensitive**: TS ≥ 0.5.

Genes with these properties are flagged when CNV changes copy number (DEL for HI, DUP for TS).

### 8.4 TAD Disruption Detection

A TAD boundary (midpoint) is disrupted if:
- **DEL**: Breakpoint(s) within `boundary_window` (10 kb default) of midpoint OR deletion spans midpoint.
- **INV/TRA/BND**: Either breakpoint within `boundary_window`.
- **DUP**: Tandem duplication spans midpoint (may re-loop boundary).

Disruption score = `min(1.0, n_boundaries × 0.4 × type_weight)` where type_weight(DEL)=1.0, TRA=0.95, INV=0.9, DUP=0.7, INS=0.3.

---

## 9. Configuration Architecture

All module-specific parameters read from `metainformant.config`, with `SV_` environment variable override support. Config keys under `structural_variants.<parameter>`.

**Key parameters:**

| Parameter | Default | Env Var | Used In |
|-----------|---------|---------|---------|
| `min_support` | 3 | `SV_MIN_SUPPORT` | SV calling, filtering |
| `min_quality` | 20.0 | `SV_MIN_QUAL` | Filtering |
| `min_sv_size` | 50 | `SV_MIN_SV_SIZE` | Filtering |
| `cnv_significance` | 0.01 | `SV_CNV_SIGNIFICANCE` | CBS segmentation |
| `max_af` | 0.01 | `SV_MAX_AF` | Frequency filter |
| `min_overlap` | 0.5 | `SV_MIN_OVERLAP` | Multi-caller merge |
| `tad_boundary_window` | 10_000 | `SV_TAD_WINDOW` | TAD disruption |

**Access:**
```python
from metainformant import config
val = config.get("structural_variants.min_support")
config.set("structural_variants.min_support", 5)
```

---

## 10. Extensibility

### Adding a New SV Type

1. Add enum member in `sv_calling.py`: `MYTYPE = "MYTYPE"`
2. Extend `classify_sv_type()` logic
3. Add to `SV_COLORS` in `plots.py`
4. Update `_classify_impact()` in `functional_impact.py`

### Adding a New Filter

```python
def filter_by_annotation(variants, annotation_db):
    # Implementation
    return filtered, FilterStats(...)

# Register in quality_filter.__init__.py __all__
# Export via filtering/__init__.py
```

### Supporting a New File Format

Add reader to `io/` (not yet formalized; current pattern: extend `core.io` for format-specific additions).

---

## 11. Performance Characteristics

### Time Complexity

| Operation | Complexity | Notes |
|-----------|------------|-------|
| CBS segmentation | **O(n)** | Per chromosome; n = number of bins |
| Split-read detection | **O(m)** | m = total reads; parse CIGAR once |
| Discordant pair detection | **O(m)** | Single pass over properly paired reads |
| Evidence clustering | **O(e log e)** | e = evidence items; sort then single-linkage |
| Gene overlap (IntervalIndex) | **O(log g + k)** | g = genes, k = overlap candidates |
| Multi-caller merge | **O(v²)** worst-case for v variants per chromosome | Early break on sorted chromosome + distance check; practical < O(v log v) |
| Population PCA | **O(n × m²)** typically use randomized SVD | n = samples, m = SVs (dimensionality) |

### Memory Layout

- **Alignments**: ~200 bytes/read. 100M reads ≈ 20 GB; streaming via `pysam.fetch()` recommended.
- **CNV**: log2 ratios: 8 bytes × n_bins. Whole-genome 1000 bp bins: ~250 KB × 8 = 2 MB.
- **Genotype matrix**: `n_svs × n_samples × 1 byte` (if stored uint8). 100K SVs × 10K samples = 1 GB.

### Scalability

- **Single sample**: 30× whole-genome BAM (~100 GB, 700M reads) processes in ~1 h (single thread) for CNV+SV calling.
- **Cohort**: Parallelize across samples (embarrassingly parallel). Population analysis bottleneck is PCA matrix factorization for >100K SVs; use randomized SVD (`sklearn.decomposition.TruncatedSVD`).

---

## 12. Design Decisions & Trade-offs

### 12.1 Dictionary-Based Read Representation

`call_structural_variants()` accepts `list[dict[str, Any]]` rather than pysam `AlignedSegment` objects directly.

**Pros:** Testable with synthetic data, no pysam dependency at function level, serializable.
**Cons:** Slight overhead creating dictionaries; consumer must transform BAM→dict.

**Pattern:** Public API supports both; transform BAM→dict at system boundary (I/O layer).

### 12.2 Environment Variable Override Pattern

All detection/annotation/filtering thresholds support `SV_` prefixed environment variables.

**Rationale:** Deploy-time configuration (Docker/K8s) without code changes; consistent with other METAINFORMANT modules.

### 12.3 No External Binary Dependencies

Although tools like SURVIVOR exist as binaries, module implements core algorithms in pure Python.

**Rationale:** Portability, reproducibility, no compilation, easier debug. Can be slower but acceptable for most cohort sizes (<10K samples). External callers can still be wrapped separately.

### 12.4 Dataclass-Only Return Types

No complex inheritance; simple dataclasses with default constructors. JSON-serializable via `asdict()`.

**Rationale:** Interoperability, transparency, avoids hidden state, easy to document.

---

## 13. Testing Philosophy

- **No mocking**: Tests use small real BAM files with known SV ground truth.
- **Golden outputs**: CBS segmentation results against R `DNAcopy` package.
- **Synthetic datasets**: `tests/fixtures/` contains simulated BAMs with injected SVs of known types/sizes.
- **End-to-end**: Pipeline integration tests cover full detection→annotation→filtering.

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
