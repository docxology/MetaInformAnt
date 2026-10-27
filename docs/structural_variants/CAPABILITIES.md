# Structural_variants Module Capabilities

## Quick Reference

### Functions by Submodule

| Function | Submodule | Purpose |
|----------|-----------|---------|
| `_classify_impact()` | `annotation.functional_impact` | Classify impact level and type. |
| `_count_coding_genes()` | `annotation.functional_impact` | Count the number of protein-coding genes in a list. |
| `_find_overlapping_genes()` | `annotation.functional_impact` | Find gene names overlapping a region. |
| `assess_dosage_sensitivity()` | `annotation.functional_impact` | Assess dosage sensitivity of a gene. |
| `predict_functional_impact()` | `annotation.functional_impact` | Predict the functional impact of a structural variant. |
| `predict_tad_disruption()` | `annotation.functional_impact` | Predict whether a structural variant disrupts TAD boundaries. |
| `score_pathogenicity()` | `annotation.functional_impact` | Compute aggregate pathogenicity score for a structural variant. |
| `_classify_relationship()` | `annotation.overlap` | Classify the spatial relationship between a variant and a feature. |
| `_compute_overlap_bp()` | `annotation.overlap` | Compute the number of overlapping base pairs between two intervals. |
| `_normalize_intervals()` | `annotation.overlap` | Normalize a list of features to GenomicInterval objects. |
| `annotate_gene_overlap()` | `annotation.overlap` | Annotate structural variants with overlapping genes. |
| `annotate_regulatory_overlap()` | `annotation.overlap` | Annotate structural variants with overlapping regulatory elements. |
| `calculate_overlap_fraction()` | `annotation.overlap` | Calculate reciprocal overlap fractions between two intervals. |
| `find_nearest_gene()` | `annotation.overlap` | Find the nearest gene to a structural variant. |
| `_cigar_ops()` | `detection.breakpoints` | Split a CIGAR into ``(operation letter, length)`` pairs. |
| `_consensus_position()` | `detection.breakpoints` | Find consensus breakpoint position from clip positions. |
| `_first_present()` | `detection.breakpoints` | Return the first value among ``keys`` that is not None (0 when none is). |
| `_inserted_between_alignments()` | `detection.breakpoints` | Bases of ``read`` between its primary and its single supplementary alignment. |
| `_junction_microhomology()` | `detection.breakpoints` | Longest microhomology among the supporting reads, each centred on its junction. |
| `_majority_strand()` | `detection.breakpoints` | Orientation of most supporting reads: '-' when more are reverse-complemented than not. |
| `_parse_clips_from_cigar()` | `detection.breakpoints` | Parse CIGAR to find soft-clip junctions. |
| `_parse_sa_alignments()` | `detection.breakpoints` | Parse an SA tag (``rname,pos,strand,CIGAR,mapQ,NM;...``) into ``(chrom, start, strand, CIGAR)``. |
| `_position_std()` | `detection.breakpoints` | Calculate standard deviation of clip positions. |
| `_query_extent()` | `detection.breakpoints` | Return ``(leading clip, aligned query length, trailing clip, reference length)`` of a CIGAR. |
| `_refine_breakpoint()` | `detection.breakpoints` | Build the refined breakpoint from the soft clips found near it. |
| `_split_read_insertion()` | `detection.breakpoints` | Most common inserted sequence among the supporting reads, or '' when not derivable. |
| `calculate_breakpoint_confidence()` | `detection.breakpoints` | Calculate confidence score for a breakpoint. |
| `cluster_breakpoints()` | `detection.breakpoints` | Cluster nearby breakpoints using single-linkage clustering. |
| `detect_microhomology()` | `detection.breakpoints` | Detect microhomology at a breakpoint junction. |
| `refine_breakpoints()` | `detection.breakpoints` | Refine breakpoint positions to base-pair resolution. |
| `_assign_cnv_state()` | `detection.cnv` | Assign a copy number state and integer CN from a log2 ratio. |
| `_calculate_segment_confidence()` | `detection.cnv` | Calculate confidence score for a CNV segment. |
| `_cbs_best_arc()` | `detection.cnv` | Find the arc whose mean differs most from the rest of the segment. |
| `_cbs_pvalue()` | `detection.cnv` | Approximate p-value for CBS test statistic. |
| `_cbs_recursive()` | `detection.cnv` | Recursive CBS implementation. |
| `_center_on_median()` | `detection.cnv` | Subtract the median of the finite values (pure-Python median centering). |
| `_check_depths()` | `detection.cnv` | Raise a ValueError naming the problem if the values cannot be read-depth counts. |
| `_gc_correct()` | `detection.cnv` | GC-content bias correction using binned median normalization. |
| `_gc_correct_py()` | `detection.cnv` | Pure-Python twin of `_gc_correct()` (binned median normalization). |
| `calculate_log2_ratio()` | `detection.cnv` | Compute log2 ratio between tumor and normal read depths. |
| `call_cnv_states()` | `detection.cnv` | Assign copy number states to segments based on log2 ratio thresholds. |
| `detect_cnv_from_depth()` | `detection.cnv` | Detect copy number variations from read depth data. |
| `merge_adjacent_segments()` | `detection.cnv` | Merge adjacent segments with the same copy number state. |
| `segment_coverage()` | `detection.cnv` | Circular Binary Segmentation (CBS) algorithm for coverage segmentation. |
| `_calculate_sv_quality()` | `detection.sv_calling` | Calculate quality score for a structural variant call. |
| `_cluster_evidence()` | `detection.sv_calling` | Cluster SV evidence by genomic proximity. |
| `_estimate_insert_size()` | `detection.sv_calling` | Estimate insert size distribution from properly paired reads. |
| `_merge_evidence()` | `detection.sv_calling` | Merge a cluster of evidence into a single SVEvidence. |
| `_parse_cigar_clips()` | `detection.sv_calling` | Parse CIGAR string to find soft-clipped regions. |
| `_parse_cigar_string()` | `detection.sv_calling` | Parse a CIGAR string into (operation, length) tuples. |
| `_parse_sa_tag()` | `detection.sv_calling` | Parse SAM SA (supplementary alignment) tag. |
| `call_structural_variants()` | `detection.sv_calling` | Call structural variants from alignment data. |
| `classify_sv_type()` | `detection.sv_calling` | Classify the type of structural variant based on evidence patterns. |
| `detect_discordant_pairs()` | `detection.sv_calling` | Detect discordant read pair evidence for structural variants. |
| `detect_split_reads()` | `detection.sv_calling` | Identify split-read evidence for structural variants. |
| `genotype_sv()` | `detection.sv_calling` | Genotype a structural variant based on supporting and reference reads. |
| `_call_quality()` | `filtering.merge` | Quality of one call (``quality``, else ``qual``); 0 when absent or None. |
| `_consensus_genotype()` | `filtering.merge` | Determine consensus genotype from multiple calls. |
| `_find()` | `filtering.merge` | Find operation for union-find with path compression. |
| `_format_genotype()` | `filtering.merge` | Format a called allele tuple (as pysam reports ``GT``) as ``'0/1'``; None when absent. |
| `_load_vcf_variants()` | `filtering.merge` | Load structural variants from a VCF file. |
| `_majority_sv_type()` | `filtering.merge` | Consensus SV type by majority vote within a merge group. |
| `_median_position()` | `filtering.merge` | Median of a list of positions (consensus breakpoint). |
| `_parse_vcf_basic()` | `filtering.merge` | Basic VCF parser without pysam dependency. |
| `_pysam_first_genotype()` | `filtering.merge` | Genotype of the first sample of a pysam record, or None when it has no ``GT``. |
| `_union()` | `filtering.merge` | Union operation for union-find with path compression. |
| `_vcf_sample_genotype()` | `filtering.merge` | Genotype of the first sample of a VCF record line (``'0/1'``), or None when there is none. |
| `calculate_reciprocal_overlap()` | `filtering.merge` | Calculate the minimum reciprocal overlap between two SVs. |
| `deduplicate_variants()` | `filtering.merge` | Remove duplicate structural variant calls. |
| `merge_callsets()` | `filtering.merge` | Merge structural variant callsets from multiple callers. |
| `survivor_merge()` | `filtering.merge` | SURVIVOR-style merging of structural variant callsets. |
| `_as_frequency()` | `filtering.quality_filter` | Return ``value`` as an allele frequency in [0, 1], or None when it is not one. |
| `_as_position()` | `filtering.quality_filter` | Return ``value`` as an int coordinate, or None when it is not an integral number. |
| `_build_population_lookup()` | `filtering.quality_filter` | Build a chromosome-keyed lookup from population database. |
| `_lookup_population_frequency()` | `filtering.quality_filter` | Look up population frequency for a variant. |
| `_parse_population_db()` | `filtering.quality_filter` | Read a population database, counting the entries that are not in a supported layout. |
| `apply_blacklist()` | `filtering.quality_filter` | Remove variants overlapping blacklist regions. |
| `filter_by_frequency()` | `filtering.quality_filter` | Filter structural variants by population allele frequency. |
| `filter_by_quality()` | `filtering.quality_filter` | Filter structural variants by quality metrics. |
| `filter_by_size()` | `filtering.quality_filter` | Filter structural variants by size. |
| `_build_consensus()` | `population.sv_population` | Build a consensus SV call from a cluster of merged calls. |
| `_genotype_by_depth()` | `population.sv_population` | Genotype an SV based on read depth ratio. |
| `_genotype_by_split()` | `population.sv_population` | Genotype an SV based on split-read and discordant-pair evidence. |
| `_linear_regression()` | `population.sv_population` | Ordinary least squares linear regression. |
| `_logistic_regression()` | `population.sv_population` | Iteratively reweighted least squares logistic regression. |
| `_normal_cdf_np()` | `population.sv_population` | Standard normal CDF approximation. |
| `_should_merge()` | `population.sv_population` | Determine if two SV calls should be merged. |
| `_total_variance()` | `population.sv_population` | Total variance: the sum of all eigenvalues of a covariance matrix. |
| `genotype_sv_population()` | `population.sv_population` | Population-scale SV genotyping across multiple samples. |
| `merge_sv_callsets()` | `population.sv_population` | Merge SV calls across samples using reciprocal overlap and breakpoint proximity. |
| `sv_allele_frequency()` | `population.sv_population` | Compute SV allele frequencies, optionally stratified by population. |
| `sv_association_test()` | `population.sv_population` | Test association between SV genotype and phenotype. |
| `sv_ld_analysis()` | `population.sv_population` | Compute LD between SVs and nearby SNPs. |
| `sv_population_structure()` | `population.sv_population` | PCA on SV genotype matrix for population structure analysis. |
| `_check_matplotlib()` | `visualization.plots` | Raise an import error if matplotlib is not available. |
| `_chrom_sort_key()` | `visualization.plots` | Sort key for chromosome names (numeric first, then alphabetic). |
| `_draw_chord()` | `visualization.plots` | Draw a chord between two angles in polar coordinates. |
| `_draw_chromosome_arcs()` | `visualization.plots` | Draw chromosome arcs around the Circos plot. |
| `_pos_to_angle()` | `visualization.plots` | Convert a genomic position to a polar angle. |
| `plot_breakpoint_detail()` | `visualization.plots` | Plot detailed view of a structural variant breakpoint. |
| `plot_circos()` | `visualization.plots` | Generate a Circos-style genome-wide structural variant plot. |
| `plot_cnv_profile()` | `visualization.plots` | Plot genome-wide CNV profile. |
| `plot_coverage_track()` | `visualization.plots` | Plot read depth coverage with structural variant overlay. |
| `plot_sv_size_distribution()` | `visualization.plots` | Plot size distribution of structural variants. |
| `plot_sv_type_summary()` | `visualization.plots` | Plot summary of SV types as bar chart and pie chart. |

---

## annotation.functional_impact

### Functions

#### `predict_functional_impact()`

**Signature**: `predict_functional_impact(variant, gene_annotations, haploinsufficiency_db, tad_boundaries, triplosensitivity_db, pli_db, loeuf_db)`

Predict the functional impact of a structural variant.

Integrates multiple lines of evidence to classify the functional
consequence of an SV:

1. Gene overlap analysis (disrupted vs. duplicated genes)
2. Dosage sensitivity of affected genes
3. TAD boundary disruption
4. Variant type-specific considerations (DEL vs DUP vs INV)
5. Aggregate pathogenicity scoring

Args:
    variant: Variant dictionary containing:
        - 'chrom': Chromosome
        - 'start': Start position
        - 'end': End position
        - 'sv_type': Type of SV (DEL, DUP, INV, TRA, INS)
        - 'id': Optional variant identifier
        - 'overlapping_genes': Optional pre-computed gene overlaps
    gene_annotations: List of gene annotation dictionaries with:
        - 'chrom': Chromosome
        - 'start': Gene start
        - 'end': Gene end
        - 'name': Gene name
        - 'is_coding': Whether gene is protein-coding (optional)
    haploinsufficiency_db: Optional dictionary mapping gene names to
        haploinsufficiency scores (0-1). Genes it marks haploinsufficient
        count as dosage-sensitive for deletions.
    tad_boundaries: Optional list of TAD boundary dictionaries with
        'chrom', 'start', 'end' keys.
    triplosensitivity_db: Optional dictionary mapping gene names to
        triplosensitivity scores (0-1). Genes it marks triplosensitive
        count as dosage-sensitive for duplications (impact type
        ``dosage_gain``); without it a duplication can never reach that
        class.
    pli_db: Optional dictionary mapping gene names to pLI scores (0-1);
        a gene with pLI >= 0.9 is haploinsufficient.
    loeuf_db: Optional dictionary mapping gene names to LOEUF values;
        a gene with LOEUF < 0.35 is haploinsufficient.

Returns:
    FunctionalImpact object with comprehensive impact assessment.

#### `assess_dosage_sensitivity()`

**Signature**: `assess_dosage_sensitivity(gene, haploinsufficiency_db, triplosensitivity_db, pli_db, loeuf_db)`

Assess dosage sensitivity of a gene.

Looks up the gene in multiple dosage sensitivity databases and
produces a composite assessment. A gene is classified as
haploinsufficient if its HI score exceeds 0.5 or pLI exceeds 0.9.

Args:
    gene: Gene name to look up.
    haploinsufficiency_db: Dictionary mapping gene names to HI scores (0-1).
    triplosensitivity_db: Dictionary mapping gene names to TS scores (0-1).
    pli_db: Dictionary mapping gene names to pLI scores (0-1).
    loeuf_db: Dictionary mapping gene names to LOEUF values.

Returns:
    DosageSensitivity object with all available information.

#### `predict_tad_disruption()`

**Signature**: `predict_tad_disruption(variant, tad_boundaries, boundary_window)`

Predict whether a structural variant disrupts TAD boundaries.

TAD (Topologically Associating Domain) boundaries are important
chromatin organization elements. Their disruption can lead to
ectopic enhancer-promoter contacts and gene misregulation.

A boundary is considered disrupted if the SV breakpoint falls
within boundary_window base pairs of the boundary midpoint.

Args:
    variant: Variant dictionary with 'chrom', 'start', 'end', 'sv_type'.
    tad_boundaries: List of TAD boundary dictionaries with:
        - 'chrom': Chromosome
        - 'start': Boundary start position
        - 'end': Boundary end position
        - 'insulation_score': Optional insulation score
        - 'genes': Optional list of genes in adjacent TADs
    boundary_window: Window size around each boundary for disruption
        detection (default 10kb).

Returns:
    TADPrediction object with disruption assessment.

#### `score_pathogenicity()`

**Signature**: `score_pathogenicity(variant, annotations)`

Compute aggregate pathogenicity score for a structural variant.

Combines multiple evidence lines into a single [0, 1] score using
a weighted logistic model:

- Size of the variant (larger = more likely pathogenic)
- Number of affected coding genes
- Dosage sensitivity of affected genes
- TAD boundary disruption
- Variant type (DEL/DUP more pathogenic than INV on average)
- Population frequency (if available)

Args:
    variant: Variant dictionary with 'chrom', 'start', 'end', 'sv_type'.
    annotations: Annotation dictionary containing:
        - 'overlapping_genes': List of overlapping gene names
        - 'dosage_sensitive_genes': List of dosage-sensitive genes
        - 'tad_disrupted': Whether TAD boundary is disrupted
        - 'impact_level': Impact level string
        - 'n_coding_genes': Number of affected coding genes
        - 'population_frequency': Allele frequency in population (optional)

Returns:
    Pathogenicity score between 0 and 1.

#### `_find_overlapping_genes()`

**Signature**: `_find_overlapping_genes(chrom, start, end, gene_annotations)`

Find gene names overlapping a region.

Args:
    chrom: Chromosome.
    start: Start position.
    end: End position.
    gene_annotations: List of gene annotation dicts.

Returns:
    List of overlapping gene names.

#### `_count_coding_genes()`

**Signature**: `_count_coding_genes(gene_names, gene_annotations)`

Count the number of protein-coding genes in a list.

Args:
    gene_names: List of gene names.
    gene_annotations: Gene annotation database.

Returns:
    Number of protein-coding genes.

#### `_classify_impact()`

**Signature**: `_classify_impact(sv_type, size, n_genes, n_dosage_genes, tad_disrupted, gene_annotations, overlapping_genes)`

Classify impact level and type.

Args:
    sv_type: Type of structural variant.
    size: Size in base pairs.
    n_genes: Number of affected genes.
    n_dosage_genes: Number of dosage-sensitive genes.
    tad_disrupted: Whether TAD boundary is disrupted.
    gene_annotations: Gene annotations for context.
    overlapping_genes: List of overlapping gene names.

Returns:
    Tuple of (impact_level, impact_type).

### Classes

#### `FunctionalImpact`

Functional impact assessment for a structural variant.

Attributes:
    variant_id: Identifier for the variant.
    impact_level: Impact level ('HIGH', 'MODERATE', 'LOW', 'MODIFIER').
    impact_type: Specific impact type (e.g., 'gene_disruption', 'gene_fusion',
        'dosage_change', 'regulatory_disruption', 'tad_disruption').
    affected_genes: List of affected gene names.
    dosage_sensitive_genes: Genes with dosage sensitivity concern.
    tad_disrupted: Whether a TAD boundary is disrupted.
    pathogenicity_score: Aggregate pathogenicity score [0, 1].
    details: Additional details about the impact.

#### `DosageSensitivity`

Dosage sensitivity information for a gene.

Attributes:
    gene_name: Gene name.
    haploinsufficiency_score: Haploinsufficiency score [0, 1].
        Higher means more likely to be haploinsufficient.
    triplosensitivity_score: Triplosensitivity score [0, 1].
        Higher means more sensitive to copy gain.
    pli_score: Probability of loss-of-function intolerance (pLI) [0, 1].
    loeuf: Loss-of-function observed/expected upper fraction (LOEUF).
        Lower values indicate more constraint.
    is_haploinsufficient: Whether gene is classified as haploinsufficient.
    is_triplosensitive: Whether gene is classified as triplosensitive.

#### `TADPrediction`

TAD boundary disruption prediction.

Attributes:
    variant_chrom: Chromosome.
    variant_start: Variant start position.
    variant_end: Variant end position.
    disrupted_boundaries: List of disrupted TAD boundary positions.
    n_boundaries_disrupted: Number of TAD boundaries disrupted.
    genes_in_affected_tads: Genes within affected TAD domains.
    disruption_score: TAD disruption severity score [0, 1].

## annotation.overlap

### Functions

#### `annotate_gene_overlap()`

**Signature**: `annotate_gene_overlap(variants, gene_db)`

Annotate structural variants with overlapping genes.

For each variant, finds all genes that overlap with the variant
interval and adds gene annotation information.

Args:
    variants: List of variant dictionaries, each containing:
        - 'chrom': Chromosome
        - 'start': Start position
        - 'end': End position
        - Additional fields are preserved in the output.
    gene_db: List of gene entries as dictionaries with 'chrom', 'start',
        'end', 'name' keys, or as GenomicInterval objects.

Returns:
    List of variant dictionaries with added fields:
        - 'overlapping_genes': List of gene names that overlap
        - 'gene_overlaps': List of OverlapResult objects
        - 'n_genes_affected': Number of genes affected

#### `annotate_regulatory_overlap()`

**Signature**: `annotate_regulatory_overlap(variants, regulatory_db)`

Annotate structural variants with overlapping regulatory elements.

For each variant, finds all regulatory elements (enhancers, promoters,
silencers, insulators, etc.) that overlap with the variant interval.

Args:
    variants: List of variant dictionaries (same format as annotate_gene_overlap).
    regulatory_db: List of regulatory element entries as dictionaries with
        'chrom', 'start', 'end', 'name', 'feature_type' keys, or as
        GenomicInterval objects.

Returns:
    List of variant dictionaries with added fields:
        - 'overlapping_regulatory': List of regulatory element names
        - 'regulatory_overlaps': List of OverlapResult objects
        - 'regulatory_types': Set of affected regulatory element types
        - 'n_regulatory_affected': Count of affected elements

#### `calculate_overlap_fraction()`

**Signature**: `calculate_overlap_fraction(interval_a, interval_b)`

Calculate reciprocal overlap fractions between two intervals.

Computes what fraction of interval A is covered by interval B and
vice versa. Both intervals must be on the same chromosome (caller
is responsible for this check).

Args:
    interval_a: First interval as dict with 'start'/'end' keys or
        (start, end) tuple.
    interval_b: Second interval in the same format.

Returns:
    Tuple of (fraction_of_a_covered, fraction_of_b_covered).
    Each value is between 0.0 and 1.0.

#### `find_nearest_gene()`

**Signature**: `find_nearest_gene(variant, gene_db, max_distance)`

Find the nearest gene to a structural variant.

If the variant overlaps a gene, the distance is 0. Otherwise, finds
the closest gene within max_distance base pairs of either breakpoint.

Args:
    variant: Variant dictionary with 'chrom', 'start', 'end' keys.
    gene_db: Gene database (see annotate_gene_overlap for format).
    max_distance: Maximum distance to search for nearest gene.

Returns:
    Dictionary with:
        - 'nearest_gene': Name of nearest gene (empty string if none found)
        - 'distance': Distance to nearest gene in bp (0 if overlapping, -1 if none found)
        - 'direction': 'upstream', 'downstream', 'overlapping', or 'none'
        - 'gene_interval': The GenomicInterval of the nearest gene, or None

#### `_normalize_intervals()`

**Signature**: `_normalize_intervals(features, feature_type)`

Normalize a list of features to GenomicInterval objects.

Args:
    features: List of dictionaries or GenomicInterval objects.
    feature_type: Default feature type if not specified in the data.

Returns:
    List of GenomicInterval objects.

#### `_compute_overlap_bp()`

**Signature**: `_compute_overlap_bp(a_start, a_end, b_start, b_end)`

Compute the number of overlapping base pairs between two intervals.

Args:
    a_start: Start of interval A.
    a_end: End of interval A.
    b_start: Start of interval B.
    b_end: End of interval B.

Returns:
    Number of overlapping base pairs (0 if no overlap).

#### `_classify_relationship()`

**Signature**: `_classify_relationship(v_start, v_end, f_start, f_end)`

Classify the spatial relationship between a variant and a feature.

Args:
    v_start: Variant start.
    v_end: Variant end.
    f_start: Feature start.
    f_end: Feature end.

Returns:
    Relationship string: 'contains', 'contained_in', or 'overlap'.

### Classes

#### `GenomicInterval`

A genomic interval.

Attributes:
    chrom: Chromosome name.
    start: Start position (0-based, inclusive).
    end: End position (0-based, exclusive).
    name: Feature name (e.g., gene name).
    feature_type: Type of feature (e.g., 'gene', 'exon', 'promoter').
    strand: Strand ('+', '-', or '.').
    info: Additional annotation fields.

#### `OverlapResult`

Result of an overlap analysis between a variant and a genomic feature.

Attributes:
    variant_chrom: Chromosome of the variant.
    variant_start: Start of the variant.
    variant_end: End of the variant.
    feature: The overlapping genomic feature.
    overlap_bp: Number of overlapping base pairs.
    overlap_fraction_variant: Fraction of the variant overlapping the feature.
    overlap_fraction_feature: Fraction of the feature overlapping the variant.
    relationship: Relationship type ('overlap', 'contained_in', 'contains').

#### `IntervalIndex`

Interval index for overlap and nearest-interval queries.

Intervals are grouped by chromosome and sorted by ``(start, end)``. Two
structures sit on that order, so a query never scans the whole chromosome:

* a segment tree of the maximum ``end`` over ranges of the sorted intervals
  reports the ``k`` intervals that overlap a region in ``O(k log n)``,
  however long the intervals are (a gene spanning most of a chromosome is
  found, and does not slow down the queries after it);
* the running maximum of the ends, together with the interval that reaches
  it, gives the nearest interval to the left of a position in ``O(1)``.

**Methods**:

| Method | Purpose |
|--------|---------|
| `__init__()` | Build an interval index. |
| `query_overlap()` | Find all intervals overlapping the query region. |
| `query_nearest()` | Find the nearest interval to a position. |

## detection.breakpoints

### Functions

#### `refine_breakpoints()`

**Signature**: `refine_breakpoints(variants, reads, window, min_clip)`

Refine breakpoint positions to base-pair resolution.

For each variant, examines split reads near the reported breakpoint
positions and determines the most likely exact breakpoint positions
by consensus of soft-clip positions.

A read is located by the reference coordinate of its soft clip, not by where
its alignment starts, so a read that spans the breakpoint from far away (a long
alignment, or a long deletion or skip in its CIGAR) counts as long as the clip
itself lies within ``window`` of the breakpoint.

The reads that clip at the consensus position are the breakpoint's supporting
reads. Their orientation, sequence and supplementary alignment give each
breakpoint's ``strand``, ``microhomology`` and ``inserted_sequence``; both
breakpoints of a variant are analysed the same way.

Args:
    variants: List of variant dictionaries, each containing:
        - 'chrom' or 'chrom1': Chromosome of breakpoint 1
        - 'start' or 'breakpoint1': Approximate position of breakpoint 1
        - 'end' or 'breakpoint2': Approximate position of breakpoint 2
        - 'pos2': Position of breakpoint 2 on 'chrom2' (translocations; it
          takes the place of 'end', which then describes breakpoint 1 only).
          ``dataclasses.asdict`` of a StructuralVariant provides it.
        - 'chrom2': Chromosome of breakpoint 2 (optional; absent, None or ''
          means the same chromosome as breakpoint 1)
        - 'sv_type': Type of structural variant (optional)
    reads: List of read alignment dictionaries with:
        - 'chrom': Chromosome
        - 'pos': Alignment start position
        - 'cigar': CIGAR string
        - 'name': Read name
        - 'is_reverse': Whether the read is reverse-complemented (optional;
          gives the breakpoint strand)
        - 'seq': Read sequence as stored in the alignment (optional; used for
          microhomology and inserted-sequence detection)
        - 'sa_tag': SA tag of a read with exactly one supplementary
          alignment (optional; with 'seq' and 'is_reverse' it gives the
          inserted sequence)
    window: Window size around each breakpoint to search for refinement.
    min_clip: Minimum soft-clip length to consider.

Returns:
    List of BreakpointPair objects with refined positions. A breakpoint with
    no soft clip in its window keeps its input position with ``support`` 0.

#### `detect_microhomology()`

**Signature**: `detect_microhomology(breakpoint, flanking_seq, max_homology)`

Detect microhomology at a breakpoint junction.

Examines the sequence around a breakpoint to find short regions of
sequence identity (microhomology), which are indicative of the
repair mechanism (e.g., NHEJ, MMEJ/alt-EJ).

The algorithm slides a window from the breakpoint position in both
directions, comparing bases for matches. The longest run of matching
bases starting from the breakpoint is reported.

Args:
    breakpoint: Breakpoint position info. Can be a dict with 'position'
        key or a Breakpoint object.
    flanking_seq: Sequence flanking the breakpoint. The breakpoint is
        assumed to be at the midpoint of this sequence.
    max_homology: Maximum microhomology length to search for.

Returns:
    Microhomology sequence string (empty if none found).

#### `cluster_breakpoints()`

**Signature**: `cluster_breakpoints(breakpoints, max_distance)`

Cluster nearby breakpoints using single-linkage clustering.

Groups breakpoints that are on the same chromosome and within
max_distance base pairs of each other. Returns clusters sorted
by cluster size (largest first).

Args:
    breakpoints: List of Breakpoint objects or dictionaries with
        'chrom' and 'position' keys.
    max_distance: Maximum distance between breakpoints to be
        assigned to the same cluster.

Returns:
    List of breakpoint clusters, each cluster is a list of Breakpoint
    objects. Clusters are sorted by size (descending).

#### `calculate_breakpoint_confidence()`

**Signature**: `calculate_breakpoint_confidence(evidence)`

Calculate confidence score for a breakpoint.

Combines multiple evidence metrics into a single confidence score:
- Number of supporting reads (more = higher confidence)
- Positional consistency (lower std = higher confidence)
- Ratio of supporting to total nearby split reads

Args:
    evidence: Dictionary or object containing:
        - 'support': Number of reads supporting the exact position (int)
        - 'total_clips': Total split reads in the region (int, optional)
        - 'position_std': Standard deviation of clip positions (float, optional)
        - 'mapq_mean': Mean mapping quality of supporting reads (float, optional)

Returns:
    Confidence score between 0 and 1.

#### `_first_present()`

**Signature**: `_first_present(variant, keys)`

Return the first value among ``keys`` that is not None (0 when none is).

#### `_refine_breakpoint()`

**Signature**: `_refine_breakpoint(chrom, approx, hits, other_chrom, other_approx, window)`

Build the refined breakpoint from the soft clips found near it.

Args:
    chrom: Chromosome of this breakpoint.
    approx: Approximate position (kept when there is no clip).
    hits: Soft clips within the window around ``approx``.
    other_chrom: Chromosome of the variant's other breakpoint.
    other_approx: Approximate position of the other breakpoint.
    window: Search window (also the tolerance for locating the other breakpoint).

Returns:
    The refined Breakpoint.

#### `_majority_strand()`

**Signature**: `_majority_strand(hits)`

Orientation of most supporting reads: '-' when more are reverse-complemented than not.

Reads that do not report ``is_reverse`` carry no orientation and are not counted;
with no orientation at all (or a tie) the breakpoint keeps the default '+'.

#### `_junction_microhomology()`

**Signature**: `_junction_microhomology(chrom, position, hits)`

Longest microhomology among the supporting reads, each centred on its junction.

`detect_microhomology()` assumes the breakpoint is at the midpoint of the
sequence it is given, so each read contributes the window of its ``seq`` that is
centred on its own clip junction (and not the whole read).

#### `_split_read_insertion()`

**Signature**: `_split_read_insertion(hits, other_chrom, other_approx, window)`

Most common inserted sequence among the supporting reads, or '' when not derivable.

Nothing but the read's own two alignments tells inserted bases from reference bases,
so a sequence is reported only for a read that has its ``seq``, its ``is_reverse``
flag and an ``sa_tag`` with exactly one supplementary alignment located at the other
breakpoint. Ties go to the longer sequence, then to the first read.

#### `_inserted_between_alignments()`

**Signature**: `_inserted_between_alignments(read, other_chrom, other_approx, window)`

Bases of ``read`` between its primary and its single supplementary alignment.

Both alignments cover an interval of the read (the CIGAR, leading clip to the end
of the aligned part). Expressed in the orientation of the stored ``seq``, the bases
between the two intervals are inserted; overlapping or abutting alignments
insert nothing. An SA entry on the opposite strand describes the reverse-complemented
read, so its interval is mirrored (its leading clip becomes its trailing clip).

Args:
    read: Read dictionary with ``seq``, ``cigar``, ``is_reverse`` and ``sa_tag``.
    other_chrom: Chromosome the supplementary alignment must lie on.
    other_approx: Position the supplementary alignment must lie within ``window`` of.
    window: Positional tolerance in bp.

Returns:
    The inserted sequence, or '' when it cannot be derived.

#### `_cigar_ops()`

**Signature**: `_cigar_ops(cigar)`

Split a CIGAR into ``(operation letter, length)`` pairs.

Args:
    cigar: A CIGAR string (``"50M30S"``) or a list of ``(BAM operation code, length)``
        tuples (``[(0, 50), (4, 30)]``).

#### `_query_extent()`

**Signature**: `_query_extent(cigar)`

Return ``(leading clip, aligned query length, trailing clip, reference length)`` of a CIGAR.

The leading and trailing clips count soft and hard clipped bases; the aligned query
length counts M, I, = and X; the reference length counts M, D, N, = and X.
Returns None for an empty or unparseable CIGAR.

#### `_parse_sa_alignments()`

**Signature**: `_parse_sa_alignments(sa_tag)`

Parse an SA tag (``rname,pos,strand,CIGAR,mapQ,NM;...``) into ``(chrom, start, strand, CIGAR)``.

``start`` is the 0-based leftmost reference position. Entries without a numeric
position or a CIGAR are skipped.

#### `_parse_clips_from_cigar()`

**Signature**: `_parse_clips_from_cigar(cigar, pos, min_clip)`

Parse CIGAR to find soft-clip junctions.

Args:
    cigar: CIGAR string, or a list of ``(BAM operation code, length)`` tuples.
    pos: Alignment start position.
    min_clip: Minimum clip length.

Returns:
    List of ``(clip_position, clip_length, query_offset)`` tuples: the reference
    coordinate of the junction between the clip and the aligned part, the clip
    length, and the offset of that junction in the read's stored sequence (hard
    clipped bases are not stored, so they are not counted).

#### `_consensus_position()`

**Signature**: `_consensus_position(clips, fallback)`

Find consensus breakpoint position from clip positions.

Uses the mode (most frequent position) as the consensus. Falls back
to the provided position if no clips are found.

Args:
    clips: List of (position, read_name) tuples.
    fallback: Fallback position if no clips are available.

Returns:
    Tuple of (consensus_position, support_count, read_names).

#### `_position_std()`

**Signature**: `_position_std(clips)`

Calculate standard deviation of clip positions.

Args:
    clips: List of (position, read_name) tuples.

Returns:
    Standard deviation of positions. Returns 0 if fewer than 2 clips.

### Classes

#### `Breakpoint`

A structural variant breakpoint.

Attributes:
    chrom: Chromosome name.
    position: Genomic position (0-based).
    strand: Strand orientation ('+' or '-'): the orientation of most of the
        reads supporting the breakpoint ('-' when more are reverse-complemented
        than not). '+' when no supporting read reports ``is_reverse``.
    confidence: Confidence score [0, 1].
    support: Number of supporting reads.
    microhomology: Microhomology sequence at the breakpoint, if any.
    inserted_sequence: Bases inserted between the two alignments of a split
        read at this breakpoint, as they read along the forward strand of the
        alignment at this breakpoint; empty when not derivable (see
        `refine_breakpoints()`).
    read_names: Names of reads supporting this breakpoint.

#### `BreakpointPair`

A pair of breakpoints defining a structural variant.

Attributes:
    bp1: First breakpoint.
    bp2: Second breakpoint.
    sv_type: Inferred structural variant type.
    confidence: Joint confidence score.

#### `_ClipHit`

One soft-clipped junction of one read (internal).

Attributes:
    position: Reference coordinate of the junction.
    name: Read name.
    clip_length: Length of the soft clip.
    query_offset: Offset of the junction in the read's ``seq`` (bases before it).
    order: ``(read index, clip index)``: the order in which the clips were read.
    read: The read dictionary the clip belongs to.

#### `_ClipIndex`

Soft-clipped junctions of a read set, sorted by reference position per chromosome.

Every read is parsed once and found by the coordinate of its clip, never by where
its alignment starts: a long alignment (or a long deletion or skip in its CIGAR)
puts the clip far from the read start, so the start position says nothing about
which reads clip near a breakpoint.

**Methods**:

| Method | Purpose |
|--------|---------|
| `near()` | Return the clips within ``window`` bp of ``position``, in the order they were read. |

## detection.cnv

### Functions

#### `detect_cnv_from_depth()`

**Signature**: `detect_cnv_from_depth(depth_data, window_size, method, significance, ploidy, min_segment_bins, merge_distance)`

Detect copy number variations from read depth data.

Analyzes per-chromosome depth arrays to identify regions of copy number
gain or loss using circular binary segmentation (CBS).

Args:
    depth_data: Dictionary mapping chromosome names to arrays of read depth
        values. Each value represents the mean depth in a genomic window and
        must be finite and non-negative (zero is a bin with no coverage).
        Can be numpy arrays or plain lists of floats. Log2 ratios are not
        depths: pass them to `segment_coverage()` instead.
    window_size: Size of each depth bin in base pairs.
    method: Detection method. Currently supports 'segmentation' (CBS).
    significance: P-value threshold for CBS change-point detection.
    ploidy: Expected ploidy of the organism (default 2 for diploid).
    min_segment_bins: Minimum number of bins for a valid segment.
    merge_distance: Maximum gap (in bp) to merge adjacent same-state segments.

Returns:
    Dictionary mapping chromosome names to CNVResult objects containing
    detected segments with copy number states.

Raises:
    ValueError: If depth_data is empty, method is unsupported, or a chromosome
        that would be analysed holds negative or non-finite depth values.

#### `_check_depths()`

**Signature**: `_check_depths(chrom, depth_array)`

Raise a ValueError naming the problem if the values cannot be read-depth counts.

Args:
    chrom: Chromosome name (for the message).
    depth_array: Depth values as a float array.

Raises:
    ValueError: If any value is negative or not finite.

#### `segment_coverage()`

**Signature**: `segment_coverage(coverage_array, significance, max_iterations, min_width)`

Circular Binary Segmentation (CBS) algorithm for coverage segmentation.

Recursively finds change-points in a coverage array by testing for
differences in means across circular permutations of the data, splitting
at the most significant change-point, and recursing on each half.

Every arc of at most half the data is scored, so the work grows with the
square of the segment length; the scan is vectorised over all arcs of one
length, which makes it about two orders of magnitude faster than scoring
the arcs one at a time (and gives exactly the same segmentation).

Args:
    coverage_array: Array of log2 ratio values to segment. Can be a numpy
        array or list of floats.
    significance: P-value threshold for accepting a change-point. Lower
        values produce fewer, more confident segments.
    max_iterations: Maximum recursion depth to prevent infinite loops.
    min_width: Minimum segment width in bins (at least 1).

Returns:
    List of (start, end, mean) tuples representing segments, where start
    is inclusive, end is exclusive, and mean is the segment mean.

Raises:
    ValueError: If min_width is less than 1.

#### `_cbs_best_arc()`

**Signature**: `_cbs_best_arc(segment_data, std_dev, min_width)`

Find the arc whose mean differs most from the rest of the segment.

Scores every arc ``[i, j)`` with ``min_width <= j - i <= n // 2`` (and at least
``min_width`` bins outside it) by the standardised difference between the mean
inside and outside the arc,
``|mean_in - mean_out| * sqrt(len_in * len_out / n) / std_dev``. A change point
is the arc that touches an end of the segment, so the single-split scan of plain
binary segmentation is the special case of this one.

The arcs of one length are scored together with array arithmetic that performs
the same floating-point operations in the same order as scoring them one at a
time, so the statistics are bit-for-bit the same. The first arc in
``(i, j)`` order wins a tie.

Args:
    segment_data: The segment's values (at least ``2 * min_width`` of them).
    std_dev: Sample standard deviation of the segment (positive).
    min_width: Minimum arc length, and minimum length outside the arc.

Returns:
    ``(statistic, i, j)`` of the best arc, or ``(0.0, -1, -1)`` if no arc scores above zero.

#### `_cbs_recursive()`

**Signature**: `_cbs_recursive(data, start, end, segments, significance, max_iterations, min_width, depth)`

Recursive CBS implementation.

Tests for a change-point in data[start:end] using the maximum circular
binary segmentation statistic. If significant, splits and recurses.

Args:
    data: Full data array.
    start: Start index (inclusive).
    end: End index (exclusive).
    segments: Accumulator list for resulting segments.
    significance: P-value threshold.
    max_iterations: Maximum recursion depth.
    min_width: Minimum segment width.
    depth: Current recursion depth.

#### `_cbs_pvalue()`

**Signature**: `_cbs_pvalue(statistic, n)`

Approximate p-value for CBS test statistic.

Uses the asymptotic approximation from Olshen et al. (2004) where the
maximum of the standardized CBS statistic converges to a Gumbel
distribution with parameters depending on the segment length.

Args:
    statistic: CBS test statistic value.
    n: Number of data points in the segment.

Returns:
    Approximate p-value.

#### `call_cnv_states()`

**Signature**: `call_cnv_states(segments, ploidy, del_threshold, dup_threshold, amp_threshold, homodel_threshold)`

Assign copy number states to segments based on log2 ratio thresholds.

Classifies each segment as homozygous deletion, deletion, neutral,
duplication, or amplification based on its mean log2 ratio value
relative to configurable thresholds.

Args:
    segments: List of segments as (start, end, mean_log2ratio) tuples
        or CNVSegment objects.
    ploidy: Expected ploidy (default 2 for diploid).
    del_threshold: Log2 ratio threshold below which a segment is called
        as a deletion (default -0.3, approximately CN=1 for diploid).
    dup_threshold: Log2 ratio threshold above which a segment is called
        as a duplication (default 0.3, approximately CN=3 for diploid).
    amp_threshold: Log2 ratio threshold for high-level amplification
        (default 1.0, approximately CN=4+ for diploid).
    homodel_threshold: Log2 ratio threshold for homozygous deletion
        (default -1.5, approximately CN=0 for diploid).

Returns:
    List of CNVSegment objects with state and cn fields populated.

#### `_assign_cnv_state()`

**Signature**: `_assign_cnv_state(mean_log2ratio, ploidy, del_threshold, dup_threshold, amp_threshold, homodel_threshold)`

Assign a copy number state and integer CN from a log2 ratio.

Args:
    mean_log2ratio: Mean log2 ratio value for the segment.
    ploidy: Expected ploidy.
    del_threshold: Deletion threshold.
    dup_threshold: Duplication threshold.
    amp_threshold: Amplification threshold.
    homodel_threshold: Homozygous deletion threshold.

Returns:
    Tuple of (state_string, integer_copy_number).

#### `_calculate_segment_confidence()`

**Signature**: `_calculate_segment_confidence(log2_ratios, segment_mean, state)`

Calculate confidence score for a CNV segment.

Confidence is based on: (1) deviation from neutral, (2) variance within
the segment (lower is better), and (3) number of supporting bins.

Args:
    log2_ratios: Array of log2 ratio values within the segment.
    segment_mean: Mean log2 ratio for the segment.
    state: Assigned CNV state.

Returns:
    Confidence score between 0 and 1.

#### `merge_adjacent_segments()`

**Signature**: `merge_adjacent_segments(segments, max_gap)`

Merge adjacent segments with the same copy number state.

Combines neighboring segments on the same chromosome that share
the same CNV state and are separated by no more than max_gap base pairs.

Args:
    segments: List of CNVSegment objects, ideally sorted by position.
    max_gap: Maximum gap in base pairs between segments to allow merging.

Returns:
    List of merged CNVSegment objects.

#### `calculate_log2_ratio()`

**Signature**: `calculate_log2_ratio(tumor_depth, normal_depth, gc_content, pseudocount)`

Compute log2 ratio between tumor and normal read depths.

Calculates the log2(tumor/normal) ratio with optional GC-content bias
correction using LOESS-like smoothing.

Args:
    tumor_depth: Array of tumor read depth values per window.
    normal_depth: Array of normal (reference) read depth values per window.
    gc_content: Optional array of GC content fractions per window (0-1)
        for GC bias correction. If provided, the ratio is corrected
        using median normalization within GC bins.
    pseudocount: Small value added to avoid division by zero.

Returns:
    Numpy array of log2 ratio values. If numpy is not available,
    returns a list of floats.

Raises:
    ValueError: If arrays have different lengths.

#### `_gc_correct()`

**Signature**: `_gc_correct(log2_ratios, gc_content, n_bins)`

GC-content bias correction using binned median normalization.

Divides the GC content range into bins, computes the median log2 ratio
within each GC bin, and subtracts the bin-specific median to remove
GC-dependent bias.

Args:
    log2_ratios: Array of log2 ratio values.
    gc_content: Array of GC content fractions (0-1).
    n_bins: Number of GC bins for correction.

Returns:
    GC-corrected log2 ratio array.

#### `_center_on_median()`

**Signature**: `_center_on_median(values)`

Subtract the median of the finite values (pure-Python median centering).

Args:
    values: Log2 ratio values.

Returns:
    Values shifted so that the median of the finite entries is zero.

#### `_gc_correct_py()`

**Signature**: `_gc_correct_py(log2_ratios, gc_content, n_bins)`

Pure-Python twin of `_gc_correct()` (binned median normalization).

Args:
    log2_ratios: Log2 ratio values.
    gc_content: GC content fractions (0-1), one per log2 ratio.
    n_bins: Number of GC bins for correction.

Returns:
    GC-corrected, median-centered log2 ratios.

### Classes

#### `CNVSegment`

A segment of constant copy number.

Attributes:
    chrom: Chromosome name.
    start: 0-based start position (inclusive).
    end: 0-based end position (exclusive).
    mean_log2ratio: Mean log2 ratio for this segment.
    n_bins: Number of bins (windows) in the segment.
    state: Copy number state ('DEL', 'DUP', 'NEUTRAL', 'AMP', 'HOMODEL').
    cn: Estimated integer copy number.
    confidence: Confidence score [0, 1].

#### `CNVResult`

Result of CNV detection.

Attributes:
    segments: List of CNV segments detected.
    log2_ratios: The log2 ratio array used for detection.
    bin_size: Size of each bin/window in base pairs.
    method: Detection method used.
    parameters: Dictionary of parameters used.

## detection.sv_calling

### Functions

#### `call_structural_variants()`

**Signature**: `call_structural_variants(alignments, min_support, insert_size_stats, min_mapq, min_sv_size)`

Call structural variants from alignment data.

Combines split-read and discordant read-pair evidence to detect
structural variants. Reads are first classified as split or discordant,
then evidence is clustered by genomic position, and variant calls are
made from clusters with sufficient support.

Args:
    alignments: List of read alignment dictionaries, each containing:
        - 'name': Read name (str)
        - 'chrom': Reference chromosome (str)
        - 'pos': Alignment start position, 0-based (int)
        - 'mapq': Mapping quality (int)
        - 'cigar': CIGAR string (str) or list of (op, length) tuples
        - 'mate_chrom': Mate chromosome (str, optional)
        - 'mate_pos': Mate start position (int, optional)
        - 'insert_size': Observed insert size (int, optional)
        - 'is_reverse': Whether read is reverse-complemented (bool, optional)
        - 'mate_is_reverse': Whether mate is reverse-complemented (bool, optional)
        - 'is_supplementary': Whether this is a supplementary alignment (bool, optional)
        - 'sa_tag': SA tag for supplementary alignments (str, optional)
    min_support: Minimum number of supporting reads for a call.
    insert_size_stats: Pre-computed insert size distribution statistics.
        If None, will be estimated from the data.
    min_mapq: Minimum mapping quality for reads to be considered.
    min_sv_size: Minimum structural variant size in base pairs. It applies
        to calls that have a measurable length (DEL, DUP, INV and
        UNKNOWN, measured as the breakpoint distance). Translocations have
        no length, and the inserted length of an insertion is not measured
        by split-read evidence, so TRA and INS calls are exempt (their
        ``size`` is 0).

Returns:
    List of StructuralVariant objects, sorted by chromosome and position.
    A translocation (TRA) carries its second breakpoint in ``chrom2`` and
    ``pos2``; its ``start``/``end`` describe the first breakpoint on
    ``chrom`` only (see `StructuralVariant`).

#### `detect_split_reads()`

**Signature**: `detect_split_reads(reads, min_clip)`

Identify split-read evidence for structural variants.

Examines CIGAR strings and supplementary alignment information to
detect reads that span structural variant breakpoints. A split read
has a soft-clipped portion of at least min_clip bases, indicating
the read crosses a breakpoint.

Args:
    reads: List of read alignment dictionaries (see call_structural_variants
        for format).
    min_clip: Minimum soft-clip length (bases) to consider as evidence
        for a breakpoint.

Returns:
    List of SVEvidence objects, one per split-read breakpoint detected.

#### `detect_discordant_pairs()`

**Signature**: `detect_discordant_pairs(pairs, insert_size_stats, n_std)`

Detect discordant read pair evidence for structural variants.

A discordant read pair has an insert size significantly different from
the expected distribution, or has reads mapping to different chromosomes,
or has unexpected orientation patterns.

Args:
    pairs: List of read alignment dictionaries (see call_structural_variants
        for format). Only reads with mate information are analyzed.
    insert_size_stats: Insert size distribution statistics.
    n_std: Number of standard deviations from the mean to classify
        a pair as discordant (default 4.0).

Returns:
    List of SVEvidence objects for discordant pairs.

#### `classify_sv_type()`

**Signature**: `classify_sv_type(evidence)`

Classify the type of structural variant based on evidence patterns.

Uses breakpoint positions, chromosome assignments, and strand
orientations to determine the SV type. Orientation signatures, stated
for the pair normalised to genomic order (upstream read first):

    - TRA: breakpoints on different chromosomes.
    - INS: split reads whose breakpoints nearly coincide (within
      ``_INS_MAX_BREAKPOINT_DISTANCE``) -- the two sides of a split
      read land at the same insertion junction.
    - INV: both reads on the same strand.
    - DUP: everted (tail-to-tail) pair -- the upstream read is on '-'
      and the downstream read on '+', so the reads point away from
      each other across the duplication junction.
    - DEL: inward-facing (head-to-head) pair -- upstream read on '+',
      downstream on '-' -- with split-read support or a span beyond a
      normal fragment (``_DEL_MIN_BREAKPOINT_DISTANCE``).

A call is only made when the strand/direction evidence supports it.
Weak or unrecognisable evidence (e.g. a normal-looking short-span FR
pair without split support, or unknown strand characters) returns
``SVType.UNKNOWN``; the classifier never manufactures a deletion from
weak evidence.

Args:
    evidence: SVEvidence object containing breakpoint and orientation info.

Returns:
    SVType enum value.

#### `genotype_sv()`

**Signature**: `genotype_sv(variant, reads, min_mapq)`

Genotype a structural variant based on supporting and reference reads.

Counts the number of reads supporting the variant (split reads and
discordant pairs) versus reads supporting the reference allele
(reads spanning the breakpoint region without split alignment).
Uses a simple allele-fraction model:
    - AF < 0.15: homozygous reference (0/0)
    - 0.15 <= AF < 0.85: heterozygous (0/1)
    - AF >= 0.85: homozygous alternate (1/1)

The reference reads are those spanning ``start - 100 .. start + 100`` (the
first breakpoint) on ``variant.chrom``.

Args:
    variant: StructuralVariant to genotype.
    reads: List of read alignment dictionaries, or a `ReadIndex`
        built from them. Pass an index when genotyping many variants
        against the same reads: it gives the same answer without scanning
        every read for every variant.
    min_mapq: Minimum mapping quality for counting reads.

Returns:
    Genotype string: '0/0', '0/1', or '1/1' ('./.' when no read supports
    either allele).

#### `_estimate_insert_size()`

**Signature**: `_estimate_insert_size(reads)`

Estimate insert size distribution from properly paired reads.

Samples insert sizes from reads that appear properly paired (same
chromosome, FR orientation) and computes robust statistics.

Args:
    reads: List of read alignment dictionaries.

Returns:
    InsertSizeStats with estimated distribution parameters.

#### `_parse_cigar_clips()`

**Signature**: `_parse_cigar_clips(cigar, pos)`

Parse CIGAR string to find soft-clipped regions.

Args:
    cigar: CIGAR string (e.g., '50M30S') or list of (operation, length) tuples.
    pos: Alignment start position.

Returns:
    List of (clip_position, clip_length, side) tuples where side is
    'left' or 'right'.

#### `_parse_cigar_string()`

**Signature**: `_parse_cigar_string(cigar)`

Parse a CIGAR string into (operation, length) tuples.

CIGAR operations: M=0, I=1, D=2, N=3, S=4, H=5, P=6, ==7, X=8

Args:
    cigar: CIGAR string (e.g., '50M30S').

Returns:
    List of (operation_code, length) tuples.

#### `_parse_sa_tag()`

**Signature**: `_parse_sa_tag(sa_tag)`

Parse SAM SA (supplementary alignment) tag.

SA tag format: rname,pos,strand,CIGAR,mapQ,NM;

Args:
    sa_tag: SA tag string.

Returns:
    List of (chromosome, position, strand) tuples.

#### `_cluster_evidence()`

**Signature**: `_cluster_evidence(evidence_list, max_distance)`

Cluster SV evidence by genomic proximity.

Groups evidence items whose breakpoints are within max_distance of each
other on the same chromosome pair, using single-linkage clustering.

Args:
    evidence_list: List of SVEvidence objects.
    max_distance: Maximum distance between breakpoints for clustering.

Returns:
    List of evidence clusters (each cluster is a list of SVEvidence).

#### `_merge_evidence()`

**Signature**: `_merge_evidence(cluster)`

Merge a cluster of evidence into a single SVEvidence.

Computes consensus breakpoint positions (median), totals support counts,
and combines read names.

Args:
    cluster: List of SVEvidence items to merge.

Returns:
    Merged SVEvidence with consensus breakpoints and combined support.

#### `_calculate_sv_quality()`

**Signature**: `_calculate_sv_quality(evidence, total_support)`

Calculate quality score for a structural variant call.

Quality is derived from: total supporting evidence, balance between
split reads and discordant pairs, and uniqueness of supporting reads.

Args:
    evidence: Merged SVEvidence.
    total_support: Total support count.

Returns:
    Quality score (Phred-scaled, capped at 999).

### Classes

#### `SVType`

Structural variant types.

#### `InsertSizeStats`

Insert size distribution statistics.

Attributes:
    mean: Mean insert size.
    std: Standard deviation of insert size.
    median: Median insert size.
    mad: Median absolute deviation.

#### `SVEvidence`

Evidence supporting a structural variant call.

Attributes:
    split_reads: Number of supporting split reads.
    discordant_pairs: Number of supporting discordant read pairs.
    depth_support: Whether read depth supports the call.
    evidence_reads: List of read names providing evidence.
    breakpoint1: Estimated position of first breakpoint.
    breakpoint2: Estimated position of second breakpoint.
    chrom1: Chromosome of first breakpoint.
    chrom2: Chromosome of second breakpoint.
    strand1: Strand orientation at breakpoint 1.
    strand2: Strand orientation at breakpoint 2.

#### `StructuralVariant`

A called structural variant.

Intra-chromosomal calls (DEL, DUP, INV, INS, UNKNOWN) occupy the interval
``chrom:start-end``. A translocation (TRA) joins two loci on different
chromosomes, so it has no interval: ``chrom``/``start`` give the first
breakpoint, ``chrom2``/``pos2`` the second, and ``end`` repeats ``start``
(a zero-length junction, like an insertion point).

Attributes:
    chrom: Primary chromosome (the chromosome of the first breakpoint).
    start: Start position (0-based). For a TRA, the first breakpoint on ``chrom``.
    end: End position (0-based). For a TRA, equal to ``start``.
    sv_type: Type of structural variant.
    size: Size of the variant in base pairs. 0 means "not measured": it is
        always 0 for TRA (no length) and for INS (split-read evidence does
        not measure the inserted length; the breakpoint distance is not it).
    quality: Quality score.
    evidence: Supporting evidence.
    genotype: Genotype string (e.g., '0/1', '1/1').
    filter_status: Filter status ('PASS' or reason).
    chrom2: Chromosome of the second breakpoint (translocations; '' otherwise).
    info: Additional info fields.
    pos2: Position (0-based) of the second breakpoint on ``chrom2``
        (translocations; None otherwise).

#### `ReadIndex`

Per-chromosome index of aligned reads by start position.

Genotyping a call only needs the reads that span the 200 bp window around its
first breakpoint. Scanning every read for every call costs
``O(calls x reads)`` (about 90% of `call_structural_variants()` at 100k
reads and 100 calls); with this index a call costs ``O(log reads)`` plus the
reads that start within one read length of the window.

Build it once from the reads and pass it to `genotype_sv()` in place of the
read list. Results are identical to scanning the list.

Args:
    reads: Read alignment dictionaries (see `call_structural_variants()`).
        ``chrom``, ``pos``, ``mapq``, ``name`` and ``read_length`` are used,
        with the same defaults as `genotype_sv()` (``read_length`` 150).

**Methods**:

| Method | Purpose |
|--------|---------|
| `count_spanning()` | Count reads that start before ``window_start`` and end after ``window_end``. |

## filtering.merge

### Functions

#### `merge_callsets()`

**Signature**: `merge_callsets(callsets, min_overlap, type_match)`

Merge structural variant callsets from multiple callers.

Uses reciprocal overlap to identify the same SV detected by different
tools. For each group of matching variants across callers, produces a
single consensus call with aggregated support: median breakpoints, the
majority SV type, the majority genotype (see `_consensus_genotype()`)
and ``consensus_quality`` = the mean ``quality`` (or ``qual``) of the merged
calls plus ``100 * log10(n_callers + 1)``.

Args:
    callsets: Dictionary mapping caller names to lists of variant
        dictionaries. Each variant should have 'chrom', 'start',
        'end', 'sv_type' keys.
    min_overlap: Minimum reciprocal overlap fraction (0-1) to consider
        two variants as the same event (default 0.5 = 50%).
    type_match: Whether SV types must match for merging (default True).

Returns:
    Tuple of (merged_variants, merge_statistics).

#### `calculate_reciprocal_overlap()`

**Signature**: `calculate_reciprocal_overlap(sv1, sv2)`

Calculate the minimum reciprocal overlap between two SVs.

Reciprocal overlap requires that the overlap constitutes at least a
certain fraction of BOTH variants. The minimum of the two fractions
is returned (i.e., the stricter criterion).

Variants on different chromosomes (or without coordinates) never
overlap and yield 0.0.

Args:
    sv1: First variant dictionary with 'chrom', 'start', 'end' keys.
    sv2: Second variant dictionary with 'chrom', 'start', 'end' keys.

Returns:
    Reciprocal overlap fraction between 0.0 and 1.0.

#### `survivor_merge()`

**Signature**: `survivor_merge(vcf_files, max_distance, min_callers, type_match, strand_match)`

SURVIVOR-style merging of structural variant callsets.

Implements the SURVIVOR merge algorithm that uses breakpoint distance
rather than reciprocal overlap for matching variants. This is more
appropriate for imprecise SV calls where exact breakpoints are uncertain.

Note: If vcf_files contains strings (file paths), this function requires
pysam or a VCF parser. If vcf_files contains dictionaries, they are
used directly.

Each merged variant combines the calls of its group as follows:

* ``genotype``: majority vote among the calls' non-missing genotypes
  (allele order ignored, so '1/0' is '0/1'). A tie goes to the genotype
  whose calls have the larger summed quality, then to the one that sorts
  first. ``'./.'`` when no call has a genotype. Genotypes are read from
  each dict's ``genotype`` and from the first sample of VCF files.
* ``consensus_quality``: the mean ``quality`` (or ``qual``; 0 when absent)
  of the merged calls, plus an agreement bonus of
  ``10 * n_callers * log10(n_callers + 1)``.

Args:
    vcf_files: Either file paths to VCF files or lists of variant
        dictionaries. If file paths, each file is treated as one caller.
        If dicts, each dict should have a '_caller' key.
    max_distance: Maximum distance between breakpoints to consider
        variants as the same event (default 1000bp).
    min_callers: Minimum number of callers required to keep a merged
        variant (default 2).
    type_match: Whether SV types must match (default True).
    strand_match: Whether strand orientations must match (default False).

Returns:
    Tuple of (merged_variants, merge_statistics).

#### `deduplicate_variants()`

**Signature**: `deduplicate_variants(variants, max_distance, min_reciprocal_overlap, type_match)`

Remove duplicate structural variant calls.

Identifies duplicate calls based on proximity and reciprocal overlap,
keeping the highest-quality call from each group of duplicates.

Args:
    variants: List of variant dictionaries with 'chrom', 'start', 'end',
        'sv_type', and optionally 'quality' keys.
    max_distance: Maximum breakpoint distance for duplicate candidates.
    min_reciprocal_overlap: Minimum reciprocal overlap to classify as
        duplicate (default 0.8 = 80%).
    type_match: Whether SV types must match (default True).

Returns:
    Tuple of (deduplicated_variants, n_removed).

#### `_union()`

**Signature**: `_union(parent, i, j)`

Union operation for union-find with path compression.

Args:
    parent: Parent array.
    i: First element.
    j: Second element.

#### `_find()`

**Signature**: `_find(parent, i)`

Find operation for union-find with path compression.

Args:
    parent: Parent array.
    i: Element to find root of.

Returns:
    Root element.

#### `_median_position()`

**Signature**: `_median_position(values)`

Median of a list of positions (consensus breakpoint).

Args:
    values: Positions from the variants in a merge group.

Returns:
    Median position as an int.

#### `_majority_sv_type()`

**Signature**: `_majority_sv_type(variant_dicts)`

Consensus SV type by majority vote within a merge group.

Args:
    variant_dicts: Variants in the merge group.

Returns:
    Most frequent SV type (enum values are unwrapped first).

#### `_call_quality()`

**Signature**: `_call_quality(variant)`

Quality of one call (``quality``, else ``qual``); 0 when absent or None.

#### `_consensus_genotype()`

**Signature**: `_consensus_genotype(genotypes, qualities)`

Determine consensus genotype from multiple calls.

Uses majority voting among non-missing genotypes, ignoring allele order.

Args:
    genotypes: List of genotype strings.
    qualities: Optional quality of each call, parallel to ``genotypes``. When
        given, a tie in the vote goes to the genotype with the larger summed
        quality, then to the one that sorts first. Without it a tie goes to the
        genotype seen first.

Returns:
    Consensus genotype string.

#### `_format_genotype()`

**Signature**: `_format_genotype(alleles)`

Format a called allele tuple (as pysam reports ``GT``) as ``'0/1'``; None when absent.

#### `_vcf_sample_genotype()`

**Signature**: `_vcf_sample_genotype(fields)`

Genotype of the first sample of a VCF record line (``'0/1'``), or None when there is none.

Phased genotypes (``0|1``) are written with ``/``; a missing call is ``'./.'``.

#### `_load_vcf_variants()`

**Signature**: `_load_vcf_variants(vcf_path)`

Load structural variants from a VCF file.

Attempts to use pysam for VCF parsing. Falls back to basic text
parsing if pysam is not available.

Args:
    vcf_path: Path to VCF file.

Returns:
    List of variant dictionaries.

#### `_pysam_first_genotype()`

**Signature**: `_pysam_first_genotype(record)`

Genotype of the first sample of a pysam record, or None when it has no ``GT``.

#### `_parse_vcf_basic()`

**Signature**: `_parse_vcf_basic(vcf_path)`

Basic VCF parser without pysam dependency.

Args:
    vcf_path: Path to VCF file.

Returns:
    List of variant dictionaries. A record with a sample column also carries
    the ``genotype`` of its first sample.

### Classes

#### `MergedVariant`

A structural variant resulting from merging multiple callsets.

Attributes:
    chrom: Chromosome.
    start: Consensus start position.
    end: Consensus end position.
    sv_type: Consensus SV type.
    size: Consensus size.
    n_callers: Number of callers that detected this variant.
    caller_names: Names of callers that detected this variant (sorted).
    support_variants: Original variant calls that were merged.
    consensus_quality: Quality based on multi-caller agreement (the rule
        differs per merge function; see `merge_callsets()` and
        `survivor_merge()`).
    genotype: Consensus genotype ('./.' when no merged call has one).

#### `MergeStats`

Statistics from a merge operation.

Attributes:
    n_input_callsets: Number of input callsets.
    n_input_variants: Total variants across all callsets.
    n_output_variants: Number of merged variants.
    n_multi_caller: Variants detected by multiple callers.
    n_single_caller: Variants detected by only one caller.
    caller_counts: Variants per caller.

## filtering.quality_filter

### Functions

#### `filter_by_quality()`

**Signature**: `filter_by_quality(variants, min_qual, min_support, min_mapq, min_gq)`

Filter structural variants by quality metrics.

Removes variants that do not meet minimum quality thresholds.
Multiple quality metrics can be applied simultaneously.

Args:
    variants: List of variant dictionaries. Expected keys:
        - 'quality' or 'qual': Quality score (float)
        - 'support' or 'n_support': Number of supporting reads (int)
        - 'mapq': Mean mapping quality of supporting reads (float, optional)
        - 'gq': Genotype quality (float, optional)
    min_qual: Minimum quality score (Phred-scaled, default 20).
    min_support: Minimum number of supporting reads (default 3).
    min_mapq: Minimum mean mapping quality (default 0, no filter).
    min_gq: Minimum genotype quality (default 0, no filter).

Returns:
    Tuple of (filtered_variants, filter_statistics).

#### `filter_by_size()`

**Signature**: `filter_by_size(variants, min_size, max_size)`

Filter structural variants by size.

Removes variants outside the specified size range. This is commonly
used to separate SVs (>=50bp) from indels (<50bp) and to exclude
extremely large artifacts.

Args:
    variants: List of variant dictionaries. Expected keys:
        - 'size': Variant size in base pairs (int)
        - Or 'start' and 'end' to compute size
        - 'sv_type': Optional; translocations (TRA, BND) are always kept
          since they don't have a meaningful "size", and so are insertions
          (INS) without a positive 'size': the inserted length of a call from
          split reads is not measured (``size`` 0), and the start-to-end span
          of such a call is its junction, not that length. An insertion with a
          positive 'size' is filtered on it like any other variant.
    min_size: Minimum size in base pairs (default 50, the standard
        SV size threshold).
    max_size: Maximum size in base pairs (default None, no upper limit).

Returns:
    Tuple of (filtered_variants, filter_statistics).

#### `filter_by_frequency()`

**Signature**: `filter_by_frequency(variants, population_db, max_af, match_window, min_reciprocal_overlap)`

Filter structural variants by population allele frequency.

Removes common variants that are likely benign based on their
frequency in population databases (e.g., gnomAD-SV, DGV).

Variant matching uses both position proximity and reciprocal overlap
to find matching database entries.

Args:
    variants: List of variant dictionaries with 'chrom', 'start', 'end',
        'sv_type' keys. May also have 'af' or 'allele_frequency' pre-computed.
    population_db: Population frequency database. Can be:
        - dict mapping 'chrom:start-end' to allele frequency
        - list of dicts with 'chrom', 'start', 'end', 'af' keys
        - None (no filtering applied, all variants pass)
    max_af: Maximum population allele frequency to keep (default 0.01 = 1%).
    match_window: Window in bp around breakpoints for matching database entries.
    min_reciprocal_overlap: Minimum reciprocal overlap fraction for matching.

Returns:
    Tuple of (filtered_variants, filter_statistics).

Raises:
    ValueError: If ``population_db`` is given but none of its entries is in a
        supported layout (so the filter would silently do nothing). Entries
        of a mostly usable database that cannot be read are skipped with a
        logged warning that counts them.

#### `apply_blacklist()`

**Signature**: `apply_blacklist(variants, blacklist_regions, min_overlap_fraction)`

Remove variants overlapping blacklist regions.

Blacklist regions are genomic intervals known to produce false-positive
SV calls, such as centromeres, telomeres, segmental duplications,
low-complexity regions, and assembly gaps.

Args:
    variants: List of variant dictionaries with 'chrom', 'start', 'end'.
    blacklist_regions: Blacklist as a list of dicts with 'chrom', 'start',
        'end' keys, or as (chrom, start, end) tuples.
    min_overlap_fraction: Minimum fraction of variant overlapping a blacklist
        region to trigger removal (default 0.0, any overlap removes the variant).
        A zero-length variant (``start == end``, such as an insertion point)
        has no bases to overlap: it counts as fully overlapping (fraction 1.0)
        when it lies strictly inside a region, and not at all on a region's
        edge, the same as a variant that only abuts a region.

Returns:
    Tuple of (filtered_variants, filter_statistics).

#### `_as_position()`

**Signature**: `_as_position(value)`

Return ``value`` as an int coordinate, or None when it is not an integral number.

#### `_as_frequency()`

**Signature**: `_as_frequency(value)`

Return ``value`` as an allele frequency in [0, 1], or None when it is not one.

#### `_parse_population_db()`

**Signature**: `_parse_population_db(population_db)`

Read a population database, counting the entries that are not in a supported layout.

Args:
    population_db: A dict ``{'chrom:start-end[-svtype]': af}`` or an iterable of dicts
        with ``chrom``, ``start``, ``end``, ``af`` (or ``allele_frequency``) and
        optionally ``sv_type``.

Returns:
    The usable entries by chromosome, plus how many entries were skipped and why.

Raises:
    ValueError: If ``population_db`` is neither a dict nor an iterable.

#### `_build_population_lookup()`

**Signature**: `_build_population_lookup(population_db)`

Build a chromosome-keyed lookup from population database.

Entries that are not in a supported layout are left out (use
`_parse_population_db()` to count them).

Args:
    population_db: Population database in either dict or list format.

Returns:
    Dictionary mapping chromosome to sorted list of (start, end, af, sv_type) tuples.

#### `_lookup_population_frequency()`

**Signature**: `_lookup_population_frequency(pop_lookup, chrom, start, end, sv_type, match_window, min_ro)`

Look up population frequency for a variant.

Args:
    pop_lookup: Population lookup structure.
    chrom: Variant chromosome.
    start: Variant start.
    end: Variant end.
    sv_type: Variant SV type.
    match_window: Position matching window.
    min_ro: Minimum reciprocal overlap.

Returns:
    Allele frequency if a match is found, None otherwise.

### Classes

#### `FilterStats`

Statistics from a filtering operation.

Attributes:
    input_count: Number of variants before filtering.
    output_count: Number of variants passing filter.
    filtered_count: Number of variants removed.
    filter_name: Name of the filter applied.
    parameters: Filter parameters used.

**Methods**:

| Method | Purpose |
|--------|---------|
| `pass_rate()` | Fraction of variants passing the filter. |

#### `_PopulationParse`

Outcome of reading a population database (internal).

Attributes:
    lookup: Chromosome -> sorted ``(start, end, af, sv_type)`` entries that could be read.
    n_records: Number of entries in the database.
    n_skipped: Entries that could not be read.
    examples: Descriptions of the first few skipped entries.

## population.sv_population

### Functions

#### `genotype_sv_population()`

**Signature**: `genotype_sv_population(sv_calls, samples, method)`

Population-scale SV genotyping across multiple samples.

For each SV, determines the genotype (0/0, 0/1, or 1/1) in each sample
using the specified method. The depth method uses read depth signals,
while the split method uses split-read and discordant-pair evidence.

Args:
    sv_calls: List of SV call dicts, each containing:
        - ``chrom`` (str): Chromosome.
        - ``start`` (int): Start position.
        - ``end`` (int): End position.
        - ``sv_type`` (str): SV type (DEL, DUP, INV, etc.).
        - ``samples`` (dict): Mapping sample name to evidence dict with
          ``depth_ratio`` (float), ``split_reads`` (int),
          ``discordant_pairs`` (int).
    samples: List of sample names to genotype.
    method: Genotyping method, ``"depth"`` or ``"split"``.

Returns:
    Dictionary with keys:
        - ``genotype_matrix``: 2D list (n_svs x n_samples) of genotype
          integers (0, 1, or 2 for number of alt alleles).
        - ``quality_scores``: 2D list (n_svs x n_samples) of genotype
          quality scores (0-100).
        - ``allele_frequencies``: List of allele frequencies per SV.

Raises:
    ValueError: If method is unrecognized.

#### `sv_allele_frequency()`

**Signature**: `sv_allele_frequency(genotype_matrix, sample_labels)`

Compute SV allele frequencies, optionally stratified by population.

Calculates minor allele frequency (MAF), site frequency spectrum (SFS),
and optionally per-population allele frequencies.

Args:
    genotype_matrix: Genotype matrix (n_svs x n_samples) with values
        0, 1, or 2 (number of alt alleles).
    sample_labels: Optional population labels for each sample (length
        n_samples). If provided, frequencies are computed per population.

Returns:
    Dictionary with keys:
        - ``frequencies``: List of allele frequencies per SV.
        - ``maf``: List of minor allele frequencies per SV.
        - ``n_polymorphic``: Number of polymorphic SVs (MAF > 0).
        - ``sfs``: Site frequency spectrum as a list of counts.
        - ``population_frequencies``: Dict mapping population label to
          list of per-SV allele frequencies (only if sample_labels given).

Raises:
    ImportError: If numpy is not available.

#### `sv_association_test()`

**Signature**: `sv_association_test(genotypes, phenotypes, covariates)`

Test association between SV genotype and phenotype.

Performs linear regression (continuous phenotype) or logistic regression
(binary phenotype) of phenotype on genotype, optionally adjusting for
covariates.

Args:
    genotypes: Genotype values per sample (0, 1, or 2).
    phenotypes: Phenotype values per sample (continuous or binary 0/1).
    covariates: Optional covariate matrix (n_samples x n_covariates).

Returns:
    Dictionary with keys:
        - ``beta``: Effect size (regression coefficient for genotype).
        - ``se``: Standard error of beta.
        - ``p_value``: Two-sided p-value.
        - ``odds_ratio``: Odds ratio (only for binary phenotype, else None).
        - ``n_samples``: Number of samples analyzed.

Raises:
    ImportError: If numpy is not available.
    ValueError: If genotypes and phenotypes have different lengths.

#### `sv_population_structure()`

**Signature**: `sv_population_structure(genotype_matrix, n_components)`

PCA on SV genotype matrix for population structure analysis.

Centers and scales the genotype matrix, then computes the top principal
components to visualize population structure driven by structural
variant genotypes.

Args:
    genotype_matrix: Genotype matrix (n_svs x n_samples) with values
        0, 1, or 2.
    n_components: Number of principal components to compute.

Returns:
    Dictionary with keys:
        - ``pcs``: 2D array (n_samples x n_components) of PC coordinates.
        - ``eigenvalues``: List of eigenvalues for each component.
        - ``variance_explained``: List of the fraction of the total variance
          (the sum of all eigenvalues of the covariance matrix, not just the
          returned components) that each component explains. It sums to less
          than 1 when ``n_components`` is smaller than the rank of the data.
        - ``loadings``: 2D array (n_svs x n_components) of SV loadings.

Raises:
    ImportError: If numpy is not available.

#### `sv_ld_analysis()`

**Signature**: `sv_ld_analysis(genotype_matrix_sv, genotype_matrix_snp, sv_positions, snp_positions)`

Compute LD between SVs and nearby SNPs.

For each SV, computes the squared Pearson correlation (r^2) with all
SNPs within a window, identifying tag SNPs (SNPs in highest LD with
each SV).

Args:
    genotype_matrix_sv: SV genotype matrix (n_svs x n_samples), values
        0/1/2.
    genotype_matrix_snp: SNP genotype matrix (n_snps x n_samples),
        values 0/1/2.
    sv_positions: Genomic position for each SV.
    snp_positions: Genomic position for each SNP.

Returns:
    Dictionary with keys:
        - ``ld_matrix``: 2D list (n_svs x n_snps) of r^2 values.
        - ``tag_snps``: List of dicts per SV with ``snp_index``,
          ``snp_position``, ``r2``.
        - ``max_r2_per_sv``: List of maximum r^2 value for each SV.

Raises:
    ImportError: If numpy is not available.

#### `merge_sv_callsets()`

**Signature**: `merge_sv_callsets(callsets, min_overlap, max_breakpoint_distance)`

Merge SV calls across samples using reciprocal overlap and breakpoint proximity.

Compares SV calls from multiple samples/callers and merges variants that
represent the same event based on reciprocal overlap fraction and
breakpoint distance criteria.

Args:
    callsets: List of callsets, each a list of SV call dicts containing:
        - ``chrom`` (str): Chromosome.
        - ``start`` (int): Start position.
        - ``end`` (int): End position.
        - ``sv_type`` (str): SV type.
        - ``sample`` (str, optional): Source sample name.
    min_overlap: Minimum reciprocal overlap fraction (0-1) for merging.
    max_breakpoint_distance: Maximum distance between breakpoints for
        merging.

Returns:
    List of merged SV dicts, each containing:
        - ``chrom``, ``start``, ``end``, ``sv_type``: Consensus call.
        - ``n_samples``: Number of samples supporting this SV.
        - ``samples``: List of sample names.
        - ``support_count``: Total number of supporting calls.
        - ``merged_calls``: List of original calls merged into this entry.

#### `_genotype_by_depth()`

**Signature**: `_genotype_by_depth(depth_ratio, sv_type)`

Genotype an SV based on read depth ratio.

For deletions, depth_ratio < 0.75 suggests het, < 0.25 suggests hom.
For duplications, depth_ratio > 1.25 suggests het, > 1.75 suggests hom.

Args:
    depth_ratio: Observed/expected read depth ratio.
    sv_type: SV type string.

Returns:
    Tuple of (genotype, quality) where genotype is 0/1/2.

#### `_genotype_by_split()`

**Signature**: `_genotype_by_split(split_reads, discordant_pairs, total_reads)`

Genotype an SV based on split-read and discordant-pair evidence.

Uses the fraction of supporting reads to estimate genotype.

Args:
    split_reads: Number of split reads supporting the SV.
    discordant_pairs: Number of discordant pairs.
    total_reads: Total reads in the region.

Returns:
    Tuple of (genotype, quality).

#### `_total_variance()`

**Signature**: `_total_variance(eigenvalues)`

Total variance: the sum of all eigenvalues of a covariance matrix.

Eigenvalues of a covariance matrix cannot be negative; the tiny negative ones
that rounding produces are counted as zero.

Args:
    eigenvalues: Every eigenvalue of the matrix, not only the leading ones.

Returns:
    The sum of the non-negative parts of the eigenvalues.

#### `_linear_regression()`

**Signature**: `_linear_regression(design, y)`

Ordinary least squares linear regression.

Args:
    design: Design matrix (n x p) including intercept.
    y: Response vector (n,).

Returns:
    Tuple of (coefficients, standard_errors, p_values).

#### `_logistic_regression()`

**Signature**: `_logistic_regression(design, y, max_iter)`

Iteratively reweighted least squares logistic regression.

Args:
    design: Design matrix (n x p) including intercept.
    y: Binary response vector (n,).
    max_iter: Maximum IRLS iterations.

Returns:
    Tuple of (coefficients, standard_errors, p_values).

#### `_normal_cdf_np()`

**Signature**: `_normal_cdf_np(x)`

Standard normal CDF approximation.

#### `_should_merge()`

**Signature**: `_should_merge(call_a, call_b, min_overlap, max_bp_dist)`

Determine if two SV calls should be merged.

Checks chromosome, SV type, reciprocal overlap, and breakpoint distance.

Args:
    call_a: First SV call dict.
    call_b: Second SV call dict.
    min_overlap: Minimum reciprocal overlap fraction.
    max_bp_dist: Maximum breakpoint distance.

Returns:
    True if calls should be merged.

#### `_build_consensus()`

**Signature**: `_build_consensus(cluster)`

Build a consensus SV call from a cluster of merged calls.

Uses median breakpoints and majority SV type.

Args:
    cluster: List of SV call dicts to merge.

Returns:
    Consensus SV call dict.

## visualization.plots

### Functions

#### `_check_matplotlib()`

**Signature**: `_check_matplotlib()`

Raise an import error if matplotlib is not available.

#### `plot_circos()`

**Signature**: `plot_circos(variants, chromosomes, output_path, title, figsize, show_labels)`

Generate a Circos-style genome-wide structural variant plot.

Draws chromosomes as arcs around a circle, with SVs represented as
colored chords connecting breakpoint pairs. Intra-chromosomal SVs
are drawn as arcs within the chromosome. Colors indicate SV type.

A chord runs from ``chrom:start`` to the second breakpoint, which is
``chrom2:pos2`` for inter-chromosomal calls (as produced by
``call_structural_variants``, e.g. via ``dataclasses.asdict``) and
``chrom:end`` for intra-chromosomal calls. A hand-built inter-chromosomal
dictionary without ``pos2`` (``chrom2`` set to another chromosome) is read
the way VCF-style callers write it: ``end`` is the position on ``chrom2``.

Args:
    variants: List of variant dictionaries with:
        - 'chrom': Chromosome
        - 'start': Start position (first breakpoint)
        - 'end': End position (second breakpoint of an intra-chromosomal SV)
        - 'sv_type': SV type string
        - 'chrom2': Chromosome of the second breakpoint (translocations, optional)
        - 'pos2': Position of the second breakpoint on 'chrom2'
          (translocations, optional; None or absent for intra-chromosomal SVs)
    chromosomes: Dictionary mapping chromosome names to sizes (in bp).
    output_path: Path to save the plot image.
    title: Plot title.
    figsize: Figure size in inches.
    show_labels: Whether to show chromosome labels.

Returns:
    matplotlib Figure object.

#### `plot_coverage_track()`

**Signature**: `plot_coverage_track(coverage, variants, region, output_path, title, figsize, bin_size)`

Plot read depth coverage with structural variant overlay.

Shows a coverage track (read depth) for a genomic region with
colored rectangles indicating structural variant positions and types.

Args:
    coverage: Array of coverage values (one per bin).
    variants: List of variant dictionaries (same format as other functions).
    region: Genomic region as dict with 'chrom', 'start', 'end' or
        tuple (chrom, start, end).
    output_path: Path to save the plot. If empty, plot is not saved.
    title: Plot title. Auto-generated from region if empty.
    figsize: Figure size in inches.
    bin_size: Size of each coverage bin in base pairs.

Returns:
    matplotlib Figure object.

#### `plot_sv_size_distribution()`

**Signature**: `plot_sv_size_distribution(variants, sv_type, output_path, title, figsize, log_scale)`

Plot size distribution of structural variants.

Creates a histogram showing the distribution of SV sizes, optionally
filtered by SV type. Uses log-scale x-axis by default since SV sizes
span several orders of magnitude.

Args:
    variants: List of variant dictionaries with 'size' or 'start'/'end',
        and 'sv_type' keys.
    sv_type: If specified, only plot variants of this type.
    output_path: Path to save the plot. If empty, plot is not saved.
    title: Plot title.
    figsize: Figure size in inches.
    log_scale: Whether to use log scale for x-axis (default True).

Returns:
    matplotlib Figure object.

#### `plot_sv_type_summary()`

**Signature**: `plot_sv_type_summary(variants, output_path, title, figsize)`

Plot summary of SV types as bar chart and pie chart.

Creates a two-panel figure showing (left) a bar chart of SV counts
by type and (right) a pie chart of the same data.

Args:
    variants: List of variant dictionaries with 'sv_type' keys.
    output_path: Path to save the plot.
    title: Plot title.
    figsize: Figure size.

Returns:
    matplotlib Figure object.

#### `plot_breakpoint_detail()`

**Signature**: `plot_breakpoint_detail(variant, reads, flanking, output_path, figsize)`

Plot detailed view of a structural variant breakpoint.

Creates a multi-panel figure showing:
- Top: Read alignment pileup around the breakpoint
- Middle: Coverage profile
- Bottom: Split/discordant read evidence

Args:
    variant: Variant dictionary with 'chrom', 'start', 'end', 'sv_type'.
    reads: List of read alignment dictionaries.
    flanking: Number of base pairs to show flanking the breakpoint.
    output_path: Path to save the plot.
    figsize: Figure size.

Returns:
    matplotlib Figure object.

#### `plot_cnv_profile()`

**Signature**: `plot_cnv_profile(segments, chromosomes, output_path, title, figsize, ploidy)`

Plot genome-wide CNV profile.

Displays log2 ratio values across the genome with colored segments
indicating copy number states. Chromosomes are laid out linearly
with alternating background shading.

Args:
    segments: List of segment dictionaries with:
        - 'chrom': Chromosome
        - 'start': Start position
        - 'end': End position
        - 'mean_log2ratio': Mean log2 ratio
        - 'state': CNV state (DEL, DUP, NEUTRAL, etc.)
    chromosomes: Dictionary mapping chromosome names to sizes.
    output_path: Path to save the plot.
    title: Plot title.
    figsize: Figure size.
    ploidy: Expected ploidy for reference line.

Returns:
    matplotlib Figure object.

#### `_chrom_sort_key()`

**Signature**: `_chrom_sort_key(chrom)`

Sort key for chromosome names (numeric first, then alphabetic).

Args:
    chrom: Chromosome name.

Returns:
    Sort key tuple.

#### `_pos_to_angle()`

**Signature**: `_pos_to_angle(chrom, pos, chrom_angles, chromosomes)`

Convert a genomic position to a polar angle.

Args:
    chrom: Chromosome.
    pos: Genomic position.
    chrom_angles: Chromosome angular ranges.
    chromosomes: Chromosome sizes.

Returns:
    Polar angle in radians, or None if chromosome not found.

#### `_draw_chromosome_arcs()`

**Signature**: `_draw_chromosome_arcs(ax, chrom_angles, outer_r, inner_r, chrom_order, show_labels)`

Draw chromosome arcs around the Circos plot.

Args:
    ax: Polar axes.
    chrom_angles: Angular ranges per chromosome.
    outer_r: Outer radius.
    inner_r: Inner radius.
    chrom_order: Ordered chromosome names.
    show_labels: Whether to show labels.

#### `_draw_chord()`

**Signature**: `_draw_chord(ax, angle1, angle2, radius, color, alpha)`

Draw a chord between two angles in polar coordinates.

Uses a quadratic Bezier curve through the center for visual clarity.

Args:
    ax: Polar axes.
    angle1: First angle in radians.
    angle2: Second angle in radians.
    radius: Radius at which the chord endpoints are placed.
    color: Color of the chord.
    alpha: Transparency.
