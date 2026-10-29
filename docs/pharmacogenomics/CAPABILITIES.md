# Pharmacogenomics Module Capabilities

This reference mirrors the checked-out implementation's signatures and docstrings. Bundled annotation tables are offline data; see the domain guides for their limits.

## Quick Reference

### Functions by Submodule

| Function | Submodule | Purpose |
|----------|-----------|---------|
| `_haplotype_allele()` | `alleles.diplotype` | Return the one allele a haplotype carries, and whether more than one was called for it. |
| `_normalize_allele_name()` | `alleles.diplotype` | Normalize an allele name for scoring-table lookup. |
| `calculate_activity_score()` | `alleles.diplotype` | Calculate activity score for an existing Diplotype. |
| `calculate_activity_score_from_alleles()` | `alleles.diplotype` | Calculate the combined activity score for two alleles. |
| `determine_diplotype()` | `alleles.diplotype` | Construct a diplotype from two star alleles with activity scoring. |
| `phased_diplotype()` | `alleles.diplotype` | Determine diplotype using phasing information. |
| `resolve_ambiguous_diplotypes()` | `alleles.diplotype` | Resolve ambiguous diplotype calls when phase is unknown. |
| `uncharacterized_alleles()` | `alleles.diplotype` | Return the alleles of a diplotype that have no activity value in the scoring table. |
| `_assess_diplotype()` | `alleles.phenotype` | Activity score, phenotype and unscorable alleles of a diplotype of a classifiable gene. |
| `_band_phenotype()` | `alleles.phenotype` | Return the phenotype of the band containing ``activity_score`` (>= 0). |
| `_canonical_gene()` | `alleles.phenotype` | Return the gene symbol as the tables spell it (stripped, upper case). |
| `_function_class()` | `alleles.phenotype` | Return the CPIC function class implied by one allele's activity value. |
| `_generate_clinical_note()` | `alleles.phenotype` | Generate a brief clinical significance note for a phenotype. |
| `_phenotype_bands()` | `alleles.phenotype` | Return a gene's score bands. |
| `_phenotype_for_alleles()` | `alleles.phenotype` | Phenotype of a diplotype from its two alleles' activity values. |
| `_unsupported_gene_error()` | `alleles.phenotype` | Build the ValueError for a gene without a CPIC phenotype definition. |
| `classify_phenotype()` | `alleles.phenotype` | Full diplotype-to-phenotype classification pipeline. |
| `get_phenotype_thresholds()` | `alleles.phenotype` | Get gene-specific phenotype threshold definitions. |
| `population_phenotype_frequencies()` | `alleles.phenotype` | Estimate expected phenotype distribution from allele frequencies. |
| `predict_metabolizer_status()` | `alleles.phenotype` | Map an activity score to a metabolizer phenotype category. |
| `supported_phenotype_genes()` | `alleles.phenotype` | Return the genes that have a CPIC phenotype definition, sorted. |
| `_load_allele_definitions_from_file()` | `alleles.star_allele` | Load allele definitions from an external file. |
| `_load_builtin_definitions()` | `alleles.star_allele` | Load allele definitions from built-in CPIC-derived tables. |
| `_parse_activity_value()` | `alleles.star_allele` | Parse the activity value of an allele definition read from a file. |
| `_select_gene_definitions()` | `alleles.star_allele` | Pick one gene's allele definitions out of a parsed JSON definition file. |
| `_validated_copy_number()` | `alleles.star_allele` | Return ``copy_number`` as an int, or raise. |
| `call_star_alleles()` | `alleles.star_allele` | Call star alleles from observed variants. |
| `detect_novel_alleles()` | `alleles.star_allele` | Detect potentially novel allele combinations. |
| `handle_cyp2d6_cnv()` | `alleles.star_allele` | Handle CYP2D6 copy number variation for star allele calling. |
| `load_allele_definitions()` | `alleles.star_allele` | Load allele definition tables for a pharmacogene. |
| `match_allele_definition()` | `alleles.star_allele` | Match observed variants against a single allele definition. |
| `_load_guidelines_from_file()` | `annotations.cpic` | Load CPIC guidelines from an external file. |
| `_normalize_phenotype_string()` | `annotations.cpic` | Normalize phenotype string for matching. |
| `get_dosing_recommendation()` | `annotations.cpic` | Get dosing recommendation for a drug based on metabolizer phenotype. |
| `list_actionable_genes()` | `annotations.cpic` | List CPIC Level A and B (actionable) gene-drug pairs. |
| `load_cpic_guidelines()` | `annotations.cpic` | Load CPIC guideline data. |
| `lookup_drug_gene()` | `annotations.cpic` | Look up CPIC guideline for a specific drug-gene pair. |
| `parse_cpic_allele_definitions()` | `annotations.cpic` | Parse CPIC allele definition tables. |
| `_contains_pgx_terms()` | `annotations.drug_labels` | Check if text contains pharmacogenomic-relevant terms. |
| `_parse_label_text()` | `annotations.drug_labels` | Parse unstructured drug label text into sections. |
| `classify_label_type()` | `annotations.drug_labels` | Classify the strength of PGx labeling for a drug. |
| `extract_biomarker_info()` | `annotations.drug_labels` | Extract biomarker testing requirements from parsed label data. |
| `parse_drug_label()` | `annotations.drug_labels` | Parse an FDA drug label file for pharmacogenomic information. |
| `search_labels_by_gene()` | `annotations.drug_labels` | Find drug labels mentioning a specific gene or biomarker. |
| `get_evidence_level()` | `annotations.pharmgkb` | Extract and interpret the evidence level from a PharmGKB annotation. |
| `get_variant_annotations()` | `annotations.pharmgkb` | Get variant-specific pharmacogenomic annotations. |
| `parse_clinical_annotations()` | `annotations.pharmgkb` | Parse PharmGKB clinical annotation data from file. |
| `query_pharmgkb_annotations()` | `annotations.pharmgkb` | Query PharmGKB clinical annotations by gene, drug, and/or variant. |
| `search_drug_pathways()` | `annotations.pharmgkb` | Get pharmacokinetic/pharmacodynamic pathway information for a drug. |
| `_abbreviate_phenotype()` | `clinical.drug_interaction` | Convert a phenotype label to its abbreviation. |
| `analyze_drug_gene_interactions()` | `clinical.drug_interaction` | Analyze drug-gene interactions for a set of drugs and genotypes. |
| `calculate_interaction_severity()` | `clinical.drug_interaction` | Calculate overall severity for a set of drug-gene interactions. |
| `check_contraindications()` | `clinical.drug_interaction` | Check if a drug is contraindicated for a given metabolizer phenotype. |
| `polypharmacy_analysis()` | `clinical.drug_interaction` | Perform comprehensive polypharmacy pharmacogenomic assessment. |
| `suggest_alternatives()` | `clinical.drug_interaction` | Suggest alternative drugs based on pharmacogenomic phenotype. |
| `_benign_classification()` | `clinical.pathogenicity` | Benign or Likely Benign call from the counts of BA, BS and BP criteria met, else None. |
| `_pathogenic_classification()` | `clinical.pathogenicity` | Pathogenic or Likely Pathogenic call from the counts of PVS, PS, PM and PP criteria met, else None. |
| `_population_frequencies()` | `clinical.pathogenicity` | Return a variant's allele frequencies by population; empty when it carries none. |
| `aggregate_evidence()` | `clinical.pathogenicity` | Combine ACMG criteria into a final 5-tier classification. |
| `apply_acmg_criteria()` | `clinical.pathogenicity` | Evaluate individual ACMG criteria from variant data. |
| `check_gnomad_frequency()` | `clinical.pathogenicity` | Check variant population frequency against gnomAD thresholds. |
| `classify_variant_acmg()` | `clinical.pathogenicity` | Classify a variant using the ACMG 5-tier system. |
| `query_clinvar()` | `clinical.pathogenicity` | Look up ClinVar classification for a variant. |
| `_canonical_phenotype()` | `clinical.reporting` | Return the MetabolizerPhenotype for a recognised label ("poor", "PM", ...), else the value as given. |
| `_format_html_report()` | `clinical.reporting` | Format report as HTML. |
| `_format_json_report()` | `clinical.reporting` | Format report as JSON string. |
| `_format_text_report()` | `clinical.reporting` | Format report as plain text. |
| `_get_abbrev()` | `clinical.reporting` | Get phenotype abbreviation from various input types. |
| `_html_escape()` | `clinical.reporting` | Basic HTML escaping for report content. |
| `add_disclaimer()` | `clinical.reporting` | Add or update the clinical disclaimer in a report. |
| `export_report()` | `clinical.reporting` | Export a clinical report to the specified format. |
| `format_recommendation()` | `clinical.reporting` | Format a single drug-gene recommendation for the clinical report. |
| `generate_clinical_report()` | `clinical.reporting` | Generate a comprehensive clinical pharmacogenomic report. |
| `generate_summary_table()` | `clinical.reporting` | Generate a summary table of all pharmacogenomic findings. |
| `_find_shared_cyp_pathways()` | `interaction.drug_interactions` | Find CYP enzymes shared by two drugs as substrates. |
| `cyp_inhibition_prediction()` | `interaction.drug_interactions` | Predict CYP enzyme inhibition and induction potential of a drug. |
| `default_interaction_database()` | `interaction.drug_interactions` | Return the built-in drug-drug interaction database. |
| `polypharmacy_risk()` | `interaction.drug_interactions` | Assess polypharmacy risk from multiple concurrent medications. |
| `predict_drug_interaction()` | `interaction.drug_interactions` | Predict drug-drug interaction between two medications. |
| `_phenotype_token()` | `metabolism.metabolizer_status` | Return the lowercase token for a phenotype (``"poor"``, ``"ultrarapid"``, ...). |
| `_split_diplotype()` | `metabolism.metabolizer_status` | Split a ``"*1/*4"`` diplotype string into its two allele names. |
| `_status_token()` | `metabolism.metabolizer_status` | Return the lowercase phenotype token for any accepted phenotype label. |
| `classify_metabolizer()` | `metabolism.metabolizer_status` | Classify metabolizer phenotype from an activity score. |
| `compute_activity_score()` | `metabolism.metabolizer_status` | Compute CPIC-style activity score from a diplotype string. |
| `default_allele_function_table()` | `metabolism.metabolizer_status` | Return the built-in allele function assignment table. |
| `dose_adjustment()` | `metabolism.metabolizer_status` | Recommend dose adjustment based on metabolizer status and drug. |
| `predict_metabolizer_status()` | `metabolism.metabolizer_status` | Predict metabolizer phenotype from genotype information. |
| `_ensure_matplotlib()` | `visualization.plots` | Raise ImportError if matplotlib is not available. |
| `_ensure_output_dir()` | `visualization.plots` | Ensure the output directory exists and return the Path. |
| `_format_bar_value()` | `visualization.plots` | Label for a bar: a whole number as an integer, anything else to four significant digits. |
| `_threshold_bands()` | `visualization.plots` | Return the phenotype score bands to mark on an activity-score plot of ``gene``. |
| `plot_acmg_criteria()` | `visualization.plots` | Plot ACMG criteria evaluation summary. |
| `plot_activity_score_distribution()` | `visualization.plots` | Plot activity score distribution as a histogram. |
| `plot_allele_frequencies()` | `visualization.plots` | Plot allele frequency distribution for a pharmacogene. |
| `plot_drug_response_heatmap()` | `visualization.plots` | Plot drug-gene response heatmap. |
| `plot_metabolizer_status()` | `visualization.plots` | Plot metabolizer status distribution as a bar chart. |
| `plot_population_comparison()` | `visualization.plots` | Plot cross-population allele frequency comparison. |

## alleles.diplotype

### Functions

#### `_haplotype_allele()`

**Signature**: `_haplotype_allele(called: list[StarAllele], gene: str, haplotype: str) -> tuple[StarAllele, bool]`

Return the one allele a haplotype carries, and whether more than one was called for it.

A haplotype is a single star allele. When the variants on it match several definitions
(for example a loss-of-function marker beside a normal-function one), the most impaired
allele is taken: a no-function variant is not offset by a normal-function marker on the same
haplotype, and the calls come back ordered by name, so taking the first would pick "*2" over
"*4". The ambiguity is logged and returned so the caller can lower its confidence.

Args:
    called: Alleles called from the haplotype's variants (:func:`call_star_alleles`)
    gene: Gene symbol
    haplotype: Label of the haplotype, for the log message

Returns:
    ``(allele, ambiguous)``. An empty ``called`` gives the reference ``*1``.

#### `_normalize_allele_name()`

**Signature**: `_normalize_allele_name(name: str, gene: str, scoring_table: Mapping[str, float] | None=None) -> str`

Normalize an allele name for scoring-table lookup.

Star allele names are prefixed with ``*`` when missing. Names already
present in the scoring table are returned unchanged, so non-star allele
designations (e.g. DPYD's "c.2846A>T" and "HapB3") are preserved.

#### `calculate_activity_score()`

**Signature**: `calculate_activity_score(diplotype: Diplotype, gene: str, scoring_table: dict[str, float] | None=None) -> float | None`

Calculate activity score for an existing Diplotype.

Args:
    diplotype: Diplotype object
    gene: Gene symbol
    scoring_table: Optional custom scoring table

Returns:
    Combined activity score, also updates the Diplotype object; ``None`` (stored on the
    diplotype too) if an allele has no activity value in the table, see
    :func:`calculate_activity_score_from_alleles`.

#### `calculate_activity_score_from_alleles()`

**Signature**: `calculate_activity_score_from_alleles(allele1_name: str, allele2_name: str, gene: str, scoring_table: dict[str, float] | None=None) -> float | None`

Calculate the combined activity score for two alleles.

Args:
    allele1_name: First allele name
    allele2_name: Second allele name
    gene: Gene symbol
    scoring_table: Optional custom scoring table

Returns:
    Combined activity score (sum of individual allele scores), or ``None`` if either allele
    has no activity value in the table. An uncharacterized allele is never assumed to be
    normal function (1.0); a warning names it, and the diplotype is indeterminate.

#### `determine_diplotype()`

**Signature**: `determine_diplotype(allele1: str | StarAllele, allele2: str | StarAllele, gene: str, scoring_table: dict[str, float] | None=None) -> Diplotype`

Construct a diplotype from two star alleles with activity scoring.

Creates a Diplotype object from two allele designations, looking up
activity values from the gene-specific scoring table to compute the
combined activity score.

Args:
    allele1: First allele name (str) or StarAllele object
    allele2: Second allele name (str) or StarAllele object
    gene: Gene symbol (e.g., "CYP2D6")
    scoring_table: Optional custom activity score mapping (allele_name -> score).
        If None, uses built-in CPIC-derived tables.

Returns:
    Diplotype object with its computed activity score. The score is ``None`` when
    either allele has no activity value in the scoring table (an uncharacterized
    allele; a warning names it): such a diplotype is indeterminate, not normal.

#### `phased_diplotype()`

**Signature**: `phased_diplotype(variants: set[str] | frozenset[str], phase_data: dict[str, int], gene: str, allele_definitions: list[StarAllele] | None=None) -> Diplotype`

Determine diplotype using phasing information.

When phase information is available (e.g., from long reads, family data,
or statistical phasing), this function assigns variants to maternal and
paternal haplotypes before calling star alleles on each independently.

A haplotype carries one allele. If its variants match more than one definition, the allele
with the lowest activity value is taken (a loss-of-function variant is not offset by a
normal-function marker on the same haplotype), a warning names the candidates, and the
diplotype's confidence is ``"moderate"`` instead of ``"high"``.

Args:
    variants: Set of all observed variant identifiers
    phase_data: Mapping of variant_id -> phase (0 for haplotype A, 1 for haplotype B).
        Variants not in phase_data are treated as unphased.
    gene: Gene symbol (e.g., "CYP2D6")
    allele_definitions: Optional pre-loaded allele definitions

Returns:
    Diplotype with phased=True and alleles called from each haplotype

#### `resolve_ambiguous_diplotypes()`

**Signature**: `resolve_ambiguous_diplotypes(possible_diplotypes: list[Diplotype]) -> Diplotype`

Resolve ambiguous diplotype calls when phase is unknown.

When genotype data is unphased, multiple diplotype configurations may be
possible. This function selects the most clinically relevant assignment
using the following priority rules (per CPIC recommendations):

1. Prefer diplotypes with known activity scores for both alleles
2. Among scored diplotypes, prefer the one with the highest activity score
   (conservative: assume best-case metabolism)

A diplotype that has not been scored (``activity_score`` is ``None``) ranks by the score
its alleles imply; if one of its alleles has no activity value it has no score and ranks
after every diplotype with as many known alleles.

Args:
    possible_diplotypes: List of possible Diplotype objects

Returns:
    The selected most likely Diplotype

Raises:
    ValueError: If possible_diplotypes is empty

#### `uncharacterized_alleles()`

**Signature**: `uncharacterized_alleles(allele1_name: str, allele2_name: str, gene: str, scoring_table: Mapping[str, float] | None=None) -> list[str]`

Return the alleles of a diplotype that have no activity value in the scoring table.

Names are normalized as in :func:`calculate_activity_score_from_alleles`, listed in the
order given and without repeats (a homozygous diplotype names its allele once).

Args:
    allele1_name: First allele name
    allele2_name: Second allele name
    gene: Gene symbol
    scoring_table: Optional custom scoring table. If None, uses the built-in table for the
        gene (an empty table for a gene without one, so every allele is uncharacterized).

Returns:
    The normalized names of the uncharacterized alleles; an empty list when both are scored.

### Classes

#### `Diplotype`

Represents a pharmacogene diplotype (pair of star alleles).

Attributes:
    allele1: First star allele (typically the lower-numbered or reference)
    allele2: Second star allele
    gene: Gene symbol
    activity_score: Combined activity score for the diplotype, or ``None`` (the default)
        when it has not been scored. ``None`` is not 0.0: it never reads as a poor
        metabolizer. Score it with :func:`calculate_activity_score`, or build the
        diplotype with :func:`determine_diplotype`, which does; ``classify_phenotype``
        computes the score from the alleles of an unscored diplotype. ``None`` also
        results when an allele has no activity value in the scoring table.
    phased: Whether the diplotype assignment is phased (True) or inferred (False)
    confidence: Confidence in the diplotype call ("high", "moderate", "low")

##### `diplotype_string`

**Signature**: `diplotype_string: str`

Return formatted diplotype string (e.g., '*1/*2').

##### `is_homozygous`

**Signature**: `is_homozygous: bool`

Return True if both alleles are the same.

##### `is_reference`

**Signature**: `is_reference: bool`

Return True if both alleles are reference (*1/*1).

## alleles.phenotype

### Functions

#### `_assess_diplotype()`

**Signature**: `_assess_diplotype(dip: Diplotype, gene: str) -> tuple[float | None, MetabolizerPhenotype, list[str]]`

Activity score, phenotype and unscorable alleles of a diplotype of a classifiable gene.

An unscored diplotype (``activity_score`` None) is scored from its alleles' values in the
built-in table; with an allele missing from the table it has no score and is indeterminate.
A scored diplotype uses the two alleles' own activity values (which the genes CPIC defines
by allele function need) when they reproduce its stored score. One scored with a custom
scoring table, or containing alleles missing from the built-in table, is classified from its
stored score alone: the score was supplied, so the alleles are not treated as uncharacterized.

Returns:
    ``(activity_score, phenotype, unknown_alleles)``

#### `_band_phenotype()`

**Signature**: `_band_phenotype(bands: tuple[tuple[float, bool, MetabolizerPhenotype], ...], activity_score: float) -> MetabolizerPhenotype`

Return the phenotype of the band containing ``activity_score`` (>= 0).

#### `_canonical_gene()`

**Signature**: `_canonical_gene(gene: str) -> str`

Return the gene symbol as the tables spell it (stripped, upper case).

#### `_function_class()`

**Signature**: `_function_class(activity: float) -> _AlleleFunction`

Return the CPIC function class implied by one allele's activity value.

#### `_generate_clinical_note()`

**Signature**: `_generate_clinical_note(phenotype: MetabolizerPhenotype, gene: str, unknown_alleles: Sequence[str]=()) -> str`

Generate a brief clinical significance note for a phenotype.

Args:
    phenotype: Metabolizer phenotype
    gene: Gene symbol
    unknown_alleles: Alleles with no activity value, named in an indeterminate note

Returns:
    Clinical significance string

#### `_phenotype_bands()`

**Signature**: `_phenotype_bands(gene_upper: str) -> tuple[tuple[float, bool, MetabolizerPhenotype], ...]`

Return a gene's score bands.

Raises:
    ValueError: If CPIC defines no phenotype for the gene. No other gene's bands stand in.

#### `_phenotype_for_alleles()`

**Signature**: `_phenotype_for_alleles(gene: str, activity1: float, activity2: float) -> MetabolizerPhenotype`

Phenotype of a diplotype from its two alleles' activity values.

Genes CPIC defines by allele function are translated from the pair of function classes;
all others from the summed activity score.

Raises:
    ValueError: If a value is negative or NaN, the gene has no CPIC phenotype definition, or
        the gene has no phenotype for that pair of allele functions (an increased-function
        SLCO1B1 or UGT1A1 allele).

#### `_unsupported_gene_error()`

**Signature**: `_unsupported_gene_error(gene_upper: str) -> ValueError`

Build the ValueError for a gene without a CPIC phenotype definition.

#### `classify_phenotype()`

**Signature**: `classify_phenotype(diplotype: Diplotype | str, gene: str) -> dict[str, Any]`

Full diplotype-to-phenotype classification pipeline.

Takes a diplotype (object or string like "*1/*4") and gene, computes the
activity score, and maps it to a metabolizer phenotype. Genes CPIC defines by
allele function (CYP2C19, SLCO1B1, UGT1A1) are classified from the function of the
two alleles, the rest from the activity score.

A ``Diplotype`` whose ``activity_score`` is ``None`` (not scored) is scored here from
its alleles; it is never read as a score of 0.0. If an allele has no activity value in
the built-in table (an uncharacterized allele), no score or phenotype can be assigned:
the result is indeterminate, with ``activity_score`` ``None`` and the allele named, as
CPIC calls such a diplotype indeterminate. It is never scored as normal function. A
``Diplotype`` that carries a score (for example one built with a custom
``scoring_table``) is classified from that score.

Args:
    diplotype: Diplotype object or string in format "*X/*Y"
    gene: Gene symbol

Returns:
    Dictionary with:
        - "gene": Gene symbol
        - "diplotype": Diplotype string
        - "activity_score": Computed activity score, or None when it cannot be computed
        - "phenotype": MetabolizerPhenotype value
        - "phenotype_abbreviation": Short form (PM, IM, NM, RM, UM, IND)
        - "clinical_significance": Brief clinical note
        - "unknown_alleles": Alleles with no activity value (empty unless the phenotype is
          indeterminate for that reason)

Raises:
    ValueError: If the gene has no CPIC phenotype definition (the message names the gene and
        lists the supported genes; CYP1A2 can be scored but not classified), the string is not
        ``*X/*Y``, or a ``Diplotype`` object belongs to a different gene.

#### `get_phenotype_thresholds()`

**Signature**: `get_phenotype_thresholds(gene: str) -> list[dict[str, Any]]`

Get gene-specific phenotype threshold definitions.

These are the activity-score bands :func:`predict_metabolizer_status` applies. For
CYP2C19, SLCO1B1 and UGT1A1, which :func:`classify_phenotype` translates from allele
function, they are the equivalent bands on the score scale.

Args:
    gene: Gene symbol (e.g., "CYP2D6")

Returns:
    List of threshold dictionaries, lowest activity first, with keys:
        - "lower_bound": Lower activity score boundary
        - "upper_bound": Upper activity score boundary (infinity for the last band)
        - "lower_inclusive": Whether a score equal to the lower bound is in the band
        - "upper_inclusive": Whether a score equal to the upper bound is in the band
        - "phenotype": MetabolizerPhenotype value
        - "abbreviation": Phenotype abbreviation (PM, IM, NM, RM, UM)
        - "gene_specific": Always True. It used to be False for the default (CYP2D6) bands
          returned for a gene with none of its own; such a gene now raises. The key is kept
          so callers written before that still work.

Raises:
    ValueError: If the gene has no CPIC phenotype definition; the message names the gene
        and lists the supported genes

#### `population_phenotype_frequencies()`

**Signature**: `population_phenotype_frequencies(allele_frequencies: dict[str, float], gene: str) -> dict[str, float]`

Estimate expected phenotype distribution from allele frequencies.

Calculates the expected frequency of each metabolizer phenotype in a
population given the allele frequency spectrum, assuming Hardy-Weinberg
equilibrium. For each possible diplotype (all combinations of alleles),
the diplotype frequency is 2*p*q for heterozygotes and p^2 for homozygotes.
Each diplotype is assigned its phenotype exactly as :func:`classify_phenotype`
does, so a diplotype with an allele that has no activity value (an uncharacterized
allele) counts as indeterminate ("IND"), never as normal.

Args:
    allele_frequencies: Mapping of allele name -> population frequency.
        Frequencies should sum to approximately 1.0. Names that normalize to the same
        allele ("4" and "*4") are added together.
    gene: Gene symbol

Returns:
    Dictionary mapping phenotype abbreviation -> expected frequency, without the
    categories that have none ("IND" appears only when an uncharacterized allele has
    a nonzero frequency).
    Example: {"PM": 0.05, "IM": 0.20, "NM": 0.65, "RM": 0.08, "UM": 0.02}

Raises:
    ValueError: If the gene has no CPIC phenotype definition (the message names the gene
        and lists the supported genes), a frequency is negative or NaN, or no frequency
        is positive.

#### `predict_metabolizer_status()`

**Signature**: `predict_metabolizer_status(activity_score: float, gene: str) -> MetabolizerPhenotype`

Map an activity score to a metabolizer phenotype category.

Uses the gene-specific bands in ``_PHENOTYPE_THRESHOLDS`` (CPIC definitions, cited
there). The bands cover every score >= 0 without gaps. For example, the CYP2D6 bands:
    AS = 0 -> PM, 0 < AS < 1.25 -> IM, 1.25 <= AS <= 2.25 -> NM, AS > 2.25 -> UM

A gene CPIC defines no phenotype for has no bands, and no other gene's bands stand in
for it: it raises (see :func:`supported_phenotype_genes`).

A score is a sum of allele values, so for CYP2C19 it cannot always identify the
alleles behind it (a decreased-function allele makes *9/*17 sum to the same 2.0 as
*1/*1). Use :func:`classify_phenotype` when the diplotype is known.

Args:
    activity_score: Combined activity score (typically sum of two allele scores)
    gene: Gene symbol (e.g., "CYP2D6")

Returns:
    MetabolizerPhenotype enum value

Raises:
    ValueError: If the activity score is negative, NaN or None (an unscored diplotype has
        no score to classify), or the gene has no CPIC phenotype definition; the message
        names the gene and lists the supported genes

#### `supported_phenotype_genes()`

**Signature**: `supported_phenotype_genes() -> list[str]`

Return the genes that have a CPIC phenotype definition, sorted.

These are the only genes the functions of this module classify; every other gene raises
``ValueError``. A gene can have activity values and still not be listed: CYP1A2 can be
scored (``diplotype.determine_diplotype``) but CPIC defines no phenotype for it.

### Classes

#### `MetabolizerPhenotype`

Standardized metabolizer phenotype categories.

Values follow CPIC terminology for consistent clinical communication.
Inherits from str for JSON serialization compatibility.

##### `abbreviation`

**Signature**: `abbreviation: str`

Return standard abbreviation (PM, IM, NM, RM, UM).

##### `from_string()`

**Signature**: `from_string(label: str) -> MetabolizerPhenotype`

Parse a phenotype label written in any of the forms used across the package.

Accepts the full CPIC name ("Poor Metabolizer"), its abbreviation ("PM"), the bare
lowercase token used by ``metabolism.metabolizer_status`` ("poor", "ultrarapid"),
variants such as "ultra-rapid", and the pre-2017 CPIC synonym "extensive" for
normal, all case-insensitively. A MetabolizerPhenotype is returned unchanged.

Raises:
    ValueError: If the label is not a recognised phenotype.

#### `_AlleleFunction`

CPIC allele function classes, ordered from least to most active.

## alleles.star_allele

### Functions

#### `_load_allele_definitions_from_file()`

**Signature**: `_load_allele_definitions_from_file(gene: str, filepath: str | Path) -> list[StarAllele]`

Load allele definitions from an external file.

Supports TSV format with columns: allele, defining_variants (semicolon-separated),
function, activity_value. Also supports JSON format with the same structure as
_BUILTIN_ALLELE_DEFINITIONS, either one gene's table or a mapping of gene symbol to table.

Every allele needs an ``activity_value`` (a number >= 0): none is assumed.

Args:
    gene: Uppercased gene symbol
    filepath: Path to the definition file

Returns:
    List of StarAllele objects

Raises:
    FileNotFoundError: If the file does not exist.
    ValueError: If an allele has no allele name or no valid activity value, or a JSON file
        keyed by gene does not contain ``gene``.

#### `_load_builtin_definitions()`

**Signature**: `_load_builtin_definitions(gene: str) -> list[StarAllele]`

Load allele definitions from built-in CPIC-derived tables.

Args:
    gene: Uppercased gene symbol

Returns:
    List of StarAllele objects

Raises:
    ValueError: If gene not found in built-in definitions

#### `_parse_activity_value()`

**Signature**: `_parse_activity_value(raw: Any, allele: str, source: str | Path) -> float`

Parse the activity value of an allele definition read from a file.

The value is required. A definition without one is not assumed to be normal function (1.0):
that would score an uncharacterized allele as normal, and a "No function" allele with no
value would read as 1.0.

Args:
    raw: The cell or JSON value
    allele: Allele name, for the error message
    source: File the definition came from, for the error message

Raises:
    ValueError: If the value is missing, blank, not a number, not finite, or negative.

#### `_select_gene_definitions()`

**Signature**: `_select_gene_definitions(data: Any, gene: str, source: str | Path) -> dict[str, Any]`

Pick one gene's allele definitions out of a parsed JSON definition file.

The file is either one gene's table (allele -> definition) or a mapping of gene symbol to
such a table. A gene-keyed file that does not contain the requested gene is an error: its
gene entries must not be read as alleles.

Raises:
    ValueError: If the JSON is not an object, or is keyed by gene but lacks ``gene``.

#### `_validated_copy_number()`

**Signature**: `_validated_copy_number(copy_number: Any) -> int`

Return ``copy_number`` as an int, or raise.

Accepts integers (including NumPy integers) and floats with an integral value (2.0).

Raises:
    ValueError: If it is negative, fractional, NaN, infinite, a bool, or not a number.

#### `call_star_alleles()`

**Signature**: `call_star_alleles(variants: set[str] | frozenset[str] | list[str], gene: str, allele_definitions: list[StarAllele] | None=None) -> list[StarAllele]`

Call star alleles from observed variants.

Matches observed variant identifiers against known allele definitions to
determine which star alleles are present. Uses a greedy algorithm that
prioritizes alleles with the most defining variants matched.

The algorithm:
1. Sort allele definitions by number of defining variants (descending)
2. For each allele, check if all defining variants are in the observed set
3. Matched alleles are collected; variants consumed by matched alleles are tracked
4. If no non-reference alleles match, the reference (*1) allele is returned

Args:
    variants: Set of observed variant identifiers (rsIDs or chr:pos:ref:alt)
    gene: Gene symbol (e.g., "CYP2D6")
    allele_definitions: Pre-loaded allele definitions. If None, loads from built-in.

Returns:
    List of matched StarAllele objects, sorted by allele name.
    Returns [*1] if no variant alleles match.

#### `detect_novel_alleles()`

**Signature**: `detect_novel_alleles(variants: set[str] | frozenset[str], known_alleles: list[StarAllele]) -> dict[str, Any]`

Detect potentially novel allele combinations.

Identifies observed variant combinations that do not exactly match any
known allele definition, which may represent novel or rare alleles.

Args:
    variants: Set of observed variant identifiers
    known_alleles: List of known StarAllele definitions

Returns:
    Dictionary with:
        - "has_novel": True if novel combinations detected
        - "unmatched_variants": Variants not explained by any known allele
        - "partial_matches": Alleles with partial but incomplete matches
        - "closest_allele": Best partial match allele name, or None
        - "closest_score": Match score of closest allele

#### `handle_cyp2d6_cnv()`

**Signature**: `handle_cyp2d6_cnv(variants: set[str] | frozenset[str], copy_number: int | float) -> dict[str, Any]`

Handle CYP2D6 copy number variation for star allele calling.

CYP2D6 is subject to gene deletion (*5), gene duplication/multiplication,
and hybrid gene arrangements. This function adjusts star allele calls
based on observed copy number.

Copy number interpretation:
    - 0: Homozygous gene deletion (*5/*5)
    - 1: Hemizygous (one allele deleted): allele/*5
    - 2: Normal diploid (standard calling)
    - 3+: Gene duplication/multiplication

For duplications (CN >= 3), the extra copies are assigned to the called allele with the
higher activity value, and that allele's activity value is multiplied by its total number
of copies (CN - 1; the other allele has one). Per CPIC guidelines, functional allele
duplications are denoted with "xN" suffix (e.g., *1x2).

Args:
    variants: Set of observed variant identifiers
    copy_number: CYP2D6 gene copy number (0, 1, 2, 3+): a non-negative whole number
        (an integer, or a float with an integral value such as 2.0)

Returns:
    Dictionary with:
        - "copy_number": Input copy number
        - "alleles": List of called StarAllele objects (adjusted for CNV)
        - "diplotype_string": Formatted diplotype string
        - "total_activity_score": Sum of all allele activity values
        - "cnv_type": "deletion", "normal", or "duplication"
        - "notes": Clinical notes about the CNV

Raises:
    ValueError: If ``copy_number`` is negative or not a whole number (a negative count is not
        read as a deletion).

#### `load_allele_definitions()`

**Signature**: `load_allele_definitions(gene: str, source: str | None=None, filepath: str | Path | None=None) -> list[StarAllele]`

Load allele definition tables for a pharmacogene.

Loads allele definitions from either the built-in CPIC-derived table or
from an external file (TSV or JSON format). The built-in tables cover
CYP2D6, CYP2C19, CYP2C9, CYP3A5, DPYD, TPMT, NUDT15, and SLCO1B1.

The source is chosen by the arguments and the combination is checked, so a
file is never silently ignored and the built-in table is never silently
substituted for one:

- ``filepath`` given, ``source`` omitted or ``"file"``: load that file.
- no ``filepath``, ``source`` omitted or ``"cpic"``: load the built-in table.
- ``source="file"`` without ``filepath``, ``source="cpic"`` with ``filepath``,
  or any other ``source``: ``ValueError``.

Args:
    gene: Gene symbol (e.g., "CYP2D6", "CYP2C19")
    source: "cpic" for the built-in table, "file" for an external file
        (case-insensitive), or None (default) to infer it from ``filepath``
    filepath: Path to an external allele definition file

Returns:
    List of StarAllele objects for the gene

Raises:
    ValueError: If ``source`` is invalid or contradicts ``filepath``, or the
        gene has no built-in definitions (the message says whether it is
        unknown or known-but-undefined), or a file allele has no valid
        ``activity_value`` (none is assumed), or a JSON file keyed by gene
        lacks ``gene``
    FileNotFoundError: If filepath specified but not found

#### `match_allele_definition()`

**Signature**: `match_allele_definition(observed_variants: set[str] | frozenset[str], allele_definition: StarAllele) -> dict[str, Any]`

Match observed variants against a single allele definition.

Performs detailed comparison between observed variants and an allele's
defining variant set, returning match statistics.

Args:
    observed_variants: Set of observed variant identifiers
    allele_definition: StarAllele definition to match against

Returns:
    Dictionary with match results:
        - "allele": Star allele name
        - "gene": Gene symbol
        - "is_match": True if all defining variants are present
        - "match_score": Fraction of defining variants matched (0.0 - 1.0)
        - "matched_variants": Set of matched variant IDs
        - "missing_variants": Set of unmatched defining variants
        - "extra_variants": Observed variants not in the definition

### Classes

#### `StarAllele`

Represents a pharmacogene star allele.

Attributes:
    name: Star allele designation (e.g., "*1", "*2", "*4")
    gene: Gene symbol (e.g., "CYP2D6", "CYP2C19")
    defining_variants: Set of variant identifiers that define this allele
        Each variant is represented as "chr:pos:ref:alt" or rsID
    function: Functional classification (e.g., "Normal function",
        "No function", "Decreased function", "Increased function")
    activity_value: Numeric activity score (e.g., 1.0 for normal, 0.0 for no function)
    clinical_significance: Clinical significance annotation
    evidence_level: Strength of evidence for this allele definition
    substrate_specificity: Optional notes on substrate-specific effects

##### `is_reference`

**Signature**: `is_reference: bool`

Return True if this is the reference (*1) allele.

##### `is_no_function`

**Signature**: `is_no_function: bool`

Return True if this allele has no function.

##### `matches_variants()`

**Signature**: `matches_variants(observed: frozenset[str] | set[str]) -> bool`

Check if observed variants match this allele's defining variants.

Args:
    observed: Set of observed variant identifiers

Returns:
    True if all defining variants are present in observed set

##### `partial_match_score()`

**Signature**: `partial_match_score(observed: frozenset[str] | set[str]) -> float`

Calculate fraction of defining variants that match.

Args:
    observed: Set of observed variant identifiers

Returns:
    Float between 0.0 and 1.0 indicating match fraction

## annotations.cpic

### Functions

#### `_load_guidelines_from_file()`

**Signature**: `_load_guidelines_from_file(filepath: str | Path) -> list[dict[str, Any]]`

Load CPIC guidelines from an external file.

Args:
    filepath: Path to the guideline file (JSON or TSV)

Returns:
    List of parsed guideline dictionaries

#### `_normalize_phenotype_string()`

**Signature**: `_normalize_phenotype_string(phenotype: str) -> str`

Normalize phenotype string for matching.

Handles every form ``MetabolizerPhenotype.from_string`` accepts: full names ("Poor
Metabolizer"), abbreviations (PM, IM, NM, RM, UM; EM is the old "extensive" for normal)
and the bare tokens ("poor", "ultra-rapid") that ``metabolism.predict_metabolizer_status``
returns. Any other text is only stripped and lowercased.

Args:
    phenotype: Raw phenotype string

Returns:
    Normalized phenotype string

#### `get_dosing_recommendation()`

**Signature**: `get_dosing_recommendation(drug: str, phenotype: str, guidelines: list[dict[str, Any]] | None=None) -> dict[str, Any] | None`

Get dosing recommendation for a drug based on metabolizer phenotype.

Searches all guidelines for the specified drug and returns the
recommendation matching the given phenotype. If the drug has guidelines
for multiple genes, returns the first matching recommendation.

Args:
    drug: Drug name (case-insensitive)
    phenotype: Metabolizer phenotype string (e.g., "Poor Metabolizer", "PM",
        "Normal Metabolizer", "NM")
    guidelines: Optional pre-loaded guidelines

Returns:
    Dictionary with recommendation details, or None if not found.
    Keys: drug, gene, cpic_level, phenotype, recommendation, classification, implications

#### `list_actionable_genes()`

**Signature**: `list_actionable_genes(guidelines: list[dict[str, Any]] | None=None, min_level: str='B') -> list[dict[str, str]]`

List CPIC Level A and B (actionable) gene-drug pairs.

Args:
    guidelines: Optional pre-loaded guidelines
    min_level: Minimum CPIC level to include ("A" for Level A only,
        "B" for Level A and B), case-insensitive

Returns:
    List of dictionaries with keys: gene, drug, cpic_level

Raises:
    ValueError: If ``min_level`` is not "A", "B", "C" or "D" (an unrecognised level is not
        read as "A").

#### `load_cpic_guidelines()`

**Signature**: `load_cpic_guidelines(filepath: str | Path | None=None) -> list[dict[str, Any]]`

Load CPIC guideline data.

Loads from either an external file (TSV or JSON) or the built-in guideline
table. The built-in table holds 9 Level A drug-gene guidelines (codeine, tramadol,
tamoxifen, clopidogrel, voriconazole, warfarin, fluorouracil, azathioprine and
simvastatin), a small subset of CPIC's: a drug-gene pair missing from it does not
mean CPIC has no guideline for it. Load the full set from a file for that.

Args:
    filepath: Optional path to external CPIC guideline file.
        Supports JSON and TSV formats. If None, uses built-in data.

Returns:
    List of guideline dictionaries with keys: drug, gene, cpic_level,
    guideline_url, recommendations

#### `lookup_drug_gene()`

**Signature**: `lookup_drug_gene(drug: str, gene: str, guidelines: list[dict[str, Any]] | None=None) -> dict[str, Any] | None`

Look up CPIC guideline for a specific drug-gene pair.

Args:
    drug: Drug name (case-insensitive)
    gene: Gene symbol (case-insensitive)
    guidelines: Optional pre-loaded guidelines. If None, uses built-in data.

Returns:
    Guideline dictionary if found, None otherwise

#### `parse_cpic_allele_definitions()`

**Signature**: `parse_cpic_allele_definitions(filepath: str | Path) -> dict[str, list[dict[str, Any]]]`

Parse CPIC allele definition tables.

Reads CPIC-formatted allele definition files (TSV format) that define the
genetic variants comprising each star allele for pharmacogenes.

Expected TSV columns: gene, allele, rsid, position, ref, alt, function, activity_value

Every TSV row needs a gene, an allele and an activity value (a number >= 0): none is assumed,
so an allele with no activity value is not read as normal function (1.0).

Args:
    filepath: Path to allele definition TSV/JSON file

Returns:
    Dictionary mapping gene -> list of allele definition dicts.
    Each allele dict contains: allele, defining_variants, function, activity_value

Raises:
    FileNotFoundError: If file does not exist
    ValueError: If a TSV row has no gene or allele, or no valid activity value

## annotations.drug_labels

### Functions

#### `_contains_pgx_terms()`

**Signature**: `_contains_pgx_terms(text: str) -> bool`

Check if text contains pharmacogenomic-relevant terms.

Terms match as whole words, so ordinary words that merely contain one ("general",
"challenge") are not pharmacogenomic content.

Args:
    text: Text to check

Returns:
    True if PGx terms found

#### `_parse_label_text()`

**Signature**: `_parse_label_text(text: str) -> dict[str, Any]`

Parse unstructured drug label text into sections.

Args:
    text: Raw drug label text

Returns:
    Dictionary with identified sections

#### `classify_label_type()`

**Signature**: `classify_label_type(label_data: dict[str, Any]) -> dict[str, Any]`

Classify the strength of PGx labeling for a drug.

Categorizes the drug label based on the strength of pharmacogenomic
testing language:
    - "required": Testing mandated before prescribing
    - "recommended": Testing strongly suggested
    - "actionable": Genotype information affects clinical decisions
    - "informational": PGx information provided for awareness

Args:
    label_data: Parsed drug label dictionary

Returns:
    Dictionary with:
        - "label_type": Classification string
        - "rationale": Explanation for classification

#### `extract_biomarker_info()`

**Signature**: `extract_biomarker_info(label_data: dict[str, Any]) -> dict[str, Any]`

Extract biomarker testing requirements from parsed label data.

Args:
    label_data: Parsed drug label dictionary

Returns:
    Dictionary with:
        - "biomarker_gene": Gene/biomarker name
        - "testing_required": Whether testing is mandated
        - "testing_recommended": Whether testing is recommended
        - "boxed_warning": Whether PGx info is in boxed warning
        - "affected_populations": Populations with specific considerations
        - "test_timing": When testing should be performed

#### `parse_drug_label()`

**Signature**: `parse_drug_label(filepath: str | Path) -> dict[str, Any]`

Parse an FDA drug label file for pharmacogenomic information.

Reads a drug label file (JSON or structured text) and extracts PGx-relevant
sections including biomarker requirements, dosing modifications, and warnings.

Args:
    filepath: Path to drug label file (JSON format preferred)

Returns:
    Dictionary with parsed label data including:
        - "drug": Drug name
        - "brand_name": Brand name
        - "gene_biomarker": Gene/biomarker mentioned
        - "sections": Dict of label sections with PGx content
        - "biomarker_info": Extracted biomarker testing information
        - "label_type": Classification of PGx labeling strength: the file's own
          ``label_type`` when it is one of the four types, otherwise classified from
          the PGx language of the sections by :func:`classify_label_type` (a label with
          no ``label_type`` is not simply "informational")

#### `search_labels_by_gene()`

**Signature**: `search_labels_by_gene(gene: str, labels_db: list[dict[str, Any]] | None=None) -> list[dict[str, Any]]`

Find drug labels mentioning a specific gene or biomarker.

Args:
    gene: Gene symbol (case-insensitive)
    labels_db: Optional pre-loaded labels database. If None, uses built-in data.

Returns:
    List of matching drug label dictionaries

## annotations.pharmgkb

### Functions

#### `get_evidence_level()`

**Signature**: `get_evidence_level(annotation: dict[str, Any]) -> dict[str, Any]`

Extract and interpret the evidence level from a PharmGKB annotation.

Args:
    annotation: Annotation dictionary with "evidence_level" key

Returns:
    Dictionary with:
        - "level": Evidence level string (e.g., "1A")
        - "description": Human-readable description
        - "strength": Qualitative strength ("Strong", "Moderate", "Weak", "Preliminary")
        - "actionable": Whether this evidence level supports clinical action

#### `get_variant_annotations()`

**Signature**: `get_variant_annotations(rsid: str) -> dict[str, Any] | None`

Get variant-specific pharmacogenomic annotations.

Args:
    rsid: Variant rsID (e.g., "rs4244285")

Returns:
    Variant annotation dictionary with gene, allele, frequency data,
    clinical significance, etc. Returns None if variant not found.

#### `parse_clinical_annotations()`

**Signature**: `parse_clinical_annotations(filepath: str | Path) -> list[dict[str, Any]]`

Parse PharmGKB clinical annotation data from file.

Reads a PharmGKB clinical annotation export file (TSV format) and
parses it into structured annotation records.

Expected TSV columns: Clinical Annotation ID, Gene, Level of Evidence,
Drug(s), Phenotype Category, PMID Count, Annotation Text

Args:
    filepath: Path to PharmGKB clinical annotations TSV file

Returns:
    List of parsed annotation dictionaries

#### `query_pharmgkb_annotations()`

**Signature**: `query_pharmgkb_annotations(gene: str | None=None, drug: str | None=None, variant: str | None=None) -> list[dict[str, Any]]`

Query PharmGKB clinical annotations by gene, drug, and/or variant.

Filters the annotation database by any combination of gene, drug, and variant.
At least one filter must be provided.

Args:
    gene: Gene symbol (case-insensitive)
    drug: Drug name (case-insensitive)
    variant: Variant rsID or identifier

Returns:
    List of matching annotation dictionaries

Raises:
    ValueError: If no filter criteria provided

#### `search_drug_pathways()`

**Signature**: `search_drug_pathways(drug: str) -> dict[str, Any] | None`

Get pharmacokinetic/pharmacodynamic pathway information for a drug.

Args:
    drug: Drug name (case-insensitive)

Returns:
    Pathway dictionary with genes involved, mechanism, active metabolite,
    etc. Returns None if drug not found.

## clinical.drug_interaction

### Functions

#### `_abbreviate_phenotype()`

**Signature**: `_abbreviate_phenotype(phenotype: str) -> str`

Convert a phenotype label to its abbreviation.

Args:
    phenotype: Phenotype in any form ``MetabolizerPhenotype.from_string`` accepts
        ("Poor Metabolizer", "poor", "PM", "Extensive Metabolizer", ...)

Returns:
    Abbreviation (PM, IM, NM, RM, UM, IND)

Raises:
    ValueError: If the label is not a recognised phenotype. No abbreviation is made up
        from its first letters.

#### `analyze_drug_gene_interactions()`

**Signature**: `analyze_drug_gene_interactions(drugs: list[str], genotype_data: dict[str, dict[str, Any]]) -> list[DrugRecommendation]`

Analyze drug-gene interactions for a set of drugs and genotypes.

For each drug, checks against all available genotype data to identify
pharmacogenomic interactions and generate clinical recommendations.

Args:
    drugs: List of drug names to analyze
    genotype_data: Dictionary mapping gene symbol to genotype information.
        Each gene entry should have:
            - "phenotype": MetabolizerPhenotype or string in any form
              ``MetabolizerPhenotype.from_string`` accepts (e.g., "Poor Metabolizer",
              "PM", "poor"; the token ``metabolism.predict_metabolizer_status`` returns)
            - "diplotype": Optional diplotype string (e.g., "*1/*4")
            - "activity_score": Optional activity score

Returns:
    List of DrugRecommendation objects for all identified interactions. A gene whose
    phenotype is Indeterminate, or a label that is not a phenotype (free text such as
    "Likely Poor Metabolizer"), gets no recommendation: a warning is logged, and the
    absence of a recommendation is not a finding of no interaction.

#### `calculate_interaction_severity()`

**Signature**: `calculate_interaction_severity(interactions: list[DrugRecommendation]) -> dict[str, Any]`

Calculate overall severity for a set of drug-gene interactions.

Assesses the combined risk from multiple pharmacogenomic interactions,
providing an overall severity classification and risk summary.

Args:
    interactions: List of DrugRecommendation objects

Returns:
    Dictionary with:
        - "overall_severity": Highest severity level across all interactions
        - "major_count": Number of major interactions
        - "moderate_count": Number of moderate interactions
        - "minor_count": Number of minor interactions
        - "total_interactions": Total interaction count
        - "risk_level": Qualitative risk assessment ("High", "Moderate", "Low")
        - "summary": Text summary of findings

#### `check_contraindications()`

**Signature**: `check_contraindications(drug: str, phenotype: str | MetabolizerPhenotype) -> dict[str, Any]`

Check if a drug is contraindicated for a given metabolizer phenotype.

Searches across all gene-drug contraindication records for the specified drug.

Args:
    drug: Drug name (case-insensitive)
    phenotype: Metabolizer phenotype: the enum, or a string in any form
        ``MetabolizerPhenotype.from_string`` accepts ("Poor Metabolizer", "PM", "poor")

Returns:
    Dictionary with:
        - "contraindicated": Whether drug is contraindicated
        - "severity": Interaction severity
        - "gene_interactions": List of specific gene-level contraindication details
        - "overall_recommendation": Summary recommendation

Raises:
    ValueError: If the phenotype is Indeterminate or not a recognised phenotype label.
        Neither can be answered with "no known contraindication", so no answer is given.

#### `polypharmacy_analysis()`

**Signature**: `polypharmacy_analysis(drug_list: list[str], genotype_data: dict[str, dict[str, Any]]) -> dict[str, Any]`

Perform comprehensive polypharmacy pharmacogenomic assessment.

Analyzes all drugs in a patient's medication list against their genotype
data to identify all PGx interactions, drug-drug-gene overlaps, and
cumulative risk.

Args:
    drug_list: List of all drugs the patient is taking
    genotype_data: Dictionary mapping gene -> genotype info (see analyze_drug_gene_interactions)

Returns:
    Dictionary with:
        - "total_drugs": Number of drugs analyzed
        - "drugs_with_interactions": Drugs that have PGx interactions
        - "interactions": Full list of DrugRecommendation objects
        - "severity_summary": Interaction severity breakdown
        - "gene_overlap": Genes affected by multiple drugs
        - "high_risk_combinations": Drug combinations with compounding PGx risk
        - "recommendations": Prioritized clinical recommendations

#### `suggest_alternatives()`

**Signature**: `suggest_alternatives(drug: str, phenotype: str | MetabolizerPhenotype, alternatives_db: dict[str, dict[str, Any]] | None=None) -> dict[str, Any]`

Suggest alternative drugs based on pharmacogenomic phenotype.

Queries the alternatives database for drugs in the same therapeutic class
that are less affected by the patient's pharmacogenomic profile.

Args:
    drug: Drug name (case-insensitive)
    phenotype: Metabolizer phenotype
    alternatives_db: Optional custom alternatives database.
        If None, uses built-in data.

Returns:
    Dictionary with:
        - "drug": Original drug name
        - "therapeutic_class": Drug's therapeutic class
        - "alternatives": List of alternative drug options with details
        - "phenotype": Patient's phenotype
        - "rationale": Why alternatives are suggested

### Classes

#### `DrugRecommendation`

Represents a pharmacogenomic drug recommendation.

Attributes:
    drug: Drug name
    gene: Gene symbol
    phenotype: Metabolizer phenotype
    recommendation: Clinical recommendation text
    evidence_level: CPIC evidence level (A, A/B, B, C, D)
    source: Source of the recommendation (e.g., "CPIC", "PharmGKB", "FDA")
    severity: Interaction severity level
    alternatives: Suggested alternative drugs

#### `InteractionSeverity`

Drug-gene interaction severity classification.

##### `requires_action`

**Signature**: `requires_action: bool`

Return True if this severity level requires clinical action.

## clinical.pathogenicity

### Functions

#### `_benign_classification()`

**Signature**: `_benign_classification(ba: int, bs: int, bp: int) -> ACMGClassification | None`

Benign or Likely Benign call from the counts of BA, BS and BP criteria met, else None.

#### `_pathogenic_classification()`

**Signature**: `_pathogenic_classification(pvs: int, ps: int, pm: int, pp: int) -> ACMGClassification | None`

Pathogenic or Likely Pathogenic call from the counts of PVS, PS, PM and PP criteria met, else None.

#### `_population_frequencies()`

**Signature**: `_population_frequencies(af_data: Any) -> dict[str, float]`

Return a variant's allele frequencies by population; empty when it carries none.

``af_data`` is a number (one "Global" frequency), a mapping of population to frequency, or
absent (``None`` or an empty string). An absent frequency is not zero: callers must not read
"no data" as "absent from controls".

#### `aggregate_evidence()`

**Signature**: `aggregate_evidence(criteria: dict[str, bool]) -> ACMGClassification`

Combine ACMG criteria into a final 5-tier classification.

Implements the combining rules from Richards et al. 2015 (Table 5):

Pathogenic:
    (i)   1 Very Strong + >=1 Strong
    (ii)  1 Very Strong + >=2 Moderate
    (iii) 1 Very Strong + 1 Moderate + 1 Supporting
    (iv)  1 Very Strong + >=2 Supporting
    (v)   >=2 Strong
    (vi)  1 Strong + >=3 Moderate
    (vii) 1 Strong + 2 Moderate + >=2 Supporting
    (viii) 1 Strong + 1 Moderate + >=4 Supporting

Likely Pathogenic:
    (i)   1 Very Strong + 1 Moderate
    (ii)  1 Strong + 1-2 Moderate
    (iii) 1 Strong + >=2 Supporting
    (iv)  >=3 Moderate
    (v)   2 Moderate + >=2 Supporting
    (vi)  1 Moderate + >=4 Supporting

Benign:
    (i) 1 Stand-alone (BA1)
    (ii) >=2 Strong

Likely Benign:
    (i) 1 Strong + 1 Supporting
    (ii) >=2 Supporting

VUS: Everything else, including contradictory evidence: when the criteria for a pathogenic
(or likely pathogenic) call and for a benign (or likely benign) call are both met, the
variant is VUS, as Table 5 says; one side does not win by being checked first.

Args:
    criteria: Dictionary of criterion name -> True/False

Returns:
    ACMGClassification enum value

#### `apply_acmg_criteria()`

**Signature**: `apply_acmg_criteria(variant_data: dict[str, Any]) -> dict[str, bool]`

Evaluate individual ACMG criteria from variant data.

Systematically evaluates each ACMG criterion based on available variant
information. Each criterion is assessed independently.

The population-frequency criteria (PM2, BA1, BS1) need an ``allele_frequency``. A variant
that carries none has not been observed to be absent from controls, so none of the three is
met: a missing frequency is not read as a frequency of 0.

Args:
    variant_data: Variant information dictionary (see classify_variant_acmg)

Returns:
    Dictionary mapping criterion name (str) to met/not-met (bool)

#### `check_gnomad_frequency()`

**Signature**: `check_gnomad_frequency(variant: dict[str, Any], population: str | None=None, threshold: float=0.05) -> dict[str, Any]`

Check variant population frequency against gnomAD thresholds.

Evaluates whether a variant exceeds frequency thresholds relevant to
ACMG benign criteria (BA1 at 5%, BS1 at disease-specific threshold).

Args:
    variant: Variant dictionary with "allele_frequency" key containing
        either a float or a dict of population -> frequency
    population: Optional specific population to check (e.g., "European",
        "East_Asian"), matched case-insensitively. If None, checks maximum
        across all populations.
    threshold: Frequency threshold for classification (default 0.05 for BA1)

Returns:
    Dictionary with:
        - "max_frequency": Maximum allele frequency observed (None if the variant carries
          no frequency: that is not a frequency of 0)
        - "population": Population with highest frequency (None if there is none)
        - "exceeds_threshold": Whether threshold is exceeded
        - "ba1_triggered": Whether BA1 benign criterion is met (AF > 5%)
        - "bs1_triggered": Whether BS1 benign criterion is met (AF > disease threshold)
        - "frequencies": All population frequencies if available

Raises:
    ValueError: If ``population`` is not one of the variant's populations. A population
        the variant has no frequency for is not read as a frequency of 0.

#### `classify_variant_acmg()`

**Signature**: `classify_variant_acmg(variant: dict[str, Any], evidence: dict[str, bool] | None=None) -> dict[str, Any]`

Classify a variant using the ACMG 5-tier system.

Takes variant information and either pre-evaluated evidence criteria or
auto-evaluates criteria from the variant data, then aggregates into a
final classification.

Args:
    variant: Variant data dictionary with keys like:
        - "rsid": rsID (e.g., "rs3892097")
        - "gene": Gene symbol
        - "consequence": Variant consequence (e.g., "missense", "frameshift")
        - "allele_frequency": Population allele frequency dict or float
        - "computational_predictions": Dict of tool -> prediction
        - "functional_data": Functional study evidence
        - "segregation_data": Family segregation data
        - "de_novo": Whether variant is de novo
    evidence: Optional pre-evaluated criteria dict mapping ACMGCriteria name
        to True/False. If None, criteria are auto-evaluated from variant data.

Returns:
    Dictionary with:
        - "classification": ACMGClassification value
        - "criteria_met": List of met ACMG criteria
        - "criteria_details": Dict of criterion -> description for met criteria
        - "score_summary": Count of criteria at each evidence level
        - "confidence": Assessment confidence ("high", "moderate", "low")

#### `query_clinvar()`

**Signature**: `query_clinvar(variant_id: str) -> dict[str, Any] | None`

Look up ClinVar classification for a variant.

Looks the variant up in the built-in table of 6 ClinVar-style records for its
classification, review status, and associated conditions. There is no ClinVar
download behind it, so ``None`` means the variant is not in that table, not that
ClinVar has no record of it.

Args:
    variant_id: Variant rsID or HGVS expression, case-insensitive (ClinVar variation
        IDs are not indexed)

Returns:
    ClinVar record dictionary, or None if not in the built-in table.
    Keys: variant_id, gene, hgvs, classification, review_status,
    stars, condition, last_evaluated, submitter_count

### Classes

#### `ACMGClassification`

ACMG 5-tier variant classification.

Per Richards et al., 2015 (Genet Med 17:405-424).

#### `ACMGCriteria`

Individual ACMG evidence criteria.

Pathogenic criteria:
    PVS1: Null variant in gene where LOF is a known disease mechanism
    PS1-4: Strong pathogenic evidence
    PM1-6: Moderate pathogenic evidence
    PP1-5: Supporting pathogenic evidence

Benign criteria:
    BA1: Stand-alone benign (allele frequency > 5%)
    BS1-4: Strong benign evidence
    BP1-7: Supporting benign evidence

## clinical.reporting

### Functions

#### `_canonical_phenotype()`

**Signature**: `_canonical_phenotype(phenotype: Any) -> Any`

Return the MetabolizerPhenotype for a recognised label ("poor", "PM", ...), else the value as given.

#### `_format_html_report()`

**Signature**: `_format_html_report(report: dict[str, Any]) -> str`

Format report as HTML.

#### `_format_json_report()`

**Signature**: `_format_json_report(report: dict[str, Any]) -> str`

Format report as JSON string.

#### `_format_text_report()`

**Signature**: `_format_text_report(report: dict[str, Any]) -> str`

Format report as plain text.

#### `_get_abbrev()`

**Signature**: `_get_abbrev(phenotype: Any) -> str`

Get phenotype abbreviation from various input types.

#### `_html_escape()`

**Signature**: `_html_escape(text: str) -> str`

Basic HTML escaping for report content.

#### `add_disclaimer()`

**Signature**: `add_disclaimer(report: dict[str, Any], custom_disclaimer: str | None=None) -> dict[str, Any]`

Add or update the clinical disclaimer in a report.

Args:
    report: Report dictionary
    custom_disclaimer: Optional custom disclaimer text.
        If None, uses the standard METAINFORMANT disclaimer.

Returns:
    Updated report dictionary with disclaimer

#### `export_report()`

**Signature**: `export_report(report: dict[str, Any], format: str='text', output_path: str | Path | None=None) -> str`

Export a clinical report to the specified format.

Args:
    report: Report dictionary (from generate_clinical_report)
    format: Output format - "text", "html", or "json"
    output_path: Optional file path to write the report to

Returns:
    Formatted report string

Raises:
    ValueError: If format is not supported

#### `format_recommendation()`

**Signature**: `format_recommendation(drug: str, gene: str, phenotype: str, guideline: dict[str, Any]) -> dict[str, Any]`

Format a single drug-gene recommendation for the clinical report.

Args:
    drug: Drug name
    gene: Gene symbol
    phenotype: Metabolizer phenotype string
    guideline: Guideline data with recommendation, evidence_level, etc.

Returns:
    Formatted recommendation dictionary

#### `generate_clinical_report()`

**Signature**: `generate_clinical_report(patient_data: dict[str, Any], genotypes: dict[str, dict[str, Any]], drugs: list[str] | None=None) -> dict[str, Any]`

Generate a comprehensive clinical pharmacogenomic report.

Analyzes patient genotype data across all pharmacogenes and generates
a structured report with phenotype classifications, drug-specific
recommendations, and clinical action items.

Args:
    patient_data: Patient demographics dictionary with optional keys:
        - "patient_id": Unique identifier
        - "name": Patient name (optional, may be de-identified)
        - "dob": Date of birth
        - "sex": Biological sex
        - "ethnicity": Self-reported ethnicity
        - "ordering_provider": Ordering provider name
    genotypes: Dictionary mapping gene symbol -> genotype info:
        - "diplotype": Diplotype string (e.g., "*1/*4"); the phenotype is classified
          from it. A diplotype with an allele that has no activity value is reported
          Indeterminate (IND, activity score ``None``), as CPIC calls it; the report
          names the allele and adds no drug recommendation for that gene.
        - "phenotype": Optional pre-determined phenotype, in any form
          ``MetabolizerPhenotype.from_string`` accepts ("Poor Metabolizer", "PM",
          "poor"); used when no diplotype is given or it cannot be classified (a gene
          CPIC defines no phenotype for, such as CYP1A2, or a malformed diplotype), in
          which case the reason is reported. With neither, the gene is Indeterminate.
        - "activity_score": Optional pre-computed activity score
    drugs: Drugs to generate specific recommendations for. If None (the default),
        every drug with a built-in CPIC guideline is analysed against the genes
        in ``genotypes``. An empty list requests no drug-specific recommendations.
        Recommendations use the phenotypes reported in "genotype_results", so a
        gene needs only a diplotype or a phenotype, not both.

Returns:
    Comprehensive report dictionary with sections:
        - "header": Report metadata
        - "patient_info": Sanitized patient demographics
        - "genotype_results": Per-gene genotype and phenotype results
        - "drug_recommendations": Per-drug clinical recommendations
        - "interaction_summary": Overall interaction risk assessment
        - "clinical_actions": Prioritized action items
        - "disclaimer": Clinical disclaimer text

#### `generate_summary_table()`

**Signature**: `generate_summary_table(results: list[dict[str, Any]]) -> list[dict[str, str]]`

Generate a summary table of all pharmacogenomic findings.

Creates a condensed tabular view of genotype results and drug
recommendations suitable for display or embedding in reports.

Args:
    results: List of genotype result or recommendation dictionaries

Returns:
    List of row dictionaries with standardized columns:
        - "Gene": Gene symbol
        - "Diplotype": Diplotype string
        - "Phenotype": Metabolizer phenotype
        - "Drug": Drug name (if recommendation)
        - "Recommendation": Brief recommendation
        - "Severity": Interaction severity

## interaction.drug_interactions

### Functions

#### `_find_shared_cyp_pathways()`

**Signature**: `_find_shared_cyp_pathways(drug_a: str, drug_b: str) -> list[str]`

Find CYP enzymes shared by two drugs as substrates.

Args:
    drug_a: Lowercased drug name.
    drug_b: Lowercased drug name.

Returns:
    List of CYP enzyme names where both drugs are substrates.

#### `cyp_inhibition_prediction()`

**Signature**: `cyp_inhibition_prediction(drug: str, cyp_profiles: dict | None=None) -> dict`

Predict CYP enzyme inhibition and induction potential of a drug.

Checks the drug against known CYP inhibitor and inducer tables to
identify which enzymes may be affected and at what potency.

Args:
    drug: Name of the drug to profile.
    cyp_profiles: Optional custom CYP profile database. If ``None``,
        the built-in inhibitor/inducer tables are used. Expected format:
        ``{"inhibitors": {cyp: [...]}, "inducers": {cyp: [...]}}``.

Returns:
    Dictionary with keys:
        - drug: Normalised drug name.
        - affected_enzymes: List of dicts, each with
          ``{enzyme, effect_type, inhibition_type, potency_estimate}``.
        - n_affected: Count of affected enzymes.
        - substrate_of: List of CYP enzymes where this drug is a substrate.
        - summary: Human-readable summary string.

#### `default_interaction_database()`

**Signature**: `default_interaction_database() -> dict`

Return the built-in drug-drug interaction database.

The database contains 39 drug pairs with severity,
mechanism, description, and evidence level. Each entry is keyed by a
frozenset-style tuple of two drug names (lowercased, sorted).

Returns:
    Dictionary mapping ``(drug_a, drug_b)`` tuples (sorted) to interaction
    records with keys: severity, mechanism, description, evidence_level,
    recommendation.

#### `polypharmacy_risk()`

**Signature**: `polypharmacy_risk(medications: list[str], interaction_db: dict | None=None) -> dict`

Assess polypharmacy risk from multiple concurrent medications.

Performs pairwise interaction checks across all medications, identifies
CYP pathway competition, and computes an overall risk score.

Args:
    medications: List of drug names currently prescribed.
    interaction_db: Optional pre-built interaction database.

Returns:
    Dictionary with keys:
        - risk_score: Float 0-1 summarising overall polypharmacy risk.
        - n_medications: Number of medications assessed.
        - interactions: List of pairwise interaction dicts (non-``"none"``).
        - competing_pathways: Dict mapping CYP enzyme to list of competing
          drugs from the input.
        - n_interactions: Count of identified interactions.
        - severity_counts: Dict mapping severity to count.
        - recommendations: List of distinct recommendation strings.

Raises:
    ValueError: If fewer than 2 medications provided.

#### `predict_drug_interaction()`

**Signature**: `predict_drug_interaction(drug_a: str, drug_b: str, interaction_db: dict | None=None) -> dict`

Predict drug-drug interaction between two medications.

Looks up the drug pair in the interaction database and returns the
interaction record including severity, mechanism, description, evidence
level, and clinical recommendation.

Args:
    drug_a: Name of the first drug.
    drug_b: Name of the second drug.
    interaction_db: Optional pre-built interaction database. If ``None``,
        the built-in database from :func:`default_interaction_database` is
        used.

Returns:
    Dictionary with keys:
        - severity: ``"major"`` | ``"moderate"`` | ``"minor"`` | ``"none"``
        - mechanism: Pharmacokinetic or pharmacodynamic mechanism.
        - description: Narrative description of the interaction.
        - evidence_level: ``"A"`` through ``"D"``.
        - recommendation: Clinical recommendation text.
        - drug_a: Normalised name of drug A.
        - drug_b: Normalised name of drug B.

## metabolism.metabolizer_status

### Functions

#### `_phenotype_token()`

**Signature**: `_phenotype_token(phenotype: MetabolizerPhenotype) -> str`

Return the lowercase token for a phenotype (``"poor"``, ``"ultrarapid"``, ...).

#### `_split_diplotype()`

**Signature**: `_split_diplotype(diplotype: str) -> tuple[str, str]`

Split a ``"*1/*4"`` diplotype string into its two allele names.

Raises:
    ValueError: If the string does not contain exactly one ``/``.

#### `_status_token()`

**Signature**: `_status_token(status: str) -> str`

Return the lowercase phenotype token for any accepted phenotype label.

``"Poor Metabolizer"``, ``"PM"``, ``"poor"`` and ``MetabolizerPhenotype.POOR`` all give
``"poor"``. A label that is not a known phenotype is returned stripped and lowercased.

#### `classify_metabolizer()`

**Signature**: `classify_metabolizer(activity_score: float, gene: str) -> str`

Classify metabolizer phenotype from an activity score.

Delegates to :func:`metainformant.pharmacogenomics.alleles.phenotype.predict_metabolizer_status`,
the package's one score-to-phenotype mapping (CPIC definitions, cited there), and returns
the phenotype as a lowercase token.

Args:
    activity_score: Numeric activity score from :func:`compute_activity_score`.
    gene: Gene symbol (e.g. ``"CYP2D6"``).

Returns:
    Phenotype string: one of ``"poor"``, ``"intermediate"``, ``"normal"``,
    ``"rapid"``, or ``"ultrarapid"``.

Raises:
    ValueError: If the activity score is negative or NaN, or the gene has no CPIC
        phenotype definition (the message names the gene and lists the supported
        genes; CYP1A2 can be scored but not classified). No other gene's bands stand in.

#### `compute_activity_score()`

**Signature**: `compute_activity_score(diplotype: str, allele_function_table: dict) -> float`

Compute CPIC-style activity score from a diplotype string.

Parses a diplotype of the form ``"*1/*4"`` and sums the function values
of each allele from the provided table.

This is the strict score: an allele with no value in the table is an error, not
normal function. :func:`predict_metabolizer_status` reports such a diplotype as
indeterminate instead.

Args:
    diplotype: Diplotype string, e.g. ``"*1/*4"``, ``"*1xN/*2"``.
    allele_function_table: Dictionary mapping allele name to numeric
        function value (e.g. ``{"*1": 1.0, "*4": 0.0}``).

Returns:
    Activity score (sum of both allele function values).

Raises:
    ValueError: If diplotype format is invalid or alleles not found.

#### `default_allele_function_table()`

**Signature**: `default_allele_function_table() -> dict`

Return the built-in allele function assignment table.

A copy of the package's single allele table,
``alleles.diplotype._ACTIVITY_SCORE_TABLES`` (CPIC allele functionality tables), so
the values here can never disagree with ``determine_diplotype`` or
``alleles.phenotype.classify_phenotype``. Change allele values there, not here.

Returns:
    Dictionary mapping gene name to a dict of
    ``{allele_name: function_value}``. Function values:
    0.0 = no function, 0.5 = decreased function, 1.0 = normal function,
    1.5 = increased function (CYP2C19 *17); CYP2D6 *10 is 0.25. Covers
    CYP2D6, CYP2C19, CYP2C9, CYP3A5, CYP1A2, DPYD, TPMT, NUDT15, SLCO1B1
    and UGT1A1. CYP1A2's values are not from CPIC (see the table in
    ``alleles.diplotype``) and CPIC defines no CYP1A2 phenotype, so it can be
    scored but not classified. The copy is independent: changing it does not
    change the built-in table.

#### `dose_adjustment()`

**Signature**: `dose_adjustment(metabolizer_status: str, drug: str, dose_guidelines: dict | None=None) -> dict`

Recommend dose adjustment based on metabolizer status and drug.

Looks up the drug in the dose guideline table and returns the
recommendation for the given metabolizer phenotype.

Args:
    metabolizer_status: Metabolizer phenotype in any common form: a token
        (``"poor"``, ``"intermediate"``, ``"normal"``, ``"rapid"``,
        ``"ultrarapid"``), a CPIC name (``"Poor Metabolizer"``), an abbreviation
        (``"PM"``) or a ``MetabolizerPhenotype``, case-insensitively.
    drug: Drug name to look up.
    dose_guidelines: Optional custom dose guidelines. If ``None``, uses
        the built-in CPIC-derived guidelines. Expected format:
        ``{drug: {phenotype: {recommendation, adjusted_dose_fraction, ...}}}``;
        the phenotype keys may be written in any of the forms above.

Returns:
    Dictionary with keys:
        - drug: Normalised drug name.
        - metabolizer_status: The phenotype token the lookup used (the stripped,
          lowercased input when it is not a known phenotype).
        - recommendation: Clinical recommendation text.
        - adjusted_dose_fraction: Multiplier for standard dose (0.0 = avoid,
          1.0 = standard, 2.0 = double); ``None`` for an ``"indeterminate"`` status
          (see below).
        - evidence_level: Evidence grade (``"A"`` through ``"D"``).
        - guideline_source: Source of the guideline.

An ``"indeterminate"`` status, which :func:`predict_metabolizer_status` returns for a
diplotype with an uncharacterized allele, is not normal: unless the guidelines hold an
entry for it, no dose is recommended (``adjusted_dose_fraction`` is ``None``, not 1.0).
An unknown drug or an unrecognised status label still returns standard dosing with
evidence level ``"D"``.

#### `predict_metabolizer_status()`

**Signature**: `predict_metabolizer_status(genotype: dict, gene: str, activity_scores: dict | None=None) -> dict`

Predict metabolizer phenotype from genotype information.

Takes a genotype dictionary containing diplotype information and
computes the metabolizer phenotype using the CPIC activity score
framework. The result is the same phenotype
:func:`metainformant.pharmacogenomics.alleles.phenotype.classify_phenotype`
gives for that diplotype (as a lowercase token), including for an allele
missing from the table: the diplotype is ``"indeterminate"`` with an activity
score of ``None`` and the allele named in ``unknown_alleles`` (CPIC's term for a
diplotype with an uncharacterized allele), never normal function.

Args:
    genotype: Dictionary with at least a ``"diplotype"`` key containing
        the diplotype string (e.g. ``"*1/*4"``), or a ``"alleles"`` key
        with a list of two allele names.
    gene: Gene symbol (e.g. ``"CYP2D6"``).
    activity_scores: Optional custom allele function table for the gene.
        If ``None``, uses the built-in table.

Returns:
    Dictionary with keys:
        - phenotype: Metabolizer classification string (``"indeterminate"`` when an
          allele has no value in the table).
        - activity_score: Computed numeric activity score, or ``None`` when an allele
          has no value in the table.
        - alleles: Tuple of two allele names.
        - gene: Gene symbol.
        - diplotype: Diplotype string.
        - unknown_alleles: The alleles with no value in the table (empty when none).
        - thresholds_used: The activity-score bands the gene uses, keyed by
          phenotype token: ``{"low", "high", "low_inclusive", "high_inclusive"}``.
          CYP2C19, SLCO1B1 and UGT1A1 are classified from allele function, which
          these bands reproduce for the built-in alleles.

Raises:
    ValueError: If genotype is missing required fields or the diplotype string is
        malformed, or if the gene has no CPIC phenotype definition (the message names
        the gene and lists the supported genes; CYP1A2 has activity values but no CPIC
        phenotype, so it can be scored with :func:`compute_activity_score` but not
        classified here).

## visualization.plots

### Functions

#### `_ensure_matplotlib()`

**Signature**: `_ensure_matplotlib() -> None`

Raise ImportError if matplotlib is not available.

#### `_ensure_output_dir()`

**Signature**: `_ensure_output_dir(output_path: str | Path) -> Path`

Ensure the output directory exists and return the Path.

#### `_format_bar_value()`

**Signature**: `_format_bar_value(value: float) -> str`

Label for a bar: a whole number as an integer, anything else to four significant digits.

Counts print as integers and frequencies keep their size (0.02 stays "0.02"; a fixed one
decimal place would print it as "0.0").

#### `_threshold_bands()`

**Signature**: `_threshold_bands(gene: str) -> list[dict[str, Any]]`

Return the phenotype score bands to mark on an activity-score plot of ``gene``.

Empty for a gene CPIC defines no phenotype for (a warning is logged): its scores are
plotted without boundaries, and no other gene's boundaries stand in.

#### `plot_acmg_criteria()`

**Signature**: `plot_acmg_criteria(criteria_results: dict[str, bool], output_path: str | Path, title: str='ACMG Variant Classification Criteria', figsize: tuple[float, float]=(14, 6)) -> Path`

Plot ACMG criteria evaluation summary.

Creates a horizontal bar chart showing which ACMG criteria are met (green)
or not met (gray), organized by evidence strength category.

Args:
    criteria_results: Dictionary mapping criterion name to met/not-met.
        Example: {"PVS1": True, "PS1": False, "PM2": True, "PP3": True, ...}
    output_path: Path to save the figure
    title: Plot title
    figsize: Figure size in inches

Returns:
    Path to saved figure

#### `plot_activity_score_distribution()`

**Signature**: `plot_activity_score_distribution(scores: list[float], gene: str, output_path: str | Path, title: str | None=None, figsize: tuple[float, float]=(10, 6)) -> Path`

Plot activity score distribution as a histogram.

Shows the distribution of activity scores across a cohort, with vertical
lines indicating phenotype threshold boundaries. Only a gene CPIC defines a
phenotype for has boundaries: for any other gene (CYP1A2 has activity values but no
CPIC phenotype) the histogram is drawn without lines and a warning is logged; no
other gene's boundaries are drawn in their place.

Args:
    scores: List of activity score values
    gene: Gene symbol
    output_path: Path to save the figure
    title: Optional custom title
    figsize: Figure size in inches

Returns:
    Path to saved figure

#### `plot_allele_frequencies()`

**Signature**: `plot_allele_frequencies(allele_freqs: dict[str, float], gene: str, output_path: str | Path, plot_type: str='bar', title: str | None=None, figsize: tuple[float, float]=(10, 6)) -> Path`

Plot allele frequency distribution for a pharmacogene.

Creates either a bar chart or pie chart showing the frequency of each
star allele in a population.

Args:
    allele_freqs: Dictionary mapping allele name to frequency.
        Example: {"*1": 0.45, "*2": 0.25, "*4": 0.15, "*10": 0.10, "*17": 0.05}
    gene: Gene symbol for the title
    output_path: Path to save the figure
    plot_type: "bar" for bar chart, "pie" for pie chart (case-insensitive)
    title: Optional custom title
    figsize: Figure size in inches

Returns:
    Path to saved figure

Raises:
    ValueError: If ``plot_type`` is neither "bar" nor "pie" (no other chart is drawn in
        its place)

#### `plot_drug_response_heatmap()`

**Signature**: `plot_drug_response_heatmap(drugs: list[str], genes: list[str], phenotypes: dict[tuple[str, str], str], output_path: str | Path, title: str='Drug-Gene Pharmacogenomic Response Matrix', figsize: tuple[float, float] | None=None) -> Path`

Plot drug-gene response heatmap.

Creates a matrix visualization showing the interaction severity or
metabolizer impact for each drug-gene combination.

Args:
    drugs: List of drug names (rows)
    genes: List of gene symbols (columns)
    phenotypes: Dict mapping (drug, gene) tuple to phenotype/severity string.
        Example: {("codeine", "CYP2D6"): "Major", ("warfarin", "CYP2C9"): "Moderate"}
    output_path: Path to save the figure
    title: Plot title
    figsize: Optional figure size (auto-calculated if None)

Returns:
    Path to saved figure

#### `plot_metabolizer_status()`

**Signature**: `plot_metabolizer_status(phenotypes: dict[str, str | int | float], output_path: str | Path, title: str='Metabolizer Status Distribution', figsize: tuple[float, float]=(10, 6)) -> Path`

Plot metabolizer status distribution as a bar chart.

Creates a bar chart showing the distribution of metabolizer phenotypes,
color-coded by category (PM=red, IM=orange, NM=green, RM=blue, UM=purple).

Args:
    phenotypes: Dictionary mapping phenotype name/abbreviation to count or frequency.
        Example: {"PM": 5, "IM": 20, "NM": 65, "RM": 8, "UM": 2}
    output_path: Path to save the figure
    title: Plot title
    figsize: Figure size in inches (width, height)

Returns:
    Path to saved figure

#### `plot_population_comparison()`

**Signature**: `plot_population_comparison(allele_freqs_by_pop: dict[str, dict[str, float]], gene: str, output_path: str | Path, title: str | None=None, figsize: tuple[float, float]=(12, 7)) -> Path`

Plot cross-population allele frequency comparison.

Creates a grouped bar chart comparing allele frequencies across
different populations (e.g., African, European, East Asian).

Args:
    allele_freqs_by_pop: Nested dict: population -> allele -> frequency.
        Example: {"European": {"*1": 0.6, "*2": 0.3}, "East Asian": {"*1": 0.4, "*10": 0.4}}
    gene: Gene symbol
    output_path: Path to save the figure
    title: Optional custom title
    figsize: Figure size in inches

Returns:
    Path to saved figure

## Data Structures

| Class | Module | Purpose | Key Attributes |
|-------|--------|---------|----------------|
| `StarAllele` | `alleles.star_allele` | A pharmacogene star allele | `.name`, `.gene`, `.defining_variants`, `.function`, `.activity_value` |
| `Diplotype` | `alleles.diplotype` | A pair of alleles with its activity score | `.allele1`, `.allele2`, `.gene`, `.activity_score`, `.phased`, `.confidence`, `.diplotype_string` |
| `MetabolizerPhenotype` | `alleles.phenotype` | Phenotype categories (a `str` enum) | `.value`, `.abbreviation`, `.from_string()` |
| `InteractionSeverity` | `clinical.drug_interaction` | Major, Moderate, Minor or None (a `str` enum) | `.value`, `.requires_action` |
| `DrugRecommendation` | `clinical.drug_interaction` | One drug-gene recommendation | `.drug`, `.gene`, `.phenotype`, `.recommendation`, `.evidence_level`, `.severity`, `.alternatives` |
| `ACMGClassification` | `clinical.pathogenicity` | ACMG 5-tier classification (a `str` enum) | `.value` |
| `ACMGCriteria` | `clinical.pathogenicity` | Individual ACMG evidence criteria (a `str` enum) | `.value` |

`classify_phenotype`, `predict_metabolizer`, `generate_clinical_report` and the `annotations`
functions return plain dictionaries; their keys are listed with each function above.
