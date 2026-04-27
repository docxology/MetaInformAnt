# Pharmacogenomics Module Capabilities

## Quick Reference

### Functions by Submodule

| Function | Submodule | Purpose |
|----------|-----------|---------|
| `calculate_activity_score()` | `alleles.diplotype` | Calculate activity score for an existing Diplotype. |
| `calculate_activity_score_from_alleles()` | `alleles.diplotype` | Calculate the combined activity score for two alleles. |
| `determine_diplotype()` | `alleles.diplotype` | Construct a diplotype from two star alleles with activity scoring. |
| `phased_diplotype()` | `alleles.diplotype` | Determine diplotype using phasing information. |
| `resolve_ambiguous_diplotypes()` | `alleles.diplotype` | Resolve ambiguous diplotype calls when phase is unknown. |
| `_generate_clinical_note()` | `alleles.phenotype` | Generate a brief clinical significance note for a phenotype. |
| `classify_phenotype()` | `alleles.phenotype` | Full diplotype-to-phenotype classification pipeline. |
| `get_phenotype_thresholds()` | `alleles.phenotype` | Get gene-specific phenotype threshold definitions. |
| `population_phenotype_frequencies()` | `alleles.phenotype` | Estimate expected phenotype distribution from allele frequencies. |
| `predict_metabolizer_status()` | `alleles.phenotype` | Map an activity score to a metabolizer phenotype category. |
| `_load_allele_definitions_from_file()` | `alleles.star_allele` | Load allele definitions from an external file. |
| `_load_builtin_definitions()` | `alleles.star_allele` | Load allele definitions from built-in CPIC-derived tables. |
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
| `_abbreviate_phenotype()` | `clinical.drug_interaction` | Convert phenotype string to abbreviation. |
| `analyze_drug_gene_interactions()` | `clinical.drug_interaction` | Analyze drug-gene interactions for a set of drugs and genotypes. |
| `calculate_interaction_severity()` | `clinical.drug_interaction` | Calculate overall severity for a set of drug-gene interactions. |
| `check_contraindications()` | `clinical.drug_interaction` | Check if a drug is contraindicated for a given metabolizer phenotype. |
| `polypharmacy_analysis()` | `clinical.drug_interaction` | Perform comprehensive polypharmacy pharmacogenomic assessment. |
| `suggest_alternatives()` | `clinical.drug_interaction` | Suggest alternative drugs based on pharmacogenomic phenotype. |
| `aggregate_evidence()` | `clinical.pathogenicity` | Combine ACMG criteria into a final 5-tier classification. |
| `apply_acmg_criteria()` | `clinical.pathogenicity` | Evaluate individual ACMG criteria from variant data. |
| `check_gnomad_frequency()` | `clinical.pathogenicity` | Check variant population frequency against gnomAD thresholds. |
| `classify_variant_acmg()` | `clinical.pathogenicity` | Classify a variant using the ACMG 5-tier system. |
| `query_clinvar()` | `clinical.pathogenicity` | Look up ClinVar classification for a variant. |
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
| `classify_metabolizer()` | `metabolism.metabolizer_status` | Classify metabolizer phenotype from an activity score. |
| `compute_activity_score()` | `metabolism.metabolizer_status` | Compute CPIC-style activity score from a diplotype string. |
| `default_allele_function_table()` | `metabolism.metabolizer_status` | Return the built-in allele function assignment table. |
| `dose_adjustment()` | `metabolism.metabolizer_status` | Recommend dose adjustment based on metabolizer status and drug. |
| `predict_metabolizer_status()` | `metabolism.metabolizer_status` | Predict metabolizer phenotype from genotype information. |
| `_ensure_matplotlib()` | `visualization.plots` | Raise ImportError if matplotlib is not available. |
| `_ensure_output_dir()` | `visualization.plots` | Ensure the output directory exists and return the Path. |
| `plot_acmg_criteria()` | `visualization.plots` | Plot ACMG criteria evaluation summary. |
| `plot_activity_score_distribution()` | `visualization.plots` | Plot activity score distribution as a histogram. |
| `plot_allele_frequencies()` | `visualization.plots` | Plot allele frequency distribution for a pharmacogene. |
| `plot_drug_response_heatmap()` | `visualization.plots` | Plot drug-gene response heatmap. |
| `plot_metabolizer_status()` | `visualization.plots` | Plot metabolizer status distribution as a bar chart. |
| `plot_population_comparison()` | `visualization.plots` | Plot cross-population allele frequency comparison. |

---

## alleles.diplotype

### Functions

#### `determine_diplotype()`

**Signature**: `determine_diplotype(allele1, allele2, gene, scoring_table)`

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
    Diplotype object with computed activity score

#### `calculate_activity_score_from_alleles()`

**Signature**: `calculate_activity_score_from_alleles(allele1_name, allele2_name, gene, scoring_table)`

Calculate the combined activity score for two alleles.

Args:
    allele1_name: First allele name
    allele2_name: Second allele name
    gene: Gene symbol
    scoring_table: Optional custom scoring table

Returns:
    Combined activity score (sum of individual allele scores)

#### `calculate_activity_score()`

**Signature**: `calculate_activity_score(diplotype, gene, scoring_table)`

Calculate activity score for an existing Diplotype.

Args:
    diplotype: Diplotype object
    gene: Gene symbol
    scoring_table: Optional custom scoring table

Returns:
    Combined activity score, also updates the Diplotype object

#### `resolve_ambiguous_diplotypes()`

**Signature**: `resolve_ambiguous_diplotypes(possible_diplotypes)`

Resolve ambiguous diplotype calls when phase is unknown.

When genotype data is unphased, multiple diplotype configurations may be
possible. This function selects the most clinically relevant assignment
using the following priority rules (per CPIC recommendations):

1. Prefer diplotypes with known activity scores for both alleles
2. Among scored diplotypes, prefer the one with the highest activity score
   (conservative: assume best-case metabolism)
3. If tied, prefer diplotypes containing more common alleles (*1, *2)

Args:
    possible_diplotypes: List of possible Diplotype objects

Returns:
    The selected most likely Diplotype

Raises:
    ValueError: If possible_diplotypes is empty

#### `phased_diplotype()`

**Signature**: `phased_diplotype(variants, phase_data, gene, allele_definitions)`

Determine diplotype using phasing information.

When phase information is available (e.g., from long reads, family data,
or statistical phasing), this function assigns variants to maternal and
paternal haplotypes before calling star alleles on each independently.

Args:
    variants: Set of all observed variant identifiers
    phase_data: Mapping of variant_id -> phase (0 for haplotype A, 1 for haplotype B).
        Variants not in phase_data are treated as unphased.
    gene: Gene symbol (e.g., "CYP2D6")
    allele_definitions: Optional pre-loaded allele definitions

Returns:
    Diplotype with phased=True and alleles called from each haplotype

### Classes

#### `Diplotype`

Represents a pharmacogene diplotype (pair of star alleles).

Attributes:
    allele1: First star allele (typically the lower-numbered or reference)
    allele2: Second star allele
    gene: Gene symbol
    activity_score: Combined activity score for the diplotype
    phased: Whether the diplotype assignment is phased (True) or inferred (False)
    confidence: Confidence in the diplotype call ("high", "moderate", "low")

**Methods**:

| Method | Purpose |
|--------|---------|
| `__post_init__()` | Normalize diplotype representation (sort alleles alphabetically). |
| `diplotype_string()` | Return formatted diplotype string (e.g., '*1/*2'). |
| `is_homozygous()` | Return True if both alleles are the same. |
| `is_reference()` | Return True if both alleles are reference (*1/*1). |

## alleles.phenotype

### Functions

#### `predict_metabolizer_status()`

**Signature**: `predict_metabolizer_status(activity_score, gene)`

Map an activity score to a metabolizer phenotype category.

Uses gene-specific threshold tables derived from CPIC guidelines. When a
gene has no defined thresholds, applies a default mapping:
    AS = 0 -> PM, 0 < AS < 1.25 -> IM, 1.25 <= AS <= 2.25 -> NM, AS > 2.25 -> UM

Args:
    activity_score: Combined activity score (typically sum of two allele scores)
    gene: Gene symbol (e.g., "CYP2D6")

Returns:
    MetabolizerPhenotype enum value

#### `get_phenotype_thresholds()`

**Signature**: `get_phenotype_thresholds(gene)`

Get gene-specific phenotype threshold definitions.

Args:
    gene: Gene symbol (e.g., "CYP2D6")

Returns:
    List of threshold dictionaries with keys:
        - "lower_bound": Lower activity score boundary (inclusive)
        - "upper_bound": Upper activity score boundary (exclusive)
        - "phenotype": MetabolizerPhenotype value
        - "abbreviation": Phenotype abbreviation (PM, IM, NM, RM, UM)

#### `classify_phenotype()`

**Signature**: `classify_phenotype(diplotype, gene)`

Full diplotype-to-phenotype classification pipeline.

Takes a diplotype (object or string like "*1/*4") and gene, computes the
activity score, and maps it to a metabolizer phenotype.

Args:
    diplotype: Diplotype object or string in format "*X/*Y"
    gene: Gene symbol

Returns:
    Dictionary with:
        - "gene": Gene symbol
        - "diplotype": Diplotype string
        - "activity_score": Computed activity score
        - "phenotype": MetabolizerPhenotype value
        - "phenotype_abbreviation": Short form (PM, IM, NM, RM, UM)
        - "clinical_significance": Brief clinical note

#### `_generate_clinical_note()`

**Signature**: `_generate_clinical_note(phenotype, gene)`

Generate a brief clinical significance note for a phenotype.

Args:
    phenotype: Metabolizer phenotype
    gene: Gene symbol

Returns:
    Clinical significance string

#### `population_phenotype_frequencies()`

**Signature**: `population_phenotype_frequencies(allele_frequencies, gene)`

Estimate expected phenotype distribution from allele frequencies.

Calculates the expected frequency of each metabolizer phenotype in a
population given the allele frequency spectrum, assuming Hardy-Weinberg
equilibrium. For each possible diplotype (all combinations of alleles),
the diplotype frequency is 2*p*q for heterozygotes and p^2 for homozygotes.

Args:
    allele_frequencies: Mapping of allele name -> population frequency.
        Frequencies should sum to approximately 1.0.
    gene: Gene symbol

Returns:
    Dictionary mapping phenotype abbreviation -> expected frequency.
    Example: {"PM": 0.05, "IM": 0.20, "NM": 0.65, "RM": 0.08, "UM": 0.02}

### Classes

#### `MetabolizerPhenotype`

Standardized metabolizer phenotype categories.

Values follow CPIC terminology for consistent clinical communication.
Inherits from str for JSON serialization compatibility.

**Methods**:

| Method | Purpose |
|--------|---------|
| `abbreviation()` | Return standard abbreviation (PM, IM, NM, RM, UM). |

## alleles.star_allele

### Functions

#### `load_allele_definitions()`

**Signature**: `load_allele_definitions(gene, source, filepath)`

Load allele definition tables for a pharmacogene.

Loads allele definitions from either the built-in CPIC-derived table or
from an external file (TSV or JSON format). The built-in tables cover
CYP2D6, CYP2C19, CYP2C9, CYP3A5, DPYD, TPMT, NUDT15, and SLCO1B1.

Args:
    gene: Gene symbol (e.g., "CYP2D6", "CYP2C19")
    source: Definition source - "cpic" for built-in, "file" for external file
    filepath: Path to external allele definition file (required if source="file")

Returns:
    List of StarAllele objects for the gene

Raises:
    ValueError: If gene not found in definitions or invalid source
    FileNotFoundError: If filepath specified but not found

#### `_load_builtin_definitions()`

**Signature**: `_load_builtin_definitions(gene)`

Load allele definitions from built-in CPIC-derived tables.

Args:
    gene: Uppercased gene symbol

Returns:
    List of StarAllele objects

Raises:
    ValueError: If gene not found in built-in definitions

#### `_load_allele_definitions_from_file()`

**Signature**: `_load_allele_definitions_from_file(gene, filepath)`

Load allele definitions from an external file.

Supports TSV format with columns: allele, defining_variants (semicolon-separated),
function, activity_value. Also supports JSON format with the same structure as
_BUILTIN_ALLELE_DEFINITIONS.

Args:
    gene: Uppercased gene symbol
    filepath: Path to the definition file

Returns:
    List of StarAllele objects

#### `call_star_alleles()`

**Signature**: `call_star_alleles(variants, gene, allele_definitions)`

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

#### `match_allele_definition()`

**Signature**: `match_allele_definition(observed_variants, allele_definition)`

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

#### `detect_novel_alleles()`

**Signature**: `detect_novel_alleles(variants, known_alleles)`

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

**Signature**: `handle_cyp2d6_cnv(variants, copy_number)`

Handle CYP2D6 copy number variation for star allele calling.

CYP2D6 is subject to gene deletion (*5), gene duplication/multiplication,
and hybrid gene arrangements. This function adjusts star allele calls
based on observed copy number.

Copy number interpretation:
    - 0: Homozygous gene deletion (*5/*5)
    - 1: Hemizygous (one allele deleted): allele/*5
    - 2: Normal diploid (standard calling)
    - 3+: Gene duplication/multiplication

For duplications (CN >= 3), the activity score of the duplicated allele
is multiplied by the additional copy count. Per CPIC guidelines, functional
allele duplications are denoted with "xN" suffix (e.g., *1x2).

Args:
    variants: Set of observed variant identifiers
    copy_number: CYP2D6 gene copy number (0, 1, 2, 3+)

Returns:
    Dictionary with:
        - "copy_number": Input copy number
        - "alleles": List of called StarAllele objects (adjusted for CNV)
        - "diplotype_string": Formatted diplotype string
        - "total_activity_score": Sum of all allele activity values
        - "cnv_type": "deletion", "normal", or "duplication"
        - "notes": Clinical notes about the CNV

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

**Methods**:

| Method | Purpose |
|--------|---------|
| `__post_init__()` | Validate and normalize star allele fields. |
| `is_reference()` | Return True if this is the reference (*1) allele. |
| `is_no_function()` | Return True if this allele has no function. |
| `matches_variants()` | Check if observed variants match this allele's defining variants. |
| `partial_match_score()` | Calculate fraction of defining variants that match. |

## annotations.cpic

### Functions

#### `load_cpic_guidelines()`

**Signature**: `load_cpic_guidelines(filepath)`

Load CPIC guideline data.

Loads from either an external file (TSV or JSON) or the built-in guideline
table. The built-in table covers major CPIC Level A/A-B guidelines.

Args:
    filepath: Optional path to external CPIC guideline file.
        Supports JSON and TSV formats. If None, uses built-in data.

Returns:
    List of guideline dictionaries with keys: drug, gene, cpic_level,
    guideline_url, recommendations

#### `_load_guidelines_from_file()`

**Signature**: `_load_guidelines_from_file(filepath)`

Load CPIC guidelines from an external file.

Args:
    filepath: Path to the guideline file (JSON or TSV)

Returns:
    List of parsed guideline dictionaries

#### `lookup_drug_gene()`

**Signature**: `lookup_drug_gene(drug, gene, guidelines)`

Look up CPIC guideline for a specific drug-gene pair.

Args:
    drug: Drug name (case-insensitive)
    gene: Gene symbol (case-insensitive)
    guidelines: Optional pre-loaded guidelines. If None, uses built-in data.

Returns:
    Guideline dictionary if found, None otherwise

#### `get_dosing_recommendation()`

**Signature**: `get_dosing_recommendation(drug, phenotype, guidelines)`

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

#### `_normalize_phenotype_string()`

**Signature**: `_normalize_phenotype_string(phenotype)`

Normalize phenotype string for matching.

Handles abbreviations (PM, IM, NM, RM, UM) and full names.

Args:
    phenotype: Raw phenotype string

Returns:
    Normalized phenotype string

#### `list_actionable_genes()`

**Signature**: `list_actionable_genes(guidelines, min_level)`

List CPIC Level A and B (actionable) gene-drug pairs.

Args:
    guidelines: Optional pre-loaded guidelines
    min_level: Minimum CPIC level to include ("A" for Level A only,
        "B" for Level A and B)

Returns:
    List of dictionaries with keys: gene, drug, cpic_level

#### `parse_cpic_allele_definitions()`

**Signature**: `parse_cpic_allele_definitions(filepath)`

Parse CPIC allele definition tables.

Reads CPIC-formatted allele definition files (TSV format) that define the
genetic variants comprising each star allele for pharmacogenes.

Expected TSV columns: gene, allele, rsid, position, ref, alt, function, activity_value

Args:
    filepath: Path to allele definition TSV/JSON file

Returns:
    Dictionary mapping gene -> list of allele definition dicts.
    Each allele dict contains: allele, defining_variants, function, activity_value

Raises:
    FileNotFoundError: If file does not exist

## annotations.drug_labels

### Functions

#### `parse_drug_label()`

**Signature**: `parse_drug_label(filepath)`

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
        - "label_type": Classification of PGx labeling strength

#### `_parse_label_text()`

**Signature**: `_parse_label_text(text)`

Parse unstructured drug label text into sections.

Args:
    text: Raw drug label text

Returns:
    Dictionary with identified sections

#### `_contains_pgx_terms()`

**Signature**: `_contains_pgx_terms(text)`

Check if text contains pharmacogenomic-relevant terms.

Args:
    text: Text to check

Returns:
    True if PGx terms found

#### `extract_biomarker_info()`

**Signature**: `extract_biomarker_info(label_data)`

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

#### `classify_label_type()`

**Signature**: `classify_label_type(label_data)`

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

#### `search_labels_by_gene()`

**Signature**: `search_labels_by_gene(gene, labels_db)`

Find drug labels mentioning a specific gene or biomarker.

Args:
    gene: Gene symbol (case-insensitive)
    labels_db: Optional pre-loaded labels database. If None, uses built-in data.

Returns:
    List of matching drug label dictionaries

## annotations.pharmgkb

### Functions

#### `query_pharmgkb_annotations()`

**Signature**: `query_pharmgkb_annotations(gene, drug, variant)`

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

#### `parse_clinical_annotations()`

**Signature**: `parse_clinical_annotations(filepath)`

Parse PharmGKB clinical annotation data from file.

Reads a PharmGKB clinical annotation export file (TSV format) and
parses it into structured annotation records.

Expected TSV columns: Clinical Annotation ID, Gene, Level of Evidence,
Drug(s), Phenotype Category, PMID Count, Annotation Text

Args:
    filepath: Path to PharmGKB clinical annotations TSV file

Returns:
    List of parsed annotation dictionaries

#### `get_evidence_level()`

**Signature**: `get_evidence_level(annotation)`

Extract and interpret the evidence level from a PharmGKB annotation.

Args:
    annotation: Annotation dictionary with "evidence_level" key

Returns:
    Dictionary with:
        - "level": Evidence level string (e.g., "1A")
        - "description": Human-readable description
        - "strength": Qualitative strength ("Strong", "Moderate", "Weak", "Preliminary")
        - "actionable": Whether this evidence level supports clinical action

#### `search_drug_pathways()`

**Signature**: `search_drug_pathways(drug)`

Get pharmacokinetic/pharmacodynamic pathway information for a drug.

Args:
    drug: Drug name (case-insensitive)

Returns:
    Pathway dictionary with genes involved, mechanism, active metabolite,
    etc. Returns None if drug not found.

#### `get_variant_annotations()`

**Signature**: `get_variant_annotations(rsid)`

Get variant-specific pharmacogenomic annotations.

Args:
    rsid: Variant rsID (e.g., "rs4244285")

Returns:
    Variant annotation dictionary with gene, allele, frequency data,
    clinical significance, etc. Returns None if variant not found.

## clinical.drug_interaction

### Functions

#### `analyze_drug_gene_interactions()`

**Signature**: `analyze_drug_gene_interactions(drugs, genotype_data)`

Analyze drug-gene interactions for a set of drugs and genotypes.

For each drug, checks against all available genotype data to identify
pharmacogenomic interactions and generate clinical recommendations.

Args:
    drugs: List of drug names to analyze
    genotype_data: Dictionary mapping gene symbol to genotype information.
        Each gene entry should have:
            - "phenotype": MetabolizerPhenotype or string (e.g., "Poor Metabolizer")
            - "diplotype": Optional diplotype string (e.g., "*1/*4")
            - "activity_score": Optional activity score

Returns:
    List of DrugRecommendation objects for all identified interactions

#### `_abbreviate_phenotype()`

**Signature**: `_abbreviate_phenotype(phenotype)`

Convert phenotype string to abbreviation.

Args:
    phenotype: Full phenotype name

Returns:
    Abbreviation (PM, IM, NM, RM, UM)

#### `check_contraindications()`

**Signature**: `check_contraindications(drug, phenotype)`

Check if a drug is contraindicated for a given metabolizer phenotype.

Searches across all gene-drug contraindication records for the specified drug.

Args:
    drug: Drug name (case-insensitive)
    phenotype: Metabolizer phenotype (string or enum)

Returns:
    Dictionary with:
        - "contraindicated": Whether drug is contraindicated
        - "severity": Interaction severity
        - "gene_interactions": List of specific gene-level contraindication details
        - "overall_recommendation": Summary recommendation

#### `calculate_interaction_severity()`

**Signature**: `calculate_interaction_severity(interactions)`

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

#### `polypharmacy_analysis()`

**Signature**: `polypharmacy_analysis(drug_list, genotype_data)`

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

**Signature**: `suggest_alternatives(drug, phenotype, alternatives_db)`

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

#### `InteractionSeverity`

Drug-gene interaction severity classification.

**Methods**:

| Method | Purpose |
|--------|---------|
| `requires_action()` | Return True if this severity level requires clinical action. |

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

## clinical.pathogenicity

### Functions

#### `classify_variant_acmg()`

**Signature**: `classify_variant_acmg(variant, evidence)`

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

#### `apply_acmg_criteria()`

**Signature**: `apply_acmg_criteria(variant_data)`

Evaluate individual ACMG criteria from variant data.

Systematically evaluates each ACMG criterion based on available variant
information. Each criterion is assessed independently.

Args:
    variant_data: Variant information dictionary (see classify_variant_acmg)

Returns:
    Dictionary mapping criterion name (str) to met/not-met (bool)

#### `aggregate_evidence()`

**Signature**: `aggregate_evidence(criteria)`

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

VUS: Everything else

Args:
    criteria: Dictionary of criterion name -> True/False

Returns:
    ACMGClassification enum value

#### `query_clinvar()`

**Signature**: `query_clinvar(variant_id)`

Look up ClinVar classification for a variant.

Queries the local ClinVar data cache for a variant's classification,
review status, and associated conditions.

Args:
    variant_id: Variant identifier (rsID, ClinVar variation ID, or HGVS)

Returns:
    ClinVar record dictionary, or None if not found.
    Keys: variant_id, gene, hgvs, classification, review_status,
    stars, condition, last_evaluated, submitter_count

#### `check_gnomad_frequency()`

**Signature**: `check_gnomad_frequency(variant, population, threshold)`

Check variant population frequency against gnomAD thresholds.

Evaluates whether a variant exceeds frequency thresholds relevant to
ACMG benign criteria (BA1 at 5%, BS1 at disease-specific threshold).

Args:
    variant: Variant dictionary with "allele_frequency" key containing
        either a float or a dict of population -> frequency
    population: Optional specific population to check (e.g., "European",
        "East_Asian"). If None, checks maximum across all populations.
    threshold: Frequency threshold for classification (default 0.05 for BA1)

Returns:
    Dictionary with:
        - "max_frequency": Maximum allele frequency observed
        - "population": Population with highest frequency
        - "exceeds_threshold": Whether threshold is exceeded
        - "ba1_triggered": Whether BA1 benign criterion is met (AF > 5%)
        - "bs1_triggered": Whether BS1 benign criterion is met (AF > disease threshold)
        - "frequencies": All population frequencies if available

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

#### `generate_clinical_report()`

**Signature**: `generate_clinical_report(patient_data, genotypes, drugs)`

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
        - "diplotype": Diplotype string (e.g., "*1/*4")
        - "phenotype": Optional pre-determined phenotype
        - "activity_score": Optional pre-computed activity score
    drugs: Optional list of drugs to generate specific recommendations for.
        If None, generates recommendations for all drugs with CPIC guidelines.

Returns:
    Comprehensive report dictionary with sections:
        - "header": Report metadata
        - "patient_info": Sanitized patient demographics
        - "genotype_results": Per-gene genotype and phenotype results
        - "drug_recommendations": Per-drug clinical recommendations
        - "interaction_summary": Overall interaction risk assessment
        - "clinical_actions": Prioritized action items
        - "disclaimer": Clinical disclaimer text

#### `_get_abbrev()`

**Signature**: `_get_abbrev(phenotype)`

Get phenotype abbreviation from various input types.

#### `format_recommendation()`

**Signature**: `format_recommendation(drug, gene, phenotype, guideline)`

Format a single drug-gene recommendation for the clinical report.

Args:
    drug: Drug name
    gene: Gene symbol
    phenotype: Metabolizer phenotype string
    guideline: Guideline data with recommendation, evidence_level, etc.

Returns:
    Formatted recommendation dictionary

#### `generate_summary_table()`

**Signature**: `generate_summary_table(results)`

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

#### `export_report()`

**Signature**: `export_report(report, format, output_path)`

Export a clinical report to the specified format.

Args:
    report: Report dictionary (from generate_clinical_report)
    format: Output format - "text", "html", or "json"
    output_path: Optional file path to write the report to

Returns:
    Formatted report string

Raises:
    ValueError: If format is not supported

#### `_format_text_report()`

**Signature**: `_format_text_report(report)`

Format report as plain text.

#### `_format_html_report()`

**Signature**: `_format_html_report(report)`

Format report as HTML.

#### `_format_json_report()`

**Signature**: `_format_json_report(report)`

Format report as JSON string.

#### `_html_escape()`

**Signature**: `_html_escape(text)`

Basic HTML escaping for report content.

#### `add_disclaimer()`

**Signature**: `add_disclaimer(report, custom_disclaimer)`

Add or update the clinical disclaimer in a report.

Args:
    report: Report dictionary
    custom_disclaimer: Optional custom disclaimer text.
        If None, uses the standard METAINFORMANT disclaimer.

Returns:
    Updated report dictionary with disclaimer

## interaction.drug_interactions

### Functions

#### `default_interaction_database()`

**Signature**: `default_interaction_database()`

Return the built-in drug-drug interaction database.

The database contains approximately 100 common drug pairs with severity,
mechanism, description, and evidence level. Each entry is keyed by a
frozenset-style tuple of two drug names (lowercased, sorted).

Returns:
    Dictionary mapping ``(drug_a, drug_b)`` tuples (sorted) to interaction
    records with keys: severity, mechanism, description, evidence_level,
    recommendation.

#### `predict_drug_interaction()`

**Signature**: `predict_drug_interaction(drug_a, drug_b, interaction_db)`

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

#### `_find_shared_cyp_pathways()`

**Signature**: `_find_shared_cyp_pathways(drug_a, drug_b)`

Find CYP enzymes shared by two drugs as substrates.

Args:
    drug_a: Lowercased drug name.
    drug_b: Lowercased drug name.

Returns:
    List of CYP enzyme names where both drugs are substrates.

#### `polypharmacy_risk()`

**Signature**: `polypharmacy_risk(medications, interaction_db)`

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

#### `cyp_inhibition_prediction()`

**Signature**: `cyp_inhibition_prediction(drug, cyp_profiles)`

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

## metabolism.metabolizer_status

### Functions

#### `default_allele_function_table()`

**Signature**: `default_allele_function_table()`

Return the built-in allele function assignment table.

Covers CYP2D6, CYP2C19, and CYP2C9 with CPIC-consensus function
assignments. Each gene maps to a dict of allele -> function_value.

Returns:
    Dictionary mapping gene name to a dict of
    ``{allele_name: function_value}``. Function values:
    0.0 = no function, 0.5 = decreased function, 1.0 = normal function,
    1.5 = increased function (CYP2C19 *17).

#### `compute_activity_score()`

**Signature**: `compute_activity_score(diplotype, allele_function_table)`

Compute CPIC-style activity score from a diplotype string.

Parses a diplotype of the form ``"*1/*4"`` and sums the function values
of each allele from the provided table.

Args:
    diplotype: Diplotype string, e.g. ``"*1/*4"``, ``"*1xN/*2"``.
    allele_function_table: Dictionary mapping allele name to numeric
        function value (e.g. ``{"*1": 1.0, "*4": 0.0}``).

Returns:
    Activity score (sum of both allele function values).

Raises:
    ValueError: If diplotype format is invalid or alleles not found.

#### `classify_metabolizer()`

**Signature**: `classify_metabolizer(activity_score, gene)`

Classify metabolizer phenotype from an activity score.

Uses gene-specific thresholds from the CPIC framework to map an
activity score to a phenotype category.

Args:
    activity_score: Numeric activity score from :func:`compute_activity_score`.
    gene: Gene symbol (e.g. ``"CYP2D6"``).

Returns:
    Phenotype string: one of ``"poor"``, ``"intermediate"``, ``"normal"``,
    ``"rapid"``, or ``"ultrarapid"``.

#### `predict_metabolizer_status()`

**Signature**: `predict_metabolizer_status(genotype, gene, activity_scores)`

Predict metabolizer phenotype from genotype information.

Takes a genotype dictionary containing diplotype information and
computes the metabolizer phenotype using the CPIC activity score
framework.

Args:
    genotype: Dictionary with at least a ``"diplotype"`` key containing
        the diplotype string (e.g. ``"*1/*4"``), or a ``"alleles"`` key
        with a list of two allele names.
    gene: Gene symbol (e.g. ``"CYP2D6"``).
    activity_scores: Optional custom allele function table for the gene.
        If ``None``, uses the built-in table.

Returns:
    Dictionary with keys:
        - phenotype: Metabolizer classification string.
        - activity_score: Computed numeric activity score.
        - alleles: Tuple of two allele names.
        - gene: Gene symbol.
        - diplotype: Diplotype string.
        - thresholds_used: The threshold ranges applied.

Raises:
    ValueError: If genotype is missing required fields.

#### `dose_adjustment()`

**Signature**: `dose_adjustment(metabolizer_status, drug, dose_guidelines)`

Recommend dose adjustment based on metabolizer status and drug.

Looks up the drug in the dose guideline table and returns the
recommendation for the given metabolizer phenotype.

Args:
    metabolizer_status: Metabolizer phenotype string (e.g. ``"poor"``,
        ``"intermediate"``, ``"normal"``, ``"ultrarapid"``).
    drug: Drug name to look up.
    dose_guidelines: Optional custom dose guidelines. If ``None``, uses
        the built-in CPIC-derived guidelines. Expected format:
        ``{drug: {phenotype: {recommendation, adjusted_dose_fraction, ...}}}``.

Returns:
    Dictionary with keys:
        - drug: Normalised drug name.
        - metabolizer_status: Input phenotype.
        - recommendation: Clinical recommendation text.
        - adjusted_dose_fraction: Multiplier for standard dose (0.0 = avoid,
          1.0 = standard, 2.0 = double).
        - evidence_level: Evidence grade (``"A"`` through ``"D"``).
        - guideline_source: Source of the guideline.

## visualization.plots

### Functions

#### `_ensure_matplotlib()`

**Signature**: `_ensure_matplotlib()`

Raise ImportError if matplotlib is not available.

#### `_ensure_output_dir()`

**Signature**: `_ensure_output_dir(output_path)`

Ensure the output directory exists and return the Path.

#### `plot_metabolizer_status()`

**Signature**: `plot_metabolizer_status(phenotypes, output_path, title, figsize)`

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

#### `plot_allele_frequencies()`

**Signature**: `plot_allele_frequencies(allele_freqs, gene, output_path, plot_type, title, figsize)`

Plot allele frequency distribution for a pharmacogene.

Creates either a bar chart or pie chart showing the frequency of each
star allele in a population.

Args:
    allele_freqs: Dictionary mapping allele name to frequency.
        Example: {"*1": 0.45, "*2": 0.25, "*4": 0.15, "*10": 0.10, "*17": 0.05}
    gene: Gene symbol for the title
    output_path: Path to save the figure
    plot_type: "bar" for bar chart, "pie" for pie chart
    title: Optional custom title
    figsize: Figure size in inches

Returns:
    Path to saved figure

#### `plot_activity_score_distribution()`

**Signature**: `plot_activity_score_distribution(scores, gene, output_path, title, figsize)`

Plot activity score distribution as a histogram.

Shows the distribution of activity scores across a cohort, with vertical
lines indicating phenotype threshold boundaries.

Args:
    scores: List of activity score values
    gene: Gene symbol
    output_path: Path to save the figure
    title: Optional custom title
    figsize: Figure size in inches

Returns:
    Path to saved figure

#### `plot_drug_response_heatmap()`

**Signature**: `plot_drug_response_heatmap(drugs, genes, phenotypes, output_path, title, figsize)`

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

#### `plot_population_comparison()`

**Signature**: `plot_population_comparison(allele_freqs_by_pop, gene, output_path, title, figsize)`

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

#### `plot_acmg_criteria()`

**Signature**: `plot_acmg_criteria(criteria_results, output_path, title, figsize)`

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

## Data Structures

| Class | Purpose | Key Attributes |
|-------|---------|----------------|
| Dataset | Container for input data | `.data`, `.metadata` |
| Result | Analysis output | `.values`, `.stats` |
| Config | Settings object | `.params`, `.paths` |
