# ACMG Variant Classification

The ACMG module implements the American College of Medical Genetics and Genomics / Association for Molecular Pathology (ACMG/AMP) 5-tier variant classification system based on the 2015 Richards et al. guidelines. It evaluates individual evidence criteria, combines them using the published combining rules, and produces a final classification.

## Key Concepts

### Five-Tier Classification

ACMG classifies variants into five categories:

1. **Pathogenic** -- Strong evidence the variant causes disease
2. **Likely Pathogenic** -- Evidence strongly suggests pathogenicity (>90% certainty)
3. **Uncertain Significance (VUS)** -- Insufficient evidence to classify
4. **Likely Benign** -- Evidence strongly suggests variant is benign
5. **Benign** -- Strong evidence the variant does not cause disease

### Evidence Criteria

Evidence criteria are organized by strength and direction:

**Pathogenic criteria:**
- PVS1: Very Strong (null variant in LOF-mechanism gene)
- PS1-PS4: Strong pathogenic evidence
- PM1-PM6: Moderate pathogenic evidence
- PP1-PP5: Supporting pathogenic evidence

**Benign criteria:**
- BA1: Stand-alone benign (allele frequency >5%)
- BS1-BS4: Strong benign evidence
- BP1-BP7: Supporting benign evidence

### Combining Rules

The module implements the full combining rules from Richards et al. 2015 Table 5. For example, Pathogenic requires: 1 Very Strong + 1 Strong, or 2 Strong, or 1 Very Strong + 2 Moderate, among other combinations.

## Data Structures

### ACMGClassification (Enum)

```python
class ACMGClassification(str, Enum):
    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely Pathogenic"
    VUS = "Uncertain Significance"
    LIKELY_BENIGN = "Likely Benign"
    BENIGN = "Benign"
```

### ACMGCriteria (Enum)

All 28 ACMG criteria as enum members: `PVS1`, `PS1`-`PS4`, `PM1`-`PM6`, `PP1`-`PP5`, `BA1`, `BS1`-`BS4`, `BP1`-`BP7`.

## Function Reference

### classify_variant_acmg

```python
def classify_variant_acmg(
    variant: dict[str, Any],
    evidence: dict[str, bool] | None = None,
) -> dict[str, Any]
```

Classify a variant using the ACMG 5-tier system. If `evidence` is provided, uses those pre-evaluated criteria; otherwise auto-evaluates from variant data. Returns `classification`, `criteria_met`, `criteria_details`, `score_summary`, and `confidence`.

### apply_acmg_criteria

```python
def apply_acmg_criteria(
    variant_data: dict[str, Any],
) -> dict[str, bool]
```

Evaluate all 28 ACMG criteria from variant data. Input keys include: `consequence`, `allele_frequency`, `computational_predictions`, `functional_data`, `segregation_data`, `de_novo`, `paternity_confirmed`, and many others. Returns a dictionary mapping criterion name to True/False.

### aggregate_evidence

```python
def aggregate_evidence(
    criteria: dict[str, bool],
) -> ACMGClassification
```

Combine met criteria into a final classification using the ACMG combining rules. Benign rules are checked first (BA1 is standalone), then Pathogenic, then Likely Pathogenic. Defaults to VUS.

### query_clinvar

```python
def query_clinvar(
    variant_id: str,
) -> dict[str, Any] | None
```

Look up ClinVar classification for a variant by rsID or HGVS notation. Returns a record with `variant_id`, `gene`, `hgvs`, `classification`, `review_status`, `stars`, `condition`, `last_evaluated`, and `submitter_count`.

### check_gnomad_frequency

```python
def check_gnomad_frequency(
    variant: dict[str, Any],
    population: str | None = None,
    threshold: float = 0.05,
) -> dict[str, Any]
```

Check variant population frequency for ACMG benign criteria. Returns `max_frequency`, `population`, `exceeds_threshold`, `ba1_triggered`, `bs1_triggered`, and `frequencies`.

## Usage Examples

```python
from metainformant.pharmacogenomics import (
    classify_variant_acmg, apply_acmg_criteria,
    aggregate_evidence, query_clinvar, check_gnomad_frequency,
)

# Full automated classification
variant = {
    "rsid": "rs3892097",
    "gene": "CYP2D6",
    "consequence": "splice_donor",
    "allele_frequency": {"European": 0.02, "East_Asian": 0.005},
    "computational_predictions": {
        "CADD": "damaging", "REVEL": "pathogenic", "SIFT": "deleterious"
    },
}
result = classify_variant_acmg(variant)
print(f"Classification: {result['classification'].value}")
print(f"Criteria met: {result['criteria_met']}")

# Manual criteria evaluation
criteria = apply_acmg_criteria(variant)
classification = aggregate_evidence(criteria)

# ClinVar lookup
clinvar = query_clinvar("rs3892097")
if clinvar:
    print(f"ClinVar: {clinvar['classification']} ({clinvar['stars']} stars)")

# Population frequency check
freq_result = check_gnomad_frequency(variant, threshold=0.01)
if freq_result["ba1_triggered"]:
    print("BA1 benign criterion met (AF > 5%)")
```

## Configuration

- **Environment prefix**: `PHARMA_`
- ClinVar data is stored locally for fast lookups (no network required)
- gnomAD frequency checks use variant-level allele frequency data from the input

## Related Modules

- `pharmacogenomics.clinical.drug_interaction` -- Uses pathogenicity to inform interaction severity
- `pharmacogenomics.clinical.reporting` -- Includes ACMG results in clinical reports
- `pharmacogenomics.annotations.cpic` -- CPIC guidelines cross-referenced with classification
