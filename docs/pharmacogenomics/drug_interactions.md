# Drug-Gene Interactions

The drug interaction module provides multi-drug pharmacogenomic interaction assessment, contraindication checking, severity scoring, polypharmacy analysis, and alternative drug suggestions. It integrates CPIC guidelines with drug-gene interaction databases to produce clinically actionable recommendations.

## Key Concepts

### Interaction Severity Levels

Drug-gene interactions are classified by severity:

- **Major**: Potentially life-threatening or causing permanent damage. Requires immediate clinical action.
- **Moderate**: May result in deterioration of clinical status. Attention recommended.
- **Minor**: Minimally clinically significant. May increase side effect frequency.
- **None**: No known pharmacogenomic interaction.

### Contraindication Database

The module maintains a contraindication database mapping (drug, gene, phenotype) tuples to interaction details including severity, clinical reason, and alternative drug suggestions. Covered drug-gene pairs include codeine/CYP2D6, clopidogrel/CYP2C19, warfarin/CYP2C9, fluorouracil/DPYD, azathioprine/TPMT, simvastatin/SLCO1B1, and others.

### Polypharmacy Risk Assessment

When patients take multiple drugs, the module identifies gene pathway overlaps (multiple drugs metabolized by the same gene), compounding risks, and high-risk drug combinations.

## Data Structures

### InteractionSeverity (Enum)

```python
class InteractionSeverity(str, Enum):
    MAJOR = "Major"
    MODERATE = "Moderate"
    MINOR = "Minor"
    NONE = "None"
```

Property: `requires_action` -- True for MAJOR and MODERATE.

### DrugRecommendation

```python
@dataclass
class DrugRecommendation:
    drug: str                                       # Drug name
    gene: str                                       # Gene symbol
    phenotype: str                                  # Metabolizer phenotype
    recommendation: str                             # Clinical recommendation
    evidence_level: str = ""                        # CPIC level (A, B, etc.)
    source: str = "CPIC"                            # Recommendation source
    severity: InteractionSeverity = InteractionSeverity.NONE
    alternatives: list[str] = field(default_factory=list)
```

## Function Reference

### analyze_drug_gene_interactions

```python
def analyze_drug_gene_interactions(
    drugs: list[str],
    genotype_data: dict[str, dict[str, Any]],
) -> list[DrugRecommendation]
```

Analyze drug-gene interactions for a set of drugs against patient genotype data. Each gene entry in `genotype_data` should have `phenotype` (string or `MetabolizerPhenotype`), and optionally `diplotype` and `activity_score`. Returns a list of `DrugRecommendation` objects.

### check_contraindications

```python
def check_contraindications(
    drug: str,
    phenotype: str | MetabolizerPhenotype,
) -> dict[str, Any]
```

Check if a drug is contraindicated for a given metabolizer phenotype. Returns `contraindicated` (bool), `severity`, `gene_interactions` (list of gene-level details), and `overall_recommendation`.

### calculate_interaction_severity

```python
def calculate_interaction_severity(
    interactions: list[DrugRecommendation],
) -> dict[str, Any]
```

Calculate overall severity for a set of interactions. Returns `overall_severity`, `major_count`, `moderate_count`, `minor_count`, `total_interactions`, `risk_level` ("Low", "Moderate", "High", "Very High"), and a text `summary`.

### polypharmacy_analysis

```python
def polypharmacy_analysis(
    drug_list: list[str],
    genotype_data: dict[str, dict[str, Any]],
) -> dict[str, Any]
```

Comprehensive polypharmacy assessment. Returns `total_drugs`, `drugs_with_interactions`, `interactions`, `severity_summary`, `gene_overlap` (genes affected by multiple drugs), `high_risk_combinations`, and prioritized `recommendations`.

### suggest_alternatives

```python
def suggest_alternatives(
    drug: str,
    phenotype: str | MetabolizerPhenotype,
    alternatives_db: dict[str, dict[str, Any]] | None = None,
) -> dict[str, Any]
```

Suggest alternative drugs in the same therapeutic class with reduced pharmacogenomic interaction risk. Returns `drug`, `therapeutic_class`, `alternatives` (list with details), `phenotype`, and `rationale`.

## Usage Examples

```python
from metainformant.pharmacogenomics import (
    analyze_drug_gene_interactions, check_contraindications,
    polypharmacy_analysis, suggest_alternatives,
)

# Patient genotype data
genotypes = {
    "CYP2D6": {"phenotype": "Poor Metabolizer", "diplotype": "*4/*4"},
    "CYP2C19": {"phenotype": "Normal Metabolizer", "diplotype": "*1/*1"},
}

# Analyze interactions for prescribed drugs
interactions = analyze_drug_gene_interactions(
    ["codeine", "clopidogrel"], genotypes
)
for rec in interactions:
    print(f"{rec.drug}/{rec.gene}: {rec.severity.value} - {rec.recommendation}")

# Check a specific contraindication
contra = check_contraindications("codeine", "Poor Metabolizer")
if contra["contraindicated"]:
    print(f"CONTRAINDICATED: {contra['overall_recommendation']}")

# Full polypharmacy assessment
poly = polypharmacy_analysis(
    ["codeine", "simvastatin", "warfarin"], genotypes
)
print(f"Risk level: {poly['severity_summary']['risk_level']}")
for rec in poly["recommendations"]:
    print(f"  {rec}")

# Get alternative drugs
alts = suggest_alternatives("codeine", "Poor Metabolizer")
for alt in alts["alternatives"]:
    print(f"  {alt['drug']}: {alt['notes']}")
```

## Configuration

- **Environment prefix**: `PHARMA_`
- Contraindication database covers major CPIC Level A drug-gene pairs
- Alternative drug database is organized by therapeutic class

## Related Modules

- `pharmacogenomics.annotations.cpic` -- CPIC guidelines used for interaction analysis
- `pharmacogenomics.alleles.phenotype` -- Provides metabolizer phenotype classification
- `pharmacogenomics.clinical.reporting` -- Includes interaction results in clinical reports
