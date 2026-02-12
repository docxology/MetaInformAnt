# CPIC Guideline Integration

The CPIC (Clinical Pharmacogenetics Implementation Consortium) module provides parsing and querying of CPIC guideline data for pharmacogenomic clinical decision support. It covers drug-gene pair lookups, phenotype-specific dosing recommendations, actionable gene lists, and allele definition table parsing.

## Key Concepts

### CPIC Evidence Levels

CPIC assigns evidence levels to drug-gene interactions:

- **Level A**: Strong genetic evidence with strong prescribing action recommended
- **Level A/B**: Strong evidence, moderate action
- **Level B**: Moderate genetic evidence with prescribing action recommended
- **Level C**: Moderate evidence, optional action
- **Level D**: Weak evidence or limited data

### Built-in Guidelines

The module ships with built-in CPIC Level A guideline data covering major drug-gene pairs including codeine/CYP2D6, clopidogrel/CYP2C19, warfarin/CYP2C9, fluorouracil/DPYD, azathioprine/TPMT, and simvastatin/SLCO1B1.

### Phenotype Normalization

The module handles both full phenotype names ("Poor Metabolizer") and abbreviations ("PM"), including legacy terminology ("Extensive Metabolizer" maps to "Normal Metabolizer").

## Function Reference

### load_cpic_guidelines

```python
def load_cpic_guidelines(
    filepath: str | Path | None = None,
) -> list[dict[str, Any]]
```

Load CPIC guideline data from built-in tables or an external JSON/TSV file. Returns a list of guideline dictionaries with keys: `drug`, `gene`, `cpic_level`, `guideline_url`, `recommendations`.

### lookup_drug_gene

```python
def lookup_drug_gene(
    drug: str,
    gene: str,
    guidelines: list[dict[str, Any]] | None = None,
) -> dict[str, Any] | None
```

Look up a CPIC guideline for a specific drug-gene pair. Case-insensitive matching. Returns the guideline dictionary or `None` if no guideline exists.

### get_dosing_recommendation

```python
def get_dosing_recommendation(
    drug: str,
    phenotype: str,
    guidelines: list[dict[str, Any]] | None = None,
) -> dict[str, Any] | None
```

Get a dosing recommendation for a drug based on metabolizer phenotype. Accepts both full names and abbreviations. Returns a dictionary with `drug`, `gene`, `cpic_level`, `phenotype`, `recommendation`, `classification`, `implications`, and `guideline_url`.

### list_actionable_genes

```python
def list_actionable_genes(
    guidelines: list[dict[str, Any]] | None = None,
    min_level: str = "B",
) -> list[dict[str, str]]
```

List CPIC-actionable gene-drug pairs at or above a specified evidence level. Returns sorted list of dictionaries with `gene`, `drug`, and `cpic_level`.

### parse_cpic_allele_definitions

```python
def parse_cpic_allele_definitions(
    filepath: str | Path,
) -> dict[str, list[dict[str, Any]]]
```

Parse CPIC allele definition TSV or JSON files. Returns a dictionary mapping gene symbol to lists of allele definitions, each containing `allele`, `defining_variants`, `function`, and `activity_value`.

## Usage Examples

```python
from metainformant.pharmacogenomics import (
    load_cpic_guidelines, lookup_drug_gene,
    get_dosing_recommendation, list_actionable_genes,
)

# Load all CPIC guidelines
guidelines = load_cpic_guidelines()

# Look up a specific drug-gene pair
guideline = lookup_drug_gene("codeine", "CYP2D6")
if guideline:
    print(f"CPIC Level: {guideline['cpic_level']}")

# Get dosing recommendation for a poor metabolizer
rec = get_dosing_recommendation("codeine", "Poor Metabolizer")
if rec:
    print(f"Recommendation: {rec['recommendation']}")
    print(f"Classification: {rec['classification']}")

# Works with abbreviations too
rec_abbrev = get_dosing_recommendation("clopidogrel", "PM")

# List all actionable gene-drug pairs
actionable = list_actionable_genes(min_level="A")
for entry in actionable:
    print(f"{entry['gene']} / {entry['drug']} (Level {entry['cpic_level']})")

# Load custom guidelines from file
custom = load_cpic_guidelines("data/custom_guidelines.json")
```

## Configuration

- **Environment prefix**: `PHARMA_`
- Guidelines can be loaded from JSON files with a `guidelines` key containing a list
- TSV files should have columns: `drug`, `gene`, `cpic_level`, `guideline_url`

## Related Modules

- `pharmacogenomics.alleles.phenotype` -- Phenotype classification (provides phenotype strings for dosing lookups)
- `pharmacogenomics.annotations.pharmgkb` -- PharmGKB annotation queries
- `pharmacogenomics.annotations.drug_labels` -- FDA drug label parsing
- `pharmacogenomics.clinical.drug_interaction` -- Drug-gene interaction analysis
