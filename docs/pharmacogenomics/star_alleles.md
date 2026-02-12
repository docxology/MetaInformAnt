# Star Allele Calling

Star allele calling identifies pharmacogene haplotype variants from observed genetic variants. Star alleles are the standard nomenclature for pharmacogene variants, where \*1 represents the reference (wild-type) allele and numbered alleles (\*2, \*3, etc.) represent variant haplotypes defined by specific SNP and indel combinations.

## Key Concepts

### Star Allele Nomenclature

Each pharmacogene has a set of defined star alleles. The \*1 allele is always the reference allele with normal function. Variant alleles carry specific defining variants (identified by rsIDs or positional notation) and have associated functional classifications: Normal function, Decreased function, No function, or Increased function.

### Activity Scores

Each star allele carries a numeric activity value (0.0 for no function, 0.5 for decreased, 1.0 for normal, 1.5 for increased). These values are summed across both alleles of a diplotype to compute the total activity score used for phenotype classification.

### Greedy Matching Algorithm

Star allele calling uses a greedy algorithm that prioritizes the most specific allele definitions first (those with the most defining variants), then progressively matches less specific alleles from the remaining unmatched variants.

## Data Structures

### StarAllele

```python
@dataclass
class StarAllele:
    name: str                              # e.g., "*1", "*4"
    gene: str                              # e.g., "CYP2D6"
    defining_variants: frozenset[str]      # e.g., frozenset({"rs3892097"})
    function: str = "Unknown function"     # Functional classification
    activity_value: float = 1.0            # Numeric activity score
    clinical_significance: str = ""
    evidence_level: str = ""
    substrate_specificity: str = ""
```

Properties: `is_reference`, `is_no_function`, `matches_variants(observed)`, `partial_match_score(observed)`.

## Function Reference

### load_allele_definitions

```python
def load_allele_definitions(
    gene: str,
    source: str = "cpic",
    filepath: str | Path | None = None,
) -> list[StarAllele]
```

Load allele definitions for a pharmacogene. Built-in tables cover CYP2D6, CYP2C19, CYP2C9, CYP3A5, DPYD, TPMT, NUDT15, and SLCO1B1. Use `source="file"` with a `filepath` for custom TSV or JSON definitions.

### call_star_alleles

```python
def call_star_alleles(
    variants: set[str] | frozenset[str] | list[str],
    gene: str,
    allele_definitions: list[StarAllele] | None = None,
) -> list[StarAllele]
```

Match observed variants against known allele definitions using a greedy algorithm. Returns matched `StarAllele` objects sorted by name. Returns `[*1]` if no variant alleles match.

### match_allele_definition

```python
def match_allele_definition(
    observed_variants: set[str] | frozenset[str],
    allele_definition: StarAllele,
) -> dict[str, Any]
```

Detailed comparison of observed variants against a single allele definition. Returns match results including `is_match`, `match_score`, `matched_variants`, `missing_variants`, and `extra_variants`.

### detect_novel_alleles

```python
def detect_novel_alleles(
    variants: set[str] | frozenset[str],
    known_alleles: list[StarAllele],
) -> dict[str, Any]
```

Identify variant combinations not matching any known allele. Returns `has_novel`, `unmatched_variants`, `partial_matches`, `closest_allele`, and `closest_score`.

### handle_cyp2d6_cnv

```python
def handle_cyp2d6_cnv(
    variants: set[str] | frozenset[str],
    copy_number: int,
) -> dict[str, Any]
```

Handle CYP2D6 copy number variation. Copy number 0 = homozygous deletion (\*5/\*5), 1 = hemizygous, 2 = normal diploid, 3+ = gene duplication. Returns adjusted alleles, diplotype string, total activity score, and clinical notes.

## Usage Examples

```python
from metainformant.pharmacogenomics import (
    call_star_alleles, load_allele_definitions,
    detect_novel_alleles, handle_cyp2d6_cnv,
)

# Load allele definitions for CYP2D6
definitions = load_allele_definitions("CYP2D6")

# Call star alleles from observed variants
observed = {"rs3892097", "rs16947"}
alleles = call_star_alleles(observed, "CYP2D6")
for allele in alleles:
    print(f"{allele.name}: {allele.function} (activity={allele.activity_value})")

# Check for novel alleles
novel = detect_novel_alleles({"rs3892097", "rs999999"}, definitions)
if novel["has_novel"]:
    print(f"Unmatched variants: {novel['unmatched_variants']}")

# Handle CYP2D6 gene duplication
cnv_result = handle_cyp2d6_cnv({"rs16947"}, copy_number=3)
print(f"Diplotype: {cnv_result['diplotype_string']}")
print(f"Activity score: {cnv_result['total_activity_score']}")
```

## Configuration

- **Environment prefix**: `PHARMA_`
- Built-in definitions are derived from CPIC allele definition tables
- Custom definitions can be loaded from TSV (columns: allele, defining_variants, function, activity_value) or JSON files

## Related Modules

- `pharmacogenomics.alleles.diplotype` -- Diplotype determination from star allele pairs
- `pharmacogenomics.alleles.phenotype` -- Metabolizer phenotype prediction from activity scores
- `pharmacogenomics.annotations.cpic` -- CPIC guideline integration
