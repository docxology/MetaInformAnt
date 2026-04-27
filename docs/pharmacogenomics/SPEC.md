# Technical Specification: Pharmacogenomics

## Core data structures

```python
@dataclass(frozen=True)
class StarAllele:
    name: str
    gene: str
    activity: float          # 0.0–2.0+; higher = more function
    defining_variants: set[str]
    function: str            # 'normal','decreased','no function','increased'

@dataclass(frozen=True)
class DiplotypeResult:
    diplotype: str
    activity_score: float
    confidence: float        # 0–1 phase certainty heuristic
    unresolved_sites: list[str] = []

@dataclass(frozen=True)
class PhenotypeResult:
    phenotype: MetabolizerPhenotype   # Enum: PM|IM|NM|RM|UM
    activity_score: float
    gene: str
    evidence: dict[str,Any]           # allele name → activity
```

## Subset-matching engine (call_star_alleles)

For gene G with alleles A₁…Aₙ, each allele has defining variant set Dᵢ.

Algorithm:
1. Filter: keep alleles where Dᵢ ⊆ observed V
2. Score: coverage = |V ∩ Dᵢ| / |Dᵢ|
3. Sort: descending coverage, then descending activity
4. Return top-k (k=1 for fast, k={2,3} for balanced, all minimal for accurate)

Runtime O(n×m) — for CYP2D6 (n≈140) still <10 ms thanks to Cython hash set ops.

## Diplotype determination rules

- Normal phase: `*A/*B` where both alleles have unique marker variants.
- Deletion handling: `*5` (no function) gets activity 0.0.
- Duplication detection: presence of CYP2D6 CNV marker rsIDs (`rs28371685`,
  `rs529216`) implies gene duplication; diplotype reflects copy number
  (`*1/*1x2`).
- Ordering: lower-activity allele first in string; alphabetical tie-breaker.

## Phenotype thresholds (CYP2D6 exemplar)

| Phenotype | Activity score range |
|-----------|----------------------|
| PM | 0.0 |
| IM | >0 – <1.0 |
| NM | 1.0 – <2.25 |
| RM | 2.25 – <3.0 |
| UM | ≥3.0 |

Other genes differ: CYP2C19 RM threshold 2.5, TPMT has no RM category.

## Input / output schema

| Function | Input | Output | Errors |
|----------|-------|--------|--------|
| `call_star_alleles(variants,gene)` | `set[str]` rsIDs | `list[StarAllele]` | `ValueError` unknown gene |
| `determine_diplotype(a1,a2,gene)` | allele names `str` | diplotype `str` | `KeyError` unknown allele |
| `classify_phenotype(diplotype,gene)` | diplotype `str` | `PhenotypeResult` | `ValueError` unknown phenotype string |
| `predict_metabolizer(variants,gene)` | variant set | `PhenotypeResult` | union of above |

## JSON custom allele schema

```json
{
  "gene": "CYP2E1",
  "star_alleles": [
    {"name":"*1","variants":["rs2031920"],"function":"normal","activity":1.0},
    {"name":"*5B","variants":["rs3813867"],"function":"decreased","activity":0.5}
  ]
}
```

Required keys: `name`, `variants` (list of dbSNP rsIDs), `activity` (float).

## Versioned external assets

| Asset | Version | Source | Update |
|-------|---------|--------|--------|
| PharmVar allele tables | 2024-01 | pharmvar.org | custom JSON override |
| CPIC guidelines | v3.0 (Jan 2024) | cpicpgx.org | `load_cpic_guidelines(path=…)` |
| PharmGKB annotations | Monthly | pharmgkb.org | cached 30d auto-refresh |

## Configuration keys (all under `pharmacogenomics.*`)

See [CONFIGURATION.md](CONFIGURATION.md) for the exhaustive list (algorithm,
chunk_size, parallel, allele_definitions_dir, cpic_guidelines_url, pharmgkb.*,
logging).
