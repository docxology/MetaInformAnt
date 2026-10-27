# Technical Specification: Pharmacogenomics

Data model, algorithms and function contracts of `metainformant.pharmacogenomics`,
as implemented. The module has no configuration layer; function arguments are the only
settings (see [CONFIGURATION.md](CONFIGURATION.md)).

## Core data structures

```python
@dataclass
class StarAllele:                        # alleles.star_allele
    name: str                            # "*4"; a missing "*" prefix is added
    gene: str                            # upper-cased
    defining_variants: frozenset[str]    # rsIDs or chr:pos:ref:alt identifiers
    function: str                        # "Normal function", "Decreased function",
                                         # "No function", "Increased function"
    activity_value: float                # 1.0 normal, 0.5 decreased, 0.0 none, 1.5 increased
    clinical_significance: str
    evidence_level: str
    substrate_specificity: str

@dataclass
class Diplotype:                         # alleles.diplotype
    allele1: str                         # the pair is stored in string order, so
    allele2: str                         # "*4/*10" is stored as "*10/*4"
    gene: str
    activity_score: float                # sum of the two allele values
    phased: bool                         # False unless built by phased_diplotype
    confidence: str                      # "high", "moderate" or "low"

class MetabolizerPhenotype(str, Enum):   # alleles.phenotype
    POOR, INTERMEDIATE, NORMAL, RAPID, ULTRARAPID, INDETERMINATE
    # .abbreviation -> PM, IM, NM, RM, UM, IND
    # .from_string(label) accepts "Poor Metabolizer", "PM", "poor", "ultra-rapid",
    #   the pre-2017 synonym "extensive" (normal), any case; ValueError otherwise
```

`classify_phenotype` returns a plain dict with `gene`, `diplotype`, `activity_score`,
`phenotype` (a `MetabolizerPhenotype`), `phenotype_abbreviation` and
`clinical_significance`. `predict_metabolizer` returns a dict whose `phenotype` is a
lowercase token (`'poor'` … `'ultrarapid'`); see the contracts table below.

## Subset-matching engine (`call_star_alleles`)

For a gene with allele definitions A₁…Aₙ, each with a defining variant set Dᵢ, and an
observed variant set V:

1. The definition named `*1` is the reference; it is never matched on its own variants,
   and a definition without defining variants is skipped.
2. The remaining definitions are sorted by `|Dᵢ|`, descending (ties keep definition
   order), so the most specific allele is tried first.
3. A definition is called when Dᵢ ⊆ V minus the variants already consumed by an earlier
   match; its variants are then consumed.
4. If nothing was called, `[*1]` is returned (`[]` when the definitions contain no `*1`).
5. The result is sorted by allele name.

Run time is O(n × m) for n definitions of m variants on average. An empty or unrecognised
variant set therefore reads as the reference allele, not as an error.

## Diplotype determination

- `determine_diplotype(allele1, allele2, gene, scoring_table=None)` looks each allele up in
  `alleles.diplotype._ACTIVITY_SCORE_TABLES[gene]` (or in `scoring_table`) and sums the two
  values. An allele missing from the table scores 1.0 and logs a warning; an unknown gene
  does not raise. Names without a `*` get one unless they are in the table (DPYD's
  `c.2846A>T` and `HapB3` are kept).
- `resolve_ambiguous_diplotypes` prefers candidates whose two alleles both have a known
  value, then the highest activity score.
- `phased_diplotype` calls alleles separately on the two haplotypes of a phase map.
- `handle_cyp2d6_cnv(variants, copy_number)` interprets copy number: 0 gives `*5/*5`, 1 gives
  the called allele with `*5`, 2 is ordinary calling, and 3 or more assigns the extra copies
  to the more active allele, multiplying its activity by `copy_number - 1` and suffixing its
  name `x<copy_number - 1>` (copy number 3 gives `*1x2/*4` for a `*1`/`*4` genotype). The
  marker identifiers `CYP2D6_DEL` and `CYP2D6_DUP` are removed before alleles are called.

## Allele values

`alleles.diplotype._ACTIVITY_SCORE_TABLES` is the package's single allele-function table.
Values come from the CPIC allele functionality tables (the source comments record the
check date, 2026-10-27): normal function 1.0, decreased 0.5, no function 0.0, CYP2C19 `*17`
(increased) 1.5. CYP2D6 carries CPIC's own activity values, including 0.25 for `*9`, `*10`
and `*41` and 0.5 for `*14`. Other values worth knowing:

| Gene | Value |
|------|-------|
| CYP2C9 | `*5`, `*8` are decreased (0.5) |
| SLCO1B1 | `*5`, `*15` are no function (0.0) |
| UGT1A1 | `*37` is decreased (0.5); `*36` is recorded as 1.0 because CPIC gives `*1/*36` the normal phenotype |
| CYP2D6 | `*1xN`, `*2xN` stand for N ≥ 2 copies and carry 2.0, which is exact for two copies only; `handle_cyp2d6_cnv` scores a real copy number |

The star-allele calling definitions in `alleles.star_allele._BUILTIN_ALLELE_DEFINITIONS`
(8 genes, 45 alleles) carry the same activity values, and
`metabolism.metabolizer_status` derives its allele table from the same dictionary.
`default_allele_function_table()` returns an independent copy for 10 genes: CYP2D6,
CYP2C19, CYP2C9, CYP3A5, DPYD, TPMT, NUDT15, SLCO1B1, UGT1A1 and CYP1A2. UGT1A1 and
CYP1A2 have values but no calling definitions (UGT1A1 alleles are defined by promoter
TA-repeat length, which a variant-identifier table cannot express); asking for their
definitions raises a `ValueError` that says so.

## Phenotype definitions

CPIC defines the phenotype of some genes from the diplotype's activity score and of others
from the function of the two alleles. Both definitions are implemented, and
`classify_phenotype` applies the one CPIC uses for the gene.

| Gene | Defined by | PM | IM | NM | RM | UM |
|------|------------|----|----|----|----|----|
| CYP2D6 | activity score (CPIC/DPWG 2019 consensus) | AS = 0 | 0 < AS < 1.25 | 1.25 ≤ AS ≤ 2.25 | none | AS > 2.25 |
| CYP2C9 | activity score | AS ≤ 0.5 | 0.5 < AS < 2 | AS ≥ 2 | none | none |
| DPYD | activity score | AS ≤ 0.5 | 0.5 < AS ≤ 1.5 | AS > 1.5 | none | none |
| CYP3A5, TPMT, NUDT15 | activity score (equal to the allele-function rule) | AS = 0 | 0 < AS < 2 | AS ≥ 2 | none | none |
| CYP2C19 | pair of allele functions | two no-function alleles | one no-function allele | `*1/*1` | `*1/*17` | `*17/*17` |
| SLCO1B1 | pair of allele functions | two reduced-function alleles | one normal and one reduced-function allele | two normal alleles | none | none |
| UGT1A1 | pair of allele functions | two reduced-function alleles | one normal and one reduced-function allele | two normal alleles | none | none |

The score bands tile every score ≥ 0 with no gap and no overlap (a test enforces this).
CYP2D6 has no rapid-metabolizer category. A score that is negative or NaN raises
`ValueError` in `alleles.phenotype.predict_metabolizer_status` and
`metabolism.metabolizer_status.classify_metabolizer` instead of being mapped to a phenotype.

**Allele-function rules.** An allele's function class is read from its activity value: 0
no function, below 1 decreased, 1 normal, above 1 increased.

- CYP2C19 follows CPIC (Moriyama 2017; Lee 2022): `*17/*17` ultrarapid, `*1/*17` rapid,
  `*1/*1` normal, a no-function allele with a normal or increased one (`*1/*2`, `*2/*17`)
  intermediate, two no-function alleles (`*2/*2`) poor. Pairs involving a decreased-function
  allele (`*9`, `*10`) map to CPIC's "likely" categories: with a normal or increased allele
  intermediate (`*1/*9`, `*9/*17`), with a decreased or no-function allele poor (`*9/*9`,
  `*2/*9`). A sum of activity values cannot separate `*9/*17` from `*1/*1` (both 2.0) or
  `*9/*9` from `*1/*2` (both 1.0); the function rules can.
- SLCO1B1 and UGT1A1 (Cooper-DeHoff 2022; Gammal 2016): two normal alleles normal, a normal
  allele with a decreased or no-function one intermediate, two reduced alleles poor.
  SLCO1B1 `*1/*5` is therefore activity score 1.0, decreased function (IM), and `*5/*5` is
  score 0, poor function (PM). A pair with an increased-function allele value (the built-in
  tables have none for these genes; a custom `activity_scores` table can supply one) has no
  category and raises `ValueError`.

**Score bands as data.** `alleles.phenotype._PHENOTYPE_THRESHOLDS[gene]` is a tuple of
`(upper_bound, upper_inclusive, phenotype)` entries. Each band starts where the previous
one ends (the first at 0) and the last runs to infinity. The CYP2C19, SLCO1B1 and UGT1A1
entries are the function rules restated on the score scale, for callers that hold only a
score. `get_phenotype_thresholds(gene)` returns them as dicts with `lower_bound`,
`upper_bound`, `lower_inclusive`, `upper_inclusive`, `phenotype`, `abbreviation` and
`gene_specific`. A gene without bands of its own (for example CYP1A2) is mapped with the
CYP2D6 bands and a logged warning; its entries have `gene_specific` False and a `note`.
`alleles.phenotype.predict_metabolizer_status(score, gene)` applies the bands to a score and
cannot apply the function rules; use `classify_phenotype` when the diplotype is known.

**One implementation.** `metabolism.metabolizer_status` holds no allele or phenotype data of
its own: it reads the allele table above and delegates classification to
`alleles.phenotype`, so `predict_metabolizer` and `classify_phenotype` return the same
phenotype for every diplotype. A diplotype scored with a custom `scoring_table`, or holding
an allele missing from the built-in table, is classified from its stored activity score.

## Input and output contracts

| Function | Input | Output | Errors |
|----------|-------|--------|--------|
| `call_star_alleles(variants, gene, allele_definitions=None)` | `set`, `frozenset` or `list` of identifiers | `list[StarAllele]`, sorted by name | `ValueError` for a gene without built-in definitions when `allele_definitions` is not given |
| `load_allele_definitions(gene, source=None, filepath=None)` | gene symbol; a built-in table or a file | `list[StarAllele]` | `ValueError` for a contradicting `source` and `filepath`, an unknown `source` or a gene without definitions; `FileNotFoundError` for a missing file |
| `determine_diplotype(allele1, allele2, gene, scoring_table=None)` | allele names or `StarAllele` objects | `Diplotype` | none for an unknown gene or allele (scores 1.0, logs a warning) |
| `classify_phenotype(diplotype, gene)` | `Diplotype` or `"*X/*Y"` | dict (see above) | `ValueError` for a string that is not `*X/*Y`; an unknown gene or allele does not raise (it scores 1.0 and the CYP2D6 bands apply, with a logged warning) |
| `alleles.phenotype.predict_metabolizer_status(activity_score, gene)` | score, gene | `MetabolizerPhenotype` | `ValueError` for a negative or NaN score |
| `predict_metabolizer(genotype, gene, activity_scores=None)` (`metabolism.metabolizer_status.predict_metabolizer_status`) | dict with `"diplotype"` or `"alleles"` | dict: `phenotype` (token), `activity_score`, `alleles`, `gene`, `diplotype`, `thresholds_used` | `ValueError` for a missing key, an allele missing from the table, a gene without a table, or an increased-function allele value for SLCO1B1 or UGT1A1 in a custom table |
| `dose_adjustment(metabolizer_status, drug, dose_guidelines=None)` | phenotype in any common form, drug | dict with `recommendation`, `adjusted_dose_fraction`, `evidence_level`, `guideline_source` | none; an unknown drug or phenotype returns standard dosing with evidence level `D` |
| `generate_clinical_report(patient_data, genotypes, drugs=None)` | patient dict, `{gene: {"diplotype": ..., "phenotype": ...}}` | report dict | none for an unclassifiable diplotype: the gene is reported as Indeterminate |

`load_allele_definitions` chooses the source from its arguments: `filepath` alone loads that
file, no arguments loads the built-in table, `source` is case-insensitive, and the
combinations `source='file'` without a path, `source='cpic'` with a path, and any other
`source` raise `ValueError`, so a file is never ignored silently.

`dose_adjustment` accepts the phenotype as a token (`"poor"`), a CPIC name
(`"Poor Metabolizer"`), an abbreviation (`"PM"`) or a `MetabolizerPhenotype`,
case-insensitively; custom guideline tables may key their phenotypes in any of those forms.

`generate_clinical_report` classifies each gene from its `diplotype` when one is given and
otherwise uses its `phenotype` (in any form `MetabolizerPhenotype.from_string` accepts). With
`drugs=None` it analyses every drug that has a built-in CPIC guideline against the genotyped
genes; `drugs=[]` requests no drug recommendations. Recommendations use the phenotypes
shown in the report's `genotype_results`. See [reporting.md](reporting.md).

## Custom allele definition files

A JSON file maps allele name to definition (optionally nested under the gene symbol); a TSV
file has the columns `allele`, `defining_variants` (semicolon-separated), `function` and
`activity_value`:

```json
{
  "*1":  {"defining_variants": [], "function": "Normal function", "activity_value": 1.0},
  "*5B": {"defining_variants": ["rs3813867"], "function": "Decreased function", "activity_value": 0.5}
}
```

`defining_variants` is the only required key for matching; `activity_value` defaults to 1.0.
`clinical_significance` and `evidence_level` are optional.

## Bundled data

All data are Python literals in the source files; nothing is downloaded or cached.

| Asset | Content | Override |
|-------|---------|----------|
| Allele calling definitions | 8 genes, 45 alleles | `load_allele_definitions(gene, filepath=...)`, passed as `allele_definitions=` |
| Allele function table | 10 genes | `scoring_table=` of `determine_diplotype`, `activity_scores=` of `predict_metabolizer` |
| CPIC guidelines | 9 guidelines, all level A | `load_cpic_guidelines(filepath=...)`, passed as `guidelines=` |
| PharmGKB annotations | 9 clinical annotations, 10 variant annotations, 7 drug pathways | `parse_clinical_annotations(filepath)` on your own TSV export |

## Configuration

None: there are no `pharmacogenomics.*` keys, environment variables or YAML sections. See
[CONFIGURATION.md](CONFIGURATION.md) for the function arguments that customise tables,
guidelines and reports.
