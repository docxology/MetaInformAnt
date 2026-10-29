# Alleles

Pharmacogenomics allele subpackage providing star allele identification from observed variants, diplotype determination with activity scoring, and metabolizer phenotype assignment per CPIC definitions. It owns the allele activity values (`diplotype._ACTIVITY_SCORE_TABLES`) and the phenotype definitions (`phenotype.py`); the `metabolism` subpackage reads both.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports star_allele, diplotype, phenotype submodules |
| `star_allele.py` | Star allele calling against CPIC allele definition tables |
| `diplotype.py` | Diplotype determination and activity score calculation |
| `phenotype.py` | Metabolizer phenotype assignment: score bands (CYP2D6, CYP2C9, DPYD, CYP3A5, TPMT, NUDT15) and allele-function rules (CYP2C19, SLCO1B1, UGT1A1) |

## Key Functions

| Function | Description |
|----------|-------------|
| `star_allele.call_star_alleles()` | Call star alleles from observed variants |
| `star_allele.load_allele_definitions()` | Load allele definition tables: the built-in table, or a file given by `filepath` (contradicting `source`/`filepath` combinations raise `ValueError`) |
| `star_allele.detect_novel_alleles()` | Detect novel alleles not in reference definitions |
| `star_allele.handle_cyp2d6_cnv()` | Handle CYP2D6 copy number variation (a negative or fractional copy number raises `ValueError`; alleles called beyond what the copies can carry are named in `notes`) |
| `diplotype.determine_diplotype()` | Determine diplotype from two star alleles; the activity score is `None` when an allele has no activity value |
| `diplotype.calculate_activity_score()` | Calculate combined activity score for a diplotype (`None` if an allele has no activity value) |
| `diplotype.uncharacterized_alleles()` | List the alleles of a diplotype that have no activity value |
| `diplotype.resolve_ambiguous_diplotypes()` | Resolve ambiguous diplotype assignments |
| `diplotype.phased_diplotype()` | Diplotype from phased variants; a haplotype that matches several alleles is called by its lowest-activity one (confidence `"moderate"`, with a warning) |
| `phenotype.classify_phenotype()` | Classify a diplotype (object or `"*X/*Y"`) into a metabolizer phenotype, returning a dict with the activity score, phenotype and any `unknown_alleles` |
| `phenotype.predict_metabolizer_status()` | Map an activity score to PM/IM/NM/RM/UM with the gene's CPIC score bands (negative, NaN or `None` scores and genes without a CPIC phenotype raise `ValueError`) |
| `phenotype.get_phenotype_thresholds()` | List a gene's score bands (a gene without a CPIC phenotype raises `ValueError`) |
| `phenotype.supported_phenotype_genes()` | The genes with a CPIC phenotype definition |
| `phenotype.population_phenotype_frequencies()` | Estimate phenotype frequencies from allele data; alleles with no activity value count as `"IND"` |
| `phenotype.MetabolizerPhenotype` | Phenotype enum; `from_string()` parses `"Poor Metabolizer"`, `"PM"` or `"poor"` |

## Usage

```python
from metainformant.pharmacogenomics.alleles import star_allele, diplotype, phenotype

alleles = star_allele.call_star_alleles(observed_variants, "CYP2D6")
names = [a.name for a in alleles] + ["*1", "*1"]  # pad: a single called allele pairs with *1
dip = diplotype.determine_diplotype(names[0], names[1], "CYP2D6")
pheno = phenotype.classify_phenotype(dip, "CYP2D6")
print(pheno["phenotype"].value, pheno["activity_score"])
```

CYP2D6 has no rapid-metabolizer category. UGT1A1 and CYP1A2 have activity values but no
calling definitions, so `call_star_alleles` raises `ValueError` for them; build their diplotypes
with `determine_diplotype` from allele names called elsewhere. See
[docs/pharmacogenomics/SPEC.md](../../../../docs/pharmacogenomics/SPEC.md#phenotype-definitions).

What has no CPIC answer is never given a normal one:

- **Genes.** Only the genes CPIC defines a phenotype for are classified
  (`phenotype.supported_phenotype_genes()`). Any other gene raises `ValueError` naming the gene
  and listing the supported ones, in `classify_phenotype`, `predict_metabolizer_status`,
  `get_phenotype_thresholds` and `population_phenotype_frequencies`; no other gene's bands stand
  in. CYP1A2 has activity values (it can be scored) but no CPIC phenotype (it cannot be
  classified), and its values have no CPIC source: see the comment in `diplotype.py`.
- **Alleles.** An allele with no activity value is uncharacterized, and a diplotype that carries
  one is indeterminate, as CPIC calls it: `classify_phenotype` returns
  `MetabolizerPhenotype.INDETERMINATE`, an activity score of `None` and the allele in
  `unknown_alleles`. Pass a `scoring_table` that lists the allele to score it.
- **Unscored diplotypes.** `Diplotype.activity_score` is `None` until the diplotype is scored (it
  is not 0.0, which would read as a poor metabolizer); `classify_phenotype` scores it from its
  alleles. A score you set yourself, including 0.0, is classified as given.
- **Definition files.** `load_allele_definitions(..., filepath=...)` needs an `activity_value` for
  every allele (none is assumed to be normal function) and refuses a JSON file keyed by gene that
  lacks the requested gene.
