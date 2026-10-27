# Pharmacogenomics Module

Predict drug response from germline genetic variation using CPIC-derived star-allele
tables and CPIC guidelines. The built-in allele tables cover CYP2D6, CYP2C19, CYP2C9,
CYP3A5, DPYD, TPMT, NUDT15 and SLCO1B1; other pharmacogenes load from your own
definition files. Clinical report generation (text, JSON, HTML) is included.

## High-level features

- **Star-allele matching** — maps observed variant identifiers (rsIDs or
  `chr:pos:ref:alt`) to star-allele names with activity values (built-in tables, each
  counting `*1`: CYP2D6 12 alleles, CYP2C19 5, CYP2C9 7, CYP3A5 4, DPYD 5, TPMT 5,
  NUDT15 3, SLCO1B1 4).
- **Diplotype & phenotype** — alleles carry CPIC activity values (normal function 1.0,
  decreased 0.5, no function 0.0; CYP2D6 `*9`, `*10` and `*41` are 0.25). The phenotype is
  assigned the way CPIC defines it for each gene: from the summed activity score for
  CYP2D6, CYP2C9, DPYD, CYP3A5, TPMT and NUDT15, and from the pair of allele functions for
  CYP2C19, SLCO1B1 and UGT1A1. Only CYP2C19 has a rapid-metabolizer category; CYP2D6 has
  none. See [Phenotype definitions](SPEC.md#phenotype-definitions).
- **Drug–gene recommendations** — 9 built-in CPIC guidelines, all level A (codeine,
  tramadol, tamoxifen / CYP2D6; clopidogrel, voriconazole / CYP2C19; warfarin / CYP2C9;
  fluorouracil / DPYD; simvastatin / SLCO1B1; azathioprine / TPMT), accessible via
  `get_dosing_recommendation(drug, phenotype)`. Further guidelines load from a file with
  `load_cpic_guidelines(filepath=...)`.
- **PharmGKB annotations** — a built-in set of 9 clinical annotations, 10 variant
  annotations and 7 drug pathways is queried offline; `parse_clinical_annotations()`
  reads your own PharmGKB TSV export. The module makes no network requests.
- **Clinical reporting** — text / JSON / HTML reports with a UTC-timestamped header,
  per-gene genotype results, drug recommendations, an interaction summary, prioritized
  action items and a disclaimer (the HTML report omits the interaction summary). Each gene
  is classified from its diplotype (or taken from a given phenotype), and with `drugs=None`
  every drug that has a built-in CPIC guideline is analysed for the genotyped genes; see
  [reporting.md](reporting.md).
- **Batch use** — there is no dedicated batch API. The per-sample functions are plain
  Python, so cohorts are processed by mapping them with
  `metainformant.core.execution.parallel` (`thread_map`, `process_map`).

## Installation

No extra is needed (there is no `pharmacogenomics` extra): the module uses only the
standard library and core dependencies (numpy and matplotlib, for the plots).

```bash
uv pip install metainformant          # or `uv pip install -e .` from a checkout
```

## Core public API

The package `__init__` files expose submodules (plus `predict_metabolizer`), so import
functions from the module that defines them.

```python
from metainformant.pharmacogenomics.alleles.star_allele import call_star_alleles   # variants → list[StarAllele]
from metainformant.pharmacogenomics.alleles.diplotype import determine_diplotype   # 2 alleles + gene → Diplotype
from metainformant.pharmacogenomics.alleles.phenotype import classify_phenotype    # Diplotype or "*1/*4" → dict
from metainformant.pharmacogenomics.annotations.cpic import (
    lookup_drug_gene,           # (drug, gene) → guideline dict, or None
    get_dosing_recommendation,  # (drug, phenotype) → recommendation dict, or None
)
from metainformant.pharmacogenomics import predict_metabolizer   # {"diplotype": "*1/*4"} + gene → dict

from metainformant.pharmacogenomics.clinical.reporting import (
    generate_clinical_report,   # (patient_data, genotypes, drugs=None) → report dict
    export_report,              # (report, format="text"|"json"|"html") → str
)
from metainformant.pharmacogenomics.clinical.pathogenicity import classify_variant_acmg

from metainformant.pharmacogenomics.annotations import cpic, pharmgkb
```

`predict_metabolizer` takes a genotype dict holding a `diplotype` string (or an `alleles`
list), not a variant set, and returns a dict whose `phenotype` is a lowercase string
(`'poor'`, `'intermediate'`, `'normal'`, `'rapid'`, `'ultrarapid'`). It is a string-and-dict
front end over the same allele table and phenotype rules as `classify_phenotype`, so the two
give the same phenotype and activity score for every diplotype (for CYP2D6 `*1/*4` both give
intermediate, activity score 1.0). They differ in the interface only: `predict_metabolizer`
raises `ValueError` for an allele that is missing from its table, whereas `classify_phenotype`
scores an unknown allele as 1.0 (normal function) and logs a warning, and `predict_metabolizer`
accepts a custom allele table through `activity_scores=`. It covers the 10 genes of
`default_allele_function_table()`, including UGT1A1 and CYP1A2 (which has no CPIC phenotype
definition, so the CYP2D6 bands apply with a logged warning).

## Typical single-patient flow

```python
variants = {'rs3892097', 'rs1065852'}          # CYP2D6 *4 and *10
alleles = call_star_alleles(variants, gene='CYP2D6')       # [*10, *4], sorted by name
dip     = determine_diplotype(alleles[0], alleles[1], gene='CYP2D6')
pheno   = classify_phenotype(dip, gene='CYP2D6')           # dict; pheno['phenotype'] is an enum
rec     = get_dosing_recommendation('codeine', pheno['phenotype'].value)
print(rec['recommendation'])   # Use codeine per standard of care. Monitor for reduced efficacy.
```

`call_star_alleles` returns every allele it matched, so a variant set that matches one
allele (or none, which returns `*1`) gives a one-element list; pad it with `*1` before
indexing a second allele. `get_dosing_recommendation` has no `gene` argument and returns
`None` when no guideline matches the drug and phenotype.

## Module layout

```
pharmacogenomics/
  alleles/             star_allele.py  diplotype.py  phenotype.py
  annotations/         cpic.py  pharmgkb.py  drug_labels.py
  clinical/            pathogenicity.py  drug_interaction.py  reporting.py
  interaction/         drug_interactions.py
  metabolism/          metabolizer_status.py
  visualization/       plots.py
```

## Data sources (bundled)

The built-in data are Python literals in the source files; nothing is downloaded or
cached.

| Resource | Built-in content | Override |
|----------|------------------|----------|
| Star-allele tables (CPIC-derived) | 8 genes, 45 alleles. UGT1A1 and CYP1A2 have activity values but no calling definitions: `load_allele_definitions` says so instead of reporting an unknown gene | `load_allele_definitions(gene, filepath=…)`, passed as `allele_definitions=` |
| CPIC guidelines | 9 guidelines, all level A | `load_cpic_guidelines(filepath=…)`, passed as `guidelines=` |
| PharmGKB annotations | 9 clinical annotations, 10 variant annotations, 7 drug pathways | `parse_clinical_annotations(filepath)` on your own TSV export |
