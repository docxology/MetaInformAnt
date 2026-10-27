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
- **Diplotype & phenotype** — summed allele activity scores plus gene-specific
  thresholds yield PM/IM/NM/RM/UM classification.
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
  action items and a disclaimer (the HTML report omits the interaction summary).
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
(`'poor'`, `'intermediate'`, `'normal'`, `'rapid'`, `'ultrarapid'`). It uses its own
allele-function table and thresholds, so its activity scores and phenotype calls can
differ from `classify_phenotype` (for CYP2D6 `*1/*4` it returns `'intermediate'`, while
`classify_phenotype` returns Normal Metabolizer).

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
| Star-allele tables (CPIC-derived) | 8 genes, 45 alleles | `load_allele_definitions(gene, source='file', filepath=…)`, passed as `allele_definitions=` |
| CPIC guidelines | 9 guidelines, all level A | `load_cpic_guidelines(filepath=…)`, passed as `guidelines=` |
| PharmGKB annotations | 9 clinical annotations, 10 variant annotations, 7 drug pathways | `parse_clinical_annotations(filepath)` on your own TSV export |
