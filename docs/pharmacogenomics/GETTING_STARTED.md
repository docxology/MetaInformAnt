# Getting Started: Pharmacogenomics

Pharmacogenomics (PGx) predicts drug response from genetic variation using star allele
conventions and CPIC clinical guidelines. The built-in allele tables cover CYP2D6,
CYP2C19, CYP2C9, CYP3A5, DPYD, TPMT, NUDT15 and SLCO1B1; other genes load from your own
definition files.

## Quick example

```python
from metainformant.pharmacogenomics.alleles.star_allele import call_star_alleles
from metainformant.pharmacogenomics.alleles.diplotype import determine_diplotype
from metainformant.pharmacogenomics.alleles.phenotype import classify_phenotype
from metainformant.pharmacogenomics.annotations.cpic import get_dosing_recommendation

variants = {'rs3892097', 'rs1065852'}  # CYP2D6 *4 and *10
alleles = call_star_alleles(variants, gene='CYP2D6')
# → [StarAllele(name='*10', ...), StarAllele(name='*4', ...)]  (every match, sorted by name)

diplotype = determine_diplotype(alleles[0], alleles[1], gene='CYP2D6')
# → Diplotype(allele1='*10', allele2='*4', gene='CYP2D6', activity_score=0.25, ...)
#   diplotype.diplotype_string == '*10/*4'  (alleles are ordered as strings)

phenotype = classify_phenotype(diplotype, gene='CYP2D6')
# → {'gene': 'CYP2D6', 'diplotype': '*10/*4', 'activity_score': 0.25,
#    'phenotype': MetabolizerPhenotype.INTERMEDIATE, 'phenotype_abbreviation': 'IM',
#    'clinical_significance': '...'}

rec = get_dosing_recommendation('codeine', phenotype['phenotype'].value)
print(rec['recommendation'])
# → "Use codeine per standard of care. Monitor for reduced efficacy."
```

`get_dosing_recommendation(drug, phenotype)` takes no gene argument; it returns `None` when
no guideline matches the drug and phenotype.

## Installation

The pharmacogenomics module ships with built-in allele definitions and CPIC guidelines.
Reports and plots work with the core dependencies, so no optional extra is needed (there
are no `pharmacogenomics` or `deconvolution` extras):

```bash
uv pip install metainformant
```

## Data sources

| Resource | What it provides | Built-in content | Update method |
|----------|------------------|------------------|---------------|
| CPIC guidelines | Drug–gene dosing recommendations per phenotype | 9 guidelines, all level A | `load_cpic_guidelines(filepath=...)`, passed as `guidelines=` |
| PharmGKB annotations | Evidence-backed variant annotations | 9 clinical annotations, 10 variant annotations, 7 drug pathways | `parse_clinical_annotations(filepath)` on your own TSV export |
| Star allele tables | Gene→star-allele definitions with activity scores | 8 genes | `load_allele_definitions(gene, filepath=...)`, passed as `allele_definitions=` |

## Core concepts

**Star alleles** – Standardized haplotype labels (e.g., `*1`, `*4`, `*2A`) that represent
specific variant combinations. The matcher returns the alleles whose defining variants
are all present among the observed identifiers, trying the most specific alleles first
and using each observed variant for at most one allele.

**Diplotype** – The unordered pair of star alleles present on the two chromosomes
(e.g., `*1/*4`). Activity score is the sum of each allele's numeric activity.

**Phenotype** – Clinical interpretation: Poor (PM), Intermediate (IM), Normal (NM),
Rapid (RM), or Ultrarapid (UM) metabolizer. Each gene is classified the way CPIC defines it:
from the activity score for CYP2D6, CYP2C9, DPYD, CYP3A5, TPMT and NUDT15, and from the pair
of allele functions for CYP2C19, SLCO1B1 and UGT1A1. CYP2D6 has no rapid category (CYP2D6
`*1/*4`, activity score 1.0, is intermediate). See
[Phenotype definitions](SPEC.md#phenotype-definitions) for every band.

**CPIC level** – Evidence strength: A (strong), B (moderate), C (optional), D (weak
evidence or limited data); the built-in guidelines are all level A.

## Batch processing

```python-snippet
import pandas as pd
from metainformant.core.execution import parallel
from metainformant.pharmacogenomics.alleles.star_allele import call_star_alleles
from metainformant.pharmacogenomics.alleles.diplotype import determine_diplotype
from metainformant.pharmacogenomics.alleles.phenotype import classify_phenotype

df = pd.read_csv('cohort.csv')  # requires 'rsids' column (semicolon-separated)

def predict(rsids):
    if pd.isna(rsids):
        return None  # missing data is not "no variants" (an empty set reads as *1/*1)
    called = [a.name for a in call_star_alleles(set(str(rsids).split(';')), 'CYP2D6')]
    called += ['*1', '*1']  # pad: a single matched allele is paired with *1
    dip = determine_diplotype(called[0], called[1], 'CYP2D6')
    return classify_phenotype(dip, 'CYP2D6')['phenotype'].value

df['CYP2D6_phenotype'] = parallel.thread_map(predict, df['rsids'], max_workers=8)
```

`thread_map` uses a thread pool (`max_workers`); the matching code is plain Python, so for
CPU-bound cohorts use `parallel.process_map` with a top-level function. For massive
cohorts switch to Dask, reusing `predict` on each partition (Dask is not a core
dependency):

```python-snippet
import pandas as pd
import dask.dataframe as dd
from dask.distributed import Client

def predict_partition(part):
    return pd.Series([predict(vs) for vs in part['rsids']], index=part.index)

client = Client(n_workers=8)
ddf = dd.read_csv('cohort_*.csv')
ddf['phenotype'] = ddf.map_partitions(predict_partition, meta=('phenotype', 'object'))
```

## Clinical reporting

```python
from metainformant.pharmacogenomics.clinical.reporting import generate_clinical_report, export_report

patient = {'patient_id': 'PT001', 'sex': 'F'}
genotypes = {
    'CYP2D6': {'diplotype': '*4/*10'},
    'CYP2C19': {'diplotype': '*2/*2'},
}
report = generate_clinical_report(patient, genotypes, drugs=['codeine', 'clopidogrel'])
print(export_report(report, format='text'))   # also 'json' and 'html'
```

The text export prints a header with the generation time, patient information, genotype
results per gene, drug recommendations with CPIC levels, prioritized clinical action
items, an interaction summary and the disclaimer.

The report is a dict with the keys `header`, `patient_info`, `genotype_results`,
`drug_recommendations`, `interaction_summary`, `clinical_actions` and `disclaimer`.
`patient_info` keeps only `patient_id`, `sex`, `ethnicity` and `ordering_provider` (names
and dates of birth are not copied). Each gene's phenotype is classified from its `diplotype`
(or taken from its `phenotype` key when there is no diplotype), and the drug
recommendations use those phenotypes, so a gene needs only one of the two. With `drugs=None`
(the default) every drug that has a built-in CPIC guideline is analysed against the genotyped
genes; a list restricts the report to those drugs, and `drugs=[]` requests none. See
[reporting.md](reporting.md).

## Configuration

The module has no configuration layer: no `pharmacogenomics.*` keys, no `PG_*`
environment variables and no YAML section. Behaviour is set through function arguments
(custom allele definitions, guidelines, scoring tables, report format); see
[CONFIGURATION.md](CONFIGURATION.md) for the list. To watch the allele matcher work,
raise the level of the standard `logging` logger:

```python
import logging
logging.getLogger('metainformant.pharmacogenomics').setLevel(logging.DEBUG)
```

## What genes are supported?

The built-in allele definitions cover 8 pharmacogenes, 45 alleles in total (counting `*1`):

- CYP2D6 (12 alleles), CYP2C19 (5), CYP2C9 (7), CYP3A5 (4)
- TPMT (5), NUDT15 (3)
- SLCO1B1 (4), DPYD (5)

UGT1A1 and CYP1A2 have activity values but no calling definitions. UGT1A1 can be classified from allele names called elsewhere; its promoter TA-repeat alleles cannot be expressed by presence-only variant identifiers. CYP1A2 can only be scored: no supported CPIC phenotype definition is bundled, so phenotype APIs raise `ValueError`. Asking for either gene's calling definitions raises a `ValueError` describing the gap.

Other genes need a definition file (see [CONFIGURATION.md](CONFIGURATION.md)). To list the
built-in genes, ask for an unknown one; the error names every gene it has:

```python
from metainformant.pharmacogenomics.alleles.star_allele import load_allele_definitions

try:
    load_allele_definitions('NOT_A_GENE')
except ValueError as exc:
    print(exc)  # "... Available genes: CYP2C19, CYP2C9, CYP2D6, CYP3A5, DPYD, NUDT15, SLCO1B1, TPMT. ..."
```

## Next steps

- **[CAPABILITIES.md](CAPABILITIES.md)** — Full class and function catalogue
- **[EXAMPLES.md](EXAMPLES.md)** — 22 detailed usage patterns
- **[PERFORMANCE.md](PERFORMANCE.md)** — Benchmarks and scaling tips
- **[TROUBLESHOOTING.md](TROUBLESHOOTING.md)** — Known failure modes and fixes
- **[CONFIGURATION.md](CONFIGURATION.md)** — Function arguments that customise tables, guidelines and reports; logging
