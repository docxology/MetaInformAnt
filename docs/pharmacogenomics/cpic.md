# CPIC Guidelines Integration

CPIC (Clinical Pharmacogenetics Implementation Consortium) publishes peer-reviewed,
evidence-based guidelines linking gene-drug pairs to prescribing guidance.
`metainformant.pharmacogenomics.annotations.cpic` bundles 9 of them as Python literals and
provides lookup functions. It makes no network requests and caches nothing.

## CPIC levels

- **A**: strong genetic evidence and a strong prescribing action recommended
- **A/B**: strong evidence, moderate action
- **B**: moderate genetic evidence and a prescribing action recommended
- **C**: moderate evidence, optional action
- **D**: weak evidence or limited data

All 9 built-in guidelines are level A.

## Built-in guidelines

| Drug | Gene | Phenotypes with a recommendation |
|------|------|----------------------------------|
| codeine, tramadol, tamoxifen | CYP2D6 | Poor, Intermediate, Normal, Ultrarapid |
| clopidogrel, voriconazole | CYP2C19 | Poor, Intermediate, Normal, Rapid, Ultrarapid |
| warfarin | CYP2C9 | Poor, Intermediate, Normal |
| fluorouracil | DPYD | Poor, Intermediate, Normal |
| azathioprine | TPMT | Poor, Intermediate, Normal |
| simvastatin | SLCO1B1 | Poor, Intermediate, Normal |

Each guideline is a dict with `drug`, `gene`, `cpic_level`, `guideline_url` and
`recommendations`, which maps a phenotype name (`'Poor Metabolizer'`) to
`{'recommendation', 'classification', 'implications'}`.

## Loading and looking up guidelines

```python
from metainformant.pharmacogenomics.annotations.cpic import (
    load_cpic_guidelines,       # the built-in list, or a list parsed from your file
    lookup_drug_gene,           # (drug, gene) -> guideline dict, or None
    get_dosing_recommendation,  # (drug, phenotype) -> recommendation dict, or None
    list_actionable_genes,      # [{'gene', 'drug', 'cpic_level'}, ...]
)

guidelines = load_cpic_guidelines()
print(len(guidelines))          # 9

guideline = lookup_drug_gene('codeine', 'CYP2D6')
print(guideline['cpic_level'])  # A
```

`lookup_drug_gene` and `get_dosing_recommendation` compare names case-insensitively.

## Single-drug lookup

`get_dosing_recommendation(drug, phenotype, guidelines=None)` takes no gene argument: it
searches every guideline for the drug and returns the first recommendation whose phenotype
matches. The phenotype is matched by its full name (`'Poor Metabolizer'`) or by a two-letter
abbreviation (`'PM'`, `'IM'`, `'NM'`, `'RM'`, `'UM'`; the old `'EM'` counts as normal), in any
case. Bare tokens such as `'poor'` are also normalized by `MetabolizerPhenotype.from_string` and match the same guideline as a full name or abbreviation.

```python
rec = get_dosing_recommendation('codeine', 'Poor Metabolizer')
print(rec['recommendation'])
# 'Avoid codeine use. Use alternative analgesics not metabolized by CYP2D6.'
print(rec['cpic_level'])        # 'A'
print(rec['classification'])    # 'Strong'
```

The returned dict has `drug`, `gene`, `cpic_level`, `phenotype`, `recommendation`,
`classification`, `implications` and `guideline_url`. The function returns `None` when no
guideline matches the drug and phenotype. To restrict the search to one gene, look the
guideline up first and pass it in:

```python
guideline = lookup_drug_gene('warfarin', 'CYP2C9')
rec = get_dosing_recommendation('warfarin', 'PM', guidelines=[guideline])
```

## Bulk queries

```python
# Every gene-drug pair at level A, sorted by gene and drug
for entry in list_actionable_genes(min_level='A'):
    print(entry['gene'], entry['drug'], entry['cpic_level'])
```

`list_actionable_genes(guidelines=None, min_level='B')` returns pairs at the given level or
better (`'B'` also includes `'A/B'`). Unrecognized `min_level` values raise `ValueError`.

## Using your own guidelines

`load_cpic_guidelines(filepath=...)` parses a JSON or TSV file and returns the list; it does
not replace the built-in data, so pass the list to the functions that take `guidelines`:

```python-snippet
guidelines = load_cpic_guidelines(filepath='cpic_guidelines_2025.json')
rec = get_dosing_recommendation('codeine', 'Poor Metabolizer', guidelines=guidelines)
```

The JSON file is a list of guideline dicts (an object with a `guidelines` list, or a single
guideline object, is also accepted). A TSV with the columns `drug`, `gene`, `cpic_level` (or
`level`) and `guideline_url` (or `url`) loads too, but its entries have no `recommendations`.
`generate_clinical_report` and `analyze_drug_gene_interactions` always use the built-in
guidelines. See [CONFIGURATION.md](CONFIGURATION.md#cpic-guidelines-override) for the schema
and the errors a missing `cpic_level` causes.

## Guideline matrix

```python-snippet
import pandas as pd

drugs = ['codeine', 'tamoxifen', 'warfarin', 'clopidogrel']
phenotypes = ['PM', 'IM', 'NM', 'RM', 'UM']
rows = []
for drug in drugs:
    for ph in phenotypes:
        rec = get_dosing_recommendation(drug, ph)
        if rec is None:
            continue  # no guidance for this combination
        rows.append({'drug': drug, 'gene': rec['gene'], 'phenotype': ph,
                     'recommendation': rec['recommendation'], 'level': rec['cpic_level']})
pd.DataFrame(rows).to_csv('cpic_matrix.csv', index=False)
```

## Integration with clinical report

`generate_clinical_report()` classifies each gene's phenotype, runs `analyze_drug_gene_interactions()`
for the requested drugs (all built-in guideline drugs when `drugs=None`), and looks each
recommendation up in these guidelines. See [reporting.md](reporting.md).

## Citations

When using CPIC data in publications, cite the appropriate CPIC guideline and acknowledge
CPIC; each built-in guideline carries its `guideline_url`.
