# CPIC Guidelines Integration

CPIC (Clinical Pharmacogenetics Implementation Consortium) publishes peer‑reviewed,
evidence‑based guidelines linking gene–drug pairs to actionable prescribing guidance.
This module bundles the guideline JSON and provides fast lookup.

## Guideline lifecycle

- **Level A** — action required (strong evidence; test recommended before prescribing)
- **Level B** — action recommended (moderate evidence)
- **Level C** — optional — may affect decision
- **Level D** — not recommended — minimal actionability

## Loading guidelines

```python
from metainformant.pharmacogenomics.annotations.cpic import (
    load_cpic_guidelines,          # built‑in snapshot
    get_dosing_recommendation,     # single lookup
    lookup_drug_gene,              # filter by level, gene, or drug
)

# Built‑in v3.0 (2024-01) loaded lazily on first call
guidelines = load_cpic_guidelines()
len(guidelines)  # ≈ 200 entries
```

## Single‑drug lookup

```python
rec = get_dosing_recommendation(
    drug='codeine',
    phenotype='Poor Metabolizer',
    gene='CYP2D6',
)
print(rec['recommendation'])
# → 'Avoid codeine; recommend alternative analgesic'
print(rec['cpic_level'])   # 'A'
print(rec['classification'])  # 'strong'
```

## Bulk queries

```python
# All Level A actionable pairs
level_a = lookup_drug_gene(level='A')
for entry in level_a:
    print(entry['drug'], entry['gene'], entry['phenotypes'].keys())
```

## Overriding / upgrading

From a downloaded newer CPIC JSON:

```python
load_cpic_guidelines(path='~/cpic_guidelines_2025.json')
```

Or set environment variable:

```bash
export PG_CPIC_GUIDELINES_URL=https://cpicpgx.org/guidelines/2025/cpic_guidelines.json
```

The module downloads once and caches under `~/.cache/metainformant/cpic/`.

## Schema validation

Each entry must contain:

```json
{
  "drug": "warfarin",
  "gene": "CYP2C9",
  "cpic_level": "A",
  "phenotypes": {
    "Poor Metabolizer": {
      "recommendation": "Reduce warfarin starting dose by …",
      "classification": "strong"
    }
  },
  "guideline_url": "https://cpicpgx.org/guidelines/…",
  "version": "3.0",
  "published": "2024-01-15"
}
```

Missing top‑level keys raise `ValueError` on load.

## Multi‑drug, multi‑phenotype matrix

```python
import pandas as pd
drugs = ['codeine','tamoxifen','warfarin','clopidogrel']
genes = ['CYP2D6','CYP2C19','CYP2C9']
rows = []
for drug in drugs:
    for gene in genes:
        for ph in ['PM','IM','NM','RM','UM']:
            try:
                rec = get_dosing_recommendation(drug, ph, gene=gene)
                rows.append({'drug':drug,'gene':gene,'phenotype':ph,
                             'recommendation':rec['recommendation'],
                             'level':rec['cpic_level']})
            except KeyError:
                continue  # no guidance for this combo
df = pd.DataFrame(rows)
df.to_csv('cpic_matrix.csv', index=False)
```

## Integration with clinical report

`generate_clinical_report()` automatically iterates over the patient's genes,
calls `get_dosing_recommendation()` for each actionable drug, and inserts a
table into the final report.

## Latest snapshot info

Current built‑in: **CPIC v3.0** published 2024‑01, covering 33 genes, 208 drug–
gene pairs. Last manual update in repo: 2024‑02. Newer versions can be loaded at
runtime without restart.

## Citations

When using CPIC data in publications, cite the appropriate guideline DOI and
acknowledge CPIC. See `cpic/CITATION.cff` in the source distribution.
