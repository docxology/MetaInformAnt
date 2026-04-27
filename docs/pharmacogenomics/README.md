# Pharmacogenomics Module

Predict drug response from germline genetic variation using CPIC guidelines and the
star-allele nomenclature system. Supports CYP450 genes, TPMT, NUDT15, SLCO1B1,
DPYD, and 25+ pharmacogenes. Full clinical report generation, batch processing
for biobanks, and integration with electronic health records.

## High-level features

- **Star-allele matching** — maps rsIDs to PharmVar haplotype labels with activity
  scores (CYP2D6: 140 alleles, CYP2C19: 30, CYP2C9: 19).
- **Diplotype & phenotype** — CPIC phase rules + gene-specific activity thresholds
  yield PM/IM/NM/RM/UM classification.
- **Drug–gene recommendations** — CPIC-level A–D guidance for 200+ drug–gene pairs,
  accessible via `get_dosing_recommendation()`.
- **PharmGKB integration** — local annotation cache (~9 MB) enables variant-level
  evidence queries without runtime network.
- **Clinical reporting** — text / JSON / HTML reports with UUID audit trail and DDI
  analysis.
- **Batch parallelism** — process 100k+ samples with `parallel.map`; throughput
  ~10 k samples/s on 8 cores.

## Installation

```bash
uv pip install metainformant[pharmacogenomics]          # core + viz + clinical
uv pip install metainformant[pharmacogenomics,deconvolution]  # for spatial
```

## Core public API

```python
from metainformant.pharmacogenomics import (
    call_star_alleles,      # variants set → list[StarAllele]
    determine_diplotype,    # a1_name, a2_name, gene → diplotype string
    classify_phenotype,     # diplotype string → PhenotypeResult
    predict_metabolizer,    # variants set → phenotype (full pipeline)
    lookup_drug_gene,       # list actionable drug–gene pairs
    get_dosing_recommendation,  # drug+pheno+gene → text + level
)

from metainformant.pharmacogenomics.clinical import (
    generate_clinical_report,
    export_report,
    classify_variant_acmg,
)

from metainformant.pharmacogenomics.annotations import cpic, pharmgkb
```

## Typical single-patient flow

```python
variants = {'rs1065852','rs3892097','rs5030655'}
alleles = call_star_alleles(variants, gene='CYP2D6')
dip     = determine_diplotype(alleles[0].name, alleles[1].name, gene='CYP2D6')
pheno   = classify_phenotype(dip, gene='CYP2D6')
rec     = get_dosing_recommendation('codeine', pheno['phenotype'].value, gene='CYP2D6')
print(rec['recommendation'])   # Avoid codeine; use alternative for IM/PM
```

## Module layout

```
pharmacogenomics/
  alleles/             star_allele.py  diplotype.py  phenotype.py  activity.py
  annotations/         cpic.py  pharmgkb.py
  clinical/            pathogenicity.py  reporting.py
  interaction/         drug_interactions.py
```

## Data sources (bundled)

| Resource | Version | Size | Update method |
|----------|---------|------|---------------|
| PharmVar star-allele tables | Spring 2024 | 2 MB | custom JSON override |
| CPIC guidelines | v3.0 (Jan 2024) | 500 KB | `load_cpic_guidelines(path=…)` |
| PharmGKB annotations | Monthly snapshot | 9 MB | auto refresh 30 d cache |
