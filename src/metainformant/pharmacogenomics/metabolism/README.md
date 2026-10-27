# Metabolism

Metabolizer status prediction from pharmacogenomic genotype data, providing CPIC activity score computation, phenotype classification and evidence-based dose adjustment recommendations. It is a string and dict front end over the `alleles` subpackage and holds no allele or phenotype data of its own: the allele table (CYP2D6, CYP2C19, CYP2C9, CYP3A5, CYP1A2, DPYD, TPMT, NUDT15, SLCO1B1 and UGT1A1) is a copy of `alleles.diplotype._ACTIVITY_SCORE_TABLES` and classification is `alleles.phenotype`'s, so a phenotype here always equals the one `alleles.phenotype.classify_phenotype` gives. Dose guidelines cover codeine, tramadol, tamoxifen, clopidogrel, warfarin, fluorouracil, azathioprine and omeprazole.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports metabolizer_status submodule |
| `metabolizer_status.py` | Activity score framework, metabolizer classification, dose adjustment |

## Key Functions

| Function | Description |
|----------|-------------|
| `metabolizer_status.default_allele_function_table()` | Get an independent copy of the built-in CPIC allele function assignments |
| `metabolizer_status.compute_activity_score()` | Compute activity score from a diplotype string and an allele table |
| `metabolizer_status.classify_metabolizer()` | Classify a gene's metabolizer phenotype token from an activity score (negative or NaN raises `ValueError`) |
| `metabolizer_status.predict_metabolizer_status()` | End-to-end metabolizer prediction from a genotype dict (also exported as `pharmacogenomics.predict_metabolizer`) |
| `metabolizer_status.dose_adjustment()` | Get the dose adjustment for a phenotype (any common form: `"poor"`, `"Poor Metabolizer"`, `"PM"`) and a drug |

## Usage

```python
from metainformant.pharmacogenomics.metabolism import metabolizer_status

table = metabolizer_status.default_allele_function_table()["CYP2D6"]
score = metabolizer_status.compute_activity_score("*1/*4", table)  # 1.0
phenotype = metabolizer_status.classify_metabolizer(score, "CYP2D6")  # "intermediate"
adjustment = metabolizer_status.dose_adjustment(phenotype, "codeine")

# or in one step from a genotype dict
result = metabolizer_status.predict_metabolizer_status({"diplotype": "*1/*4"}, "CYP2D6")
```
