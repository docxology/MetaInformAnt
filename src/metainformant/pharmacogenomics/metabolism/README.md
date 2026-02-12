# Metabolism

Metabolizer status prediction from pharmacogenomic genotype data, providing CPIC activity score computation and evidence-based dose adjustment recommendations for CYP2D6, CYP2C19, and CYP2C9.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports metabolizer_status submodule |
| `metabolizer_status.py` | Activity score framework, metabolizer classification, dose adjustment |

## Key Functions

| Function | Description |
|----------|-------------|
| `metabolizer_status.default_allele_function_table()` | Get built-in CPIC allele function assignments |
| `metabolizer_status.compute_activity_score()` | Compute activity score from diplotype alleles |
| `metabolizer_status.classify_metabolizer()` | Classify metabolizer phenotype from activity score |
| `metabolizer_status.predict_metabolizer_status()` | End-to-end metabolizer prediction from genotype |
| `metabolizer_status.dose_adjustment()` | Get dose adjustment recommendation for a drug/phenotype |

## Usage

```python
from metainformant.pharmacogenomics.metabolism import metabolizer_status

score = metabolizer_status.compute_activity_score("*1", "*4", gene="CYP2D6")
phenotype = metabolizer_status.classify_metabolizer(score, gene="CYP2D6")
adjustment = metabolizer_status.dose_adjustment("codeine", phenotype)
```
