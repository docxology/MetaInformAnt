# Alleles

Pharmacogenomics allele calling subpackage providing star allele identification from observed variants, diplotype determination with activity scoring, and metabolizer phenotype prediction per CPIC guidelines.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports star_allele, diplotype, phenotype submodules |
| `star_allele.py` | Star allele calling against CPIC allele definition tables |
| `diplotype.py` | Diplotype determination and activity score calculation |
| `phenotype.py` | Metabolizer phenotype prediction from activity scores |

## Key Functions

| Function | Description |
|----------|-------------|
| `star_allele.call_star_alleles()` | Call star alleles from observed variants |
| `star_allele.load_allele_definitions()` | Load allele definition tables (built-in or file) |
| `star_allele.detect_novel_alleles()` | Detect novel alleles not in reference definitions |
| `star_allele.handle_cyp2d6_cnv()` | Handle CYP2D6 copy number variation |
| `diplotype.determine_diplotype()` | Determine diplotype from two star alleles |
| `diplotype.calculate_activity_score()` | Calculate combined activity score for a diplotype |
| `diplotype.resolve_ambiguous_diplotypes()` | Resolve ambiguous diplotype assignments |
| `phenotype.classify_phenotype()` | Map activity score to metabolizer phenotype |
| `phenotype.predict_metabolizer_status()` | Predict PM/IM/NM/RM/UM from genotype data |
| `phenotype.population_phenotype_frequencies()` | Estimate phenotype frequencies from allele data |

## Usage

```python
from metainformant.pharmacogenomics.alleles import star_allele, diplotype, phenotype

alleles = star_allele.call_star_alleles("CYP2D6", observed_variants)
dip = diplotype.determine_diplotype("CYP2D6", alleles)
pheno = phenotype.classify_phenotype(dip.activity_score, "CYP2D6")
```
