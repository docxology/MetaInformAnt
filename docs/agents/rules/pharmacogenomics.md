# Pharmacogenomics Module Rules

## Purpose
Clinical pharmacogenomics analysis with star allele calling, diplotype determination, metabolizer phenotype prediction, CPIC/PharmGKB/ACMG annotation, drug-gene interactions, and clinical reporting.

## Source Structure
```
src/metainformant/pharmacogenomics/
├── alleles/
│   ├── star_allele.py   # Star allele calling, definition matching, CYP2D6 CNV
│   ├── diplotype.py     # Diplotype determination, activity score calculation
│   └── phenotype.py     # Metabolizer status prediction, population frequencies
├── annotations/
│   ├── cpic.py          # CPIC guideline loading, drug-gene lookup, dosing
│   ├── pharmgkb.py      # PharmGKB queries, clinical annotations, evidence levels
│   └── drug_labels.py   # FDA drug label parsing, biomarker extraction
├── clinical/
│   ├── pathogenicity.py    # ACMG variant classification, ClinVar/gnomAD queries
│   ├── drug_interaction.py # Drug-gene interactions, contraindications, polypharmacy
│   └── reporting.py        # Clinical report generation, disclaimers
└── visualization/
    └── plots.py          # Metabolizer status, allele frequencies, ACMG criteria charts
```

## Dependencies
- **Required**: numpy, pandas
- **Optional**: ClinVar API, gnomAD API, CPIC database, PharmGKB API

## Import Patterns
```python
from metainformant.pharmacogenomics.alleles import star_allele, diplotype, phenotype
from metainformant.pharmacogenomics.annotations import cpic, pharmgkb, drug_labels
from metainformant.pharmacogenomics.clinical import pathogenicity, drug_interaction, reporting
```

## Configuration
- Environment prefix: `PHARMA_` (e.g., `PHARMA_THREADS`, `PHARMA_DB_PATH`)
- Output path: `output/pharmacogenomics/<analysis_type>/` (e.g., `alleles/`, `clinical/`, `reports/`)

## Integration
- **Pharmacogenomics → GWAS**: Variant data from association studies
- **Pharmacogenomics → DNA**: Genomic coordinates and variant calling
- **Pharmacogenomics → Phenotype**: Clinical phenotype data

## Testing
- Use `@pytest.mark.network` for tests requiring ClinVar/PharmGKB APIs
- Generate test VCF and allele data programmatically
- All test outputs to `tmp_path`
