# Pharmacogenomics Module

Clinical pharmacogenomics analysis with star allele calling, drug-gene interactions, and clinical reporting.

## Submodules

### Alleles (`alleles/`)
- **star_allele.py** - Star allele calling, definition matching, novel allele detection, CYP2D6 CNV handling
- **diplotype.py** - Diplotype determination, activity score calculation, ambiguity resolution
- **phenotype.py** - Metabolizer status prediction, phenotype thresholds, population frequencies

### Annotations (`annotations/`)
- **cpic.py** - CPIC guideline loading, drug-gene lookup, dosing recommendations
- **pharmgkb.py** - PharmGKB queries, clinical annotations, evidence levels, drug pathways
- **drug_labels.py** - FDA drug label parsing, biomarker extraction, label classification

### Clinical (`clinical/`)
- **pathogenicity.py** - ACMG variant classification, ClinVar/gnomAD queries, evidence aggregation
- **drug_interaction.py** - Drug-gene interaction analysis, contraindications, polypharmacy, alternatives
- **reporting.py** - Clinical report generation, recommendation formatting, disclaimers

### Visualization (`visualization/`)
- **plots.py** - Metabolizer status charts, allele frequencies, activity scores, drug response heatmaps, ACMG criteria

## Usage

```python
from metainformant.pharmacogenomics.alleles import star_allele, diplotype, phenotype
from metainformant.pharmacogenomics.annotations import cpic, pharmgkb
from metainformant.pharmacogenomics.clinical import pathogenicity, drug_interaction, reporting
```
