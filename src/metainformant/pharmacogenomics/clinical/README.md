# Clinical

Clinical pharmacogenomics subpackage providing ACMG variant classification, drug-gene interaction analysis with severity scoring, and clinical report generation in text, HTML, and JSON formats.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports pathogenicity, drug_interaction, reporting submodules |
| `pathogenicity.py` | ACMG/AMP 5-tier variant classification (PVS, PS, PM, PP, BA, BS, BP) |
| `drug_interaction.py` | Drug-gene interaction assessment and contraindication checking |
| `reporting.py` | Clinical pharmacogenomic report generation |

## Key Functions

| Function | Description |
|----------|-------------|
| `pathogenicity.classify_variant_acmg()` | Classify variant using ACMG 5-tier system |
| `pathogenicity.apply_acmg_criteria()` | Apply individual ACMG evidence criteria; PM2, BA1 and BS1 need an `allele_frequency` (a missing one is not read as 0) |
| `pathogenicity.aggregate_evidence()` | Combine criteria into final classification; pathogenic and benign criteria that both reach a call give VUS |
| `pathogenicity.query_clinvar()` | Look a variant up (rsID or HGVS) in a built-in table of 6 ClinVar-style records; `None` means "not in that table" |
| `pathogenicity.check_gnomad_frequency()` | Check BA1/BS1 against population frequencies; no frequency gives `None`, and a population the variant lacks raises `ValueError` |
| `drug_interaction.analyze_drug_gene_interactions()` | Assess drug-gene interactions for a patient; the phenotype may be written in any form (`"poor"`, `"PM"`, `"Poor Metabolizer"`), and an Indeterminate one gets no recommendation (a warning is logged) |
| `drug_interaction.check_contraindications()` | Check for genotype-based contraindications; an Indeterminate or unrecognised phenotype raises `ValueError` instead of reporting "no known contraindication" |
| `drug_interaction.calculate_interaction_severity()` | Score interaction severity (Major/Moderate/Minor) |
| `drug_interaction.suggest_alternatives()` | Suggest alternative drugs based on genotype |
| `reporting.generate_clinical_report()` | Generate a comprehensive clinical PGx report dict; `drugs=None` analyses every drug with a built-in CPIC guideline, `drugs=[]` none |
| `reporting.export_report()` | Export a report as text, HTML, or JSON (`format=`), optionally writing it to `output_path` |

## Usage

```python
from metainformant.pharmacogenomics.clinical import pathogenicity, drug_interaction, reporting

classification = pathogenicity.classify_variant_acmg(variant, evidence)
interactions = drug_interaction.analyze_drug_gene_interactions(drugs, genotypes)
report = reporting.generate_clinical_report(patient_data, genotypes, drugs)
reporting.export_report(report, format="html", output_path="output/pgx_report.html")
```
