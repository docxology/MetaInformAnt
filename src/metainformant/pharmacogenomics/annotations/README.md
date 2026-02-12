# Annotations

Pharmacogenomics annotation subpackage providing CPIC guideline lookups, PharmGKB clinical annotation queries, and FDA drug label parsing for pharmacogenomic biomarker information.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports cpic, pharmgkb, drug_labels submodules |
| `cpic.py` | CPIC guideline data, drug-gene lookups, dosing recommendations |
| `pharmgkb.py` | PharmGKB clinical annotations, evidence levels, variant annotations |
| `drug_labels.py` | FDA drug label parsing, biomarker extraction, PGx classification |

## Key Functions

| Function | Description |
|----------|-------------|
| `cpic.lookup_drug_gene()` | Look up CPIC drug-gene pair information |
| `cpic.get_dosing_recommendation()` | Get dosing recommendation for drug/gene/phenotype |
| `cpic.list_actionable_genes()` | List all CPIC-actionable pharmacogenes |
| `cpic.parse_cpic_allele_definitions()` | Parse CPIC allele definition tables |
| `pharmgkb.query_pharmgkb_annotations()` | Query PharmGKB clinical annotations by gene or drug |
| `pharmgkb.get_evidence_level()` | Get PharmGKB evidence level for a variant-drug pair |
| `pharmgkb.get_variant_annotations()` | Get variant-level annotations from PharmGKB |
| `drug_labels.parse_drug_label()` | Parse FDA drug label for PGx information |
| `drug_labels.extract_biomarker_info()` | Extract biomarker details from label text |
| `drug_labels.search_labels_by_gene()` | Search drug labels by pharmacogene |

## Usage

```python
from metainformant.pharmacogenomics.annotations import cpic, pharmgkb, drug_labels

rec = cpic.get_dosing_recommendation("codeine", "CYP2D6", "Poor Metabolizer")
annots = pharmgkb.query_pharmgkb_annotations(gene="CYP2D6")
labels = drug_labels.search_labels_by_gene("HLA-B")
```
