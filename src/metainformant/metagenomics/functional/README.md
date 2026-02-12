# Functional

Functional metagenomics analysis providing ORF prediction, HMM-based gene annotation, gene family classification, and metabolic pathway reconstruction from metagenomic sequences.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports annotation and pathways submodules |
| `annotation.py` | ORF prediction, gene annotation, gene family classification |
| `pathways.py` | Metabolic pathway reconstruction and completeness scoring |

## Key Functions

| Function | Description |
|----------|-------------|
| `annotation.predict_orfs()` | Predict open reading frames from assembled contigs |
| `annotation.annotate_genes()` | HMM-based functional annotation against COG/KEGG/Pfam |
| `annotation.classify_gene_families()` | Classify annotated genes into gene family categories |
| `pathways.reconstruct_pathways()` | Map annotations to KEGG/MetaCyc metabolic pathways |
| `pathways.calculate_pathway_completeness()` | Score pathway completeness based on required reactions |
| `pathways.compare_pathway_profiles()` | Compare pathway profiles across samples |
| `pathways.find_differential_pathways()` | Find differentially abundant pathways between groups |

## Usage

```python
from metainformant.metagenomics.functional import annotation, pathways

orfs = annotation.predict_orfs(contigs, min_length=100)
annotations = annotation.annotate_genes(orfs, database="KEGG")
pathway_results = pathways.reconstruct_pathways(annotations)
completeness = pathways.calculate_pathway_completeness(pathway_results)
```
