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
| `annotation.predict_orfs()` | Predict open reading frames in a nucleotide sequence (six frames, NCBI `translation_table`) |
| `annotation.annotate_genes()` | Score proteins against user-supplied profile matrices (COG/KEGG/Pfam by profile ID prefix); no E-values are computed. `profile_coverage` and `consensus_identity` describe the best hit |
| `annotation.classify_gene_families()` | Assign genes to COG/KEGG/Pfam families from `reference_profiles` (no family without them) |
| `pathways.reconstruct_pathways()` | Map KO/EC annotations to metabolic pathways (8 built-in KEGG pathways or custom definitions) |
| `pathways.calculate_pathway_completeness()` | Score pathway completeness based on required reactions |
| `pathways.compare_pathway_profiles()` | Compare pathway profiles across samples |
| `pathways.find_differential_pathways()` | Report pathways whose mean completeness differs between groups by at least `min_diff` (descriptive; no statistical test) |

## Usage

```python
from metainformant.metagenomics.functional import annotation, pathways

orfs = annotation.predict_orfs(contig_sequence, min_length=100)
annotations = annotation.annotate_genes(protein_sequences, hmm_db=profiles)
pathway_results = pathways.reconstruct_pathways(gene_to_kos)
completeness = pathways.calculate_pathway_completeness(pathway_def, observed_kos)
```
