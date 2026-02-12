# Protein Sequence

Core protein sequence processing: FASTA I/O, physicochemical property calculation, pairwise and multiple sequence alignment, and proteome-level analysis.

## Contents

| File | Purpose |
|------|---------|
| `sequences.py` | FASTA reading/writing, amino acid composition, molecular weight, pI, hydropathy |
| `alignment.py` | Needleman-Wunsch and Smith-Waterman alignment, identity scoring, MSA |
| `proteomes.py` | Taxonomy ID handling, proteome download from UniProt, cross-proteome comparison |

## Key Functions

| Function | Description |
|----------|-------------|
| `read_fasta()` | Parse protein sequences from FASTA file into dict |
| `write_fasta()` | Write sequence dict to FASTA file |
| `molecular_weight()` | Calculate protein molecular weight from sequence |
| `isoelectric_point()` | Estimate isoelectric point (pI) |
| `hydropathy_score()` | Kyte-Doolittle hydropathy windowed scores |
| `amino_acid_composition()` | Fractional amino acid frequencies |
| `instability_index()` | Guruprasad instability index |
| `global_align()` | Needleman-Wunsch global alignment of two sequences |
| `local_align()` | Smith-Waterman local alignment |
| `multi_sequence_alignment()` | Progressive MSA for protein families |
| `read_taxon_ids()` | Load taxonomy IDs from file |
| `download_proteome_fasta()` | Download proteome FASTA from UniProt by taxon |
| `compare_proteomes()` | Statistical comparison of two proteome FASTA files |

## Usage

```python
from metainformant.protein.sequence.sequences import read_fasta, molecular_weight
from metainformant.protein.sequence.alignment import global_align

seqs = read_fasta("proteins.fasta")
mw = molecular_weight("MKLVLSDEL")
result = global_align("MKLVL", "MKLAL")
```
