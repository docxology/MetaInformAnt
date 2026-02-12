# Alignment

Pairwise and multiple sequence alignment for DNA sequences, with evolutionary distance calculations for phylogenetic analysis.

## Contents

| File | Purpose |
|------|---------|
| `distances.py` | Evolutionary distance measures (Jukes-Cantor, Kimura, p-distance, k-mer) |
| `msa.py` | Multiple sequence alignment via MUSCLE/MAFFT/ClustalW with fallback |
| `pairwise.py` | Needleman-Wunsch (global) and Smith-Waterman (local) alignment |

## Key Functions

| Function | Description |
|----------|-------------|
| `global_align()` | Needleman-Wunsch global alignment with configurable scoring |
| `local_align()` | Smith-Waterman local alignment for subsequence matching |
| `progressive_alignment()` | MSA via external tools with simple pairwise fallback |
| `generate_consensus_from_alignment()` | Derive consensus sequence from aligned sequences |
| `calculate_alignment_quality()` | Compute gap fraction, identity, and conservation scores |
| `jukes_cantor_distance()` | JC69 evolutionary distance assuming equal mutation rates |
| `kimura_distance()` | K2P distance accounting for transition/transversion bias |
| `distance_matrix()` | Build pairwise distance matrix from multiple sequences |
| `kmer_distance_matrix()` | Alignment-free distance matrix using k-mer profiles |
| `calculate_alignment_identity()` | Fraction of identical positions in an alignment |
| `find_conserved_regions()` | Detect contiguous identical regions in aligned sequences |

## Usage

```python
from metainformant.dna.alignment.pairwise import global_align, local_align
from metainformant.dna.alignment.distances import distance_matrix, jukes_cantor_distance
from metainformant.dna.alignment.msa import progressive_alignment

result = global_align("ATCGATCG", "ATCAATCG")
dm = distance_matrix({"s1": "ATCG", "s2": "ATCA", "s3": "ATGG"}, method="kimura")
aligned = progressive_alignment({"s1": "ATCG", "s2": "ATCGATCG"}, method="muscle")
```
