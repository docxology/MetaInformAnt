### DNA: Distances

Functions: `p_distance`, `jc69_distance`, `kmer_distance`, `kmer_distance_matrix`

```mermaid
flowchart LR
  Asequence1[Sequence 1] & Bsequence2[Sequence 2] --> CpDistance[p_distance]
  A & B --> Djc69Distance[jc69_distance] 
  A & B --> EkmerDistance[kmer_distance]
  FmultipleSequences[Multiple Sequences] --> GkmerDistanceMatrix[kmer_distance_matrix]
```

Example

```python
from metainformant.dna import distances

# Simple proportion distance
p_dist = distances.p_distance("ATCG", "ATGG")  # 0.25

# Jukes-Cantor 69 model distance  
jc_dist = distances.jc69_distance("ATCG", "ATGG")

# K-mer based distance
kmer_dist = distances.kmer_distance("ATCGATCG", "GCTAGCTA", k=3, metric="cosine")

# Distance matrix for multiple sequences
id_to_seq = {"A": "ATCG", "B": "GCTA", "C": "AAAA"}
matrix = distances.kmer_distance_matrix(id_to_seq, k=2)
```

Distance methods:
- **p_distance**: Simple proportion of differing sites (Hamming distance / length)
- **jc69_distance**: Jukes-Cantor 69 evolutionary model correction
- **kmer_distance**: K-mer frequency vector distance (cosine or euclidean)
- **kmer_distance_matrix**: Pairwise k-mer distances for multiple sequences
