### DNA: Mutations

Functions: `apply_point_mutations`, `hamming_distance`, `random_point_mutations`

```mermaid
flowchart TD
  AoriginalSequence[Original Sequence] --> BapplyPointMutations[apply_point_mutations]
  CmutationMap[Mutation Map] --> B
  B --> DmutatedSequence[Mutated Sequence]
  
  Esequence1[Sequence 1] & Fsequence2[Sequence 2] --> GhammingDistance[hamming_distance]
  
  HoriginalSequence[Original Sequence] --> IrandomPointMutations[random_point_mutations]
  Jnumber&Seed[Number & Seed] --> I
  I --> KrandomlyMutatedSequence[Randomly Mutated Sequence]
```

Example

```python
from metainformant.dna import mutations

# Apply specific point mutations
seq = "ATCGATCG"
changes = {0: "G", 3: "A"}  # Position -> new base
mutated = mutations.apply_point_mutations(seq, changes)  # "GTCAATCG"

# Calculate Hamming distance
seq1 = "ATCG"
seq2 = "ATGG" 
dist = mutations.hamming_distance(seq1, seq2)  # 1

# Generate random mutations
seq = "ATCGATCGATCG"
random_mutated = mutations.random_point_mutations(
    seq, 
    num_mutations=3, 
    seed=42  # For reproducibility
)
```

Features:
- **Precise mutations**: Apply specific changes at defined positions
- **Distance calculation**: Hamming distance between sequences
- **Random mutagenesis**: Introduce controlled random mutations
- **Case preservation**: Maintains original case when possible
- **Unique positions**: Random mutations avoid duplicate positions
- **Reproducible**: Seeded random number generation

Mutation types:
- **Point mutations**: Single nucleotide changes
- **Substitutions**: Replace existing base with different base
- **Random selection**: Avoids changing to same base when possible

Related: Used in evolutionary analyses and [phylogeny](./phylogeny.md) simulations.
