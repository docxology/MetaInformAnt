# Sequence Alignment

Tools for aligning DNA and protein sequences using dynamic programming and external aligners.

## Purpose

This module provides functions for:
- Pairwise alignment (global and local)
- Multiple sequence alignment (MSA) via external tools
- Calculation of evolutionary distances for phylogenetic analysis

## Key Components

| File | Description |
|------|-------------|
| [pairwise.py](pairwise.py) | Needleman-Wunsch (global) and Smith-Waterman (local) alignment |
| [msa.py](msa.py) | Interfaces for MUSCLE, Clustal, MAFFT |
| [distances.py](distances.py) | Jukes-Cantor, Kimura, and p-distance calculations |

## Usage

```python
from metainformant.dna.alignment import global_align, local_align

result = global_align("ACGT", "ACGT")
print(result.score)
```

## Related Documentation

- **Parent**: [src/metainformant/dna/README.md](../README.md)
- **SPEC**: [SPEC.md](SPEC.md)
- **AGENTS**: [AGENTS.md](AGENTS.md)
- **Phylogeny Module**: [src/metainformant/dna/phylogeny/README.md](../phylogeny/README.md)
