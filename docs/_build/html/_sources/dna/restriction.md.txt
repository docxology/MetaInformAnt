### DNA: Restriction Sites

Function: `find_restriction_sites`

```mermaid
flowchart LR
  AdnaSequence[DNA Sequence] & Benzyme->motifDict[Enzyme->Motif Dict] --> CfindRestrictionSites[find_restriction_sites]
  C --> Denzyme->positionsDict[Enzyme->Positions Dict]
```

Example

```python
from metainformant.dna import restriction

# Define restriction enzymes and their recognition motifs
enzymes = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC", 
    "HindIII": "AAGCTT",
    "NotI": "GCGGCCGC"
}

seq = "ATCGAATTCGATGGATCCTAG"
sites = restriction.find_restriction_sites(seq, enzymes)
# Returns: {"EcoRI": [4], "BamHI": [13], "HindIII": [], "NotI": []}
```

Features:
- **Multiple enzymes**: Process many restriction enzymes simultaneously
- **IUPAC motifs**: Supports ambiguous nucleotide codes in recognition sequences
- **Forward strand only**: Searches only the provided sequence orientation
- **Zero-based positions**: Returns 0-based start positions of recognition sites
- **Complete results**: Returns empty list for enzymes with no sites found

Common restriction enzymes:
- **Type II**: Most common, recognize palindromic sequences
- **Rare cutters**: 8+ base recognition (NotI, SfiI, AscI)
- **Frequent cutters**: 4-6 base recognition (MspI, HaeIII, TaqI)

Depends on: [motifs](./motifs.md) for IUPAC pattern matching.
