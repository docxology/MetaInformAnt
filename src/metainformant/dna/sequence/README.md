# DNA Sequence Analysis

Core tools for fundamental DNA sequence operations.

## Purpose

This module provides:
- Basic sequence string operations
- Composition analysis (GC content, skews)
- Motif discovery
- Restriction enzyme simulation
- Consensus sequence generation

## Key Components

| File | Description |
|------|-------------|
| [core.py](core.py) | `reverse_complement`, `validate_dna_sequence` |
| [composition.py](composition.py) | `gc_content`, `at_skew`, `gc_skew` |
| [motifs.py](motifs.py) | `find_motifs`, `find_repeats` |
| [restriction.py](restriction.py) | `cut_with_enzyme`, `find_restriction_sites` |
| [consensus.py](consensus.py) | `generate_consensus` from MSA |

## Usage

```python
from metainformant.dna.sequence import gc_content, reverse_complement

rc = reverse_complement("ACGT")
gc = gc_content("ACGTAGCT")
```

## Related Documentation

- **Parent**: [src/metainformant/dna/README.md](../README.md)
- **SPEC**: [SPEC.md](SPEC.md)
- **AGENTS**: [AGENTS.md](AGENTS.md)
- **Alignment Module**: [../alignment/README.md](../alignment/README.md)
