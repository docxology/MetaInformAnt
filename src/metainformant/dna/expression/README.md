# Gene Expression Analysis

Tools for modeling and analyzing gene expression processes from DNA sequences.

## Purpose

This module provides:
- Transcription simulation
- Translation with customizable genetic codes
- Codon usage bias analysis

## Key Components

| File | Description |
|------|-------------|
| [transcription.py](transcription.py) | DNA→RNA transcription models |
| [translation.py](translation.py) | RNA→Protein translation with genetic codes |
| [codon.py](codon.py) | Codon usage statistics and bias analysis |

## Usage

```python
from metainformant.dna.expression import translate_sequence, codon_usage

protein = translate_sequence(mrna_seq, genetic_code="standard")
usage = codon_usage(dna_seq)
```

## Related Documentation

- **Parent**: [src/metainformant/dna/README.md](../README.md)
- **SPEC**: [SPEC.md](SPEC.md)
- **AGENTS**: [AGENTS.md](AGENTS.md)
- **Sequence Module**: [../sequence/README.md](../sequence/README.md)
