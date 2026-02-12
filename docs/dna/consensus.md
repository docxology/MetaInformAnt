### DNA: Consensus Sequences

The `metainformant.dna.sequence.consensus` module provides tools for generating consensus
sequences from multiple aligned DNA sequences. It supports majority-rule voting,
IUPAC ambiguity codes, quality-weighted consensus, bootstrap confidence estimation,
and alignment conservation analysis.

---

## Functions Overview

| Function                      | Description                                                 |
|-------------------------------|-------------------------------------------------------------|
| `generate_consensus`          | Threshold-based consensus (configurable frequency cutoff)   |
| `consensus_with_ambiguity`    | Consensus using full IUPAC ambiguity codes                  |
| `majority_consensus`          | Simple majority vote (threshold=0.5)                        |
| `strict_consensus`            | 100% agreement required (ambiguity codes for disagreements) |
| `quality_weighted_consensus`  | Phred quality-weighted base calling                         |
| `consensus_from_alignment`    | Comprehensive analysis returning consensus + statistics     |
| `consensus_statistics`        | Conservation metrics for an alignment                       |
| `find_consensus_breaks`       | Sliding window detection of low-conservation regions        |
| `bootstrap_consensus`         | Bootstrap resampling with confidence score                  |

---

## Basic Consensus Generation

### Majority Rule Consensus

The simplest approach: at each position, the base appearing in more than 50% of
sequences wins. Gaps are ignored during voting.

```python
from metainformant.dna.sequence.consensus import majority_consensus

seqs = ["ATCGATCG", "ATCGATCG", "AGCGATCG"]
result = majority_consensus(seqs)
# "ATCGATCG" -- position 2: T appears 2/3 times, wins majority
```

### Threshold-Based Consensus

Customize the frequency threshold. When no base meets the threshold, an IUPAC
ambiguity code is used instead.

```python
from metainformant.dna.sequence.consensus import generate_consensus

seqs = ["ATCG", "AGCG", "ACCG"]
result = generate_consensus(seqs, threshold=0.7)
# Position 2: T(1/3), G(1/3), C(1/3) -- none meets 0.7 threshold
# Uses IUPAC ambiguity code instead
```

**Parameters:**

| Parameter    | Type         | Default | Description                              |
|--------------|--------------|---------|------------------------------------------|
| `sequences`  | `List[str]`  | required| List of aligned, equal-length sequences  |
| `threshold`  | `float`      | `0.5`   | Minimum frequency for consensus base     |

### Strict Consensus with Ambiguity Codes

Uses IUPAC codes to represent all observed bases at each position.

```python
from metainformant.dna.sequence.consensus import strict_consensus

seqs = ["ATCG", "AGCG"]
result = strict_consensus(seqs)
# Position 2: T and G observed -> "K" (IUPAC code for G/T)
```

### IUPAC Ambiguity Codes

The module handles the full IUPAC degenerate base alphabet:

| Code | Bases    | Code | Bases     |
|------|----------|------|-----------|
| R    | A, G     | V    | A, C, G   |
| Y    | C, T     | H    | A, C, T   |
| S    | G, C     | D    | A, G, T   |
| W    | A, T     | B    | C, G, T   |
| K    | G, T     | N    | A, C, G, T|
| M    | A, C     |      |           |

---

## Quality-Weighted Consensus

When Phred quality scores are available, bases are weighted by their quality to
produce a more reliable consensus.

```python
from metainformant.dna.sequence.consensus import quality_weighted_consensus

seqs = ["ATCG", "AGCG"]
quals = [
    [30, 30, 30, 30],  # High quality for all bases in seq 1
    [30, 10, 30, 30],  # Low quality at position 2 in seq 2
]
result = quality_weighted_consensus(seqs, quals)
# Position 2: T (high quality) wins over G (low quality)
```

**Parameters:**

| Parameter    | Type              | Default | Description                              |
|--------------|-------------------|---------|------------------------------------------|
| `sequences`  | `List[str]`       | required| Aligned sequences                        |
| `qualities`  | `List[List[int]]` | required| Phred quality scores (one list per seq)  |

---

## Comprehensive Alignment Analysis

### `consensus_from_alignment`

Generates both majority and strict consensus, calculates statistics, and returns
a structured analysis dictionary.

```python
from metainformant.dna.sequence.consensus import consensus_from_alignment

alignment = ["ATCG-T", "ATCG-T", "ATCGGT", "AGCG-T"]
consensus, analysis = consensus_from_alignment(alignment)

# analysis contains:
# {
#     "majority_consensus": "ATCG-T",
#     "strict_consensus": "...",
#     "final_consensus": "...",      # Majority if conservation > 0.8, else strict
#     "statistics": { ... },
#     "alignment_length": 6,
#     "sequence_count": 4,
# }
```

### `consensus_statistics`

Calculate conservation metrics across the alignment.

```python
from metainformant.dna.sequence.consensus import consensus_statistics

seqs = ["ATCGATCG", "ATCGATCG", "AGCAATCG"]
stats = consensus_statistics(seqs)
# {
#     "total_positions": 8,
#     "conserved_positions": 6,    # All sequences agree
#     "conservation": 0.75,        # 6/8
#     "ambiguous_positions": 2,    # Disagreement at positions 2, 4
#     "gap_positions": 0,
#     "sequences_count": 3,
# }
```

---

## Advanced Analysis

### Finding Low-Conservation Regions

Scan the alignment with a sliding window to identify regions where consensus
breaks down.

```python
from metainformant.dna.sequence.consensus import find_consensus_breaks

seqs = ["ATCGATCGATCG", "ATCGATCGATCG", "GCTAGCTAGCTA"]
breaks = find_consensus_breaks(seqs, window_size=4)
# Returns [(position, conservation_score), ...]
# Low scores indicate regions of disagreement
```

### Bootstrap Confidence

Estimate consensus confidence through bootstrap resampling of the input sequences.

```python
from metainformant.dna.sequence.consensus import bootstrap_consensus

seqs = ["ATCG", "ATCG", "ATCG", "AGCG"]
consensus, confidence = bootstrap_consensus(seqs, n_bootstraps=100)
# confidence is the fraction of bootstrap replicates producing the same consensus
```

---

## Requirements

- All input sequences must be pre-aligned (equal length).
- Gap characters (`-`) are ignored when counting bases for consensus voting.
- Empty sequence lists raise `ValueError`.

---

## See Also

- **[Accessions](accessions.md)** -- NCBI accession validation
- `metainformant.dna.sequence.consensus` -- Source module
