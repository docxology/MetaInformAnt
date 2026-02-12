### Simulation: Sequences

The `metainformant.simulation.models.sequences` module provides functions for generating
random biological sequences (DNA and protein), introducing mutations, simulating
evolutionary divergence over generations, and analyzing sequence family relationships.

All functions accept an optional `rng` parameter (`random.Random` instance) for
reproducible results.

---

## Constants

| Constant       | Value                        | Description                         |
|----------------|------------------------------|-------------------------------------|
| `DNA_BASES`    | `"ATCG"`                    | Standard DNA nucleotides            |
| `RNA_BASES`    | `"AUCG"`                    | Standard RNA nucleotides            |
| `AMINO_ACIDS`  | `"ACDEFGHIKLMNPQRSTVWY"`    | Standard 20 amino acids             |
| `GENETIC_CODE` | `Dict[str, str]` (64 codons)| Standard codon-to-amino-acid mapping|

---

## Sequence Generation

### `generate_random_dna`

Generate a random DNA sequence with controllable GC content.

```python
from metainformant.simulation import generate_random_dna

seq = generate_random_dna(100)                          # Default 50% GC
seq = generate_random_dna(100, gc_content=0.65)         # AT-rich suppressed
seq = generate_random_dna(100, gc_content=0.3, rng=rng) # Reproducible
```

**Parameters:**

| Parameter    | Type    | Default | Description                                     |
|--------------|---------|---------|-------------------------------------------------|
| `length`     | `int`   | required| Sequence length (must be >= 1)                  |
| `gc_content` | `float` | `0.5`   | Target GC fraction (0.0 to 1.0)                |
| `rng`        | `Random`| `None`  | Random number generator for reproducibility     |

### `generate_random_protein`

Generate a random protein sequence from the 20 standard amino acids with uniform
frequency.

```python
from metainformant.simulation import generate_random_protein

protein = generate_random_protein(50)
```

### `generate_coding_sequence`

Generate a random coding DNA sequence and its protein translation simultaneously.
The length must be divisible by 3.

```python
from metainformant.simulation import generate_coding_sequence

dna_seq, protein_seq = generate_coding_sequence(300, gc_content=0.5)
# dna_seq is 300 bases, protein_seq is the translated amino acid string
```

---

## Mutation and Evolution

### `mutate_sequence`

Introduce a specified number of point mutations at random positions. Works for both
DNA and protein sequences (auto-detected from character set).

```python
from metainformant.simulation import mutate_sequence

original = "ATCGATCGATCG"
mutated = mutate_sequence(original, n_mut=3)
# Exactly 3 positions will differ from the original
```

**Parameters:**

| Parameter | Type    | Default | Description                               |
|-----------|---------|---------|-------------------------------------------|
| `seq`     | `str`   | required| Input sequence (DNA or protein)           |
| `n_mut`   | `int`   | required| Number of mutations (0 to len(seq))       |
| `rng`     | `Random`| `None`  | Random number generator                   |

### `evolve_sequence`

Simulate sequence evolution over multiple generations. The number of mutations per
generation is drawn from a Poisson distribution with mean = `len(seq) * mutation_rate`.

```python
from metainformant.simulation import evolve_sequence

ancestor = generate_random_dna(500)
descendant = evolve_sequence(ancestor, generations=1000, mutation_rate=0.001)
```

**Parameters:**

| Parameter       | Type    | Default | Description                                  |
|-----------------|---------|---------|----------------------------------------------|
| `sequence`      | `str`   | required| Starting sequence                            |
| `generations`   | `int`   | required| Number of generations (>= 0)                |
| `mutation_rate`  | `float` | `0.001` | Per-base per-generation mutation probability |
| `rng`           | `Random`| `None`  | Random number generator                      |

---

## Translation

### `translate_dna_to_protein`

Translate a DNA sequence to a protein string using the standard genetic code. Stops
at the first stop codon encountered.

```python
from metainformant.simulation import translate_dna_to_protein

protein = translate_dna_to_protein("ATGAAAGCGTGA")
# "MKA*" -- M(Met), K(Lys), A(Ala), *(Stop)
```

The `frame` parameter selects the reading frame (0, 1, or 2).

### `reverse_transcribe_protein_to_dna`

Reverse translate a protein to one possible DNA sequence (random codon usage).

```python
from metainformant.simulation import reverse_transcribe_protein_to_dna

dna = reverse_transcribe_protein_to_dna("MKA")
# One of many possible codon combinations encoding Met-Lys-Ala
```

---

## Sequence Families and Divergence

### `generate_sequence_family`

Generate a family of related sequences from a common ancestor, simulating independent
evolutionary lineages.

```python
from metainformant.simulation import generate_sequence_family

ancestor = generate_random_dna(200)
family = generate_sequence_family(
    ancestor, n_descendants=5, generations=500, mutation_rate=0.002,
)
# Returns [ancestor, descendant_1, ..., descendant_5]
```

### `analyze_sequence_divergence`

Compute pairwise divergence statistics for a set of equal-length sequences.

```python
from metainformant.simulation import analyze_sequence_divergence

stats = analyze_sequence_divergence(family)
# {
#     "num_sequences": 6,
#     "sequence_length": 200,
#     "mean_similarity": 0.87,
#     "mean_divergence": 0.13,
#     "variable_positions": 42,
#     "variable_fraction": 0.21,
#     "pairwise_similarities": [...],
#     "pairwise_divergences": [...],
# }
```

### `simulate_gene_duplication`

Simulate gene duplication followed by independent divergence of each copy.

```python
from metainformant.simulation import simulate_gene_duplication

copies = simulate_gene_duplication(
    original_gene=generate_random_dna(600),
    n_copies=3,
    divergence_time=10000,
    mutation_rate=1e-8,
)
# Returns 3 diverged copies of the original gene
```

### `calculate_sequence_similarity`

Calculate the fraction of identical positions between two equal-length sequences.

```python
from metainformant.simulation import calculate_sequence_similarity

sim = calculate_sequence_similarity("ATCGATCG", "ATCAATCG")
# 0.875 (7 of 8 positions identical)
```

---

## See Also

- **[RNA Count Simulation](rna_counts.md)** -- Simulate RNA-seq count matrices
- **[Simulation Overview](index.md)** -- Full module architecture
- `metainformant.simulation.models.sequences` -- Source module
