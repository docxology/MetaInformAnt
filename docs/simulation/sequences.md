### Simulation: Sequences

Functions: `generate_random_dna`, `mutate_sequence`, `generate_random_protein`.

```python
from metainformant.simulation import generate_random_dna, mutate_sequence, generate_random_protein

dna = generate_random_dna(50)
mut = mutate_sequence(dna, n_mut=5)
prot = generate_random_protein(50)
```
