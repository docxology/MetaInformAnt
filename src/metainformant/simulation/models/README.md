# Simulation Models

Synthetic data generators for agent-based ecology, population genetics, RNA-seq counts, and sequence evolution simulations.

## Contents

| File | Purpose |
|------|---------|
| `agents.py` | Agent-based ecosystem simulation with population dynamics |
| `popgen.py` | Population genetics: genotype matrices, bottleneck, admixture, SFS |
| `rna.py` | RNA-seq count simulation: bulk, single-cell, differential expression |
| `sequences.py` | DNA/protein sequence generation, mutation, evolution, divergence |

## Key Functions

| Function | Description |
|----------|-------------|
| `create_ecosystem()` | Initialize agent-based ecosystem with species parameters |
| `run_simulation()` | Run multi-step agent simulation with logging |
| `generate_population_sequences()` | Generate polymorphic DNA sequences for a population |
| `generate_genotype_matrix()` | SNP genotype matrix with allele frequency control |
| `simulate_bottleneck_population()` | Population passing through bottleneck event |
| `simulate_rnaseq_counts()` | Simulated RNA-seq count matrix with biological variation |
| `simulate_differential_expression()` | Counts with known DE genes for benchmarking |
| `generate_random_dna()` | Random DNA sequence with configurable GC content |
| `evolve_sequence()` | Simulate sequence divergence along a phylogeny |
| `generate_sequence_family()` | Generate related sequence family with known divergence |

## Usage

```python
from metainformant.simulation.models.rna import simulate_rnaseq_counts
from metainformant.simulation.models.sequences import generate_random_dna

counts = simulate_rnaseq_counts(n_genes=1000, n_samples=6)
seq = generate_random_dna(length=500, gc_content=0.45)
```
