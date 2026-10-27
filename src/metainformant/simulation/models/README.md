# Simulation Models

Synthetic data generators for agent-based ecology, population genetics, RNA-seq counts, and sequence evolution simulations.

## Contents

| File | Purpose |
|------|---------|
| `agents.py` | Agent-based ecosystem simulation with population dynamics |
| `popgen.py` | Population genetics: neutral coalescent engine and sequence generators (bottleneck, expansion), genotype matrices, two-population Fst, linkage disequilibrium, SFS, admixture, selection |
| `rna.py` | RNA-seq count simulation: bulk, single-cell, differential expression |
| `sequences.py` | DNA/protein sequence generation, mutation, evolution, divergence |

## Key Functions

| Function | Description |
|----------|-------------|
| `create_ecosystem()` | Initialize agent-based ecosystem with species parameters |
| `run_simulation()` | Run multi-step agent simulation with logging |
| `generate_population_sequences()` | Aligned DNA sequences as one neutral coalescent sample: one genealogy, infinite sites, no recombination; per-site θ from `nucleotide_diversity`, `wattersons_theta`, `theta` or `mutation_rate` (first given wins, default 0.01) |
| `generate_genotype_matrix()` | SNP genotype matrix with allele frequency control |
| `generate_two_populations()` | Two populations at a target Fst (Balding-Nichols allele frequencies) |
| `simulate_bottleneck_population()` | Coalescent sample drawn at the end of a bottleneck with `n_sequences`; a heuristic size trajectory without it |
| `simulate_population_expansion()` | Coalescent sample after exponential growth with `n_sequences`; a heuristic size trajectory without it |
| `generate_linkage_disequilibrium_data()` | 0/1 alleles with a chosen adjacent-site r² (first-order Markov chain) |
| `PopulationEpoch` | One epoch of a piecewise-constant or exponential demographic history |
| `simulate_coalescent_spectra()` | Site frequency spectra of many replicates: Poisson at θ, or exactly `segregating_sites` mutations (`ms -s`) |
| `simulate_coalescent_haplotypes()` | One coalescent sample as a 0/1 haplotype matrix |
| `haplotypes_to_sequences()` | Write a haplotype matrix as DNA with seeded positions and bases |
| `simulate_rnaseq_counts()` | Simulated RNA-seq count matrix with biological variation |
| `simulate_differential_expression()` | Two-group counts with known DE genes for benchmarking; each fold change is a signed ratio on the group-1 mean (positive multiplies, negative divides by its magnitude, so `-4.0` gives about 0.25x) |
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
