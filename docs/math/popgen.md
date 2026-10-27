### Math: Population Genetics

Functions: `hardy_weinberg_genotype_freqs`, `selection_update`, `mutation_update`, `fixation_probability`, `watterson_theta`, `heterozygosity_decay`, `inbreeding_coefficient`, `equilibrium_heterozygosity_infinite_alleles`, `island_model_update`, `mutation_selection_balance_recessive`, `mutation_selection_balance_dominant`.

They live in the submodules of `metainformant.math.population_genetics`:

| Submodule | Functions |
|-----------|-----------|
| `core` | `hardy_weinberg_genotype_freqs`, `heterozygosity_decay`, `inbreeding_coefficient`, `mutation_selection_balance_recessive`, `mutation_selection_balance_dominant` |
| `selection` | `selection_update`, `mutation_update` (and its own copies of the two mutation-selection balance functions) |
| `statistics` | `fixation_probability`, `equilibrium_heterozygosity_infinite_alleles` |
| `demography` | `island_model_update` |
| `coalescent` | `watterson_theta` |

`metainformant.math.popgen` is only a thin compatibility module: it offers `hardy_weinberg_genotype_freqs`, `hardy_weinberg_allele_freqs` and `fst_from_freqs`.

Example

```python
from metainformant.math.population_genetics import coalescent, core, demography, selection, statistics

p2, two_pq, q2 = core.hardy_weinberg_genotype_freqs(0.3)  # (AA, Aa, aa) frequencies
pn = selection.selection_update(0.3, fitness_AA=1.1, fitness_Aa=1.0, fitness_aa=0.9)
pmut = selection.mutation_update(0.3, mu=1e-6, nu=1e-6)
u = statistics.fixation_probability(0.01, Ne=1000, selection_coefficient=0.001)

H_t = core.heterozygosity_decay(0.5, effective_size=1000, generations=10)
F_t = core.inbreeding_coefficient(1000, generations=10)
Heq = statistics.equilibrium_heterozygosity_infinite_alleles(theta=4 * 1000 * 1e-5)  # theta = 4 Ne mu

pnext = demography.island_model_update(0.2, migration_rate=0.1, migrant_allele_frequency=0.8)
q_rec = core.mutation_selection_balance_recessive(mutation_rate=1e-6, selection_coefficient=1e-2)
q_dom = core.mutation_selection_balance_dominant(mutation_rate=1e-6, selection_coefficient=1e-2)

theta_w = coalescent.watterson_theta(25, 10, sequence_length=1000)  # S = 25 sites, n = 10 sequences, per site
```

Notes:
- `selection_update` takes either explicit genotype fitnesses (`fitness_AA`, `fitness_Aa`, `fitness_aa`) or `selection_coefficient` with `dominance_coefficient` (default 0.5, additive).
- `mutation_update` needs both rates: `mutation_rate_forward` (A to a, alias `mu`) and `mutation_rate_backward` (a to A, alias `nu`).
- `fixation_probability(initial_frequency, population_size, selection_coefficient)` uses Kimura's diffusion approximation; `Ne` and `p0` are keyword aliases for `population_size` and `initial_frequency`, and a neutral allele (`selection_coefficient=0`) fixes with probability equal to its frequency.
- `equilibrium_heterozygosity_infinite_alleles(theta)` returns θ / (1 + θ). When `mutation_rate` is also given, the first argument is read as Ne instead and θ = 4 · Ne · `mutation_rate`.
- `watterson_theta(S, n, sequence_length=1)` returns S / (a_n · `sequence_length`) with a_n = Σ 1/i for i < n, so it is per site only when a sequence length is passed.
