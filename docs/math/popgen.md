### Math: Population Genetics

Functions: `hardy_weinberg_genotype_freqs`, `selection_update`, `mutation_update`, `fixation_probability`, `watterson_theta`, `heterozygosity_decay`, `inbreeding_coefficient`, `equilibrium_heterozygosity_infinite_alleles`, `island_model_update`, `mutation_selection_balance_recessive`, `mutation_selection_balance_dominant`.

Example

```python
from metainformant.math import popgen as pg

p2, two_pq, q2 = pg.hardy_weinberg_genotype_freqs(0.3)
pn = pg.selection_update(0.3, fitness_AA=1.1, fitness_Aa=1.0, fitness_aa=0.9)
pmut = pg.mutation_update(0.3, mu=1e-6, nu=1e-6)
u = pg.fixation_probability(0.01, effective_population_size=1000, selection_coefficient=0.001)

H_t = pg.heterozygosity_decay(0.5, effective_population_size=1000, generations=10)
F_t = pg.inbreeding_coefficient(1000, generations=10)
Heq = pg.equilibrium_heterozygosity_infinite_alleles(Ne=1000, mutation_rate=1e-5)

pnext = pg.island_model_update(0.2, migration_rate=0.1, migrant_pool_frequency=0.8)
q_rec = pg.mutation_selection_balance_recessive(mutation_rate=1e-6, selection_coefficient=1e-2)
q_dom = pg.mutation_selection_balance_dominant(mutation_rate=1e-6, selection_coefficient=1e-2)
```
