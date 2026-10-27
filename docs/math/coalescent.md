### Math: Coalescent Helpers

Functions: `expected_time_to_mrca`, `expected_total_branch_length`, `expected_pairwise_diversity`, `watterson_theta`, `expected_segregating_sites`, `tajima_constants`, `tajimas_D`, `expected_sfs_counts`, `expected_coalescent_waiting_times`, `simulate_coalescent`.

Example

```python
from metainformant.math.population_genetics import coalescent as coal

Tmrca = coal.expected_time_to_mrca(10, 1e5)
L = coal.expected_total_branch_length(10, 1e5)
pi = coal.expected_pairwise_diversity(1e5, 1e-8)

thetaW = coal.watterson_theta(42, 10)
Es = coal.expected_segregating_sites(n_samples=10, theta=0.01)  # theta * a_n

const = coal.tajima_constants(10)
D = coal.tajimas_D(k, 42, 10)  # k: mean pairwise differences per sequence, 42 segregating sites, n = 10

sfs = coal.expected_sfs_counts(n_samples=10, theta=0.01)  # theta / i for i = 1..n-1
Tk = coal.expected_coalescent_waiting_times(10, 1e5)
```

# Math: Coalescent theory helpers

Neutral coalescent expectations with diploid scaling conventions.

Functions

- expected_time_to_mrca(n, Ne): E[T_MRCA] = 4Ne ∑_{k=2}^n 1/(k(k−1))
- expected_total_branch_length(n, Ne): E[L] = 4Ne H_{n−1}
- expected_coalescent_waiting_times(n, Ne): list of E[T_k] = 4Ne/(k(k−1)) for k=n..2
- expected_pairwise_diversity(Ne, μ): π ≈ 4Neμ (per site)
- expected_pairwise_diversity_from_theta(θ): E[π] = θ (per site)
- watterson_theta(S, n, sequence_length=None): θ_W = S/(a1) or per-site with length
- expected_segregating_sites(n, θ, sequence_length=None): E[S] = a1 θ (or ×L)
- expected_sfs_counts(n, θ, sequence_length=None): E[X_i] = θ/i (or θL/i)
- tajima_constants(n): a1,a2,b1,b2,c1,c2,e1,e2
- tajimas_D(π, S, n): standard normalized D returning 0 when undefined; π is the mean number of pairwise differences per sequence, not per site
- simulate_coalescent(n_samples, effective_size=10000, mutation_rate=1e-8, *, random_state=None): one simulated genealogy (see Simulation below)

Examples

```python
from metainformant.math.population_genetics.coalescent import (
    expected_time_to_mrca, expected_total_branch_length,
    expected_sfs_counts, watterson_theta, expected_segregating_sites,
)

tmrca = expected_time_to_mrca(10, 1_000.0)
L = expected_total_branch_length(10, 1_000.0)
sfs = expected_sfs_counts(10, theta=0.01)
theta_hat = watterson_theta(num_segregating_sites=123, sample_size=20, sequence_length=10_000)
E_S = expected_segregating_sites(20, theta=0.005, sequence_length=10_000)
```

Simulation

`simulate_coalescent` draws one genealogy. While k lineages remain, the wait to the next coalescence is exponential with rate k(k−1)/(4Ne) generations and two randomly chosen lineages merge. The number of segregating sites is Poisson with mean `mutation_rate` × total branch length. Averaged over many runs, `total_time` matches `expected_time_to_mrca`.

```python
from metainformant.math.population_genetics.coalescent import simulate_coalescent

genealogy = simulate_coalescent(10, effective_size=1_000.0, mutation_rate=1e-3, random_state=42)
genealogy["coalescence_times"]    # n_samples - 1 cumulative times; the last one is the time to the MRCA
genealogy["total_time"]           # time to the most recent common ancestor
genealogy["total_branch_length"]  # sum over intervals of lineages × interval length
genealogy["n_segregating_sites"]  # Poisson(mutation_rate × total_branch_length)
genealogy["theta_watterson"]      # watterson_theta(n_segregating_sites, n_samples)
```

The dictionary also repeats `n_samples`, `effective_size` and `mutation_rate`. `ValueError` is raised for fewer than 2 samples, a non-positive `effective_size` or a negative `mutation_rate`.

`random_state` is keyword-only and takes an `int` seed, a `numpy.random.Generator` or a `numpy.random.RandomState`. An int or a generator gives a reproducible genealogy that does not depend on numpy's global random state, and a seeded run leaves that global state unchanged. With the default `None` the draws come from numpy's global generator, which `np.random.seed(...)` controls.

Notes

- These helpers assume standard neutral coalescent and per-site mutation rates when θ is used.
- For outputs generated in examples or notebooks, write files under `output/` by default.
