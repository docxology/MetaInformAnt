### Math: Coalescent Helpers

Functions: `expected_time_to_mrca`, `expected_total_branch_length`, `expected_pairwise_diversity`, `watterson_theta`, `expected_segregating_sites`, `tajima_constants`, `tajimas_D`, `expected_sfs_counts`, `expected_coalescent_waiting_times`.

Example

```python
from metainformant.math import coalescent as coal

Tmrca = coal.expected_time_to_mrca(10, 1e5)
L = coal.expected_total_branch_length(10, 1e5)
pi = coal.expected_pairwise_diversity(1e5, 1e-8)

thetaW = coal.watterson_theta(42, 10)
Es = coal.expected_segregating_sites(theta=0.01, sample_size=10)

const = coal.tajima_constants(10)
D = coal.tajimas_D(num_segregating_sites=42, pairwise_diversity=pi, sample_size=10)

sfs = coal.expected_sfs_counts(sample_size=10, theta=0.01)
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
- tajimas_D(S, π, n): standard normalized D returning 0 when undefined

Examples

```python
from metainformant.math.coalescent import (
    expected_time_to_mrca, expected_total_branch_length,
    expected_sfs_counts, watterson_theta, expected_segregating_sites,
)

tmrca = expected_time_to_mrca(10, 1_000.0)
L = expected_total_branch_length(10, 1_000.0)
sfs = expected_sfs_counts(10, theta=0.01)
theta_hat = watterson_theta(num_segregating_sites=123, sample_size=20, sequence_length=10_000)
E_S = expected_segregating_sites(20, theta=0.005, sequence_length=10_000)
```

Notes

- These helpers assume standard neutral coalescent and per-site mutation rates when θ is used.
- For outputs generated in examples or notebooks, write files under `output/` by default.


