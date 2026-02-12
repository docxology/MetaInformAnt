# Information Geometry

Fisher-Rao geodesic distance, natural gradient, information projection, alpha-divergences, exponential family entropy, Hellinger distance, entropy power inequality, and information dimension (Grassberger-Procaccia).

## Key Concepts

**Fisher-Rao distance** is the geodesic distance on the statistical manifold equipped with the Fisher information metric: d_FR(p, q) = 2 * arccos(sum(sqrt(p_i * q_i))). It is a true metric satisfying the triangle inequality. Equals 0 for identical distributions and pi for distributions with disjoint support.

**Natural gradient** rescales the ordinary Euclidean gradient by the inverse Fisher information matrix: nat_grad = F^{-1} @ grad. This makes parameter updates invariant to reparameterisation, yielding faster convergence in optimisation on statistical manifolds.

**Information projection** (m-projection) finds the distribution in a constraint set closest to a reference distribution in KL divergence: q* = argmin_{q in C} D_KL(q || p). Computed via iterative proportional fitting (Sinkhorn scaling).

**Alpha-divergence** is a one-parameter family interpolating between forward KL (alpha->1), reverse KL (alpha->0), and the squared Hellinger distance (alpha=0.5). General formula: D_alpha = (4/(1-alpha^2)) * (1 - sum(p^{(1+alpha)/2} * q^{(1-alpha)/2})).

**Exponential family entropy** exploits the Legendre transform: H[p] = -eta^T E[T(x)] + A(eta), where eta is the natural parameter, E[T(x)] is the expected sufficient statistic, and A(eta) is the log-partition function. Returns entropy in nats.

**Hellinger distance** H(p,q) = (1/sqrt(2)) * sqrt(sum((sqrt(p_i) - sqrt(q_i))^2)) takes values in [0,1] and is related to the Bhattacharyya coefficient: H = sqrt(1 - BC).

**Entropy power inequality** states N(X+Y) >= N(X) + N(Y) for independent random variables, where N(X) = (1/(2*pi*e)) * exp(2*h(X)) is the entropy power. Tight for Gaussians.

**Information dimension** (Grassberger-Procaccia) estimates the correlation dimension d_2 from point-cloud data by fitting log C(r) vs log r, where C(r) is the correlation integral.

## Function Reference

### `fisher_rao_distance(p, q) -> float`

Compute the Fisher-Rao geodesic distance between two discrete distributions.

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `p` | `Sequence[float]` | First distribution (sums to 1) |
| `q` | `Sequence[float]` | Second distribution (sums to 1, same length) |

**Returns** non-negative distance. Returns 0.0 for identical distributions, pi for disjoint support.

### `natural_gradient(loss_gradient, fisher_info_matrix) -> np.ndarray`

Compute the natural gradient F^{-1} @ grad. Solves the linear system directly for numerical stability; falls back to pseudo-inverse for singular Fisher matrices.

### `information_projection(p, constraint_set, method="iterative_scaling", max_iter=100, tol=1e-8) -> Dict`

Project distribution p onto a constraint set via iterative proportional fitting.

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `p` | `Sequence[float]` | Reference distribution |
| `constraint_set` | `List[Tuple[List[int], List[float]]]` | (index_set, target_marginal) pairs |
| `max_iter` | `int` | Maximum iterations |
| `tol` | `float` | Convergence tolerance |

**Returns** dict with `projected` (list), `kl_divergence` (float), `n_iterations` (int).

### `statistical_divergence(p, q, alpha=0.5) -> float`

Compute alpha-divergence between two distributions. Falls back to KL divergence for alpha near 0 or 1.

### `exponential_family_entropy(natural_params, sufficient_stats_expectation, log_partition) -> float`

Compute entropy of an exponential family distribution via the Legendre transform. Returns entropy in nats.

### `hellinger_distance(p, q) -> float`

Compute Hellinger distance in [0, 1] between two discrete distributions.

### `entropy_power_inequality(variances) -> Dict`

Compute entropy powers and verify the EPI for independent Gaussian random variables.

**Returns** dict with `entropy_powers` (list), `sum_entropy_power` (float), `epi_bound` (float), `epi_satisfied` (bool).

### `information_dimension(samples, r_values=None, method="correlation") -> Dict`

Estimate the correlation dimension d_2 from point-cloud data using the Grassberger-Procaccia algorithm.

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `samples` | `np.ndarray` | Point cloud, shape (n_samples, n_dims) or 1-D |
| `r_values` | `np.ndarray` or `None` | Radii for correlation integral; auto-chosen if None |
| `method` | `str` | Only `"correlation"` is currently supported |

**Returns** dict with `dimension` (float), `r_values` (array), `correlation_integral` (array), `r_squared` (float, quality of log-log fit).

## Usage Examples

```python
import numpy as np
from metainformant.information.metrics.advanced.geometry import (
    fisher_rao_distance, natural_gradient, information_projection,
    statistical_divergence, exponential_family_entropy,
    hellinger_distance, entropy_power_inequality, information_dimension,
)

# Fisher-Rao distance between nucleotide frequency profiles
genome_a = [0.30, 0.20, 0.20, 0.30]
genome_b = [0.25, 0.25, 0.25, 0.25]
d_fr = fisher_rao_distance(genome_a, genome_b)
print(f"Fisher-Rao distance: {d_fr:.4f}")

# Hellinger distance
d_h = hellinger_distance(genome_a, genome_b)
print(f"Hellinger distance: {d_h:.4f}")

# Alpha-divergence (alpha=0.5 relates to Hellinger)
d_alpha = statistical_divergence(genome_a, genome_b, alpha=0.5)
print(f"Alpha-divergence (alpha=0.5): {d_alpha:.4f}")

# Natural gradient for optimisation on a statistical manifold
grad = np.array([1.0, -0.5, 0.3])
fim = np.array([[2.0, 0.1, 0.0], [0.1, 1.5, 0.2], [0.0, 0.2, 1.0]])
nat_grad = natural_gradient(grad, fim)
print(f"Natural gradient: {nat_grad}")

# Information projection with marginal constraints
p_ref = [0.25, 0.25, 0.25, 0.25]
constraints = [([0, 1], [0.6]),  # indices 0,1 should have mass 0.6
               ([2, 3], [0.4])]  # indices 2,3 should have mass 0.4
result = information_projection(p_ref, constraints)
print(f"Projected: {result['projected']}")
print(f"KL divergence: {result['kl_divergence']:.6f}")

# Exponential family entropy (Gaussian example)
# For N(mu, sigma^2): eta = [mu/sigma^2, -1/(2*sigma^2)]
# E[T] = [mu, mu^2 + sigma^2], A(eta) = -eta1^2/(4*eta2) - 0.5*log(-2*eta2)
eta = [0.0, -0.5]
e_t = [0.0, 1.0]  # E[x] = 0, E[x^2] = 1 for N(0,1)
log_A = 0.5 * np.log(2 * np.pi)
h = exponential_family_entropy(eta, e_t, log_A)
print(f"Gaussian entropy: {h:.4f} nats")

# Entropy power inequality for three independent Gaussians
epi = entropy_power_inequality([1.0, 2.0, 0.5])
print(f"EPI satisfied: {epi['epi_satisfied']}")

# Correlation dimension of gene expression point cloud
samples = np.random.randn(200, 3)
dim = information_dimension(samples)
print(f"Estimated dimension: {dim['dimension']:.2f} (R^2={dim['r_squared']:.3f})")
```

## Configuration

Environment variable prefix: `INFO_`

Requires `numpy` and `scipy.spatial.distance.pdist` for the information dimension estimator.

## Related Modules

- `metainformant.information.metrics.core.syntactic` -- Shannon entropy, KL divergence
- `metainformant.information.metrics.advanced.channel` -- channel capacity, rate-distortion
- `metainformant.information.metrics.core.estimation` -- bias-corrected estimation
- `metainformant.information.metrics.advanced.decomposition` -- multivariate information
