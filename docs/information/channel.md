# Channel Capacity and Rate-Distortion

Blahut-Arimoto channel capacity computation, rate-distortion function, information bottleneck method, channel mutual information, and analytical capacity formulas for standard noisy channel models.

## Key Concepts

**Channel capacity** C = max_{p(x)} I(X;Y) is the maximum rate at which information can be reliably transmitted through a noisy channel defined by transition matrix P(Y|X). The Blahut-Arimoto algorithm iteratively refines the input distribution to achieve this maximum.

**Rate-distortion theory** characterises the minimum rate R(D) required to represent a source within average distortion D. The R(D) curve is computed by sweeping a Lagrange multiplier (beta) that trades off rate against distortion.

**Information bottleneck** (Tishby et al. 2000) finds a compressed representation T of input X that preserves maximal information about a relevant variable Y. Optimises L = I(T;X) - beta * I(T;Y), balancing compression against relevance. Applied in biological sequence clustering and neural decoding.

**Noisy channel models** have closed-form capacity expressions: Binary Symmetric Channel C = 1 - H(p), Binary Erasure Channel C = 1 - epsilon, AWGN Channel C = 0.5 * log2(1 + SNR).

## Function Reference

### `channel_capacity(transition_matrix, max_iter=100, tol=1e-8) -> Dict`

Compute channel capacity using the Blahut-Arimoto algorithm.

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `transition_matrix` | `np.ndarray` | P(Y\|X) matrix, shape (n_x, n_y), rows sum to 1 |
| `max_iter` | `int` | Maximum iterations |
| `tol` | `float` | Convergence tolerance |

**Returns** dict with `capacity` (bits), `optimal_input` (list of floats), `n_iterations`.

### `rate_distortion(source_dist, distortion_matrix, max_rate=None, n_points=50) -> Dict`

Compute the rate-distortion function R(D) for a discrete source using the Blahut-Arimoto algorithm for rate-distortion.

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `source_dist` | `Sequence[float]` | Source probability distribution P(X) |
| `distortion_matrix` | `np.ndarray` | d(x, x_hat) matrix, shape (n_x, n_x_hat) |
| `max_rate` | `float` or `None` | Upper bound on rate; defaults to source entropy |
| `n_points` | `int` | Number of points on the R(D) curve |

**Returns** dict with `rates` (list), `distortions` (list), `rd_curve` (list of (rate, distortion) tuples).

### `information_bottleneck(p_xy, beta=1.0, n_clusters=2, max_iter=100, tol=1e-6, seed=42) -> Dict`

Compute the information bottleneck compression of X with respect to Y.

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `p_xy` | `np.ndarray` | Joint distribution P(X,Y), shape (n_x, n_y), sums to 1 |
| `beta` | `float` | Tradeoff: higher preserves more info about Y |
| `n_clusters` | `int` | Number of compressed states T |
| `max_iter` | `int` | Maximum iterations |
| `tol` | `float` | Convergence tolerance |
| `seed` | `int` | Random seed for initialisation |

**Returns** dict with `assignments` (cluster per X value), `I_T_X` (compression cost), `I_T_Y` (relevance preserved), `beta`.

### `channel_mutual_information(transition_matrix, input_dist) -> float`

Compute I(X;Y) = H(Y) - H(Y|X) for a channel with a given input distribution. Useful for evaluating sub-optimal input distributions.

### `noisy_channel_capacity(noise_level, channel_type="binary_symmetric") -> float`

Analytical channel capacity for standard models.

**Channel Types:**
| Type | Formula | Noise Parameter |
|------|---------|----------------|
| `binary_symmetric` | C = 1 - H(p) | Crossover probability p in [0,1] |
| `binary_erasure` | C = 1 - epsilon | Erasure probability epsilon in [0,1] |
| `awgn` | C = 0.5*log2(1+SNR) | Signal-to-noise ratio (linear) |

## Usage Examples

```python
import numpy as np
from metainformant.information.metrics.advanced.channel import (
    channel_capacity, rate_distortion, information_bottleneck,
    channel_mutual_information, noisy_channel_capacity,
)

# Channel capacity of a binary symmetric channel (transition matrix)
p = 0.1  # crossover probability
bsc = np.array([[1 - p, p], [p, 1 - p]])
result = channel_capacity(bsc)
print(f"BSC capacity: {result['capacity']:.4f} bits")
print(f"Optimal input: {result['optimal_input']}")
print(f"Converged in {result['n_iterations']} iterations")

# Analytical capacity (closed-form)
c_bsc = noisy_channel_capacity(0.1, "binary_symmetric")
c_bec = noisy_channel_capacity(0.3, "binary_erasure")
c_awgn = noisy_channel_capacity(10.0, "awgn")  # SNR=10
print(f"BSC: {c_bsc:.4f}, BEC: {c_bec:.4f}, AWGN: {c_awgn:.4f} bits")

# Rate-distortion for a binary source with Hamming distortion
source = [0.6, 0.4]
hamming = np.array([[0, 1], [1, 0]], dtype=float)
rd = rate_distortion(source, hamming, n_points=30)
for r, d in rd["rd_curve"][:5]:
    print(f"R={r:.3f} bits, D={d:.3f}")

# Information bottleneck: compress gene expression patterns
# Joint distribution over 4 expression levels x 3 phenotypes
p_xy = np.random.dirichlet(np.ones(3), size=4)
p_xy /= p_xy.sum()  # normalise to valid joint distribution
ib = information_bottleneck(p_xy, beta=2.0, n_clusters=2)
print(f"Compression: I(T;X) = {ib['I_T_X']:.4f}")
print(f"Relevance:   I(T;Y) = {ib['I_T_Y']:.4f}")
print(f"Assignments: {ib['assignments']}")

# Mutual information for a specific input distribution
uniform_input = [0.5, 0.5]
mi = channel_mutual_information(bsc, uniform_input)
print(f"I(X;Y) with uniform input: {mi:.4f} bits")
```

## Configuration

Environment variable prefix: `INFO_`

## Related Modules

- `metainformant.information.metrics.core.syntactic` -- Shannon entropy, MI
- `metainformant.information.metrics.advanced.geometry` -- Fisher-Rao, information geometry
- `metainformant.information.metrics.core.estimation` -- bias-corrected estimation
- `metainformant.information.metrics.advanced.decomposition` -- PID
