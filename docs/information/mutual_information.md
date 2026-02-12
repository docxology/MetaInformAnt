# Mutual Information and Divergence Measures

Mutual information, conditional mutual information, KL divergence, Jensen-Shannon divergence, total correlation, transfer entropy, normalized mutual information, and information coefficient.

## Key Concepts

**Mutual information** I(X;Y) = H(X) + H(Y) - H(X,Y) quantifies the amount of information one variable provides about another. It is symmetric, non-negative, and equals zero if and only if the variables are independent.

**Conditional mutual information** I(X;Y|Z) measures the information X provides about Y when Z is already known. Computed as I(X;Y|Z) = H(X,Z) + H(Y,Z) - H(Z) - H(X,Y,Z). Central to Markov blanket discovery and causal inference in gene regulatory networks.

**KL divergence** D_KL(p||q) = sum(p_i * log(p_i/q_i)) measures the information lost when approximating distribution p with q. Asymmetric, non-negative, zero iff p = q. Returns infinity when q assigns zero probability to an event with nonzero probability under p.

**Jensen-Shannon divergence** JSD(p||q) = 0.5 * D_KL(p||m) + 0.5 * D_KL(q||m) where m = (p+q)/2. Symmetric, bounded in [0, 1] (base 2), always finite. The square root of JSD is a true metric.

**Total correlation** TC(X1,...,Xn) = sum(H(Xi)) - H(X1,...,Xn) measures the total multivariate dependence among n variables. Equals mutual information for n=2.

**Transfer entropy** T(X->Y) = I(Y_{t+1}; X_t | Y_t) measures directed information flow from time series X to Y, capturing the additional predictive power X provides about Y's future beyond Y's own past.

**Normalized mutual information** scales MI to [0,1] using one of several normalisation schemes: arithmetic mean, geometric mean, max, or min of the marginal entropies.

## Function Reference

### `mutual_information(x, y, base=2.0) -> float`

Calculate mutual information I(X;Y) from paired observations.

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `x` | `Sequence[Any]` | First variable observations |
| `y` | `Sequence[Any]` | Second variable observations (same length) |
| `base` | `float` | Logarithm base |

**Returns** mutual information I(X;Y).

### `conditional_mutual_information(x, y, z, base=2.0) -> float`

Calculate conditional mutual information I(X;Y|Z). All three sequences must have the same length.

### `kl_divergence(p, q, base=2.0) -> float`

Calculate Kullback-Leibler divergence D_KL(p||q). Both distributions must sum to 1.0 and have equal length. Returns `float("inf")` when q has zero probability where p does not.

### `jensen_shannon_divergence(p, q, base=2.0) -> float`

Calculate Jensen-Shannon divergence JSD(p||q). Always finite, symmetric, and bounded.

### `total_correlation(variables, base=2.0) -> float`

Calculate total correlation (multi-information) for a list of variable sequences. All sequences must have the same length.

### `transfer_entropy(x, y, lag=1, base=2.0) -> float`

Calculate transfer entropy from X to Y with a given time lag. Measures directed information flow: T(X->Y) = I(Y_{t+1}; X_t | Y_t). Requires sequences longer than lag+1.

### `normalized_mutual_information(x, y, method="arithmetic", base=2.0) -> float`

Calculate normalized mutual information (NMI) in [0,1].

**Methods:**
| Method | Formula |
|--------|---------|
| `arithmetic` | 2*I(X;Y) / (H(X) + H(Y)) |
| `geometric` | I(X;Y) / sqrt(H(X) * H(Y)) |
| `max` | I(X;Y) / max(H(X), H(Y)) |
| `min` | I(X;Y) / min(H(X), H(Y)) |

### `information_coefficient(x, y, base=2.0) -> float`

Calculate the information coefficient IC = I(X;Y) / (H(X) + H(Y) - I(X;Y)). Values range from 0 (independent) to 1 (deterministic relationship). The denominator equals the joint entropy H(X,Y).

## Usage Examples

```python
from metainformant.information.metrics.core.syntactic import (
    mutual_information, conditional_mutual_information,
    kl_divergence, jensen_shannon_divergence,
    total_correlation, transfer_entropy,
    normalized_mutual_information, information_coefficient,
)

# Mutual information between gene expression and phenotype
expression = [0, 0, 1, 1, 2, 2, 0, 1, 2, 0]
phenotype  = [0, 0, 1, 1, 1, 1, 0, 1, 1, 0]
mi = mutual_information(expression, phenotype)
print(f"I(expression; phenotype) = {mi:.3f} bits")

# Conditional MI: does gene X inform phenotype beyond gene Y?
gene_x = [1, 0, 1, 0, 1, 0, 1, 0]
gene_y = [1, 1, 0, 0, 1, 1, 0, 0]
target = [1, 0, 0, 0, 1, 1, 0, 0]
cmi = conditional_mutual_information(target, gene_x, gene_y)
print(f"I(target; gene_x | gene_y) = {cmi:.3f} bits")

# KL divergence between observed and expected nucleotide frequencies
observed = [0.31, 0.19, 0.21, 0.29]
expected = [0.25, 0.25, 0.25, 0.25]
kl = kl_divergence(observed, expected)
print(f"D_KL(observed || uniform) = {kl:.4f} bits")

# Jensen-Shannon divergence (symmetric)
jsd = jensen_shannon_divergence(observed, expected)
print(f"JSD = {jsd:.4f} bits")

# Total correlation among three regulatory genes
reg_a = [0, 1, 0, 1, 1, 0, 1, 0]
reg_b = [0, 1, 1, 1, 1, 0, 0, 0]
reg_c = [1, 0, 0, 1, 1, 1, 0, 0]
tc = total_correlation([reg_a, reg_b, reg_c])
print(f"Total correlation: {tc:.3f} bits")

# Transfer entropy: information flow in a signalling cascade
signal   = [0, 1, 1, 0, 0, 1, 1, 0, 1, 0]
response = [0, 0, 1, 1, 0, 0, 1, 1, 0, 1]
te = transfer_entropy(signal, response, lag=1)
print(f"T(signal -> response) = {te:.3f} bits")

# Normalized MI for clustering comparison
nmi = normalized_mutual_information(expression, phenotype, method="arithmetic")
print(f"NMI = {nmi:.3f}")

# Information coefficient
ic = information_coefficient(expression, phenotype)
print(f"IC = {ic:.3f}")
```

## Configuration

Environment variable prefix: `INFO_`

## Related Modules

- `metainformant.information.metrics.core.syntactic` -- entropy measures
- `metainformant.information.metrics.core.estimation` -- bias-corrected MI estimation
- `metainformant.information.metrics.advanced.decomposition` -- PID, co-information
- `metainformant.information.metrics.advanced.channel` -- channel capacity
