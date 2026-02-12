# Entropy Measures

Shannon entropy, Renyi entropy, Tsallis entropy, joint entropy, conditional entropy, cross-entropy, and entropy-from-counts for discrete probability distributions and observed sequences.

## Key Concepts

**Shannon entropy** quantifies the average surprise or uncertainty in a discrete probability distribution: H(X) = -sum(p_i * log(p_i)). Measured in bits (base 2) or nats (base e). Maximum entropy occurs for the uniform distribution.

**Joint entropy** H(X,Y) measures the total uncertainty of two variables observed together. It is bounded by max(H(X), H(Y)) <= H(X,Y) <= H(X) + H(Y), with the upper bound achieved when X and Y are independent.

**Conditional entropy** H(X|Y) = H(X,Y) - H(Y) measures the remaining uncertainty in X after observing Y. The chain rule decomposes joint entropy as H(X,Y) = H(Y) + H(X|Y).

**Cross-entropy** H(p, q) = -sum(p_i * log(q_i)) measures the average number of bits needed to encode data from distribution p using a code optimised for distribution q. Equals H(p) + D_KL(p||q).

**Renyi entropy** generalises Shannon entropy with a parameter alpha: H_alpha = (1/(1-alpha)) * log(sum(p_i^alpha)). Converges to Shannon entropy as alpha approaches 1. Alpha=0 gives Hartley entropy (log of support size), alpha=inf gives min-entropy.

**Tsallis entropy** is a non-extensive generalisation: S_q = (1 - sum(p_i^q)) / (q - 1). Unlike Renyi entropy, Tsallis entropy is non-additive for independent systems, making it suitable for modelling long-range correlations in biological sequences.

## Function Reference

### `shannon_entropy(probs, base=2.0) -> float`

Calculate Shannon entropy from a probability distribution.

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `probs` | `Sequence[float]` | Probability distribution (must sum to 1.0) |
| `base` | `float` | Logarithm base (2.0 for bits, e for nats) |

**Returns** entropy value as a float.

### `shannon_entropy_from_counts(counts) -> float`

Calculate Shannon entropy in bits from frequency counts. Accepts either a list of integers or a dictionary mapping items to counts.

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `counts` | `Sequence[int]` or `Dict[Any, int]` | Frequency counts |

**Returns** entropy in bits (base 2).

### `joint_entropy(x, y, base=2.0) -> float`

Calculate joint entropy H(X,Y) from paired observations. Counts joint occurrences of (x_i, y_i) pairs and computes entropy of the joint distribution.

### `conditional_entropy(x, y, base=2.0) -> float`

Calculate conditional entropy H(X|Y) = H(X,Y) - H(Y). Measures remaining uncertainty in X after observing Y.

### `cross_entropy(p, q, base=2.0) -> float`

Calculate cross-entropy H(p, q) = -sum(p_i * log(q_i)). Returns infinity if q assigns zero probability where p is nonzero.

### `renyi_entropy(probs, alpha=2.0, base=2.0) -> float`

Calculate Renyi entropy of order alpha. Supports special cases: alpha=0 (Hartley entropy), alpha=inf (min-entropy). Raises ValueError if alpha=1 (use Shannon entropy instead, since Renyi converges to Shannon at alpha=1).

### `tsallis_entropy(probs, q=2.0, base=2.0) -> float`

Calculate Tsallis entropy of order q. Non-extensive generalisation: S_q = (1 - sum(p_i^q)) / (q - 1). Raises ValueError if q=1.

## Usage Examples

```python
from metainformant.information.metrics.core.syntactic import (
    shannon_entropy, shannon_entropy_from_counts,
    joint_entropy, conditional_entropy, cross_entropy,
    renyi_entropy, tsallis_entropy,
)

# Shannon entropy from a probability distribution
probs = [0.25, 0.25, 0.25, 0.25]
h = shannon_entropy(probs, base=2.0)
print(f"Shannon entropy (uniform over 4): {h:.3f} bits")  # 2.000 bits

# Shannon entropy from observed counts
counts = {"A": 100, "C": 80, "G": 90, "T": 70}
h_counts = shannon_entropy_from_counts(counts)
print(f"Entropy from nucleotide counts: {h_counts:.3f} bits")

# Joint and conditional entropy
codons = [1, 1, 2, 2, 3, 3, 1, 2]
amino_acids = [1, 1, 2, 2, 3, 3, 1, 2]
h_joint = joint_entropy(codons, amino_acids)
h_cond = conditional_entropy(codons, amino_acids)
print(f"H(codon, aa) = {h_joint:.3f}, H(codon|aa) = {h_cond:.3f}")

# Cross-entropy between true and predicted distributions
p_true = [0.7, 0.2, 0.1]
q_pred = [0.5, 0.3, 0.2]
ce = cross_entropy(p_true, q_pred)
print(f"Cross-entropy: {ce:.3f} bits")

# Renyi entropy (order 2 = collision entropy)
h_renyi = renyi_entropy([0.5, 0.3, 0.2], alpha=2.0)
print(f"Renyi entropy (alpha=2): {h_renyi:.3f} bits")

# Tsallis entropy
h_tsallis = tsallis_entropy([0.5, 0.3, 0.2], q=2.0)
print(f"Tsallis entropy (q=2): {h_tsallis:.3f}")
```

## Configuration

Environment variable prefix: `INFO_`

## Related Modules

- `metainformant.information.metrics.core.estimation` -- bias-corrected entropy estimation
- `metainformant.information.metrics.core.syntactic` -- mutual information, KL divergence
- `metainformant.information.metrics.advanced.decomposition` -- multivariate decomposition
- `metainformant.information.metrics.advanced.geometry` -- information geometry
