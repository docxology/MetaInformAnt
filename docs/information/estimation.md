# Entropy Estimation

Bias-corrected entropy estimation methods including plugin, Miller-Madow, Chao-Shen, and jackknife estimators, plus mutual information estimation, KL divergence estimation, bootstrap confidence intervals, entropy rate estimation, and advanced bias corrections (Panzeri-Treves, effective sample size).

## Key Concepts

**Plugin estimator** (maximum likelihood) computes entropy directly from sample frequencies: H_hat = -sum(f_i/n * log(f_i/n)). This is negatively biased for finite samples, underestimating true entropy by approximately (k-1)/(2n) where k is the number of categories.

**Miller-Madow correction** subtracts the leading bias term: H_mm = H_plugin + (k-1)/(2n). Simple and effective for moderate sample sizes.

**Chao-Shen estimator** addresses sparse data by using a Good-Turing coverage correction. It upweights contributions from rare species/categories, making it particularly suitable for biodiversity and microbiome data where many species are observed only once.

**Jackknife estimator** reduces bias by computing leave-one-out entropy estimates and combining them: H_jack = k*H_plugin - (k-1)*mean(H_{-i}). Provides first-order bias reduction.

**Panzeri-Treves correction** uses response-specific frequencies for more accurate small-sample bias correction, originally developed for neural coding analysis.

**Bootstrap confidence intervals** resample the observed data with replacement to estimate the sampling distribution of the entropy estimator, yielding confidence intervals without distributional assumptions.

**Entropy rate** estimates H(X_{n+1}|X_1^n) for sequential data using block entropies: h = H(blocks of length n+1) - H(blocks of length n). Characterises the per-symbol uncertainty in a stationary process.

## Function Reference

### `entropy_estimator(counts, method="plugin", bias_correction=True) -> float`

Estimate Shannon entropy with selectable method and optional bias correction.

**Parameters:**
| Parameter | Type | Description |
|-----------|------|-------------|
| `counts` | `Dict[Any, int]` or `List[int]` | Frequency counts |
| `method` | `str` | `"plugin"`, `"miller_madow"`, `"chao_shen"`, or `"jackknife"` |
| `bias_correction` | `bool` | Apply bias correction (plugin method only) |

**Returns** entropy estimate in bits.

### `mutual_information_estimator(x, y, method="plugin", bias_correction=True) -> float`

Estimate mutual information I(X;Y) = H(X) + H(Y) - H(X,Y) with bias-corrected entropy estimates. Ensures non-negative result.

### `kl_divergence_estimator(p, q, method="plugin", bias_correction=True) -> float`

Estimate KL divergence D_KL(P||Q) from sample lists. Converts to empirical probability distributions and computes divergence. Returns infinity when Q has zero probability for an event observed in P.

### `bias_correction(entropy, sample_size, alphabet_size) -> float`

Apply general bias correction: H_corrected = H - (d-1)/(2n). Standalone function for custom pipelines.

### `entropy_bootstrap_confidence(counts, method="plugin", n_bootstraps=1000, confidence_level=0.95, random_state=None) -> Dict`

Calculate bootstrap confidence interval for entropy.

**Returns** dict with:
- `entropy`: Point estimate
- `ci_lower`: Lower bound of confidence interval
- `ci_upper`: Upper bound of confidence interval
- `confidence_level`: The confidence level used
- `n_bootstraps`: Number of bootstrap samples

### `effective_sample_size_correction(entropy, sample_size, alphabet_size) -> float`

Apply heuristic effective sample size correction for dependencies in the data.

### `panzeri_treves_bias_correction(entropy, sample_size, alphabet_size, response_frequencies=None) -> float`

Apply the Panzeri-Treves (1996) bias correction. More accurate than Miller-Madow for small samples. Accepts optional response-specific frequencies; falls back to uniform assumption.

### `entropy_rate_estimator(sequence, order=1, method="plugin") -> float`

Estimate the entropy rate of a stationary process from a sequence. Uses block entropy differences: h = H(n+1-blocks) - H(n-blocks) for given Markov order.

## Usage Examples

```python
from metainformant.information.metrics.core.estimation import (
    entropy_estimator, mutual_information_estimator,
    kl_divergence_estimator, bias_correction,
    entropy_bootstrap_confidence, entropy_rate_estimator,
    panzeri_treves_bias_correction,
)

# Compare estimation methods on sparse data
counts = {"A": 50, "C": 30, "G": 15, "T": 3, "N": 1, "X": 1}
for method in ["plugin", "miller_madow", "chao_shen", "jackknife"]:
    h = entropy_estimator(counts, method=method)
    print(f"{method:15s}: {h:.4f} bits")

# Bias-corrected mutual information
gene_a = [0, 1, 0, 1, 1, 0, 1, 0, 0, 1]
gene_b = [0, 1, 1, 1, 1, 0, 0, 0, 0, 1]
mi = mutual_information_estimator(gene_a, gene_b, method="miller_madow")
print(f"MI (Miller-Madow): {mi:.4f} bits")

# Bootstrap confidence interval
ci = entropy_bootstrap_confidence(
    counts, method="plugin",
    n_bootstraps=2000, confidence_level=0.95, random_state=42
)
print(f"H = {ci['entropy']:.4f} [{ci['ci_lower']:.4f}, {ci['ci_upper']:.4f}]")

# Entropy rate of a DNA sequence
import numpy as np
sequence = list("ATCGATCGATCGATCG" * 10)
h_rate = entropy_rate_estimator(sequence, order=2, method="plugin")
print(f"Entropy rate (order 2): {h_rate:.4f} bits/symbol")

# Panzeri-Treves correction for neural/expression data
import numpy as np
h_raw = 2.5  # raw plugin estimate
h_pt = panzeri_treves_bias_correction(
    h_raw, sample_size=100, alphabet_size=10,
    response_frequencies=np.array([15, 12, 10, 10, 10, 10, 10, 8, 8, 7])
)
print(f"Panzeri-Treves corrected: {h_pt:.4f} bits")

# Manual bias correction
h_plugin = 2.3
h_corrected = bias_correction(h_plugin, sample_size=200, alphabet_size=20)
print(f"Corrected entropy: {h_corrected:.4f} bits")
```

## Configuration

Environment variable prefix: `INFO_`

## Related Modules

- `metainformant.information.metrics.core.syntactic` -- exact entropy from known distributions
- `metainformant.information.metrics.advanced.channel` -- channel capacity
- `metainformant.information.metrics.advanced.geometry` -- information geometry
- `metainformant.information.metrics.advanced.decomposition` -- multivariate decomposition
