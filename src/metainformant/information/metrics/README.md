# Information Metrics

Information theory library: core measures (entropy, MI, KL divergence), advanced methods (channel capacity, geometry, decomposition), and sequence analysis.

## Contents

| Directory | Purpose |
|-----------|---------|
| `core/` | Discrete and continuous entropy, MI, KL divergence, estimation |
| `advanced/` | Channel capacity, information geometry, PID, hypothesis testing |
| `analysis/` | Sequence information profiles, complexity, dataset comparison |

### core/

| File | Key Functions |
|------|---------------|
| `syntactic.py` | `shannon_entropy()`, `mutual_information()`, `kl_divergence()`, `jensen_shannon_divergence()`, `transfer_entropy()`, `renyi_entropy()` |
| `continuous.py` | `differential_entropy()`, `mutual_information_continuous()`, `copula_entropy()`, `information_flow_network()` |
| `estimation.py` | `entropy_estimator()` (plugin, Miller-Madow, Chao-Shen, jackknife), `entropy_bootstrap_confidence()` |

### advanced/

| File | Key Functions |
|------|---------------|
| `channel.py` | `channel_capacity()`, `rate_distortion()`, `information_bottleneck()` |
| `decomposition.py` | `partial_information_decomposition()`, `co_information()`, `o_information()` |
| `hypothesis.py` | `mi_permutation_test()`, `independence_test()`, `entropy_confidence_interval()` |
| `semantic.py` | `semantic_similarity()`, `information_content()`, `semantic_entropy()` |
| `geometry.py` | Re-exports from `fisher_rao.py` and `information_projection.py` |

### analysis/

| File | Key Functions |
|------|---------------|
| `analysis.py` | `information_profile()`, `information_signature()`, `analyze_sequence_information()` |
| `advanced_analysis.py` | `fisher_information()`, `variation_of_information()` |

## Usage

```python
from metainformant.information.metrics.core.syntactic import shannon_entropy, mutual_information
from metainformant.information.metrics.advanced.channel import channel_capacity

h = shannon_entropy([0.5, 0.25, 0.25])
mi = mutual_information(x_data, y_data)
cap = channel_capacity(transition_matrix)
```
