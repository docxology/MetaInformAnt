# Specification: core

## 🎯 Scope
Core information theory metrics.

## 🧱 Architecture
- **Dependency Level**: Core
- **Component Type**: Source Code

## 💾 Data Structures
- **Modules**: 4 Python modules
- **Key Concepts**: Refer to Pydantic models in source.

## 🔌 API Definition
### Exports
- `__init__.py`
- `continuous.py`
- `estimation.py`
- `syntactic.py`

### Estimator contracts

**Units**: `continuous.py` returns **nats**; `estimation.py` returns **bits**.

**`continuous.differential_entropy(samples, method, bins)`** — 1D input is a
single variable; 2D `(n, d)` input is joint samples and the JOINT differential
entropy is estimated in the full d-dimensional space (inputs are never
flattened).
- `"histogram"`: equipartition histogram, `H = -sum_i p_i log(p_i / V)` with
  bin hypervolume `V` (Sturges bins when `bins is None`).
- `"kde"`: leave-one-out Gaussian KDE, Scott per-dimension bandwidth.
- `"knn"`: Kozachenko-Leonenko (KSG-style) k-NN estimator,
  `H = psi(n) - psi(k) + log(c_d) + (d/n) sum_i log(eps_i)`, `k=3`,
  `c_d = pi^(d/2) / Gamma(d/2 + 1)`, `eps_i` = distance to the k-th NN.
  1D private helpers (`_differential_entropy_{histogram,kde,knn}`) reject
  `ndim > 1`; joint data must go through the public dispatcher.

`mutual_information_continuous`, `conditional_entropy_continuous`,
`conditional_entropy_continuous_3d` and `transfer_entropy_continuous` estimate
all joint entropies in the true joint space.

**`estimation.entropy_estimator(counts, method, bias_correction)`** — plugin
(`+ (k-1)/(2n*ln2)` Miller-Madow when `bias_correction`), `miller_madow`,
`chao_shen` (Chao & Shen 2003: coverage `C = 1 - f1/n`, adjusted probs
`pa = C*p`, inclusion probs `la = 1-(1-pa)^n`, `H = -sum pa*log2(pa)/la`;
returns 0 when `C <= 0`), `jackknife`.

**`estimation.mutual_information_estimator(x, y, ...)`** — first-order MI bias
correction uses the POSSIBLE joint alphabet `|X|*|Y|`:
`I_corrected = I_plugin + (|X||Y| - |X| - |Y| + 1)/(2n*ln2)` bits; result
clipped at 0.

**`estimation.bias_correction(entropy, n, d)`** — canonical additive
Miller-Madow helper; `effective_sample_size_correction` is a
backward-compatible alias of it.

**`estimation.panzeri_treves_bias_correction(entropy, n, d, freqs)`** — real
PT (1996) correction: Bayesian support estimate `R` (occupied responses plus
estimated unobserved quasi-responses, faithful to pyentropy's
`pt_bayescount`), `H_PT = H_plugin + (R-1)/(2n*ln2)` bits (ADDED; reduces to
Miller-Madow when all alphabet responses are observed).
