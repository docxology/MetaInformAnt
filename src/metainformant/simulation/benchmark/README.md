# Simulation Benchmark

Synthetic benchmark dataset generators for evaluating bioinformatics methods. Produces datasets with known ground truth for classification, regression, clustering, differential expression, and GWAS benchmarks, with built-in evaluation and method comparison utilities.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Re-exports `generators` module |
| `generators.py` | Dataset generators, evaluation functions, and benchmark suite runner |

## Key Functions

| Function | Description |
|----------|-------------|
| `generate_benchmark_dataset()` | Generate synthetic dataset for classification/regression/clustering |
| `generate_synthetic_variants()` | Create synthetic GWAS variant data with known associations |
| `generate_synthetic_expression()` | Generate synthetic RNA-seq expression with DE genes |
| `evaluate_benchmark()` | Evaluate method predictions against known ground truth |
| `benchmark_suite()` | Run multiple methods on a benchmark and compare results |

## Usage

```python
from metainformant.simulation.benchmark import generators

data = generators.generate_benchmark_dataset(task="classification", n_samples=500, seed=42)
expression = generators.generate_synthetic_expression(n_genes=1000, n_samples=20)
report = generators.benchmark_suite(dataset, methods=[method_a, method_b])
```
