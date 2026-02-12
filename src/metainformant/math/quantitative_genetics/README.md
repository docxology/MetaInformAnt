# Quantitative Genetics

Heritability estimation, selection response, and Price equation decomposition for quantitative trait analysis.

## Contents

| File | Purpose |
|------|---------|
| `core.py` | Heritability calculations and Lande's multivariate response equation |
| `price.py` | Price equation decomposition and weighted statistical functions |

## Key Functions

| Function | Description |
|----------|-------------|
| `narrow_sense_heritability()` | Calculate h-squared from additive genetic and phenotypic variance |
| `realized_heritability()` | Heritability from selection experiment response |
| `lande_equation_response()` | Multivariate response to selection via Lande equation |
| `price_equation()` | Decompose evolutionary change into selection and transmission |
| `delta_mean_trait()` | Change in mean trait value under selection |
| `weighted_variance()` | Variance weighted by fitness values |
| `weighted_covariance()` | Fitness-weighted covariance between traits |

## Usage

```python
from metainformant.math.quantitative_genetics.core import narrow_sense_heritability
from metainformant.math.quantitative_genetics.price import price_equation

h2 = narrow_sense_heritability(additive_genetic_variance=2.0, phenotypic_variance=5.0)
cov_term, trans_term, total = price_equation(
    fitness=[1.2, 0.8, 1.0],
    parent=[5.0, 3.0, 4.0],
    offspring=[5.5, 2.8, 4.1],
)
```
