# Quantitative Genetics

This submodule focuses on the evolution of continuous traits.

## Components

- **core.py**: Heritability and response to selection (Breeder's equation support).
- **price.py**: The Price equation for decomposing evolutionary change.

## Usage

```python
from metainformant.math.quantitative_genetics import price_equation

# Decompose evolutionary change
cov, trans, total = price_equation(fitness, parent_traits, offspring_traits)
```
