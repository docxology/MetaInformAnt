### Math: Price Equation

Functions: `expectation`, `covariance`, `price_equation`.

```mermaid
graph TD
  A[fitness w] & B[trait z] --> C[covariance]
  D[offspring trait z'] --> E[transmission]
  C & E --> F[total change]
```

Example

```python
from metainformant.math import price

cov, trans, total = price.price_equation([1.0, 1.2, 0.9], [0.2, 0.4, 0.1], [0.25, 0.35, 0.15])
```


