### Math: Selection Models

Functions: `kin_selection_response`, `multilevel_selection_decomposition`.

```python
from metainformant.math import selection

r_b_minus_c = selection.kin_selection_response(relatedness=0.5, benefit=0.4, cost=0.1)
between, within, total = selection.multilevel_selection_decomposition(
    group_means=[1.0, 1.1, 0.9],
    individual_deviations=[0.2, -0.1, 0.0],
    selection_strength_group=0.5,
    selection_strength_individual=0.3,
)
```
