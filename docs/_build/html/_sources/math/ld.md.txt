### Math: Linkage Disequilibrium and Maps

Functions: `ld_coefficients`, `r_squared`, `ld_decay_r2`, `haldane_d_to_c`, `haldane_c_to_d`, `kosambi_d_to_c`, `kosambi_c_to_d`, `expected_r2_from_Ne_c`.

Example

```python
from metainformant.math import ld_coefficients, r_squared
from metainformant.math import haldane_d_to_c, kosambi_c_to_d, expected_r2_from_Ne_c

D, Dp = ld_coefficients(0.6, 0.4, 0.7, 0.3, 0.5)
r2 = r_squared(0.6, 0.4, 0.7, 0.3, 0.5)

c = haldane_d_to_c(0.01)  # map distance to recombination fraction
d = kosambi_c_to_d(c)     # inverse mapping under Kosambi

E_r2 = expected_r2_from_Ne_c(Ne=1000, recombination_fraction=0.01)
```
