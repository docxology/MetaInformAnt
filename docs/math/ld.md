### Math: Linkage Disequilibrium and Maps

Functions: `ld_coefficients`, `r_squared`, `ld_decay_r2`, `haldane_d_to_c`, `haldane_c_to_d`, `kosambi_d_to_c`, `kosambi_c_to_d`, `expected_r2_from_Ne_c`.

Example

```python
from metainformant.math.population_genetics.ld import ld_coefficients
from metainformant.math.core.utilities import r_squared
from metainformant.math.population_genetics.ld import haldane_d_to_c, kosambi_c_to_d
from metainformant.math.population_genetics.statistics import expected_r2_from_Ne_c

D, Dp = ld_coefficients(0.6, 0.4, 0.7, 0.3, 0.5)
r2 = r_squared(0.6, 0.4, 0.7, 0.3, 0.5)

c = haldane_d_to_c(0.01)  # map distance to recombination fraction
d = kosambi_c_to_d(c)     # inverse mapping under Kosambi

E_r2 = expected_r2_from_Ne_c(recombination_rate=0.01, Ne=1000)  # Sved (1971): 1 / (1 + 4 Ne c)
```

#### LD decay with distance: `ld_decay_r2`

`ld_decay_r2(distances, r_squared_values, max_distance=None)` summarises how r² falls with the distance between SNP pairs. It fits

    r² = r²₀ · exp(−k · d)

to the (distance, r²) pairs by bounded nonlinear least squares (r²₀ ≥ 0 and k ≥ 0, using SciPy). Every pair enters the fit with equal weight and nothing is binned. Distances may be in any unit (bp, kb, Mb): the rate comes back per that unit. `max_distance` drops the pairs farther apart than that.

The fit returns a dictionary:

| Key | Meaning |
|-----|---------|
| `decay_rate` | k, the exponential rate constant per unit of distance |
| `half_decay_distance` | ln 2 / k, the distance at which the fitted r² is half its value at distance 0 (`inf` when there is no decay) |
| `initial_r2` | r²₀, the fitted r² at distance 0 |

```python
import numpy as np
from metainformant.math.population_genetics.ld import ld_decay_r2

distances = np.linspace(1, 5000, 400)  # distance between each SNP pair, in bp
r2 = 0.8 * np.exp(-0.002 * distances)  # an exact exponential decay
fit = ld_decay_r2(list(distances), list(r2), max_distance=4000)
# fit["decay_rate"] = 0.002, fit["half_decay_distance"] = 346.57 bp, fit["initial_r2"] = 0.8 (to fitting precision)
```

Cases without a fitted decay:

- Flat or rising data (the exponential fits no better than the mean r²): `decay_rate` 0.0, `half_decay_distance` `inf` and `initial_r2` the mean r².
- Every r² equal to zero: the same three keys, with `initial_r2` 0.0.
- No pair left after `max_distance`, or fewer than three distinct distances: only `decay_rate` (0.0) and `half_decay_distance` (`inf`) are returned.

`ld_decay_r2` raises `ValueError` when the two lists differ in length, hold a non-finite value, or a distance is negative, and `RuntimeError` when the optimiser does not converge.

`k` is a descriptive decay constant. The model has no finite-sample floor (a sample of n haplotypes has E[r²] near 1/n at long range) and it is not the Hill and Weir (1988) expectation, so k does not estimate 4Nₑc; use `expected_r2_from_Ne_c` for the Sved (1971) expectation.

With a scalar first argument the function is a closed form instead: `ld_decay_r2(r2_start, recombination_rate=c, generations=t)` returns `r2_start * (1 - c) ** (2 * t)`, and both keyword arguments are required.
