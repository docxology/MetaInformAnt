### Math: Epidemiology

Functions: `sir_step`, `seir_step`, `sis_step`, `basic_reproduction_number`, `effective_reproduction_number`, `herd_immunity_threshold`.

Example

```python
from metainformant.math import sir_step, seir_step, sis_step
from metainformant.math import basic_reproduction_number, effective_reproduction_number, herd_immunity_threshold

Sn, In, Rn = sir_step(0.99, 0.01, 0.0, beta=0.5, gamma=0.25, dt=0.1)
Sn, En, In, Rn = seir_step(0.99, 0.0, 0.01, 0.0, beta=0.5, sigma=0.2, gamma=0.25, dt=0.1)
Sn, In = sis_step(0.99, 0.01, beta=0.5, gamma=0.25, dt=0.1)

R0 = basic_reproduction_number(0.5, 0.25)
Re = effective_reproduction_number(R0, susceptible_fraction=0.9)
HIT = herd_immunity_threshold(R0)
```


