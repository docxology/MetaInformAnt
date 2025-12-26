# Math: Overview

Analytical helpers for evolutionary and behavioral models.

- Price equation
- Kin and multilevel selection
- Drift–Diffusion Model (DDM)
- Population dynamics (logistic growth, Lotka-Volterra)
- Coalescent theory helpers (π, Watterson's θ, Tajima's D, SFS)
- Linkage disequilibrium and genetic maps (Haldane, Kosambi)
- Population genetics (drift, mutation–selection balance, migration)
- Epidemiology (SIR/SEIR/SIS, R0, herd immunity)

See: [Price](./price.md), [Selection](./selection.md), [DDM](./ddm.md), [Dynamics](./dynamics.md), [LD](./ld.md), [Coalescent](./coalescent.md), [Epidemiology](./epidemiology.md), [Population genetics](./popgen.md).

Run examples with uv:

```bash
uv run python -c "from metainformant.math import coalescent as c; print(c.expected_time_to_mrca(10,1e3))"
```