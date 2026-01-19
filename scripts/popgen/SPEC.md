# SPEC: PopGen Scripts

Simulation and analysis of population genetics data.

## Core Workflows

- `simulate_populations.py`: Large-scale population simulations using coalescent or forward-time models.
- `calculate_fst_spectrum.py`: Computes Fst across genomic windows.

## Standards

- **Theory Integration**: Pulls core mathematical implementations from `metainformant.math.population_genetics`.
- **Performance**: Window-based calculations are parallelized across chromosomes.
