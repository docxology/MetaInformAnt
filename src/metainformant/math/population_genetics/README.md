# Population Genetics

This submodule contains mathematical models and statistical tools for population genetics.

## Components

- **core.py**: Fundamental Hardy-Weinberg and mutation-selection balance models.
- **statistics.py**: Summary statistics (diversity, Tajima's D, etc.).
- **coalescent.py**: Coalescent theory simulations and theory.
- **demography.py**: Demographic models (growth, bottlenecks).
- **effective_size.py**: Effective population size calculations.
- **fst.py**: Fixation index (Fst) calculations.
- **ld.py**: Linkage disequilibrium analysis.
- **selection.py**: Selection coefficients and response to selection.

## Usage

```python
from metainformant.math.population_genetics import tajimas_D, expected_pairwise_diversity

# Calculate Tajima's D
d = tajimas_D(segregating_sites=[1, 0, 1])
```
