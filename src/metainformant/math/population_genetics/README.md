# Population Genetics

Mathematical foundations for population genetics: coalescent theory, demography, F-statistics, linkage disequilibrium, selection, and summary statistics.

## Contents

| File | Purpose |
|------|---------|
| `coalescent.py` | Coalescent theory: TMRCA, Watterson theta, Tajima's D, SFS |
| `demography.py` | Growth models: exponential, logistic, age structure, bottleneck |
| `effective_size.py` | Effective population size under sex ratio and demographic scenarios |
| `fst.py` | F-statistics: FST from allele frequencies, Weir's FST, pairwise matrices |
| `ld.py` | Linkage disequilibrium coefficients, LD decay, Haldane/Kosambi maps |
| `selection.py` | Selection theory: breeder's equation, Hamilton's rule, mutation-selection balance |
| `statistics.py` | Summary statistics: diversity, segregating sites, fixation probability |

## Key Functions

| Function | Description |
|----------|-------------|
| `expected_time_to_mrca()` | Expected coalescent time to most recent common ancestor |
| `watterson_theta()` | Watterson's estimator of population mutation rate |
| `tajimas_D()` | Tajima's D test for neutrality |
| `exponential_growth_model()` | Simulate exponential population growth trajectory |
| `fst_from_allele_freqs()` | FST between two populations from allele frequencies |
| `ld_coefficients()` | D, D-prime, and r-squared from allele or genotype data |
| `breeders_equation_response()` | Response to selection via breeder's equation |
| `kin_selection_response()` | Hamilton's rule for kin selection |
| `fixation_probability()` | Probability of allele fixation under selection |

## Usage

```python
from metainformant.math.population_genetics.coalescent import expected_time_to_mrca
from metainformant.math.population_genetics.fst import fst_from_allele_freqs
from metainformant.math.population_genetics.selection import breeders_equation_response

tmrca = expected_time_to_mrca(n_samples=10, effective_size=10000)
fst = fst_from_allele_freqs([0.6, 0.4, 0.8], [0.3, 0.7, 0.2])
response = breeders_equation_response(heritability=0.5, selection_differential=2.0)
```
