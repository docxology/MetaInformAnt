# Population

Population genetics analysis including diversity statistics, neutrality tests, F_ST estimation, linkage disequilibrium, and specialized visualization.

## Contents

| File | Purpose |
|------|---------|
| `core.py` | Core population statistics: allele frequencies, heterozygosity, Tajima's D, F_ST, LD |
| `analysis.py` | Extended analyses: neutrality test suites with normal-approximation or simulated p-values, selection detection, population structure |
| `visualization.py` | Population genetics plots: F_ST heatmaps, SFS, PCA, LD decay, diversity comparisons |

## Key Functions

| Function | Description |
|----------|-------------|
| `allele_frequencies()` | Compute allele frequencies from a genotype matrix |
| `nucleotide_diversity()` | Average pairwise differences (pi) across sequences |
| `tajimas_d()` | Tajima's D statistic for testing neutral evolution |
| `wattersons_theta()` | Watterson's estimator of population mutation rate |
| `hudson_fst()` | F_ST between two populations using Hudson's method |
| `observed_heterozygosity()` | Fraction of heterozygous genotypes |
| `linkage_disequilibrium()` | D' coefficient between two polymorphic sites |
| `calculate_fst()` | Pairwise F_ST from raw sequence data |
| `fu_and_li_d_star_from_sequences()`, `fu_and_li_f_star_from_sequences()` | Fu and Li's D* and F* (no outgroup) |
| `normalized_fay_wu_h()` | Zeng et al.'s normalized Fay and Wu's H from an outgroup-polarized spectrum |
| `detect_selection()` | Run Tajima's D, Fu and Li's D*, or the McDonald-Kreitman test (normal-approximation p-values) |
| `calculate_tajima_d()`, `calculate_fu_li_d()`, `calculate_fu_li_f()`, `calculate_fay_wu_h()` | A statistic and its two-sided p-value; `p_value_method="normal"` (default) or `"simulation"`, with `n_simulations` and `random_state` |
| `neutrality_test_suite()` | Tajima's D, Fu and Li's D*/F* and Fay and Wu's H with p-values and interpretation; takes the same `p_value_method`, `n_simulations` and `random_state` options |
| `detect_population_structure()` | Infer K subpopulations from sequence data |
| `plot_fst_matrix()` | Heatmap of pairwise F_ST values |
| `plot_allele_frequency_spectrum()` | Site frequency spectrum bar chart |
| `plot_pca_results()` | PCA scatter plot colored by population |
| `plot_ld_decay()` | LD (r-squared) as a function of physical distance |

## Usage

```python
from metainformant.dna.population.core import nucleotide_diversity, tajimas_d, hudson_fst
from metainformant.dna.population.analysis import neutrality_test_suite

pi = nucleotide_diversity(sequences)
d_stat = tajimas_d(sequences)
fst = hudson_fst(pop1_sequences, pop2_sequences)
results = neutrality_test_suite(sequences)

# Monte Carlo p-values against a simulated neutral coalescent null (seeded)
simulated = neutrality_test_suite(sequences, p_value_method="simulation", n_simulations=2000, random_state=1)
```

## Neutrality p-values

With `p_value_method="normal"` (default) a p-value is `erfc(abs(z) / sqrt(2))` of the standardized statistic *z*. With `"simulation"` it is a two-sided Monte Carlo p-value against `n_simulations` (default 2000) samples from the standard neutral coalescent (constant size, infinite sites, no recombination) with the same number of sequences, each conditioned on the observed number of segregating sites (the `ms -s` approach): twice the smaller one-sided tail, `(1 + #{replicates at least as extreme}) / (n_simulations + 1)`, capped at 1. `random_state` (an integer, `random.Random` or `numpy.random.Generator`) seeds the replicates; `None` draws fresh entropy.

The null is a constant-size neutral population, so a small p-value can reflect demography as well as selection; the Fay and Wu null assumes the outgroup polarizes every site correctly. `calculate_summary_statistics()` and `detect_selection()` use the normal approximation.

The replicates are drawn by `dna/integration/neutrality_null.py`, which `analysis.py` imports lazily, so `import metainformant.dna` does not load the simulation package. See [Neutrality p-values](../../../../docs/dna/population.md#neutrality-p-values-analysis).
