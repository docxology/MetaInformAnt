# Population

Population genetics analysis including diversity statistics, neutrality tests, F_ST estimation, linkage disequilibrium, and specialized visualization.

## Contents

| File | Purpose |
|------|---------|
| `core.py` | Core population statistics: allele frequencies, heterozygosity, Tajima's D, F_ST, LD |
| `analysis.py` | Extended analyses: neutrality test suites, selection detection, population structure |
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
| `detect_selection()` | Run Tajima's D, Fu-Li, or Fay-Wu tests for selection |
| `neutrality_test_suite()` | Combined D, F, and H tests with interpretation |
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
```
