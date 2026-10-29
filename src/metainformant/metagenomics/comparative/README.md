# Comparative

Comparative metagenomics analysis providing differential abundance testing, indicator species analysis, effect size ranking, and ML-based biomarker discovery for microbiome studies.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Exports differential_abundance submodule |
| `differential_abundance.py` | ALDEx2, ANCOM, IndVal, LEfSe, and random forest methods |

## Key Functions

| Function | Description |
|----------|-------------|
| `differential_abundance.differential_abundance()` | Run differential abundance testing between groups |
| `differential_abundance.clr_transform()` | Centered log-ratio transformation for compositional data |
| `differential_abundance.indicator_species()` | IndVal analysis for group-associated taxa |
| `differential_abundance.effect_size_analysis()` | LEfSe-inspired univariate effect sizes with a Kruskal-Wallis test (not a multivariate LDA) |
| `differential_abundance.biomarker_discovery()` | ML-based biomarker discovery using feature importance; `cv_accuracy` is `None` (with `cv_unavailable_reason`) when cross-validation cannot run |

## Usage

```python
from metainformant.metagenomics.comparative import differential_abundance

results = differential_abundance.differential_abundance(
    counts, groups, taxa_names, method="aldex2_like"
)
indicators = differential_abundance.indicator_species(counts, groups, taxa_names, seed=42)
biomarkers = differential_abundance.biomarker_discovery(counts, groups, taxa_names, seed=42)
```

`counts` may be a list of lists, a NumPy array or a pandas DataFrame. The randomised procedures (`aldex2_like` Monte Carlo, permutations, random forest) are unseeded by default (`seed=None`); pass `seed` for reproducible results. Without SciPy the Kruskal-Wallis test is computed in pure Python for any number of groups, and Welch's t-test always uses the exact Student's t tail (a regularized incomplete beta function). A test that cannot be computed reports a NaN p-value (never 1.0), which is left out of the Benjamini-Hochberg correction; `ancom_like` reports the W statistic (`w_statistic`, `detected`) and no p-value.
