# Comparative Metagenomics

Differential abundance testing, compositional data analysis, indicator species analysis, effect size ranking (LEfSe-style), and machine-learning biomarker discovery for metagenomic data.

## Key Concepts

**Compositional data analysis** addresses the constraint that microbiome abundances are relative (sum to a constant). The Centered Log-Ratio (CLR) transform projects counts into an unconstrained space suitable for standard statistics: CLR(x_i) = log(x_i) - mean(log(x)).

**ALDEx2-like analysis** draws `n_monte_carlo` Dirichlet Monte Carlo instances per sample (pseudocount 0.5 prior), averages their CLR transforms, and applies Welch's t-test and Cohen's d, accounting for the compositional nature and the count-sampling noise of sequencing data.

**ANCOM-like analysis** uses pairwise log-ratio comparisons to identify taxa differentially abundant relative to most other taxa, providing robust results under varying reference frames.

**LEfSe-style effect size** computes a Kruskal-Wallis p-value for every taxon (no significance filter is applied) and ranks taxa by a simplified LDA-like score on log10 scale: the standardized difference of group means on CLR-transformed data, as log10(1 + score). With more than two groups the score is instead the largest pairwise difference of group means, unstandardized.

**Biomarker discovery** uses random forest classification with cross-validation to rank taxa by feature importance, identifying the most discriminative taxa between groups.

## Function Reference

All functions live in `metainformant.metagenomics.comparative.differential_abundance` (the `metainformant.metagenomics.comparative` package exposes that module, not the functions). Inputs share one layout: `counts` is a samples x taxa list of lists of counts (not a numpy array: most functions raise `ValueError` for one, so pass `array.tolist()`), `groups` holds one integer label per sample, and `taxa_names` holds one name per column. Mismatched lengths raise `ValueError`.

### `differential_abundance(counts, groups, taxa_names, method="aldex2_like", n_monte_carlo=128, seed=None) -> List[Dict]`

Test differential abundance between two groups. `groups` must contain exactly two distinct labels; the lower label is "group 1". `n_monte_carlo` (at least 1) and `seed` apply to `aldex2_like` only; without `seed` its p-values vary between runs.

**Methods:**
| Method | Description |
|--------|-------------|
| `aldex2_like` | Dirichlet Monte Carlo CLR (point-estimate CLR if numpy is missing) + Welch's t-test + Cohen's d effect size |
| `ancom_like` | Pairwise log-ratio comparisons (Welch's t-test, p < 0.05 per taxon pair); `effect_size` is the W fraction of rejected comparisons and `p_value` is `1 - W` |
| `simple_deseq` | Median-of-ratios normalization + Welch's t-test and Cohen's d on the normalized counts (no negative binomial model is fitted) |

**Returns** list of dicts sorted by adjusted p-value, each with: `taxon`, `log2fc` (log2((mean_group1 + 1) / (mean_group2 + 1)); positive when higher in group 1; computed on normalized means for `simple_deseq`), `p_value`, `adjusted_p` (BH-corrected), `effect_size`, `mean_group1`, `mean_group2` (means of the raw counts).

### `clr_transform(counts, pseudocount=0.5) -> List[List[float]]`

Apply centered log-ratio transformation. Adds pseudocount before log to handle zeros.

### `indicator_species(counts, groups, taxa_names, n_permutations=999, seed=None) -> List[Dict]`

IndVal indicator species analysis. Combines specificity (a taxon's mean abundance in a group relative to the sum of its group means) and fidelity (fraction of the group's samples where it occurs). Significance via permutation test; works with more than two groups.

**Returns** per-taxon dicts, sorted by `indicator_value` (descending), with `taxon`, `indicator_value` (0--1), `p_value`, `associated_group` (the group label with the highest IndVal).

### `effect_size_analysis(counts, groups, taxa_names) -> List[Dict]`

LEfSe-style analysis: Kruskal-Wallis p-value plus a simplified LDA-like effect size on CLR-transformed data. All taxa are scored and returned; the p-values are not used for filtering and are not multiple-testing corrected.

**Returns** per-taxon dicts, sorted by `lda_score` (descending), with `taxon`, `lda_score` (log10 scale), `p_value`, `direction` (`"enriched_in_<label>"`, for example `"enriched_in_0"`).

### `biomarker_discovery(counts, groups, taxa_names, method="random_forest", n_estimators=100, cv_folds=5, seed=None) -> Dict`

Machine-learning biomarker identification using random forest on CLR-transformed data (`method="random_forest"` is the only accepted value). Falls back to Cohen's d ranking when scikit-learn is unavailable.

**Returns** dict with:
- `selected_taxa`: Top-ranked taxa (importance > mean).
- `importances`: Per-taxon importance scores (dict of taxon -> score).
- `cv_accuracy`: Cross-validated classification accuracy (0.0 in the Cohen's d fallback).
- `model_summary`: For the random forest: `method`, `n_estimators`, `cv_folds`, `cv_scores`, `n_features`, `n_selected`. `cv_folds` is the number of folds actually used: the requested value capped by the number of distinct groups and by the size of the smallest group, and at least 2. The fallback reports `method="effect_size_ranking"`, a `note`, `n_features`, and `n_selected`.

## Usage Examples

```python
from metainformant.metagenomics.comparative.differential_abundance import (
    differential_abundance, clr_transform,
    indicator_species, effect_size_analysis,
    biomarker_discovery,
)

# Count matrix: 6 samples x 5 taxa
counts = [
    [100, 20, 5, 0, 50],
    [90, 25, 8, 2, 45],
    [110, 18, 3, 1, 55],
    [10, 80, 60, 40, 5],
    [15, 75, 55, 35, 8],
    [12, 85, 65, 45, 3],
]
groups = [0, 0, 0, 1, 1, 1]
taxa = ["Bacteroides", "Prevotella", "Ruminococcus", "Faecalibacterium", "Akkermansia"]

# ALDEx2-like differential abundance
da_results = differential_abundance(counts, groups, taxa, method="aldex2_like", seed=42)
for r in da_results:
    if r["adjusted_p"] < 0.05:
        print(f"{r['taxon']}: log2FC={r['log2fc']:.2f}, q={r['adjusted_p']:.4f}")

# CLR transformation
clr_data = clr_transform(counts)

# Indicator species
indicators = indicator_species(counts, groups, taxa, n_permutations=999, seed=42)
for ind in indicators[:3]:
    print(f"{ind['taxon']}: IndVal={ind['indicator_value']:.3f}, p={ind['p_value']:.4f}")

# LEfSe-style effect sizes
lefse = effect_size_analysis(counts, groups, taxa)
for r in lefse[:3]:
    print(f"{r['taxon']}: LDA={r['lda_score']:.3f}, {r['direction']}")

# Random forest biomarkers
biomarkers = biomarker_discovery(counts, groups, taxa, seed=42)
print(f"Selected biomarkers: {biomarkers['selected_taxa']}")
print(f"CV accuracy: {biomarkers['cv_accuracy']:.3f}")
```

## Statistical Methods

### Benjamini-Hochberg FDR Correction

The p-values from `differential_abundance` (all three methods) are adjusted using the Benjamini-Hochberg procedure to control the false discovery rate and returned as `adjusted_p`. The p-values from `indicator_species` and `effect_size_analysis` are raw.

### Welch's t-test

Used for group comparisons with unequal variances. Falls back to a pure Python implementation (normal approximation) when scipy is unavailable.

### Kruskal-Wallis test

Used by `effect_size_analysis` (`scipy.stats.kruskal`). Without scipy, a rough rank-based approximation (`exp(-H/2)`) is used instead.

### Cohen's d Effect Size

Standardized mean difference: d = (mean_A - mean_B) / pooled_std. Used for ranking and as a fallback importance metric.

## Configuration

This module reads no environment variables; the `META_` prefix used elsewhere in the package does not apply here.

## Optional Dependencies

All three are installed with METAINFORMANT, but the code falls back when they are missing:

- `scipy` -- Welch's t-test, Kruskal-Wallis test (pure-Python approximations otherwise)
- `numpy` -- Dirichlet Monte Carlo draws in `aldex2_like` (single-point CLR otherwise)
- `scikit-learn` -- random forest and cross-validation in `biomarker_discovery` (Cohen's d ranking otherwise)

## Related Modules

- `metainformant.metagenomics.diversity` -- alpha/beta diversity
- `metainformant.metagenomics.amplicon` -- OTU/ASV input data
- `metainformant.metagenomics.shotgun` -- community profiling
- `metainformant.ecology.analysis.indicators` -- equivalent ecological tests
