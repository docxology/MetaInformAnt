# Comparative Metagenomics

Differential abundance testing, compositional data analysis, indicator species analysis, effect size ranking (LEfSe-style), and machine-learning biomarker discovery for metagenomic data.

## Key Concepts

**Compositional data analysis** addresses the constraint that microbiome abundances are relative (sum to a constant). The Centered Log-Ratio (CLR) transform projects counts into an unconstrained space suitable for standard statistics: CLR(x_i) = log(x_i) - mean(log(x)).

**ALDEx2-like analysis** combines CLR transformation with Welch's t-test and effect size estimation for differential abundance, accounting for the compositional nature of sequencing data.

**ANCOM-like analysis** uses pairwise log-ratio comparisons to identify taxa differentially abundant relative to most other taxa, providing robust results under varying reference frames.

**LEfSe-style effect size** first identifies significantly different taxa (Kruskal-Wallis), then ranks them by a simplified Linear Discriminant Analysis score on log10 scale.

**Biomarker discovery** uses random forest classification with cross-validation to rank taxa by feature importance, identifying the most discriminative taxa between groups.

## Function Reference

### `differential_abundance(counts, groups, taxa_names, method="aldex2_like", n_monte_carlo=128) -> List[Dict]`

Test differential abundance between two groups.

**Methods:**
| Method | Description |
|--------|-------------|
| `aldex2_like` | CLR + Welch's t-test + Cohen's d effect size |
| `ancom_like` | Pairwise log-ratio comparisons, W statistic |
| `simple_deseq` | Median-of-ratios normalization + Wald-like test |

**Returns** list of dicts sorted by adjusted p-value, each with: `taxon`, `log2fc`, `p_value`, `adjusted_p` (BH-corrected), `effect_size`, `mean_group1`, `mean_group2`.

### `clr_transform(counts, pseudocount=0.5) -> List[List[float]]`

Apply centered log-ratio transformation. Adds pseudocount before log to handle zeros.

### `indicator_species(counts, groups, taxa_names, n_permutations=999, seed=None) -> List[Dict]`

IndVal indicator species analysis. Combines specificity (how concentrated a taxon is in a group) and fidelity (how consistently it occurs). Significance via permutation test.

**Returns** per-taxon dicts with `indicator_value` (0--1), `p_value`, `associated_group`.

### `effect_size_analysis(counts, groups, taxa_names) -> List[Dict]`

LEfSe-style analysis: Kruskal-Wallis test followed by simplified LDA effect size on CLR-transformed data.

**Returns** per-taxon dicts with `lda_score` (log10 scale), `p_value`, `direction`.

### `biomarker_discovery(counts, groups, taxa_names, method="random_forest", n_estimators=100, cv_folds=5) -> Dict`

Machine-learning biomarker identification using random forest on CLR-transformed data. Falls back to Cohen's d ranking when scikit-learn is unavailable.

**Returns** dict with:
- `selected_taxa`: Top-ranked taxa (importance > mean).
- `importances`: Per-taxon importance scores.
- `cv_accuracy`: Cross-validated classification accuracy.
- `model_summary`: Method parameters and metrics.

## Usage Examples

```python
from metainformant.metagenomics import (
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
da_results = differential_abundance(counts, groups, taxa, method="aldex2_like")
for r in da_results:
    if r["adjusted_p"] < 0.05:
        print(f"{r['taxon']}: log2FC={r['log2fc']:.2f}, q={r['adjusted_p']:.4f}")

# CLR transformation
clr_data = clr_transform(counts)

# Indicator species
indicators = indicator_species(counts, groups, taxa, n_permutations=999)
for ind in indicators[:3]:
    print(f"{ind['taxon']}: IndVal={ind['indicator_value']:.3f}, p={ind['p_value']:.4f}")

# LEfSe-style effect sizes
lefse = effect_size_analysis(counts, groups, taxa)
for r in lefse[:3]:
    print(f"{r['taxon']}: LDA={r['lda_score']:.3f}, {r['direction']}")

# Random forest biomarkers
biomarkers = biomarker_discovery(counts, groups, taxa)
print(f"Selected biomarkers: {biomarkers['selected_taxa']}")
print(f"CV accuracy: {biomarkers['cv_accuracy']:.3f}")
```

## Statistical Methods

### Benjamini-Hochberg FDR Correction

All p-values from differential abundance tests are adjusted using the Benjamini-Hochberg procedure to control the false discovery rate.

### Welch's t-test

Used for group comparisons with unequal variances. Falls back to a pure Python implementation (normal approximation) when scipy is unavailable.

### Cohen's d Effect Size

Standardized mean difference: d = (mean_A - mean_B) / pooled_std. Used for ranking and as a fallback importance metric.

## Configuration

Environment variable prefix: `META_`

## Optional Dependencies

- `scipy` -- Welch's t-test, Kruskal-Wallis test
- `numpy` -- array operations
- `scikit-learn` -- random forest biomarker discovery (falls back to effect size ranking)

## Related Modules

- `metainformant.metagenomics.diversity` -- alpha/beta diversity
- `metainformant.metagenomics.amplicon` -- OTU/ASV input data
- `metainformant.metagenomics.shotgun` -- community profiling
- `metainformant.ecology.indicators` -- equivalent ecological tests
