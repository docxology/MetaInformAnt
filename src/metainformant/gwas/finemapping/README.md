# GWAS Fine-Mapping

Statistical fine-mapping methods for identifying causal variants from GWAS signals, including credible set analysis, SuSiE regression, colocalization, and eQTL integration.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Package exports for `credible_sets`, `colocalization`, `eqtl` |
| `credible_sets.py` | ABF credible sets, SuSiE regression, Bayes factors, conditional analysis |
| `colocalization.py` | Multi-trait colocalization, GWAS-eQTL coloc, CLPP, regional coloc |
| `eqtl.py` | cis/trans eQTL scanning, conditional eQTL, effect size estimation |

## Key Functions

| Function | Description |
|----------|-------------|
| `compute_credible_set()` | ABF credible set from Z-scores with PIP calculation |
| `susie_regression()` | Sum of Single Effects regression for multi-causal fine-mapping |
| `compute_bayes_factors()` | Wakefield approximate Bayes factors |
| `colocalization()` | Pairwise trait colocalization (coloc-style) |
| `conditional_analysis()` | Stepwise conditional association analysis |
| `annotate_credible_set()` | Functional annotation enrichment of credible sets |
| `multi_trait_coloc()` | Multi-trait colocalization extending pairwise to N traits |
| `eqtl_coloc()` | GWAS-eQTL colocalization with gene-level summaries |
| `compute_clpp()` | Colocalization posterior probability from individual PIPs |
| `cis_eqtl_scan()` | Scan for cis-eQTLs within a window around each gene |
| `trans_eqtl_scan()` | Genome-wide trans-eQTL association testing |
| `eqtl_effect_sizes()` | Estimate eQTL effect sizes and confidence intervals |

## Usage

```python
from metainformant.gwas.finemapping import credible_sets, colocalization, eqtl

cs = credible_sets.compute_credible_set(z_scores, coverage=0.95)
susie = credible_sets.susie_regression(z_scores, ld_matrix, n_effects=5)
coloc = colocalization.multi_trait_coloc({"gwas": z1, "eqtl": z2})
results = eqtl.cis_eqtl_scan(expression, genotypes, gene_pos, var_pos)
```
