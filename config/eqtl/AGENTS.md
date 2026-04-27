# Agent Directives: config/eqtl

## Role

eQTL integration pipeline configuration linking GWAS variants with RNA-seq expression data.

## Contents

| File | Description |
| :--- | :--- |
| `eqtl_amellifera.yaml` | A. mellifera eQTL analysis config |

## Configuration Structure

```yaml
# eQTL-specific configuration
gwas_results: output/gwas/association_results.tsv
expression_matrix: output/amalgkit/amellifera/work/curate/curated_matrix.tsv
species: apis_mellifera
significance_threshold: 5e-8
```

## Rules

- Validate with schema before committing new configs
- Follow NO MOCKING policy — tests use real config files
- Use `uv` for dependency management
- Environment overrides use `AK_` prefix
