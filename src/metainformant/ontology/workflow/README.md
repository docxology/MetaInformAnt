# workflow

Sub-package of `metainformant.ontology`: the stage-10 GWAS ontology pipeline. See the module's
`AGENTS.md`/`README.md` for the domain overview.

## Contents

| File | Purpose |
|------|---------|
| `__init__.py` | Package exports for `run_ontology` |
| `run_ontology.py` | `run_ontology_analysis(config, phenotype, model, project_root)`: GO/HPO enrichment of GWAS hits |

## Running

`metainformant ontology run --input <config.yaml> --phenotype <label> --model <label>` loads the YAML configuration
(it must be a mapping) and calls `run_ontology_analysis` with the working directory as the project root. The pipeline
reads `<paths.results_dir>/<phenotype>/<model>/summary_statistics.tsv` (the default `results_dir` is `output`) and,
when present, `post_gwas/post_gwas_results.json` in that same directory for gene annotations. It writes its results to
`<paths.results_dir>/<phenotype>/<model>/ontology/`: `ora_results.tsv`, `gsea_results.tsv`, `enrichment_dotplot.png`,
`pathway_network.png`, `hpo_mapping.json` and `ontology_summary.json`. The GO gene sets come from the QuickGO REST API
and the HPO terms from the EBI OLS API, so a run needs network access.

## Configuration

`run_ontology_analysis` reads the `ontology` section of the configuration. Every key is optional:

| Key | Default | Meaning |
|-----|---------|---------|
| `taxon_id` | `7460` | NCBI taxon whose GO annotations are fetched |
| `top_n_hits` | `20` | Number of top GWAS hits that supply the hit gene list (`top_n` of `gwas_hits_to_genes`) |
| `go_aspects` | all three GO aspects | GO aspects to fetch (`biological_process`, `molecular_function`, `cellular_component`) |
| `gsea_permutations` | `500` | Permutations per gene set; the GSEA step uses at most 200 of them |
| `min_go_set_size`, `max_go_set_size` | `3`, `500` | Size range of the GO gene sets tested by ORA and GSEA |
| `hpo_phenotype` | the `--phenotype` label | Text mapped to HPO terms |
| `fdr_threshold` | `0.05` | Adjusted p-value below which an ORA term counts as significant in the summary |
| `jaccard_threshold` | `0.3` | Minimum Jaccard similarity for an edge of the pathway network |
| `rate_limit_s` | `0.2` | Pause between API calls, in seconds |
| `max_go_pages` | `25` | Maximum number of QuickGO result pages fetched for the taxon-wide gene sets |
| `random_state` | unset | Seed of the GSEA permutations (see below) |

`random_state` is passed to `gsea(..., random_state=...)`. With the same value, the same ranked genes and the same GO
gene sets, the GSEA p-values, NES and FDR are reproducible, and the global `random` state is neither seeded nor
consumed. Leave it unset (or `null`) for unseeded permutations. The value goes through `int()`, so an integer or a
numeric string such as `"11"` (as a YAML or environment value may arrive) is accepted. The ORA step draws no random
numbers, and the pathway network plot uses the default layout seed of `pathway_network_plot` (`layout_seed=42`).

```yaml
ontology:
  taxon_id: 7460
  gsea_permutations: 500
  random_state: 42   # seed of the GSEA permutations; omit for unseeded runs
paths:
  results_dir: output
```
