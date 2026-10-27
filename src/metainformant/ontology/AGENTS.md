# Agent Directives: ontology

**Context**: Gene ontology and functional annotation module for METAINFORMANT.

## Capabilities

This module provides functionality organized into the following structure:

## Subpackages

- `annotation/` — exports: `build_background_from_vcf_genes`, `genes_to_go_annotations`, `gwas_hits_to_genes`, `rank_genes_by_pvalue`
- `core/` — exports: `go`, `obo`, `types`
- `pathway_enrichment/` — exports: `enrichment`
- `query/` — exports: `query`, `serialize`
- `visualization/` — exports: `enrichment_dotplot`, `enrichment_score_plot`, `pathway_network_plot`
- `workflow/` — exports: `run_ontology`

## Reproducibility

- `pathway_enrichment.enrichment.gsea(..., random_state=...)` draws its permutation gene sets from a generator of its
  own (an integer seed or a `random.Random`); the global `random` state is neither seeded nor consumed. The stage-10
  workflow reads the `ontology.random_state` configuration key and passes it to `gsea` (see `workflow/README.md`).
- Layout seeds: `visualization.visualization.plot_go_dag`, `plot_ontology_network` and `create_interactive_go_network`
  take `layout_seed` (default `None`: the layout is unseeded and draws from NumPy's global generator, as networkx does).
  `visualization.plots.pathway_network_plot` takes `layout_seed` too, with the default `42`, so its layout is
  reproducible out of the box. The defaults differ by design: `visualization.py` stays unseeded unless a seed is given,
  `plots.py` is seeded unless another seed is given. A given seed never touches the global random state.

## Rules

- Use `metainformant.core.utils.logging` for all logging
- Use `metainformant.core.io` for domain data file I/O. Direct stdlib parsing is allowed in core, protocol adapters, subprocess/CLI glue, and narrow parser internals when covered by tests.
- Follow REAL IMPLEMENTATION policy — all tests must use real implementations
- Use `uv` for dependency management

## Related Documentation

- **Module guide**: [../../../docs/ontology/](../../../docs/ontology/) — In-depth usage, architecture, and examples
- **API reference**: [SPEC.md](SPEC.md) — Type signatures, data structures, error codes
- **Core infrastructure**: [../core/AGENTS.md](../core/AGENTS.md) — Shared utilities (logging, config, I/O)
- **Full module index**: [../../../docs/index.md](../../../docs/index.md) — Overview of all METAINFORMANT modules
- **Phenotype module**: [../phenotype/AGENTS.md](../phenotype/AGENTS.md)
- **Epigenome module**: [../epigenome/AGENTS.md](../epigenome/AGENTS.md)
- **Visualization module**: [../visualization/AGENTS.md](../visualization/AGENTS.md)
