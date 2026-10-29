# AGENTS.md — `MetaInformAnt/tests/ontology/`

tests for the ontology module.
Files (verified 2026-10-28): `__init__.py`, `test_get_subontology_hierarchy.py`,
`test_go_enrichment_alpha_and_fallbacks.py`, `test_go_enrichment_background.py`, `test_gsea_permutation_null.py`,
`test_obo_synonyms.py`, `test_ontology_annotation_regression.py`, `test_ontology_api.py`,
`test_ontology_comprehensive.py`, `test_ontology_enrichment.py`, `test_ontology_go_basic.py`,
`test_ontology_json_roundtrip.py`, `test_ontology_layout_and_workflow_seeds.py`, `test_ontology_obo_definitions.py`,
`test_ontology_obo_parser.py`, `test_ontology_query.py`, `test_ontology_query_exports.py`,
`test_ontology_query_semantics.py`, `test_ontology_query_subsets_and_cache.py`,
`test_ontology_relation_roots_leaves.py`, `test_ontology_seeded_randomness.py`, `test_ontology_serialization.py`,
`test_ontology_serialize.py`, `test_ontology_types.py`, `test_ontology_visualization.py`,
`test_ora_correction_names.py`, `test_ora_exact_log_factorials.py` (26 test files plus `__init__.py`).

## Coverage notes
- `test_ontology_obo_definitions.py`: `def:` parsing (quoted text, dbxref list, escapes, trailing modifiers, malformed values), `metadata["definition_xrefs"]`, `[Typedef]` definitions, and OBO save/load round trips of definitions and dbxrefs.
- `test_obo_synonyms.py`: `parse_obo_synonym` (text, scope, dbxrefs, escapes, a `!` inside the quotes, synonym types, malformed values), synonyms read whole from `[Term]` stanzas, and OBO save/load round trips of synonyms written with `quote_obo_string`.
- `test_ontology_query_exports.py`: every public function of `query.py` and `serialize.py` is re-exported by the `metainformant.ontology.query` package and listed in `__all__`; named imports, star imports and the submodules stay importable.
- `test_ontology_query_semantics.py`: `ancestors`, `descendants` and `common_ancestors` exclude the queried term itself, for `is_a` and other relation types; the docstring examples run and give the documented results.
- `test_ontology_query_subsets_and_cache.py`: `filter_by_namespace` and `subgraph` carry the induced parent/child links (roots, parents and ancestors checked against an independent `networkx` digraph); the `ancestors`/`descendants` cache is tied to its ontology by weak reference and structure stamp, refreshes on `add_term`/`add_relationship` and disappears with the ontology.
- `test_get_subontology_hierarchy.py`: `get_subontology` keeps only the hierarchy inside the sub-ontology (a seed whose parents were left out is a root, `get_ancestors` stays inside, no dangling or duplicated `is_a` relationships) and accepts a generator for `root_terms`.
- `test_ontology_json_roundtrip.py`: `load_ontology(format="json")` restores OBO-parsed ontologies (the `relationships` metadata key) and does not duplicate `is_a` links; the loaded ontology equals the saved one.
- `test_go_enrichment_background.py`: `enrich_genes` counts only the annotated genes that are in the background, so `"fisher"` and `"hypergeometric"` agree (oracles `scipy.stats.fisher_exact` and `hypergeom.sf`).
- `test_go_enrichment_alpha_and_fallbacks.py`: `alpha` sets the boolean `significant` column and every tested term is returned; without scipy (the real module runs in a subprocess that cannot import it) the p-value is the exact upper tail; the `semantic_similarity` score ranges (Resnik in bits, Lin and Jiang-Conrath in [0, 1]).
- `test_ora_exact_log_factorials.py`: the log-factorials behind `over_representation_analysis` equal `log(n!)` (oracle `math.factorial`) and the p-values match `scipy.stats.hypergeom.sf`.
- `test_gsea_permutation_null.py`: `gsea` uses the preranked gene-set permutation null (exact-rational ES oracle, replayed seeded draws for p-value and NES, calibration of p-values over random gene sets).
- `test_ontology_seeded_randomness.py`: `gsea(random_state=...)` and `pathway_network_plot(layout_seed=...)` reproduce with the same seed and leave the global random state alone.
- `test_ontology_layout_and_workflow_seeds.py`: `layout_seed` of `plot_go_dag`, `plot_ontology_network` and `create_interactive_go_network` (checked through the node positions matplotlib draws; the Plotly function through the source), and the `ontology.random_state` configuration key of `run_ontology_analysis` (read by `_gsea_random_state` and handed to `gsea`; the workflow itself needs QuickGO and is not run).

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/ontology`.
  `test_ontology_api.py` probes internet connectivity when it is imported and queries public services, so it was left
  out of the check on 2026-10-27: `.venv/bin/python -m pytest -q tests/ontology --ignore=tests/ontology/test_ontology_api.py
  -p no:cacheprovider -p no:rerunfailures` gave 373 passed, 1 skipped (the skip needs plotly).
Repo-wide policy: see the repository-root `AGENTS.md`.
