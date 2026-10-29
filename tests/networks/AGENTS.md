# AGENTS.md — `MetaInformAnt/tests/networks/`

tests for the networks module.
Files (verified 2026-10-28): `__init__.py`, `test_community_fallback_kwargs.py`,
`test_community_seeds_and_comparison.py`, `test_graph_algorithms.py`, `test_graph_core_directed_edges.py`,
`test_grn_hypergeometric_exact.py`, `test_grn_motifs_self_loops.py`, `test_grn_tf_list.py`,
`test_grn_validation_ties.py`, `test_information_flow.py`, `test_louvain_global_state.py`,
`test_motif_binomial_exact.py`, `test_network_summary_fallback.py`, `test_networks_community.py`,
`test_networks_community_detection.py`, `test_networks_comprehensive.py`, `test_networks_graph.py`,
`test_networks_graph_algorithms_depth.py`, `test_networks_grn_inference_depth.py`, `test_networks_pathway.py`,
`test_networks_pathway_depth.py`, `test_networks_ppi.py`, `test_networks_ppi_impl_depth.py`,
`test_networks_regulatory.py`, `test_networks_seeded_randomness.py`, `test_pathway_correction_names.py`,
`test_pathway_cosine_empty_sets.py`, `test_pathway_enrichment_fallbacks.py`, `test_pathway_fdr_correction.py`,
`test_pathway_hypergeometric.py`, `test_ppi_centrality_and_json.py`, `test_ppi_functional_enrichment_population.py`,
`test_ppi_network_enrichment_exactness.py`, `test_regulatory_feed_forward_loops.py`,
`test_regulatory_mutual_info_seed.py` (34 test files plus `__init__.py`).

## Coverage notes
- `test_pathway_enrichment_fallbacks.py`: `pathway_enrichment` computes the exact upper-tail p-value (oracles
  `scipy.stats.hypergeom.sf` and `fisher_exact`) and applies `correction="fdr"` when scipy cannot be imported (the
  real module runs in a subprocess without scipy); an unknown `method` raises `ValueError`; `create_pathway_network`
  keeps pathway IDs that contain underscores.
- `test_pathway_cosine_empty_sets.py`: a pathway without genes has similarity 0.0 with every pathway under the
  Jaccard, overlap and cosine methods, in `pathway_similarity_analysis` and `create_pathway_network`; the cosine
  oracle is NumPy over 0/1 indicator vectors.
- `test_community_fallback_kwargs.py`: the Leiden-to-Louvain fallback of `detect_communities` keeps only the options
  Louvain accepts (a signature-level check, plus an end-to-end run where python-louvain is installed and leidenalg is
  not).
- `test_louvain_global_state.py`: `louvain_communities` hands its seed to python-louvain and does not seed the global
  `random` / `numpy.random` generators (a static check of the source, plus an end-to-end run where python-louvain is
  installed).
- `test_networks_seeded_randomness.py`: `compute_network_motifs` and `analyze_regulatory_dynamics` take
  `random_state`, the same seed reproduces the output, and the global random states are neither seeded nor consumed.
- `test_ppi_centrality_and_json.py`: the centrality section of `ppi_network_analysis` holds the ten most central
  proteins, most central first; `save_ppi_network(format="json")` output loads through `load_ppi_network` and
  `ProteinNetwork.from_file` (including the `links` key of older networkx releases).
- `test_regulatory_feed_forward_loops.py`: a feed-forward loop needs the direct TF1 -> target edge in
  `regulatory_motifs` and `GeneRegulatoryNetwork.regulatory_motifs`; the method lists in the docstrings of `infer_grn`
  and `infer_regulatory_network_from_expression` name only methods that run.
- `test_regulatory_mutual_info_seed.py`: `random_state` of the `'mutual_info'` method of `infer_grn` and
  `infer_regulatory_network_from_expression` (one generator shared by all gene pairs), checked against scikit-learn's
  `mutual_info_regression`.

## Conventions
- Real implementations with small deterministic data; the lexical no-mocks gate
  applies (no `MagicMock`/`unittest.mock`).
- Run: `env -u VIRTUAL_ENV uv run pytest -q tests/networks`.
  Checked 2026-10-27 with `.venv/bin/python -m pytest -q tests/networks -p no:cacheprovider -p no:rerunfailures`:
  697 passed, 10 skipped (the skips need python-louvain, or leidenalg and igraph, which this environment lacks).
Repo-wide policy: see the repository-root `AGENTS.md`.
