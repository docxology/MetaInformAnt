# MetaInformAnt TODO

> **Last Updated**: 2026-10-26
>
> This TODO contains only active unfinished work. Completed work is recorded in
> `CHANGELOG.md`, the review and benchmark documents it references, and git
> history, not in this list.

The Hymenoptera Amalgkit project has a separate, data-root-scoped backlog at
[`projects/hymenoptera_amalgkit/TODO.md`](projects/hymenoptera_amalgkit/TODO.md).
Keep that work there rather than mixing a live external-volume transcriptome
campaign with the BeeWAS release backlog below.

## Release engineering

Remaining work is gated in order: engineering checks; nested publication and
reachable gitlinks; operational manifests, receipts, locks, and campaign
snapshots; then owner-reviewed scientific promotion. A green engineering gate
never establishes cohort completion or biological inference. The live RNA
campaign remains a moving operational snapshot.

- [ ] Raise `COVERAGE_MINIMUM` (currently 70.00%, one point under the hosted
  full-suite measurement of 71.02% in Test Suite run 40240714837) in measured
  increments as maintained-domain tests expand; each increase requires a hosted
  full-suite receipt and a corresponding test/documentation change.

## Methods

- [ ] Replace the heuristic sequence generators in `simulation.models.popgen`
  with coalescent simulation. `generate_population_sequences` mutates every
  sequence independently from one ancestor (a star genealogy), so a "neutral"
  sample of 20 gives Tajima's D near -2.6. The bottleneck and expansion
  scenarios resample founders, which leaves the bottleneck's D near zero and
  pushes expansion's D positive, the opposite of the textbook signature.
- [ ] Make `generate_two_populations(fst=...)` reach its target: Hudson's
  F_ST of its output is near 0.5 for targets 0.05, 0.2 and 0.5.
- [ ] Offer simulated null distributions for the neutrality tests. The
  current p-values use the normal approximation, which is rough for small
  samples.
- [ ] Make `visualization.genomics.trees.plot_phylo_tree` draw the trees that
  `dna.phylogeny.neighbor_joining_tree` returns. That builder represents
  leaves as `None` in its nested dict, and the plot raises `AttributeError`.
- [ ] Port `dna.external.genomes.list_genome_assemblies`,
  `get_genome_metadata` and `download_reference_genome` to the current NCBI
  Datasets API. A documentation review saw the v2 endpoints they call return
  404, so they come back empty; `download_genome_package` works. Verify with
  network access.

- [ ] Pass genotypes to `write_summary_statistics` from the GWAS workflow
  (`workflow_execution.py`) so pipeline outputs carry MAF instead of `NA`.
- [ ] Accept the `chrom` key in `gwas/visualization/genomic/genome.py` and
  `statistical/effects.py`, as `interactive_manhattan` now does, and review
  the calibration of the `multi_trait_coloc` shared-signal heuristic.

## Known defects from the documentation review

Found while making the documentation executable; each was reproduced but is
not yet fixed.

- [ ] `simulation.models.rna` treats negative fold changes as magnitudes,
  and `generate_linkage_disequilibrium_data` ignores `recombination_rate`
  (its `r_squared_target` acts as a copy probability).
- [ ] `math.population_genetics`: `effective_size_from_family_size_variance`
  asserts on an integer size, and `effective_size_sex_ratio` jumps from 400
  to 200 at a sex ratio of exactly 0.5.
- [ ] `core.io.download`: `resume=False` over an existing file double-counts
  bytes, a missing destination directory is treated as a file name, and a
  complete file with `resume=True` fails with HTTP 416.
- [ ] `dna.variation.mutations.generate_point_mutations` uses the global
  random state and ignores `mutation_rate`.
- [ ] `networks`: pathway IDs containing underscores are split, an unknown
  enrichment method runs as hypergeometric, Leiden-only arguments reach the
  Louvain fallback, `granger`/`grnboost2` are listed but raise, feed-forward
  loops ignore the direct edge, and the PPI "top 10" centrality is the first
  ten nodes.
- [ ] `ml`: `compare_validation_strategies`' default strategy always errors,
  and `OllamaConfig(max_retries=0)` crashes `generate()`.
- [ ] `ontology`: `filter_by_namespace` and `subgraph` lose the parent/child
  maps, `enrich_genes` ignores `alpha`, and Resnik similarity is documented as
  0-1 but returns bits.
- [ ] Visualization: duplicate `cmap`/`marker` and unsupported
  `color`/`alpha` arguments raise in `analysis/information.py`,
  `analysis/timeseries.py` and `plot_circular_barplot`; Sankey node indices,
  per-chromosome marker sizes, hierarchical network layout, sequence-logo
  heights and single-panel `multi_panel` are wrong.
- [ ] `dna.population` plots ignore `output_path` in
  `plot_population_diversity`, return `None` from `plot_pca_results` for a
  bare array, and drop summary and neutrality entries beyond six and four.
- [ ] GWAS plots: genome-wide x positions are inflated (one scatter per
  variant), effect-plot positions are a running sum, the interactive suite
  skips TSV inputs and ignores `vcf_file`, regional plots compare `chrom` by
  exact type, and QQ "expected" values are random and fail at p = 0.
- [ ] `pharmacogenomics`: the CYP2D6, CYP2C19 and SLCO1B1 phenotype bands
  contradict CPIC (`*1/*4` reads Normal, CYP2C19 `*1/*1` Rapid), the
  metabolizer module uses a second, inconsistent table, and
  `generate_clinical_report(drugs=None)` returns no recommendations.
- [ ] `longread`: the default adapter table holds two 8 bp entries that match
  random sequence (about 190 bp trimmed from adapter-free 3 kb reads),
  `haplotag_reads` tags every read HP=0 for `phase_reads` blocks, the
  "bayesian" methylation model is a capped linear map, and merged SV calls
  lose mate and inserted-sequence fields.
- [ ] `structural_variants`: translocation calls take start and end from two
  chromosomes, insertion calls never reach the default minimum size,
  breakpoint refinement drops `chrom2` and strand, gene-overlap queries miss
  long genes, and `dosage_gain` is unreachable.
- [ ] `metagenomics`: bin quality does no marker search (any 3 kb bin scores
  complete), the OTU k-mer prefilter rejects pairs above the identity
  threshold, ASV denoising inverts its Poisson test, Fisher's alpha and
  pathway completeness are wrong, and functional annotation invents labels.
- [ ] `metabolomics`: `read_csv` turns unparseable cells into 0.0, unknown
  `ion_mode` selects the negative adducts, enrichment p-values fall to 0.0 on
  error, and ComBat correction depends on batch order.
- [ ] `epigenome`: peak files do not survive a save/load round trip,
  overlapping-peak detection misses pairs, `find_tf_binding_sites` is
  simulated, and the ChIP/ATAC integration step is a placeholder.
- [ ] `singlecell`: `leiden_clustering` silently runs k-means on plain
  arrays, several steps mutate their input, the `run_*` I/O helpers only touch
  output files, and diffusion pseudotime is effectively Euclidean distance.
- [ ] `spatial`: the `radius` graph method always raises, mutual kNN weights
  double, Xenium and MERFISH scatter plots swap axes, and loaders fall back to
  zero coordinates when the cell file is missing.
- [ ] Smaller contract bugs: `biological_embedding("temporal")` raises for
  any `EventSequence`, `event_importance("frequency")` ranks by earliness,
  `apply_qc_filters` raises `IndexError` when variant and genotype counts
  differ, `regional_coloc` reports success when every window errors,
  `cross_validate_biological(method="ensemble", random_state=...)` raises,
  and `load_ontology(format="json")` fails on OBO-parsed ontologies.
- [ ] Packaging: the `louvain` extra names the wrong package, and
  `np.trapezoid` in metabolomics plots needs NumPy 2 while `numpy>=1.24` is
  allowed. `canonical_correlation` defaults to no regularization, which gives
  correlations of exactly 1 when features outnumber samples.

## Documentation

- [ ] Audit the module cards and runtime figures in
  `docs/comparisons/methods_matrix.md`, and the source trees in
  `docs/agents/rules/*.md` (only the GWAS tree has been corrected). Imports
  and attribute chains in active documentation are checked by
  `tests/quality/test_doc_imports.py`; call signatures and prose are not.
- [ ] Bring call signatures and claims up to date in pages the executable
  sweep did not cover: `docs/cloud/SPEC.md`; `docs/visualization/`
  animations, dimred, expression, trees, plots, ontology, simulation, SPEC and
  AGENTS; `docs/ml/features.md` and `VALIDATION_REPORT.md`;
  `docs/dna/distances.md`, `translation.md` and `sequences.md`;
  `docs/ontology/query.md`; `docs/ecology/` ordination and visualization;
  `docs/gwas/` workflow, structure, index and ASSESSMENT; and the
  capability, performance, troubleshooting and topic pages the omics sweep
  left untouched in `docs/{metabolomics,metagenomics,singlecell,spatial,
  pharmacogenomics,structural_variants,longread,epigenome}/`.

---

## Current Scope

BeeWAS four-population statistical GWAS reporting is operational for the
combined cohort and population scopes `C`, `I`, `M`, and `R`. Manuscript
rendering, visual QA, artifact indexing, guarded dashboard generation, GO
curation-atlas outputs, statistical deep-dive summaries, and release packaging
are generated and validated. The active work below is evidence-gated: it
requires human-approved metadata, human-reviewed GO mapping decisions, or
independent replication evidence before the repository can mark the item
complete.

Use the package-module sweep to refresh status after any evidence is supplied:

```bash
cd projects/apis_gwas
uv run python -m beewas.gwas.review_release build
```
The full `apis_gwas` tree lives in the nested checkout on the producer
machine; this checkout carries only a README/AGENTS stub, so this command
runs where that tree exists.

---

## Minor Updates
> **Note**: `projects/apis_gwas/...` paths in the items below are references
> to the full BeeWAS tree in the nested `apis_gwas` checkout on the producer
> machine (external absolute-path targets). This checkout holds only a
> README/AGENTS stub; the TODO items themselves remain active.

- [ ] **Finalize BeeWAS manuscript metadata** - Complete
  `projects/apis_gwas/docs/manuscript/renderable/final_metadata_intake.tsv`
  (to be created; the apis_gwas tree is currently a README/AGENTS stub)
  with final title, authors, affiliations, corresponding author, funding,
  contributions, conflicts, data availability, code availability, reviewer
  name, review date, and final status values. Current status is tracked in
  `projects/apis_gwas/output/pdf/beewas_gwas/manuscript_metadata_status.tsv`.

## Medium Improvements

- [ ] **Review BeeWAS GO mapping evidence** - Review
  `projects/apis_gwas/doc/operations/reviews/gwas_reviewed_go_mapping_intake_template.tsv`,
  set justified rows to `review_status=reviewed`, retain provenance fields, and
  ingest the reviewed TSV with `uv run python -m beewas.gwas.reviewed_go ingest
  path/to/reviewed_go.tsv`. Current intake status is tracked in
  `projects/apis_gwas/results/gwas_79cram_clean_reviewed_fullgemma_guarded_20260624_summary/tables/gwas_reviewed_go_mapping_intake_status.tsv`.

## Major Improvements

- [ ] **Supply independent replication evidence** - Complete
  `projects/apis_gwas/doc/operations/reviews/gwas_replication_review_scaffold.tsv`
  with external cohort details, comparable trait definitions, model
  compatibility, evidence URIs, reviewer decisions, reviewer names, and review
  dates before any replication claim is made.

---

---

## Operating Rules

- Use package-module commands for BeeWAS work, for example
  `uv run python -m beewas.gwas.guarded ...`.
- Keep runtime outputs under `output/`, `results/`, or documented data/output
  directories.
- For the live Hymenoptera cloud campaign, reattach an interrupted controller
  with `--skip-bundle-upload` and the original state/source/input IDs; never
  rebuild or overwrite an immutable source object under an existing run ID.
- Do not convert curation evidence into biological, causal, gene-function,
  replication, breeding, or final heritability claims without a reviewed
  decision artifact.
- Keep this TODO forward-looking; move resolved evidence and dated snapshots to
  the appropriate status or release document.
- Do not edit nested repository contents from the parent. Reconcile nested
  gitlinks only after their own branches/PRs contain validated, remotely
  reachable commits.
- Do not touch `/Volumes/external_drive/Data/amalgkit` (AMALGKIT_DATA_ROOT) except for explicitly read-only
  monitoring.
- Do not reset, discard, recreate the active environment, or claim
  campaign/scientific completion from a moving operational snapshot.
- Before each release claim, recapture parent/nested commits, reachability,
  dirty paths, PR state, campaign status, current manifest, provenance, and
  reproducible evidence-bundle status.
