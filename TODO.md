# MetaInformAnt TODO

> **Last Updated**: 2026-10-27
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

- [ ] Raise `COVERAGE_MINIMUM` (currently 71.80%, about one point under the
  hosted full-suite measurement of 72.84% in Test Suite run 40446465878) in measured
  increments as maintained-domain tests expand; each increase requires a hosted
  full-suite receipt and a corresponding test/documentation change.

## Methods

- [ ] Make `visualization.genomics.trees.plot_phylo_tree` draw the trees that
  `dna.phylogeny.neighbor_joining_tree` returns. That builder represents
  leaves as `None` in its nested dict, and the plot raises `AttributeError`.
- [ ] Port `dna.external.genomes.list_genome_assemblies`,
  `get_genome_metadata` and `download_reference_genome` to the current NCBI
  Datasets API. A documentation review saw the v2 endpoints they call return
  404, so they come back empty; `download_genome_package` works. Verify with
  network access.
- [ ] Review the calibration of the `multi_trait_coloc` shared-signal heuristic.
- [ ] Choose a default for `multiomics.canonical_correlation` regularization:
  with none, correlations are exactly 1 whenever features outnumber samples.
- [ ] Add recombination to the coalescent simulator (each sample has a single
  genealogy), so `generate_linkage_disequilibrium_data` can honour
  `recombination_rate`; replace the heuristic expansion and bottleneck models
  of `generate_site_frequency_spectrum` and the statistics-dict forms of the
  demographic simulators; give `detect_selection` and
  `calculate_summary_statistics` the simulated null the neutrality tests have.
- [ ] Use a first-passage-time series (Navarro & Fuss 2009) in
  `ddm_log_likelihood` and `fit_ddm_parameters`, which approximate the
  decision-time density by a normal distribution.

## Known defects

Reproduced while fixing the previous defect backlog and not yet fixed.

- [ ] `pharmacogenomics`: unknown genes fall back to the CYP2D6 bands (decide
  whether to raise); unknown alleles score 1.0; an unscored `Diplotype` reads
  as Poor Metabolizer; the CYP1A2 activity values have no CPIC source.
- [ ] `metabolomics` and `multiomics`: `read_mgf` drops a final block without
  `END IONS`; `extract_chromatogram` and `prepare_spectrum_plot` accept any
  `ppm_tolerance`; `active_module_detection` hard-codes permutation seed 42;
  `detect_batch_effects` uses `1 - cdf` where `sf` is exact.
- [ ] `simulation` and `math`: `simulate_sequence_evolution` applies one
  mutation per generation whatever `mutation_rate` is; the list form of
  `effective_size_from_family_size_variance` (k/(V-1)) is not a standard
  estimator.
- [ ] `metagenomics`: the SciPy-free Welch test is a normal approximation;
  `biomarker_discovery` reports `cv_accuracy=0.0` when cross-validation cannot
  run; NMDS has no rank step and reports `1 - stress` as variance explained;
  `refine_bins` filters split bins on completeness only; PCoA ordination plots
  use power iteration; annotation confidence is `score/100`; scaffolds do not
  infer contig orientation.
- [ ] `structural_variants`: `merge_callsets`, `survivor_merge`,
  `deduplicate_variants` and `filter_by_frequency` ignore a translocation's
  `chrom2`/`pos2`, and `predict_tad_disruption` checks only its first
  breakpoint; population PCA eigenvalues are rescaled when samples are fewer
  than SVs; `find_nearest_gene` measures from the midpoint and rebuilds its
  index per call; genotyping ignores the caller's `min_mapq`;
  `detect_microhomology` is a mirror-match heuristic.
- [ ] `longread`: `PacBio_barcoded` in the adapter table is the Illumina
  Nextera read-1 sequence, and full-length LSK110 adapters need verified
  sequences; nested adapters can claim a shorter adapter's site; MM-tag
  decoding does not re-orient `-` strand or reverse-aligned reads, misaligns
  multi-code entries such as `C+mh`, and skips `.` flags; float read lengths
  are skipped; `variant_haplotypes` keys do not survive JSON.
- [ ] `epigenome`: CpG-island edges resolve only to the window size; per-CpG
  site lists and the methylation-versus-peak association do not scale to
  whole-genome bisulfite data; `differential_peaks` pairs each peak with its
  first overlap only; `learn_chromatin_states` fixes its seed.
- [ ] `singlecell` and `spatial`: `compute_diffusion_map` doubles mutual kNN
  weights and `diffusion_map_reduction` takes `eigh` of a non-symmetric
  matrix; Leiden ignores `use_weights`; PAGA and Slingshot are approximations;
  `filter_cells` leaves `uns["neighbors"]` stale for consumers other than DPT;
  `fetch_atlas_datasets` returns an empty directory for ArrayExpress; the
  Visium loader keeps expression rows with no coordinates when no barcodes
  match; spatial Leiden/Louvain fall back to greedy modularity invisibly.
- [ ] Decide whether `core.utils.errors.ValidationError` should also subclass
  `ValueError`: `singlecell` raises a bare `ValueError` for some invalid
  arguments and `ValidationError` for others, so `except ValueError` catches
  only some (`tests/singlecell/test_singlecell_docstring_exceptions.py` pins the
  current hierarchy).
- [ ] `networks` and `information`: `integration.networks._von_neumann_entropy`
  uses the row-normalised adjacency, which is not a density matrix; the GRN
  hypergeometric p-value falls back to a normal approximation without SciPy.
- [ ] `gwas`, `protein` and `eqtl`: the summary-statistics branch of
  `download_variant_data` needs network verification (numeric p-value
  formatting, dict-valued variant column, silent non-200 responses, 1000-row
  truncation); regional plots raise `KeyError` for records with no chromosome
  key; `genome.manhattan_plot` and `effect_direction_plot` accept every
  chromosome spelling but only `POS`/`pos` and `p_value`/`P`, unlike the regional
  plots; `suite.py` names regional files from upcast floats; the protein
  secondary-structure, SASA, domain and similarity functions remain simplified.

## Documentation

- [ ] Audit the module cards and runtime figures in
  `docs/comparisons/methods_matrix.md` and `docs/CAPABILITY_MATRIX.md`, and
  the source trees in `docs/agents/rules/*.md` (only the GWAS tree has been
  corrected). Imports, attribute chains, call signatures and relative links in
  active documentation are checked by `tests/quality/test_doc_imports.py` and
  `tests/quality/test_doc_calls_and_links.py`; prose claims are not.
- [ ] Review prose claims (described return values, quoted outputs,
  capability statements) in pages no domain sweep has covered:
  `docs/cloud/SPEC.md`; `docs/visualization/` animations, dimred, expression,
  trees and plots; `docs/ml/features.md`; `docs/dna/sequences.md`;
  `docs/ecology/` ordination and visualization; and `docs/gwas/` structure,
  index and ASSESSMENT.

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
