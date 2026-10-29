# MetaInformAnt TODO

> **Last Updated**: 2026-10-29
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

- [ ] Raise `COVERAGE_MINIMUM` (currently 74.90%, about one point under the
  hosted full-suite measurement of 75.91% in Test Suite run 40594851520) in measured
  increments as maintained-domain tests expand; each increase requires a hosted
  full-suite receipt and a corresponding test/documentation change.

## Methods

- [ ] Recheck live NCBI paging consistency and contact-configured sequence
  reports. Genome/annotation downloads now have live data checks; duplicate
  assembly pages are rejected, and upstream inconsistencies remain explicit.
- [ ] Validate multi-trait colocalization calibration and prior sensitivity
  beyond its independent-study, single-causal-variant model.

## Remaining limitations

Resolved defect entries have moved to the changelog. The following work
still needs implementation, source evidence or an explicit API decision.

- [ ] Trace the built-in CYP1A2 activity values to primary sources. Phenotype
  calls for unsupported genes now raise; unknown alleles are Indeterminate.
- [ ] Replace the short LSK110 adapter motifs with verified full-length
  sequences before advertising full adapter coverage.
- [ ] Assess whole-genome scaling of per-CpG site lists and windowed
  methylation entropy. Peak association now uses an interval index.
- [ ] Separate BEDgraph methylation levels from measured read counts before
  using reconstructed counts for inferential DMR tests.
- [ ] Implement an ArrayExpress atlas downloader; unsupported accessions now
  raise instead of returning an empty dataset.
- [ ] Replace or explicitly scope the simplified Slingshot-style trajectories,
  heuristic trajectory pseudotime, and sequence-based secondary-structure
  predictions. Structure-based secondary-element assignment and structural
  similarity now have explicit algorithms and parameters.
- [ ] Review pharmacogenomic CNV and ambiguous-diplotype decision rules,
  default dose adjustments for unknown drugs, and built-in annotation-table
  provenance before clinical use.
- [ ] Review remaining fixed RNG seeds, degenerate-input conventions and
  approximate statistical fallbacks against their documented contracts.

## Documentation

- [ ] Audit the module cards and runtime figures in
  `docs/comparisons/methods_matrix.md` and `docs/CAPABILITY_MATRIX.md`, and
  the source trees in `docs/agents/rules/*.md` (only the GWAS tree has been
  corrected). Imports, attribute chains, call signatures and relative links in
  active documentation are checked by `tests/quality/test_doc_imports.py` and
  `tests/quality/test_doc_calls_and_links.py`; prose claims are not.
- [ ] Review prose claims (described return values, quoted outputs,
  capability statements) in pages no domain sweep has covered:
  `docs/visualization/` animations, dimred, expression, trees and plots;
  `docs/ml/features.md`;
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
