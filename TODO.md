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

## Workflowz review — 2026-11-10

- [ ] **WFZ-20261110-01 — P1: Preserve normalized frequencies in `replicator_step`.** Observable outcome: updating a valid non-normalized frequency vector uses the same normalized state as the derivative calculation and returns the correct replicator update instead of rescaling the derivative by the input sum. Evidence: `src/metainformant/math/evolutionary_dynamics/core.py:100-103` normalizes only a local NumPy array, while `src/metainformant/math/evolutionary_dynamics/core.py:126-132` advances the original list. Acceptance probe from repository root: `uv run python -c 'from metainformant.math.evolutionary_dynamics.core import replicator_step; import numpy as np; actual=replicator_step([2.0,2.0],np.array([[1.0,0.0],[0.0,2.0]]),dt=0.1); expected=[0.45,0.55]; assert np.allclose(actual,expected),(actual,expected)'`. Dependencies/constraints: none; preserve the existing normalization convention and public signature.
- [ ] **WFZ-20261110-02 — P1: Reject an undefined all-zero replicator state.** Observable outcome: an all-zero frequency vector raises a clear `ValueError` rather than producing NaN derivatives and a NaN-filled next state. Evidence: `src/metainformant/math/evolutionary_dynamics/core.py:100-103` divides by the zero sum, and `src/metainformant/math/evolutionary_dynamics/core.py:129-133` leaves the resulting non-finite values unchanged when normalization is skipped. Acceptance probe from repository root: `uv run python -c 'from metainformant.math.evolutionary_dynamics.core import replicator_derivative,replicator_step; import numpy as np; A=np.zeros((2,2));
for fn in (replicator_derivative,replicator_step):
 try: fn([0.0,0.0],A)
 except ValueError: pass
 else: raise AssertionError(f"{fn.__name__} accepted zero total frequency")'`. Dependencies/constraints: coordinate with WFZ-20261110-01 because both change the same functions; retain normalization for valid nonzero inputs and avoid changing the documented return type.
- [ ] **WFZ-20261110-03 — P2: Enforce the dominant mutation-selection equilibrium domain.** Observable outcome: nonpositive selection, negative mutation/selection inputs, and non-finite results are rejected explicitly; the zero-selection case no longer escapes as a raw `ZeroDivisionError`, and estimates above frequency 1 are not returned as valid probabilities. Evidence: `src/metainformant/math/population_genetics/core.py:92-106` documents q=μ/s but returns the division without validation, unlike the explicit domain checks at `src/metainformant/math/population_genetics/core.py:57-62` and `src/metainformant/math/population_genetics/core.py:79-82`. Acceptance probe from repository root: `uv run python -c 'from metainformant.math.population_genetics.core import mutation_selection_balance_dominant as f; import math; assert f(0.01,0.1)==0.1;
for args in ((0.01,0.0),(0.01,-0.1),(-0.01,0.1),(0.2,0.1),(float("nan"),0.1)):
 try: f(*args)
 except ValueError: pass
 else: raise AssertionError(f"accepted invalid equilibrium inputs {args}")'`. Dependencies/constraints: this tightens degenerate-input behavior already called out for review in `TODO.md:57-58`; decide whether μ/s>1 is invalid input or a documented saturated boundary before implementing, and keep valid finite values unchanged.
- [ ] **WFZ-20261110-04 — P2: Enforce the recessive mutation-selection equilibrium domain.** Observable outcome: negative inputs, nonpositive selection, NaN, and equilibrium frequencies above 1 produce deliberate `ValueError`s rather than raw division/domain errors or invalid probabilities; valid calls still return sqrt(μ/s). Evidence: `src/metainformant/math/population_genetics/core.py:109-123` directly evaluates `sqrt(μ/s)` without validating its biological domain, despite nearby explicit validation at `src/metainformant/math/population_genetics/core.py:57-62`. Acceptance probe from repository root: `uv run python -c 'from metainformant.math.population_genetics.core import mutation_selection_balance_recessive as f; import math; assert math.isclose(f(0.01,0.04),0.5);
for args in ((0.01,0.0),(0.01,-0.1),(-0.01,0.1),(0.2,0.1),(float("nan"),0.1)):
 try: f(*args)
 except ValueError: pass
 else: raise AssertionError(f"accepted invalid equilibrium inputs {args}")'`. Dependencies/constraints: align the boundary policy and exception messages with WFZ-20261110-03; do not silently clamp results, and update the docstring to state the selected domain.
